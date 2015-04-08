#!/usr/bin/env python
"""
TODO:
-similar case with other parameters
-it would save disk space if we do not add go categories with 0 genes
to all assoc_tables, but only add them in the end. (I think the most recent version supports
this, try it...)
-Package functions so that we can run enrichment on candidate gene list directly (without a rod_s).
-Package into classes?
"""
from __future__ import print_function 
import os
import pandas as pd
import numpy as np
import gc, logging
from hs_vervet.tools import hs_pandas as hp
eu = os.path.expanduser
jn = os.path.join

def make_df(df_or_array,new_columns):
    try:
        old_columns = tuple(df_or_array.columns)
        df = pd.DataFrame(df_or_array.values, columns=new_columns)
    except AttributeError:
        old_columns = (None for _ in new_columns)
        df = pd.DataFrame(df_or_array, columns=new_columns)
    return df, old_columns


def init_rank_table(assoc):
    rt = pd.DataFrame({"value":assoc.values,"rank":0,"out_of":0},index=assoc.index)
    rt.index.name = 'category'
    return rt

def shift_rod(rod_df, rnd, mode = "grid"):
    """
    shift reference ordered data across the whole genome

    Input:
    rod_df ... pandas dataframe or series with mulitiindex (chrom, pos)

    modes ...
        'grid' ... just rotate the index of the rod data frame
                    this means that the positions stay the same only the
                    value for each position becomes different
                    Faster, but means that you only hit the same grid-point
                    this should make it conservative on large grids. Large
                    grids are problematic if the fraction of top windows 
                    considered becomes large 
        'continuous' ... add the random shift to each index value.
                         NOT IMPLEMENTED
    """
    if mode == "grid":
        new_start_i = int(len(rod_df)*rnd)
        rotate_data = np.concatenate((rod_df.iloc[new_start_i:].values,rod_df.iloc[:new_start_i].values))
        if  isinstance(rod_df,pd.core.series.Series):
            r = pd.Series(rotate_data,index=rod_df.index)
            return r
        elif isinstance(rod_df,pd.core.frame.DataFrame):
            r = pd.DataFrame(rotate_data,index=rod_df.index,columns=rod_df.columns)
            return r
    else:
        raise UserException("Only mode grid supported.")


class CandidateEnrichment(object):
    """
    Test enrichment of a candidate feature list
    against a feature to category mapping.
    In a common example, features will be genes 
    and categories will gene ontology (GO) categories.
    
    Input:
    candidate_features ... list of feature names to be tested
    feature_to_category ... data frame with arbitrary index
    """

    def __init__(self, candidate_features, feature_to_category):

        self.candidate_features = np.unique(candidate_features)
        if len(self.candidate_features) < len(candidate_features):
            logging.warning("There are duplicates in candidate_features, "
                                                    "going to be removed...")

        (ftc,(self.feature_name, self.category_name)) = make_df(feature_to_category,
                                                                ['feature', 'category'])

        n_candidates = len(self.candidate_features)
        self.all_features = ftc['feature'].unique()
        self.candidate_features = np.intersect1d(self.candidate_features,
                                                 self.all_features, assume_unique=True)
        assert len(self.candidate_features) > 0, ("No candidate feature found in the "
                                   "first column of feature_to_category.")
        if len(self.candidate_features) < n_candidates:
            logging.warning("Only {} of the {} candidates are present "
                            "in feature_to_category mapping. "
                            "Only those are going to be used.".format(len(self.candidate_features), n_candidates))
        self.feature_to_category = ftc
        self.make_real_assoc_rank_table(self.candidate_features)


    def get_association(self,candidate_features):
        """
        Get series with number of candidate features associated with each
        category in feature_to_category.
        """
        assoc = self.feature_to_category.set_index("feature").ix[candidate_features].groupby("category").apply(len)
        assoc.name = "n_features"
        assoc.index.name = "category"
        return assoc

    def make_real_assoc_rank_table(self,candidate_features):
        real_assoc = self.get_association(candidate_features)
        rt = init_rank_table(real_assoc)
        self.rank_table = rt

    def permute(self,n_permut):
        """
        Update the supplied rank table (rt) with the 
        """
        rt = self.rank_table
        for i in xrange(n_permut):
            permut_features = np.random.choice(self.all_features,
                                            size=len(self.candidate_features), replace=False)
            assoc = self.get_association(permut_features)
            rt["rank"] += (rt["n_features"] > assoc.reindex(rt.index).fillna(0))
            rt["out_of"] += 1
            gc.collect() #needed for low memory profile
        rt.sort('rank',ascending=False,inplace=True)
        #return rt






class SummaryEnrichment(CandidateEnrichment):
    """
    summary ... name of the function of the groupby object
                         to apply to the data (e.g. 'mean', 'max',...)
    max_dist ... is not implemented yet!!!


    two scenarios: 
    -- summary across features + summary across categories
    -- summary across features -> take top features 
        (can be modelled by more complicated summary function, eg, lambda x: mean(x)>thresh)? 
         -> count features per cat (summary with sum)
    """
    def __init__(self, value_s, feature_df,  feature_to_category, feature_summary=None, feature_summary_fun=None,
                 category_summary=None,category_summary_fun=None, max_dist=0):
                #

        assert (feature_summary is not None) != (feature_summary_fun is not None), \
                                "Specify either feature_summary OR feature_summary_fun."
        assert (category_summary is not None) != (category_summary_fun is not None), \
                                "Specify either category_summary OR category_summary_fun."
        self.value_s = value_s.copy()
        self.feature_df = feature_df.copy()
        (self.feature_to_category, (self.feature_name, self.category_name)) = \
                                                         make_df(feature_to_category,
                                                                ['feature', 'category'])
        self.feature_summary = feature_summary
        self.category_summary = category_summary
        self.feature_summary_fun = feature_summary_fun
        self.category_summary_fun = category_summary_fun

    def real_assoc(self):
        self.make_real_assoc_rank_table(self.value_s)


    def get_association(self, value_s):
        summary_per_feature = self.get_summary_per_feature(value_s)
        assoc = self.get_summary_per_category(summary_per_feature)
        return assoc

    def get_summary_per_feature(self,value_s):
        values_per_feature = hp.data_per_feature(value_s,self.feature_df)
        #print(values_per_feature)
        if self.feature_summary is not None:
            groups = values_per_feature.groupby(lambda idx: idx[1])
            return getattr(groups, self.feature_summary)()
        elif self.feature_summary_fun is not None:
            return self.feature_summary_fun(values_per_feature)
        #summary_per_feature = hp.apply_to_feature(values_per_feature,
        #                                          groupby_func_name=self.feature_summary,
        #                                         function=self.feature_summary_fun)
        #return summary_per_feature


    def get_summary_per_category(self,value_per_feature):
        """
        Calculates summary (e.g. mean) of values for the features
        in each of the given categories.

        Returns:
        series
        """
        value_to_category = self.feature_to_category.copy()
        values_per_feature_to_category = value_per_feature.ix[value_to_category['feature'].values].values
        del value_to_category['feature']
        value_to_category['value'] =  values_per_feature_to_category
        if self.category_summary is not None:
            groups = value_to_category.groupby('category')
            print(getattr(groups, self.category_summary)())
            return getattr(groups, self.category_summary)()['value']
        elif self.category_summary_fun is not None:
            return self.category_summary_fun(value_to_category)['value']
        #summary_per_category = getattr(value_to_category.groupby('category'), self.category_summary)()
        #return summary_per_category['value']



    def permute(self, n_permut):
        rank_table = self.rank_table
        for rnd in np.random.rand(n_permut):
            s = shift_rod(self.value_s, rnd)
            assoc = self.get_association(s)
            rank_table["rank"] += (rank_table["value"] > \
                                            assoc.reindex(rank_table.index).fillna(0))
            rank_table["out_of"] += 1
            gc.collect() #needed for low memory profile
        rank_table.sort('rank',ascending=False,inplace=True)


class TopScoresEnrichment(SummaryEnrichment):
    """
    Test for enrichment by shifting the
    input reference ordered scores.

    First, we search for features
    close to the top values of the reference 
    ordered input scores, and get associations
    between features and categories.

    Then, we repeat the following step many times:
    Randomly shift the input scores,
    against the (chrom, pos) index and
    repeat the first step for the shifted data.

    Repeating the second step many times,
    gives a null distribution of
    associations between features and categories,
    against which we compare the real assocations
    from step one.
    """
    def __init__(self,value_s, feature_df,  feature_to_category,
                    top_type, top, ascending=False):
        top_types = ['count','threshold','quantile']
        assert top_type in top_types, "top_type must be one of {}".format(top_types)

        self.value_s = value_s.copy()
        self.feature_df = feature_df.copy()
        (self.feature_to_category, (self.feature_name, self.category_name)) = \
                                                         make_df(feature_to_category,
                                                                ['feature', 'category'])
        if top_type == 'count':
            top_n = top
        elif top_type == 'quantile':
            top_n = int(len(self.value_s)*top)
        elif top_type == 'threshold':
            value_s_s = value_s.sort(ascending=ascending, inplace=False)
            if ascending:
                top_n = np.argmax(value_s_s.values>top)
            else:
                top_n = np.argmax(value_s_s.values<top)
            del value_s_s
        self.top_n = top_n
        self.ascending = ascending


    def get_association(self, value_s):
        value_s = value_s.sort(ascending=self.ascending, inplace=False)
        top_s = value_s.iloc[:self.top_n]
        del value_s
        cand_genes = get_genes(top_s, gene_df, max_dist=args.max_dist)
        
        return assoc

   
#    def __init__(self, value_s, feature_df,  feature_to_category, feature_summary, category_summary, ascending, max_dist, top_n):
#        super(TopScoresEnrichment, self).__init__()
#
#
#
#            value_s = value_s.sort(ascending=args.ascending, inplace=False)
#            top_s = value_s.iloc[:top_n]
#            del value_s
#            cand_genes = get_genes(top_s, gene_df, max_dist=args.max_dist)
#
#            #save peak info for candidate genes:
#            peaks_per_gene = get_peaks(cand_genes,gene_df,top_s,args.max_dist)
#            peaks_per_gene.to_csv(args.peaks_per_gene_fn,sep=ppg_sep)
#
#            #save list of top peaks
#            peak_info = get_peak_info(top_s,peaks_per_gene)
#            peak_info.to_csv(args.top_peaks_fn,sep=tp_sep)
#
#            assoc = get_go_assoc(cand_genes, gene_to_cat)
#            assoc.to_csv(args.assoc_fn,sep=assoc_sep,header=True)
#
#
#    def permut_assoc(rod_s, rnd, gene_df, gene_to_go, top_n, max_dist, ascending):
#        """
#        This is the main function.
#        Use with rnd = 0 to get the real assoc.
#        """
#        s = shift_rod(rod_s, rnd)
#        s.sort(ascending=ascending, inplace=True)
#        top_s = s.iloc[:top_n].copy()
#        del s
#        gc.collect() #otherwise memory seems to leak
#        #print(gc.garbage)
#        cand_genes = get_genes(top_s, gene_df, max_dist=max_dist)
#        assoc = get_go_assoc(cand_genes, gene_to_go)
#        return assoc
#
#    def multiple_permut_assoc_low_mem(rod_s, init_rank_table, gene_df, gene_to_go, top_n, max_dist, n_runs, ascending, rnd_seed=None):
#        """
#        todo: Document!!!
#        """
#        if rnd_seed is not None:
#            np.random.seed(rnd_seed)
#        if not rod_s.index.is_monotonic:
#            rod_s = rod_s.sort_index()
#        for rnd in np.random.rand(n_runs):
#            assoc = permut_assoc(rod_s, rnd, gene_df, gene_to_go, 
#                                           top_n, max_dist, ascending)
#            init_rank_table["rank"] += (init_rank_table["n_genes"] > \
#                                        assoc.reindex(init_rank_table.index).fillna(0))
#            init_rank_table["out_of"] += 1
#        return init_rank_table
#
#
#        rank_table = multiple_permut_assoc_low_mem(value_s, rank_table,
#                                            gene_df,
#                                            gene_to_cat, 
#                                            top_n,
#                                            args.max_dist,
#                                            args.n_runs,
#                                            ascending=args.ascending,
#                                            rnd_seed=None)


#general l I/O functions

def get_sep(fn):
    ext = os.path.splitext(fn)[-1]
    if ext == ".tsv":
        sep = "\t"
    elif ext == ".csv":
        sep = ","
    else:
        raise Exception("File ending must be .tsv or .csv but file is {}".format(fn))
    return sep

def read_table(file_handle,sep=None,**kwargs):
    if sep is None:
        sep = get_sep(file_handle.name)
    try:
        return pd.read_csv(file_handle,sep=sep,**kwargs)
    except pd.parser.CParserError:
        print(file_handle)
        raise
    


#these functions are needed in all modes


def get_genes(peak_s, gene_df, max_dist):
    """
    take the input series and gets 
    names of genes nearby

    Input:
    peak_s ... pandas series with (chrom, pos) index and value of
                the statistic ('peak height'). Series should be named.
    gene_df ... data frame with gene info 
    """
    all_genes = []
    if not gene_df.index.is_monotonic:
        gene_df = gene_df.sort_index()
    tot_hit_df = pd.DataFrame()
    for chrom in peak_s.index.droplevel(1).unique():
        loc_gene_df = gene_df.ix[chrom]
        #loc_gene_df = loc_gene_df.append(pd.DataFrame(np.nan,index=[np.inf],columns=loc_gene_df.columns))
        pos_rel_to_start = np.searchsorted(loc_gene_df.index.values-max_dist,peak_s.ix[chrom].index.values)
        pos_rel_to_end = np.searchsorted(loc_gene_df["end"].values+max_dist,peak_s.ix[chrom].index.values)
        genes = list(set(loc_gene_df["gene_id"].iloc[np.hstack([range(a,b) for a,b in zip(pos_rel_to_end,pos_rel_to_start)])]))
        all_genes += genes
    return all_genes

def get_top_n(value_s,top_q):
    return int(len(value_s)*top_q)

def permut_assoc(rod_s, rnd, gene_df, gene_to_go, top_n, max_dist, ascending):
    """
    This is the main function.
    Use with rnd = 0 to get the real assoc.
    """
    s = shift_rod(rod_s, rnd)
    s.sort(ascending=ascending, inplace=True)
    top_s = s.iloc[:top_n].copy()
    del s
    gc.collect() #otherwise memory seems to leak
    #print(gc.garbage)
    cand_genes = get_genes(top_s, gene_df, max_dist=max_dist)
    assoc = get_go_assoc(cand_genes, gene_to_go)
    return assoc

def get_go_assoc(gene_ls, gene_to_go):
    """
    Get series with number of genes associated with each
    category in gene_to_go
    """
    s = gene_to_go.set_index("gene_symbol").ix[gene_ls].groupby("go_identifier").apply(len)
    s.name = "n_genes"
    s.index.name = "category"
    return s




def multiple_permut_assoc_low_mem(rod_s, init_rank_table, gene_df, gene_to_go, top_n, max_dist, n_runs, ascending, rnd_seed=None):
    """
    todo: Document!!!
    """
    if rnd_seed is not None:
        np.random.seed(rnd_seed)
    if not rod_s.index.is_monotonic:
        rod_s = rod_s.sort_index()
    for rnd in np.random.rand(n_runs):
        assoc = permut_assoc(rod_s, rnd, gene_df, gene_to_go, 
                                       top_n, max_dist, ascending)
        init_rank_table["rank"] += (init_rank_table["n_genes"] > \
                                    assoc.reindex(init_rank_table.index).fillna(0))
        init_rank_table["out_of"] += 1
    return init_rank_table

def save_permut_assoc_table(assoc_table,fn):
    assoc_table.to_csv(fn)

#these functions are only needed in the reduce mode
def get_gene_info(gene_ls,gene_df):
    """
    for a list of gene ids,
    get a data frame with their 
    position
    """
    gi = gene_df[gene_df["gene_id"].apply(lambda x: x in gene_ls)]
    return gi

def get_peaks(gene_ls,gene_df,top_s,max_dist):
    """
    For each gene in gene_info get the
    peaks within max_dist in top_s. This 
    is basically reverse engineering to get
    the peak info for each gene that was found 
    to be associated with a peak. 
    
    Input:
    gene_info ... data frame with index ('chrom','start')
                and columns 'gene_id' and 'end'
    top_s ... series of peak positions with index (chrom, pos)
                and values peak height
    max_dist ... maximum distance between gene and peak
    """
    gene_info = get_gene_info(gene_ls, gene_df)
    def get_dist(df,gene_pos):
        """
        calculate distance
        """
        s = pd.Series(df.index.droplevel(0).values - gene_pos.ix[df.index[0][0]],
                                                  index=df.index.droplevel(0).values)
        #df = pd.DataFrame(df.index.droplevel(0).values - gene_pos.ix[df.index[0][0]],
        #                                          index=df.index.droplevel(0).values)
        return s
    tot_gene_peaks_df = pd.DataFrame()
    if not top_s.index.is_monotonic:
        top_s = top_s.sort_index()
    for chrom in gene_info.index.droplevel(1).unique():
        loc_top_s = top_s.ix[chrom]
        start = np.searchsorted(loc_top_s.index.values+max_dist,gene_info.ix[chrom].index.values)
        end = np.searchsorted(loc_top_s.index.values-max_dist,gene_info.ix[chrom]["end"].values)
        
        x = pd.concat([loc_top_s.iloc[st:ed] for st,ed in zip(start,end)],
                      keys=gene_info.ix[chrom]["gene_id"].values)
        x.name = "peak_height"
        


        dist_start = x.groupby(lambda i: i[0]).\
                    apply(lambda df: get_dist(df,
                                              gene_info.ix[chrom].reset_index().set_index("gene_id")["pos"]))
        dist_start.name = "dist_start"
        dist_end = x.groupby(lambda i: i[0]).\
                    apply(lambda df: get_dist(df,
                                              gene_info.ix[chrom].set_index("gene_id")["end"]))
        dist_end.name = "dist_end"
        gene_peaks_df = pd.concat([x,dist_start,dist_end],axis=1)
        gene_peaks_df.index = pd.MultiIndex.from_arrays([gene_peaks_df.index.droplevel(1),
                                         [chrom]*len(x),
                                         gene_peaks_df.index.droplevel(0)])
        tot_gene_peaks_df = pd.concat([tot_gene_peaks_df, gene_peaks_df],axis=0)
            

    tot_gene_peaks_df.index.names = ["gene_id","chrom","peak_pos"]
    return tot_gene_peaks_df


def get_peak_info(top_s,peaks_per_gene):
    peak_height_name = top_s.name
    gene_list_peak_pos = peaks_per_gene.reset_index([0])["gene_id"].groupby(lambda x: x).apply(list)
    gene_list_peak_pos.name = "genes"
    gene_list_peak_pos.index = pd.MultiIndex.from_tuples(gene_list_peak_pos.index)
    peak_info = pd.concat([top_s,gene_list_peak_pos],axis=1)
    peak_info.sort(peak_height_name,ascending=False,inplace=True)
    peak_info.index.names = ["chrom","pos"]
    return peak_info

def get_initial_rank_table(real_assoc):
    df = pd.DataFrame({"n_genes":real_assoc.values,"rank":0,"out_of":0},index=real_assoc.index)
    df.index.name = "category"
    return df

def get_genes_per_go(gene_ls, gene_to_go):
    s = gene_to_go.set_index("gene_symbol").ix[gene_ls].groupby("go_identifier").apply(lambda x: list(x.index))
    s.name = "genes"
    return s

def get_p_val(rank_table):
    """
    Input:
    
    """
    r =  1-rank_table["rank"]*1./(rank_table["out_of"]+1)
    r.sort()
    r.name = "p_value"
    return r

def update_rank(rank_table,assoc_table):
    """
    UNUSED, remove????
    todo: describe what assoc_table is, isn't this a series?
    """
    r = assoc_table.apply(lambda row: get_rank(row,rank_table),axis=1)
    r.index.name = "category"
    return r



def total_rank(rank_table, permut_fns):
    """
    Unused, remove?????
    """
    rt = rank_table.copy()
    for f in permut_fns:
        rt = update_rank(rt,f)
    return rt


def empirical_rank(value,dist):
    """
    get empirical p value of
    value with respect to list of 
    values in dist
    """
    array = np.append(value,dist)
    temp = array.argsort()
    ranks = np.empty(len(array), int)
    ranks[temp] = np.arange(len(array))
    return ranks[0]

def get_rank(series,rank_df):
    go = series.name
    try:
        go_s = rank_df.ix[go]
    except KeyError:
        go_s = pd.Series({"n_genes":0,"rank":0,"out_of":0})
    real_val = go_s["n_genes"]
    old_rank = go_s["rank"]
    old_out_of = go_s["out_of"]
    rank = empirical_rank(real_val,series.values)
    new_rank = old_rank + rank
    new_out_of = old_out_of + len(series)

    return pd.Series({"n_genes":real_val,"rank":new_rank,"out_of":new_out_of})






if __name__ == "__main__":
    import argparse
    import pdb

    parser = argparse.ArgumentParser(description=
                                      "Test enrichment of genes in certain categories "
                                      "compared to the empirical distribution\n "
                                      "of genes per category, produced by \n "
                                      "shifting (rotating) the input data across "
                                                                        "the genome.")
    parser.add_argument("--gene_to_cat_fn", type=argparse.FileType('r'),
                                     help="Filename for file that links genes and"
                                          " categories. E.g. go associations. "
                                          "Index should be unique.",
                                                                        required=True)
    parser.add_argument("-v","--verbose",
                            type=int,choices=[0,1],default=1)
    subparsers = parser.add_subparsers(dest='mode',
                                        help='Mode "real_assoc", "permutation" or "reduce"')


    #######mode real assoc##########
    parser_0 = subparsers.add_parser('real_assoc', 
                                      help='Get assocications for the real original'
                                                          ' data set and save them.')
    parser_0.add_argument("--peaks_per_gene_fn",type=argparse.FileType('w'),
                            required=True,
                                   help="File to write peak info for each gene.")
    parser_0.add_argument("--top_peaks_fn",type=argparse.FileType('w'),required=True,
                                                  help="File to write top peak info.")
    parser_0.add_argument("--assoc_fn",type=argparse.FileType('w'),
                            required=True,help="File to write associations"
                                                                   "per gene.")
    
    parser_0.add_argument("--log_fn",type=argparse.FileType('w'),
                                    help="File to write output log to.")
    ######mode permutation##########
    parser_a = subparsers.add_parser('permutation', help='Do random shifts of input rod and write \n'
                                                         'the results to a file.')
    parser_a.add_argument("--n_runs", type=int, required=True)
    parser_a.add_argument("--out_fn",type=argparse.FileType('w'),required=True)
    parser_a.add_argument("--real_assoc_fn",type=argparse.FileType('r'),required=True,
                                                                help="File to read real associations"
                                                                                     "from.")
    parser_a.add_argument("--log_fn",type=argparse.FileType('a'),help="File to append output log to.")
    #####shared by real assoc and permutation#####
    
    for p in [parser_0, parser_a]:
        p.add_argument("--in_rod",type=argparse.FileType('r'),
                                     help="Input reference ordered data containing the values per position.",
                                        required=True)
        p.add_argument("--cols",type=str, nargs=3, help="Column labels or positions of [chrom pos value]. "
                                                        "Expects three "
                                                        " integers or three strings.",
                                                                                   required=True)
        group = p.add_mutually_exclusive_group(required=True)
        group.add_argument("--top_n",type=int, help="Number of top values of the input rod to consider.")
        group.add_argument('--top_q', type=float,help="Input quantile to consider, e.g. 0.01 for top 1%.")
        group.add_argument("--thresh",type=float,help="Set a threshold for peak height. "
                                                      "Enrichment is tested for the peaks with "
                                                      "value>thresh (or < if acending)")
        p.add_argument("--max_dist",type=int, help="Maximum distance (in bp) between gene and rod feature"
                                                    "in order to consider a gene.", required=True)
  
        p.add_argument("--gene_df_fn", type=argparse.FileType('r'), 
                                        help="Filename of the gene info dataframe "
                                             " with index (chrom,start) and columns 'gene_id', 'end'",
                                                                required=True)
        g2 = p.add_mutually_exclusive_group(required=False)
        g2.add_argument("--ascending",dest="ascending",action="store_true",help="Sort ascending, (e.g., for p-values).")
        g2.add_argument("--descending",dest="ascending",action="store_false",help="Sort descending, (e.g., for scores).")

        g2.set_defaults(ascending=False)

    ######mode reduce##########
    parser_b = subparsers.add_parser('reduce', help='Get rank and p values for real data \n'
                                                    'compared to shifted data runs.')
    parser_b.add_argument("--permut_fns",nargs='*',default=None,type=str)#type=argparse.FileType('r'))
    parser_b.add_argument("--cat_to_name_fn",default=None,type=argparse.FileType('r'),help="Table that maps categories"
                                                                                      "to category descriptions.")
    parser_b.add_argument("--pval_out",type=argparse.FileType('w'),help="File to write enrichment p_values to.")
    
    parser_b.add_argument("--pval_sign_out",type=argparse.FileType('w'),
                                                        help="File to write significant enrichment p_values to.")
    parser_b.add_argument("--pval_sign_thresh",type=float,
                                                        help="P value threhold for pval_sign_out.")
    parser_b.add_argument("--peaks_per_gene_fn",type=argparse.FileType('r'),help="File path of gene info as produced in"
                                                                            " real assoc mode. must have a column gene_id.")

    parser_b.add_argument("--log_fn",type=argparse.FileType('w'),help="File to append output log to.")

    args = parser.parse_args()


    def printv(*text,**kwa):
        if args.log_fn is not None:
            print(*text,file=args.log_fn,**kwa)
        if args.verbose > 0:
            print(*text,**kwa)

    if args.mode != "permutation":
        printv("Running",os.path.basename(__file__),"in mode",args.mode+".")
        printv("Interpreted input arguments as:")
        printv(args)

    gene_to_cat = read_table(args.gene_to_cat_fn,index_col=[0])


    if args.mode == "real_assoc" or args.mode == "permutation":
        #get data
        try:
            cols = [int(i) for i in args.cols]
        except ValueError:
            cols = args.cols
        value_s = read_table(args.in_rod,index_col=[0,1],
                                    usecols=cols, squeeze=True)
        value_s.index.names = ["chrom","pos"]
        assert isinstance(value_s,pd.core.series.Series)

        gene_df = read_table(args.gene_df_fn,index_col=[0,1])

        if args.top_n is not None:
            top_n = args.top_n

        elif args.top_q is not None:
            top_n = int(len(value_s)*args.top_q)
        else:
            value_s_s = value_s.sort(ascending=args.ascending, inplace=False)
            if args.ascending:
                top_n = np.argmax(value_s_s.values>args.thresh)
            else:
                top_n = np.argmax(value_s_s.values<args.thresh)
            del value_s_s
        

        if args.mode == "real_assoc":

            printv("Using the top", top_n, "peaks.")

            ppg_sep = get_sep(args.peaks_per_gene_fn.name)
            tp_sep = get_sep(args.top_peaks_fn.name)
            assoc_sep = get_sep(args.assoc_fn.name)
           
            value_s = value_s.sort(ascending=args.ascending, inplace=False)
            top_s = value_s.iloc[:top_n]
            del value_s
            cand_genes = get_genes(top_s, gene_df, max_dist=args.max_dist)

            #save peak info for candidate genes:
            peaks_per_gene = get_peaks(cand_genes,gene_df,top_s,args.max_dist)
            peaks_per_gene.to_csv(args.peaks_per_gene_fn,sep=ppg_sep)

            #save list of top peaks
            peak_info = get_peak_info(top_s,peaks_per_gene)
            peak_info.to_csv(args.top_peaks_fn,sep=tp_sep)

            assoc = get_go_assoc(cand_genes, gene_to_cat)
            assoc.to_csv(args.assoc_fn,sep=assoc_sep,header=True)

        elif args.mode == "permutation":
    
            out_sep = get_sep(args.out_fn.name)
            real_assoc = read_table(args.real_assoc_fn,squeeze=True,index_col=0)
            
            rank_table = get_initial_rank_table(real_assoc)

            rank_table = multiple_permut_assoc_low_mem(value_s, rank_table,
                                                gene_df,
                                                gene_to_cat, 
                                                top_n,
                                                args.max_dist,
                                                args.n_runs,
                                                ascending=args.ascending,
                                                rnd_seed=None)

            del value_s
            #rank_table = update_rank(rank_table,assoc_table)
            rank_table = rank_table[["n_genes","rank","out_of"]]
            rank_table.to_csv(args.out_fn,sep=out_sep)

    elif args.mode == "reduce":
        
        if args.pval_sign_out is not None:
            assert args.pval_sign_thresh is not None, "specify argument pval_sign_thresh or remove argument pval_sign_out "
        if args.permut_fns is None:
            args.permut_fns = []
        
        permut_fhs = []
        for fn in args.permut_fns:
            try:
                if os.stat(fn).st_size>0:
                    permut_fhs.append(open(fn))
                else:
                    printv(fn, "seems to be empty. Skipping it.")
            except Exception, e:
                printv("Can't open file, skipping it.")
                printv(str(e))

        out_sep = get_sep(args.pval_out.name)
        sign_out_sep = get_sep(args.pval_sign_out.name)
        tot_rank = read_table(permut_fhs[0],index_col=0).dropna()
        tot_rank["index"] = tot_rank.index
        tot_rank.drop_duplicates(subset="index",inplace=True)
        del tot_rank["index"]
        
        for fh in permut_fhs[1:]:
            rank_table = read_table(fh,index_col=0).dropna()
            rank_table["index"] = rank_table.index
            rank_table.drop_duplicates(subset="index",inplace=True)
            try:
                tot_rank["rank"] = tot_rank["rank"].add(rank_table["rank"],fill_value=0)
                tot_rank["out_of"] = tot_rank["out_of"].add(rank_table["out_of"],fill_value=0)
            except:
                raise

        p_vals = get_p_val(tot_rank)

        p_val_df = tot_rank.join(p_vals)
        p_val_df.index.name = "category"
        if args.cat_to_name_fn is not None:
            cat_to_name = read_table(args.cat_to_name_fn,index_col=[0])
            cat_to_name = cat_to_name.set_index("go_identifier",drop=True)
            p_val_df = p_val_df.join(cat_to_name)

        # get a list of genes per go category
        try:
            cand_genes = np.unique(read_table(args.peaks_per_gene_fn,usecols=["gene_id"]).values)
            gene_per_go_s = get_genes_per_go(cand_genes, gene_to_cat)
            p_val_df = pd.concat([p_val_df, gene_per_go_s], axis=1)
            def check_len(el):
                try:
                    return(len(el))
                except TypeError:
                    return 0
            if not (p_val_df["genes"].apply(check_len) == p_val_df["n_genes"]).all():
                assert_df = p_val_df[["n_genes","genes"]]
                assert_df["len_genes"] = assert_df["genes"].apply(check_len)
                assert_df = assert_df[["n_genes","len_genes","genes"]]
                printv("Genes per category from",args.peaks_per_gene_fn,
                        "inconsistent with n_genes reported in",permut_fhs[0].name,":",
                        assert_df[assert_df["n_genes"]!=assert_df["len_genes"]])
        except NameError, e:
            printv("not adding gene names, no peaks_per_gene_file")
            printv(str(e))

        p_val_df[["n_genes","rank","out_of"]] = p_val_df[["n_genes","rank","out_of"]].astype(int)
        p_val_df.sort(["p_value", "n_genes"],inplace=True,ascending=[True, False])
        p_val_df["benjamini_hochberg"] = p_val_df["p_value"] * \
                                        len(gene_to_cat["go_identifier"].unique()) /\
                                        np.arange(1,len(p_val_df)+1)
        c = list(p_val_df.columns.values)
        p_val_df = p_val_df[c[:-3] + [c[-1]] + c[-3:-1]]
        p_val_df.index.name = "category"
        if args.pval_out is not None:
            p_val_df.to_csv(args.pval_out, sep=out_sep)
        if args.pval_sign_out is not None:
            sign_p_val = p_val_df[p_val_df["p_value"]<args.pval_sign_thresh]
            sign_p_val.to_csv(args.pval_sign_out, sep=sign_out_sep)
            printv("Found",len(sign_p_val),"associations above p_value of",args.pval_sign_thresh)







