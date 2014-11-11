#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
eu = os.path.expanduser
jn = os.path.join

print pd.__version__

#general I/O functions

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
    return pd.read_csv(file_handle,sep=sep,**kwargs)


#these functions are needed in all modes

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
        #print chrom, genes
        all_genes += genes
    return all_genes

def permut_assoc(rod_s, rnd, gene_df, gene_to_go, top_n, max_dist):
    """
    This is the main function.
    Use with rnd = 0 to get the real assoc.
    """
    s = shift_rod(rod_s, rnd)
    s.sort(ascending=False, inplace=True)
    top_s = s.iloc[:top_n]
    cand_genes = get_genes(top_s, gene_df, max_dist=max_dist)
    assoc = get_go_assoc(cand_genes, gene_to_go)
    return assoc

def get_go_assoc(gene_ls, gene_to_go):
    """
    Get series with number of genes associated with each
    category in gene_to_go
    """
    s = gene_to_go.set_index("gene_symbol").ix[gene_ls].groupby("go_identifier").apply(len)
    return s

def multiple_permut_assoc(rod_s, gene_df, gene_to_go, top_n, max_dist, n_runs, rnd_seed=None):
    if rnd_seed is not None:
        np.random.seed(rnd_seed)
    assoc_table = pd.concat([permut_assoc(rod_s, rnd, gene_df, gene_to_go, top_n, max_dist) for rnd in np.random.rand(n_runs)],axis=1)
    assoc_table = assoc_table.fillna(0).astype(int)
    return assoc_table

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

def get_peaks(gene_info,top_s,max_dist):
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
    #print "gene_info"
    #print gene_info
    #print 'top_s'
    #print top_s
    #print 'max_dist'
    #print max_dist
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
        #dist_start.columns = ["dist_start"]
        dist_end = x.groupby(lambda i: i[0]).\
                    apply(lambda df: get_dist(df,
                                              gene_info.ix[chrom].set_index("gene_id")["end"]))
        dist_end.name = "dist_end"
        #dist_end.columns = "dist_end"
        #the following has problems if there is only a single gene
        gene_peaks_df = pd.concat([x,dist_start,dist_end],axis=1)
        #gene_peaks_df = x.join(dist_start).join(dist_end)
        gene_peaks_df.index = pd.MultiIndex.from_arrays([gene_peaks_df.index.droplevel(1),
                                         [chrom]*len(x),
                                         gene_peaks_df.index.droplevel(0)])
        try:
            tot_gene_peaks_df = pd.concat([tot_gene_peaks_df, gene_peaks_df],axis=0)
        except:
            print "loc_top_s"
            print loc_top_s
            print "start"
            print start
            print "end"
            print end
            print "x"
            print x
            print "dist_start"
            print dist_start
            print "dist_end"
            print dist_end
            print "gene_peaks_df"
            print gene_peaks_df
            print "gene_info"
            print gene_info
            print 'gene_info.ix[chrom].reset_index().set_index("gene_id")["pos"]'
            print gene_info.ix[chrom].reset_index().set_index("gene_id")["pos"]
            

    tot_gene_peaks_df.index.names = ["gene_id","chrom","peak_pos"]
    return tot_gene_peaks_df

def get_initial_rank_table(real_assoc):
    return pd.DataFrame({"n_genes":real_assoc.values,"rank":0,"out_of":0},index=real_assoc.index)

def get_p_val(rank_table):
    """
    Input:
    
    """
    print "get_p_val input rank table"
    print rank_table
    try:
        r =  1-rank_table["rank"]*1./(rank_table["out_of"]+1)
    except:
        print "rank_table within get_p_val 2"
        print rank_table
        raise
    r.sort()
    r.name = "p_value"
    return r

def update_rank(rank_table,permut_assoc_fh):
    
    assoc_table = read_table(permut_assoc_fh,index_col=0)
    r = assoc_table.apply(lambda row: get_rank(row,rank_table),axis=1)
    print "assoc_table"
    print assoc_table
    print "rank_table"
    print rank_table
    print "new_rank_table"
    print r
    return r
    
def total_rank(rank_table, permut_fns):
    rt = rank_table.copy()
    for f in permut_fns:
        rt = update_rank(rt,f)
    return rt
    
#def save_p_value_total(rank_table, permut_fns):
    
    
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


    parser = argparse.ArgumentParser(description="Test enrichment of genes in certain categories "
                                                 "compared to the empirical distribution\n "
                                                 "of genes per category, produced by \n "
                                                 "shifting (rotating) the input data across the genome.")
    subparsers = parser.add_subparsers(dest='mode',help='Mode "permutation" or "reduce"')
    
    parser_a = subparsers.add_parser('permutation', help='Do random shifts of input rod and write \n'
                                                         'the results to a file.')
    parser_a.add_argument("--n_runs", type=int, required=True)
    parser_a.add_argument("--out_fn",type=argparse.FileType('w'),required=True)

    parser_b = subparsers.add_parser('reduce', help='Get rank and p values for real data \n'
                                                    'compared to shifted data runs.')
    parser_b.add_argument("--permut_fns",nargs='*',default=None,type=argparse.FileType('r'))
    parser_b.add_argument("--cat_to_name_fn",default=None,type=argparse.FileType('r'),help="Table that maps categories"
                                                                                      "to category descriptions.")
    parser_b.add_argument("--pval_out",type=argparse.FileType('w'),required=True,help="File to write enrichment p_values to.")
    parser_b.add_argument("--peaks_per_gene_fn",type=argparse.FileType('w'),required=True,
                                                                help="File to write peak info for each gene.")



    parser.add_argument("--in_rod",type=argparse.FileType('r'),
                                     help="Input reference ordered data containing the values per position.",
                                        required=True)
    parser.add_argument("--col",type=str,help="Column label of the input table that contains the values"
                                                    " by which the rod positions should be ranked",
                                                                                    required=True)
    parser.add_argument("--top_n",type=int,help="number of top values of the input rod to consider.",
                                                                        required=True)
    parser.add_argument("--max_dist",type=int, help="Maximum distance (in bp) between gene and rod feature"
                                                    "in order to consider a gene.", required=True)
    parser.add_argument("--gene_df_fn", type=argparse.FileType('r'), 
                                        help="Filename of the gene info dataframe "
                                             " with index (chrom,start) and columns 'gene_id', 'end'",
                                                                required=True)
    parser.add_argument("--gene_to_cat_fn",type=argparse.FileType('r'),
                                            help="Filename for file that links genes and"
                                                 " categories. E.g. go associations. Index should be unique.",
                                                    required=True)


    args = parser.parse_args()

    value_s = read_table(args.in_rod,index_col=[0,1],
                            usecols=["chrom", "pos", args.col])
    gene_df = read_table(args.gene_df_fn,index_col=[0,1])
    gene_to_cat = read_table(args.gene_to_cat_fn,index_col=[0])



    if args.mode == "reduce":
        if args.permut_fns is None:
            args.permut_fns = []

        out_sep = get_sep(args.pval_out.name)
        ppg_sep = get_sep(args.peaks_per_gene_fn.name)
        #print "value_s"
        #print value_s
        value_s.sort(args.col,ascending=False, inplace=True)
        #print "value_s"
        #print value_s
        top_s = value_s.iloc[:args.top_n]
        #print "top_s"
        #print top_s
        del value_s
        cand_genes = get_genes(top_s, gene_df, max_dist=args.max_dist)

        print "gene_df"
        print "gene_df"
        #save peak info for candidate genes:
        gene_info = get_gene_info(cand_genes, gene_df)
        peaks_per_gene = get_peaks(gene_info,top_s,args.max_dist)
        peaks_per_gene.to_csv(args.peaks_per_gene_fn,sep=ppg_sep)


        assoc = get_go_assoc(cand_genes, gene_to_cat)

        #print "real_assoc:"
        #print real_assoc
        rank_table = get_initial_rank_table(assoc)
        #if no files are provided we just get the intital rank table
        tot_rank = total_rank(rank_table, args.permut_fns)
        p_vals = get_p_val(tot_rank)

        p_val_df = tot_rank.join(p_vals)
        p_val_df.sort("p_value",inplace=True)
        p_val_df.index.name = "category"
        if args.cat_to_name_fn is not None:
            cat_to_name = read_table(args.cat_to_name_fn,index_col=[0])
            cat_to_name = cat_to_name.set_index("go_identifier",drop=True)
            p_val_df = p_val_df.join(cat_to_name)
        p_val_df.to_csv(args.pval_out, sep=out_sep)
        


    elif args.mode == "permutation":
        out_sep = get_sep(args.out_fn.name)

        assoc_table = multiple_permut_assoc(value_s, 
                                            gene_df,
                                            gene_to_cat, 
                                            args.top_n,
                                            args.max_dist,
                                            args.n_runs,
                                            rnd_seed=None)
        assoc_table.to_csv(args.out_fn,sep=out_sep)



