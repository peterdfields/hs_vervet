import numpy as np
import pandas as pd

def index_rolling(s,window,func,overlap=0,min_n_values=0,*args,**kwargs):
    """
    Apply function in rolling windows, where the window
    size is defined with respect to the index values.
    This means that different windows can comprise different
    numbers of elements.
    
    s ... pandas Series
    window ... window size in units of the index values
    func ... function to apply to the series values within each
                window
    overlap ... oberlap size of windows
    min_n_values ... minimum number of values (e.g. SNPs) per window
                    None will be returned for windows that have less values.

    args, kwarg ... additional arguments for func
    
    Example: index_rolling(pd.Series(range(20),index=np.random.rand(20))),0.1,max)
    """
    assert isinstance(s, pd.Series)
    #note that basis must be sorted in order for this to work properly
    windows_min = s.index.min()
    windows_max = s.index.max()
    window_starts = np.arange(windows_min, windows_max-window, window-overlap)
    window_starts = pd.Series(window_starts, index = window_starts+window/2)
    def applyToWindow(val):
        # using slice_indexer rather that what.loc [val:val+window] allows
        # window limits that are not specifically in the index
        try:
            indexer = s.index.slice_indexer(val,val+window,1)
        except IndexError:
            print val, val+window
            print s
            raise
        chunk = s.iloc[indexer]
        try:
            if len(chunk) < min_n_values:
                #print indexer, "too few snps"
                return None
        except TypeError:
            return None
        try:
            return func(chunk,*args,**kwargs)
        except ValueError, e:
            if "empty sequence" in str(e):
                #print indexer, chunk
                return None
            else:
                raise
    rolled = window_starts.apply(applyToWindow)
    return rolled

def data_per_feature(rod,feature_df):
    """
    Get the entires in rod which lie within a feature
    (e.g. gene) in feature_df.
    Input:
    rod (reference ordered data)... pandas series or data frame with multiindex (chrom, pos)
                                    such as SNP genotypes
    feature_df (gene annotation data frame)... index must be (chrom,feature_name), must have columns 'start', 'end'
    """
    rod = pd.DataFrame(rod)
    feature_hit_df = pd.DataFrame()
    chrpos_names = rod.index.names
    feature_name = feature_df.index.names[1]
    for chrom in rod.index.droplevel(1).unique():
        pos_rel_to_start = feature_df.ix[chrom]['start'].searchsorted(rod.ix[chrom].index)
        pos_rel_to_end = np.searchsorted(feature_df.ix[chrom]["end"].values,rod.ix[chrom].index.values)
        in_feature = (pos_rel_to_start - pos_rel_to_end) == 1
        feature_id = feature_df.ix[chrom].iloc[pos_rel_to_end[in_feature]].index
        snp_df = rod.ix[chrom][in_feature]
        snp_df[chrpos_names[0]] = chrom
        snp_df[feature_name] = feature_id
        snp_df.set_index([chrpos_names[0],feature_name],append=True,inplace=True)
        snp_df = snp_df.reorder_levels([chrpos_names[0], feature_name,chrpos_names[1]])
        feature_hit_df = feature_hit_df.append(snp_df)
    return feature_hit_df

def get_features(peak_s, feature_df, max_dist):
    """
    take the input series and gets.
    names of features nearby

    Input:
    peak_s ... pandas series with (chrom, pos) index and value of
                the statistic ('peak height'). Series should be named.
    feature_df ... data frame with feature info.
    """
    all_features = []
    if not feature_df.index.is_monotonic:
        feature_df = feature_df.sort_index()
    tot_hit_df = pd.DataFrame()
    for chrom in peak_s.index.droplevel(1).unique():
        loc_feature_df = feature_df.ix[chrom]
        #loc_feature_df = loc_feature_df.append(pd.DataFrame(np.nan,index=[np.inf],columns=loc_feature_df.columns))
        pos_rel_to_start = np.searchsorted(loc_feature_df['start']-max_dist,peak_s.ix[chrom].index.values)
        pos_rel_to_end = np.searchsorted(loc_feature_df["end"].values+max_dist,peak_s.ix[chrom].index.values)
        features = list(set(loc_feature_df["feature_id"].iloc[np.hstack([range(a,b) for a,b in zip(pos_rel_to_end,pos_rel_to_start)])]))
        all_features += features
    return all_features


def apply_to_feature(feature_df,groupby_func_name=None,function=None):
    """
    Apply a function to the entries for each feature.
    feature_df ... dataframe with index (chrom, feature_name, pos)
                   (Such as the output of data_per_feature())
    groupby_func_name ... name of the function of the groupby object
                         to apply to the data
                          This is faster than applying a function object.
    function ... alternatively: function object to apply
    """
    groups = feature_df.groupby(lambda idx: idx[1])
    if groupby_func_name is not None:
        return getattr(groups,groupby_func_name)()
    elif function is not None:
        return groups.apply(function)
    else:
        raise ValueError("Either groupby_func_name or function have to be given.")
