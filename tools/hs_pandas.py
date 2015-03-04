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

def coding_variants_per_gene(gen_df,gene_df):
    """
    Get the variants in gen_df which lie within a feature
    (e.g. gene) in gene_df.
    Input:
    gen_df (genotype data frame)... rod such as SNP genotypes
    gene_df (gene annotation data frame)... index must be (chrom,feature_name), must have columns 'start', 'end'
    """
    gene_hit_df = pd.DataFrame()
    for chrom in gen_df.index.droplevel(1).unique():
        pos_rel_to_start = gene_df.ix[chrom]['start'].searchsorted(gen_df.ix[chrom].index)
        pos_rel_to_end = np.searchsorted(gene_df.ix[chrom]["end"].values,gen_df.ix[chrom].index.values)
        in_gene = (pos_rel_to_start - pos_rel_to_end) == 1
        gene_id = gene_df.ix[chrom].iloc[pos_rel_to_end[in_gene]].index
        snp_df = gen_df.ix[chrom][in_gene]
        snp_df['chrom'] = chrom
        snp_df['gene_id'] = gene_id
        snp_df.set_index(['chrom','gene_id'],append=True,inplace=True)
        snp_df = snp_df.reorder_levels(['chrom', 'gene_id','pos'])
        gene_hit_df = gene_hit_df.append(snp_df)
    return gene_hit_df
