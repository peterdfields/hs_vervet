#!/usr/bin/env python
"""
convert angsd maf.gz output to
SweepFinder format.
"""
import os
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Convert angsd maf.gz output to\n"
                                             "SweepFinder format.")
parser.add_argument("in_maf_gz",type = argparse.FileType('r'), 
                    default = '-',
                     help="Input angsd maf.gz. "
                          "Assumes that column anc is available.")
parser.add_argument("out_fn",
                    type = argparse.FileType('w'), 
                    default = '-', 
                    help="Output sweepfinder tsv.")

parser.add_argument("--sfs_fn",
                    type = argparse.FileType('w'),default = False,
                    help="Output filename for site frequenct spectrum.")
args = parser.parse_args()




def angsd_to_sweepfinder(angsd):
    #remove fixed sites
    ac = (angsd["knownEM"]*2*angsd["nInd"]).round()
    angsd = angsd[(ac>0)*
                        (ac<2*angsd["nInd"])]
    ac = (angsd["knownEM"]*2*angsd["nInd"]).round()
    fold = (angsd["anc"] != angsd["major"]) * \
            (angsd["anc"] != angsd["minor"])
    
    x = (angsd["anc"] == angsd["major"]) * ac + \
        (angsd["anc"] == angsd["minor"]) * (2*angsd["nInd"]-ac) + \
        fold * ac.where(ac<2*angsd["nInd"]-ac,2*angsd["nInd"]-ac) 
    sf_df = pd.DataFrame({'x': x.astype(int),
                          'n':2*angsd["nInd"],
                          'folded':fold.astype(int)})
    sf_df.index = sf_df.index.droplevel()
    sf_df.index.name = "position"
    return  sf_df.reindex_axis(["x","n","folded"], axis=1)

def sfs(sweepfinder_df):
    """
    use sweepfinder df to calculate the site
    frequency spectrum of the segregating sites with ancestral
    state information
    input:
        df as produced by angsd_to_sweepfinder(angsd)
    returns:
        series n:number of sites
    """
    hist, bins = np.histogram(sweepfinder_df["x"][sweepfinder_df["folded"]==0],
                                               bins=np.arange(0.5,n_chrom))
    index =  (bins[:-1]+0.5).astype(int)
    return pd.Series(hist,index=index)

angsd_df = pd.read_csv(args.in_maf_gz,
                       sep="\t",compression="gzip",index_col=(0,1))


sf_df = angsd_to_sweepfinder(angsd_df)

sf_df.to_csv(args.out_fn,sep="\t")

if args.sfs_fn:
    sfs(sf_df).to_csv(args.sfs_fn,header=None,sep="\t")

