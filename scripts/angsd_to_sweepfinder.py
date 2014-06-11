#!/usr/bin/env python
"""
convert angsd maf.gz output to
SweepFinder format.
"""
import os
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
    sf_df = pd.DataFrame({'x': x,
                          'n':2*angsd["nInd"],
                          'folded':fold.astype(int)})
    sf_df.index = sf_df.index.droplevel()
    sf_df.index.name = "location"
    return  sf_df.reindex_axis(["x","n","folded"], axis=1)



angsd_df = pd.read_csv(args.in_maf_gz,
                       sep="\t",compression="gzip",index_col=(0,1))


angsd_to_sweepfinder(angsd_df).to_csv(args.out_fn,sep="\t")
