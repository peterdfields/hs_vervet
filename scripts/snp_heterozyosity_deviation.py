#!/usr/bin/env python
"""
Calculates the relative difference between expected
and observed heterozygosity.
"""
import sys, os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description=   "Calculates the relative difference between\n"
                                                "expected and observed heterozygosity.")
parser.add_argument("in_012_tsv",type = argparse.FileType('r'), default = '-', help="Input genotype-file filename.")
parser.add_argument("out_fn",type = argparse.FileType('w'), default = '-', help="Output tsv filename.")

args = parser.parse_args()

genotype_df = pd.read_csv(args.in_012_tsv,sep="\t",index_col=[0,1],na_values=['N'])


af = genotype_df.mean(axis=1)*1./2
Ho = (genotype_df==1).sum(axis=1)*1./pd.notnull(genotype_df).sum(axis=1)
del genotype_df
He = 2*af*(1-af)
del af

deviation = (Ho-He)*1./He
del He
del Ho
deviation.name = "(Ho-He)/He"
deviation.to_csv(args.out_fn,sep="\t",header=True,index_label=("chrom","pos"))


