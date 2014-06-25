#!/usr/bin/env python
"""
This module provides the analysis tools
used (and partly developed and defined)
in notebooks/20140607_vervetpopgen_manuscript.ipynb.
Try to keep the functions and classed here up-to-date
with those in the notebook.
"""

import os, sys
import numpy as np
import pandas as pd
eu = os.path.expanduser
jn = os.path.join
#meta_dir = eu("~/vervet_project/metadata")
#var_ana_dir = eu("~/vervet_project/analyses/20140403_UnifiedGenotyper_ref3500_non_VRC/_data")
#manuscript_dir = eu("~/vervet_project/manuscript/163_pop_manuscript")
#pops = ["aet","cyn","pyn","pys","sab","tan"]
#colors = ["blue","magenta","cyan","green","orange","red"]
#meta_df = pd.read_csv(jn(meta_dir,"163_population_ucla_id_taxon.csv"),index_col=0)
#meta_df["index"] = meta_df.index
#meta_df.drop_duplicates(subset='index', take_last=True, inplace=True)
#del meta_df["index"]
#ucla_ids = meta_df.index.unique()
#chrom_length = pd.read_csv(eu("~/vervet_project/metadata/ref3500.tsv"),sep="\t",index_col=0,squeeze=True,header=None)
#autosomes = ["CAE" + str(i) for i in range(1,30)]
#chromosomes = ["CAE" + str(i) for i in range(1,30)+["X","Y"]]

#----------------------------------------------------------------
#General
#very general

#general functions
def df_byline(file,start_line,end_line,**kwargs):
    header = pd.read_csv(file,nrows=1,**kwargs)
    columns = header.columns
    return pd.read_csv(file,
                       skiprows=start_line,
                       nrows=end_line-start_line,
                       header=None,names=columns,**kwargs)

def get_gen_df(file,start_line=0,end_line=sys.maxint):
    """
    get 0,1,2 genotype data frame
    """
    df = df_byline(file,start_line,end_line,index_col=[0,1],
                    na_values="N",sep="\t")
    return df.astype(np.float16)


def get_anc_df(chrom):
    """
    get ref, alt, ancestral data frame
    """
    fn = jn(var_ana_dir,"GATK_UG_163_ref3500_{}_filtered_ancestral_allele.tsv".format(chrom))
    return pd.read_csv(fn,index_col=[0,1],sep="\t",na_values="NA")
    

def anc_der_df(gen_df,anc_df):
    """
    returns genotype df poliarised with ancestral state
    sites where ancestral state is not ref or alt are omitted
    """
    ref_anc = (anc_df["REF"] == anc_df["AA"])
    alt_anc = (anc_df["ALT"] == anc_df["AA"])
    return pd.concat([gen_df[ref_anc],(gen_df[alt_anc]-2).abs()]).sort()

#-------------------------------------------------------------
#Global diversity
#Site frequency spectrum
#from genotypes

def get_ac(gen_df,mode="sum"):
    """
    mode = "sum" ... just sums up the entries, this means that nans are treated as 0
    mode = "mean" ... calculates average ignoring nans -> the spectrum is pulled towards
                        intermediate frequencies
    """
    n = gen_df.shape[1]
    if mode == "sum":
        return gen_df.sum(axis=1)
    elif mode == "mean":
        af = gen_df.mean(axis=1)/2.
        return (af*2*n).round()
    else:
        raise

def get_acs(gen_df):
    """
    get allele count spectrum (histogram)
    """
    n = gen_df.shape[1]
    ac = get_ac(gen_df)
    bins = np.arange(-0.5,2*n+0.501)
    h, bins = np.histogram(ac,bins=bins)
    return h, bins

def acs_from_chrom(chrom,anc_der=True):
    """
    allele count spectrum (histogram)
    for either ancestral-derived (anc_der=True)
    or ref-alt alleles
    """
    gen_df = get_gen_df(chrom)
    if anc_der:
        anc_df = get_anc_df(chrom)
        df = anc_der_df(gen_df,anc_df)
        del anc_df
    else:
        df = gen_df
    del gen_df
    return get_acs(df)

#-------------------------------------------------
#pairwise difference matrix

def pairwise_diff_numpy(gen_arr):
    """Squared pairwise distances between all 
    columns of 0,1,2 genotype array arr.
    This matrix based function is at least 10 
    times faster than iterating over columns.
    """
    gen_arr = gen_arr.astype(np.float64)-1
    #compare heterozygous with hom alt
    mat1 = np.where(gen_arr==0,-1,gen_arr)
    mat1[np.isnan(gen_arr)]=0
    #and hom ref
    mat2 = np.where(gen_arr==0,1,gen_arr)
    mat2[np.isnan(gen_arr)]=0
    #account for heterozygous comparisons
    mat3 = np.where(gen_arr==0,1,0)
    #don't count nan comparisons
    n = np.dot((~np.isnan(gen_arr)*1.).T,~np.isnan(gen_arr)*1.)#gen_arr.shape[0]
    B = (np.dot(mat1.T,mat1)+np.dot(mat2.T,mat2))/2.
    het = np.dot(mat3.T,mat3)
    return ((n- B)+het)/2.

def pairwise_diff_mat(df):
    """
    Calculate pairwise difference data frame.
    For Genotype 0,1,2 data frame.
    Uses numpy matrix multiplication.
    """
    diff = pairwise_diff_numpy(df.values)
    return pd.DataFrame(diff,index=df.columns,columns=df.columns)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Perform different popgen analysis."
                                                    "Mostly for chunks of the genome.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(help='Type of analysis.',dest="analysis")
    #GENERAL REDUCE
    reduce = subparsers.add_parser("reduce",help="reduce tsvs/ data frames")
    reduce.add_argument("mode",choices=["add","cat"])
    reduce.add_argument('files', nargs='+',type=argparse.FileType('r'),help="Tsv files to reduce.")
    reduce.add_argument("--index_col",nargs="+",type=int,default=[0,1],help="Index columns in the tsv")
    reduce.add_argument("-o",type=argparse.FileType('w'), default = "sys.stdout",help="Output filename.")
    #PAIRWISE DIFF
    pairwise_diff = subparsers.add_parser("pairwise_diff",help="Calculate pairwise distance matrix.")
 

    pairwise_diff.add_argument("file", type=argparse.FileType('r'), default = sys.stdin,
                                                help="Input 0,1,2 genotype matrix filename.")
    pairwise_diff.add_argument("-s",type=int,default=1,help="Start line to parse.")
    pairwise_diff.add_argument("-e",type=int,default=sys.maxint,help="End line to parse.")
    pairwise_diff.add_argument("-o",type=argparse.FileType('w'), default = "sys.stdout",help="Output filename.")
    


    args = parser.parse_args()

    if args.analysis == "reduce":
        if args.mode == "add":
            i=0
            total_df = pd.read_csv(args.files[0],sep="\t",index_col=args.index_col)
            for f in args.files[1:]:
                i += 1
                df = pd.read_csv(f,sep="\t",index_col=args.index_col)
                total_df += df
            total_df.to_csv(args.o,sep="\t")
        elif args.mode == "cat":
            raise ValueError("Mode 'cat' not implemented.")
    elif args.analysis == "pairwise_diff":
        assert args.s >= 1
        #get genotype data frame
        gen_df = get_gen_df(args.file,start_line=args.s,end_line=args.e)
        #print gen_df
        #calculate pairwise difference data frame
        pw_diff_df = pairwise_diff_mat(gen_df)
        #print pw_dif_df
        pw_diff_df.to_csv(args.o,sep="\t")





