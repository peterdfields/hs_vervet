#!/usr/bin/env python
import sys, os
import pandas as pd
import numpy as np


#define parse function
def filtered_missing_stats(reader):
    snp = ["A","G","C","T"]
    novar = [None,"X",'.']
    sites_dic = {"total":0,"pass_nosnp":0,"pass_snp":0,"filter_nosnp":0,"filter_snp":0}
    prev_pos = 0
    N_df = pd.DataFrame(0,columns=sites_dic.keys(),index=reader.samples)
    Nxy = np.zeros((len(reader.samples),len(reader.samples)))

    for rec in reader:
        while rec.POS == prev_pos:  #jump over indel entries or duplicated entries
            rec = reader.next()
        sites_dic["total"] += 1
        ns = np.array([1 if s.data.GT is None else 0 for s in rec.samples]) #vector of Ns
        N_df["total"] +=  ns

        Nxy += np.outer(ns,ns)
        try:
            alt = rec.ALT[0].sequence
        except AttributeError:
            alt = rec.ALT[0] 
        is_variant = not alt in novar
        if is_variant:
            if not rec.FILTER:
                category = "pass_snp"

            else:
                category = "filter_snp"
        else:
            if not rec.FILTER or (len(rec.FILTER)==1 and ("LowQual" in rec.FILTER) and rec.QUAL>5):
                category = "pass_nosnp"
            else:
                category = "filter_nosnp"
        sites_dic[category] += 1
        N_df[category] += ns
        prev_pos = rec.POS
    
    Nxy = pd.DataFrame(Nxy,index=reader.samples,columns=reader.samples)
    return (sites_dic, N_df, Nxy)


def reduce_stats(stats_ls):
    sites_dic = {k:sum([d[0][k] for d in stats_ls]) for k in stats_ls[0][0].keys()}
    N_df = sum([j[1] for j in stats_ls])
    nNxy = sum([j[2] for j in stats_ls])
    #correlation
    corr = (Nxy-1./sites_dic["total"]*np.outer(N_df["total"],N_df["total"]))/\
            np.sqrt(np.outer(N_df["total"]*(1-1./sites_dic["total"]*N_df["total"]),N_df["total"]*(1-1./sites_dic["total"]*N_df["total"])))
    return sites_dic, N_df, corr

if __name__ == "__main__":
    import argparse
    import json
    from hs_vervet.tools import bioparallel

    parser = argparse.ArgumentParser(description="Produce statistics on filtered sites "
                                                 "and missing genotypes in a VCF.\n"
                                                 "Outputs to different files:\n"
                                                 "(1) the number of filtered/nonfiltered "
                                                 "SNPs and non-SNPs. \n"
                                                 "(2) for each individuals the number of Ns "
                                                 "in each of these categories\n"
                                                 "(3) matrix of the correlation of Ns across individuals.")
    
    parser.add_argument("in_vcf",help="Input vcf filename.")
    parser.add_argument("out_filter_count",type=argparse.FileType('w'),help="Output filter count filename.")
    parser.add_argument("out_N_count",type=argparse.FileType('w'),help="Output N count filename.")
    parser.add_argument("out_N_corr",type=argparse.FileType('w'),help="Output filter count filename.")

    #parallel specific arguments
    parser.add_argument("--chrom",nargs='+',help="Chromosome")
    parser.add_argument("--chrom_len",type=int,default=None,nargs='*',
                                                    help="Chromosome length")
    parser.add_argument("--ncpus",type=int,default=None,
                            help="Number of processed to spawn.")
    parser.add_argument("--tmp_dir",default=".",
                        help="Directory to write temporary files to.")

    args = parser.parse_args()

    parser =  bioparallel.VCFParser(args.in_vcf,filtered_missing_stats,
                                    chromosomes=args.chrom,
                                    chrom_len=args.chrom_len,
                                    mode="count",
                                    reduce_fun=reduce_stats,
                                    tmp_dir=args.tmp_dir)
    parser.run(ncpus=args.ncpus)

    #save results to file
    json.dump(parser.out[0],args.out_filter_count)
    parser.out[1].to_csv(args.out_N_count,sep="\t")
    parser.out[2].to_csv(args.out_N_corr,sep="\t")


