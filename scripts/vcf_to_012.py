#!/usr/bin/env python
"""
Creates 0,1,2 genotype matrix from VCF
"""
#import sys, os



def get_012(gt_str):
    if gt_str[:3] in ["0/0","0|0"]:
        return "0"
    elif gt_str[:3] in ["1/1","1|1"]:
        return "2"
    elif gt_str[:3] in ["0/1","0|1","1|0"]:
        return "1"
    elif gt_str[:3] == "./.":
        return "N"
    else:
        raise ValueError("Unsupported genotype " + gt_str)

def vcf_to_012(vcf_fh, tsv_fh):
    for line in vcf_fh:
        if line[0] == "#":
            if line[1:6] == "CHROM":
                tsv_fh.write("chrom\tpos\t"+line.split("\t",9)[-1])
            continue
        d = line.split("\t")
        gt = map(get_012,d[9:])
        tsv_fh.write("\t".join(d[:2])+"\t"+"\t".join(gt)+"\n")

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description="Convert a vcf into a tsv with 0,1,2 genotypes.\n"
                                                   "0..reference,1..het,2..alt.")
    parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-', help="Input vcf filename.")
    parser.add_argument("out_fn",type = argparse.FileType('w'), default = '-', help="Output tsv filename")
    #parser.add_argument("--ancestral_derived",action = 'store_true', help=  "0,1,2 corresponds to ancestral vs\n"
    #                                                                        "derived instead of ref vs alt.\n"
    #                                                                        "Using AA tag from VCF info.\n"
    #                                                                        "Sites for which this tag is "
    #                                                                        "unavailable or neither ref \n"
    #                                                                        "nor alt are ignored.")

    args = parser.parse_args()
    vcf_to_012(args.in_vcf, args.out_fn)
