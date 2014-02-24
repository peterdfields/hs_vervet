#!/usr/bin/env python
import os, sys
import vcf
import argparse

parser = argparse.ArgumentParser(description="Creates a bed file with the intervals of all indels in the input vcf.")
parser.add_argument("chrom",help="Chromosome name.")
parser.add_argument('in_vcf', type = argparse.FileType('r'), default = '-',help="INDEL vcf filename.")
parser.add_argument("-o","--out-fn",type = argparse.FileType('w'), default = '-',help="name of the outgroup")

args = parser.parse_args()

reader = vcf.Reader(args.in_vcf)


#attention bed has 0-based coordinates while VCF has 1-based coordinates
for record in reader:
    if record.is_indel:
        start = record.POS - 1
        end = start + len(record.REF)
        args.out_fn.write("{}\t{}\t{}\n".format(args.chrom,start,end))
    


