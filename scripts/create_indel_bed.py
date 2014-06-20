#!/usr/bin/env python
import os, sys
import vcf
import argparse

parser = argparse.ArgumentParser(description="Creates a bed file with the intervals of all indels in the input vcf.")
parser.add_argument('in_vcf', type = argparse.FileType('r'), default = '-',help="INDEL vcf filename.")
parser.add_argument("-o","--out-fn",type = argparse.FileType('w'), default = '-',help="name of the output bed.")
parser.add_argument("-e","--extend_interval",default=0,type=int,help="extend the interval by n bases to left and right")

args = parser.parse_args()

reader = vcf.Reader(args.in_vcf)


#attention bed has 0-based coordinates while VCF has 1-based coordinates
for record in reader:
    if record.is_indel:
        start = max(0, record.POS - 1 - args.extend_interval)
        # I think in the following line we to not add -1 to the contig length
        # because bed files report right open intervals
        end = min(reader.contigs[record.CHROM].length ,start + len(record.REF) + args.extend_interval)
        args.out_fn.write("{}\t{}\t{}\n".format(record.CHROM,start,end))
    


