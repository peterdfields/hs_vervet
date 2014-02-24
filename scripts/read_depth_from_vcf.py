#!/usr/bin/env python
import os, sys
import vcf
import argparse

parser = argparse.ArgumentParser(description="Count the total depth across all SNPs in a VCF.")
parser.add_argument("in_vcf",help="Input vcf filename.")

args = parser.parse_args()

vcf_reader = vcf.Reader(filename=os.path.abspath(args.in_vcf))
tot_depth = 0

for (i,record) in enumerate(vcf_reader):
            tot_depth += int(record.INFO['DP'])

sys.stdout.write('{}\t{}\n'.format(i+1,tot_depth))
