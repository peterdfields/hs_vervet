#!/usr/bin/env python
import os, sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Convert fastq quality scores between illumina <1.8 (phred64) and sanger (phred33) format.")
parser.add_argument("in_fastq",help="Input fastq file.")
parser.add_argument("out_fastq",help="Output fstq file.")
parser.add_argument("-c","--conversion_direction",default="64to33",choices=['64to33','33to64'],help="Direction of the conversion.")

args = parser.parse_args()

if args.conversion_direction == "64to33":
    in_form = "fastq-illumina"
    out_form = "fastq-sanger"
else:
    in_form = "fastq-sanger"
    out_form = "fastq-illumina"

count = SeqIO.convert(args.in_fastq, in_form, args.out_fastq, out_form)
print "Converted {} records from {} to {}".format(count,in_form,out_form)
