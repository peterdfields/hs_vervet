#!/usr/bin/env python
"""
Creates 0,1,2 genotype matrix from VCF
"""
import sys, os
import vcf
import argparse

parser = argparse.ArgumentParser(description="Convert a vcf into a tsv with 0,1,2 genotypes.\n"
                                               "0..reference,1..het,2..alt.")
parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-', help="Input vcf filename.")
parser.add_argument("out_fn",type = argparse.FileType('w'), default = '-', help="Output tsv filename")

args = parser.parse_args()

vcf_reader = vcf.Reader(args.in_vcf)
tf = args.out_fn

tf.write("chrom"+"\t"+"pos"+"\t"+"\t".join(vcf_reader.samples)+"\n")
for record in vcf_reader:
    tf.write(record.CHROM+"\t"+str(record.POS))
    for sample in record.samples:
        tf.write("\t")
        if sample["GT"] is None:
            tf.write("N")
        elif sample["GT"] == "0/0":
            tf.write("0")
        elif sample["GT"] == "0/1":
            tf.write("1")
        elif sample["GT"] == "1/1":
            tf.write("2")
        else:
            raise ValueError("Unsupported genotype " + sample["GT"])
    tf.write("\n")

