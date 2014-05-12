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
parser.add_argument("--ancestral_derived",action = 'store_true', help=  "0,1,2 corresponds to ancestral vs\n"
                                                                        "derived instead of ref vs alt.\n"
                                                                        "Using AA tag from VCF info.\n"
                                                                        "Sites for which this tag is "
                                                                        "unavailable or neither ref \n"
                                                                        "nor alt are ignored.")

args = parser.parse_args()

vcf_reader = vcf.Reader(args.in_vcf)
tf = args.out_fn

allele0 = "0/0"
allele1 = "1/1"


tf.write("chrom"+"\t"+"pos"+"\t"+"\t".join(vcf_reader.samples)+"\n")
for record in vcf_reader:

    if args.ancestral_derived:
        try:
            aa = record.INFO['AA']
        except KeyError:
            continue
        if aa == record.REF:
            allele0 = "0/0"
            allele1 = "1/1"
        elif aa == record.ALT[0]:
            allele0 = "1/1"
            allele1 = "0/0"
        else:
            continue

    tf.write(record.CHROM+"\t"+str(record.POS))

    for sample in record.samples:
        tf.write("\t")
        if sample["GT"] is None:
            tf.write("N")
        elif sample["GT"] == allele0:
            tf.write("0")
        elif sample["GT"] == "0/1":
            tf.write("1")
        elif sample["GT"] == allele1:
            tf.write("2")
        else:
            raise ValueError("Unsupported genotype " + sample["GT"])
    tf.write("\n")

