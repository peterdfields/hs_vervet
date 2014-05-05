#!/usr/bin/env python
"""
Calculates the relative difference between expected
and observed heterozygosity.
"""
import sys, os, math
import argparse

parser = argparse.ArgumentParser(description=   "Calculates the relative difference between\n"
                                                "expected and observed heterozygosity.")
parser.add_argument("in_012_tsv",type = argparse.FileType('r'), default = '-', help="Input genotype-file filename.")
parser.add_argument("out_fn",type = argparse.FileType('w'), default = '-', help="Output tsv filename.")

args = parser.parse_args()

args.out_fn.write("chrom\tpos\t(Ho-He)/He\n")

args.in_012_tsv.readline()

i = 0
for line in args.in_012_tsv:
    dat = line.strip().split("\t")
    chrom = dat[0]
    pos = dat[1]
    genotypes = dat[2:]
    nans = genotypes.count("N")
    hom_ref = genotypes.count("0")
    het = genotypes.count("1")
    hom_alt = genotypes.count("2")
    #print "info",nans, hom_ref,het, hom_alt
    assert (nans+hom_ref+het+hom_alt) == len(genotypes), genotypes
    af = (hom_alt+het/2.)*1./(hom_alt+het+hom_ref) if (hom_alt+het+hom_ref) > 0 else float("NaN")
    Ho = het*1./(hom_alt+het+hom_ref) if (hom_alt+het+hom_ref) > 0 else float("NaN")
    He = 2.*af*(1-af)
    deviation = (Ho-He)*1./He if He > 0 else float("NaN")
    if not math.isnan(deviation):
        args.out_fn.write("{}\t{}\t{}\n".format(chrom,pos,deviation))
        #"entry written " + "{}\t{}\t{}\n".format(chrom,pos,deviation)
        print chrom,pos,af,He,Ho,deviation
    else:
        print "NaN:",chrom,pos,af,He,Ho,deviation
    i+=1
    if i>50:
        break




