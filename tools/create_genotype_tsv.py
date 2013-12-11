#!/usr/bin/env python
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import vcf
import pandas as pd


#def 

#reader = vcf.Reader(open(vcf_fname, 'r'))



if __name__ == '__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Parse VCF into a 0/1 genotype tsv file.")
    parser.add_argument('--input', type = argparse.FileType('r'), default = '-')
    #parser.add_argument("input_vcf",help="VCF filename to parse.")
    #parser.add_argument("output_tsv",help="Filename of tsv to output.")
    parser.parse_args()
    print parser.input
