#!/usr/bin/env python
import os, sys
#import vcf
import argparse

parser = argparse.ArgumentParser(description="Add or update filter information in vcf header.")
parser.add_argument("in_vcf",help="Input vcf filename.")
parser.add_argument("out_vcf",help="Output vcf filename.")
parser.add_argument("-f","--filter-expression",default=None, help="name of the outgroup",nargs="+")
parser.add_argument("-d","--filter-description",default=None,help="name of the outgroup",nargs="+")

args = parser.parse_args()

def get_absolute_path(fn):
    if fn[0]=='/':
        return fn
    else:
        return os.path.join(os.getcwd(),fn)

print args


