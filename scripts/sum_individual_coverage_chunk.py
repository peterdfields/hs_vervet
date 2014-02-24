#!/usr/bin/env python
import pandas as pd
import sys, os, gc
import argparse

parser = argparse.ArgumentParser(usage='%(prog)s [-h] [-o OUT] coverage_file coverage_file [coverage_file ...]',
                            description="Read several files tab-delimited files, "\
                                 "where the columns are chrom,pos,coverage. "\
                                 "Output a file with the total coverage for each site, "\
                                 "summed over all files. The input can be produced "
                                 "by 'samtools depth'")
parser.add_argument('file1', nargs=1, metavar='coverage_file',help="Files as produced by samtools depth.")
parser.add_argument('file2', nargs='+', metavar='coverage_file', help=argparse.SUPPRESS)
parser.add_argument("-o","--out",default=sys.stdout,help="Output filename. (default=stdout)")

args = parser.parse_args()


coverage_file1 = os.path.abspath(args.file1[0]) 


coverage_files = [os.path.abspath(fn) for fn in args.file2]


tp = pd.read_csv(coverage_file1,sep='\t',index_col=[0,1],header=None,iterator=True, chunksize=1000)
total = pd.concat(tp, ignore_index=True)

for file in coverage_files:
    tp1 = pd.read_csv(file,sep='\t',index_col=[0,1],header=None,iterator=True, chunksize=1000)
    series = pd.concat(tp1, ignore_index=True)
    total = total.add(series,fill_value=0).astype(int)
    del series
    gc.collect()


total.to_csv(os.path.abspath(args.out),sep='\t',header=None)

sys.exit(0)
