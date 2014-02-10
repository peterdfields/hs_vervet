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

def get_absolute_path(fn):
    if fn[0]=='/':
        return fn
    else:
        return os.path.join(os.getcwd(),fn)


coverage_file1 = get_absolute_path(args.file1[0]) 


coverage_files = [get_absolute_path(fn) for fn in args.file2]


total = pd.read_csv(coverage_file1,sep='\t',index_col=[0,1],header=None,names=['depth'],squeeze=True)

for file in coverage_files:
    series = pd.read_csv(file,sep='\t',index_col=[0,1],header=None,names=['depth'],squeeze=True)
    total = total.add(series,fill_value=0).astype(int)
    #del series
    gc.collect()


total.to_csv(args.out,sep='\t',header=None)

