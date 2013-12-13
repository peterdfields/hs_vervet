#!/usr/bin/env python
"""
this script adds a Filter to the filter column of a VCF files if a variant is adjacent to another variant
"""
import sys, os
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

def get_absolute_path(fn):
    if fn[0]=='/':
        return fn
    else:
        return os.path.join(os.getcwd(),fn)


in_vcf_fn = get_absolute_path(sys.argv[1])
out_vcf_fn = get_absolute_path(sys.argv[2])
#minimum distance between variants, I first use 1bp 
min_dist = int(sys.argv[3])

def add_filter(line):
    if "AdjacentSNP" not in line[6]:
        if line[6] == "." or line[6] == "PASS":
            line[6] = "AdjacentSNP"
        else:
            line[6] += ";AdjacentSNP"


with open(in_vcf_fn,'r') as in_vcf:
    with open(out_vcf_fn,'w') as out_vcf:
        last_line=''
        current_line=''
        for line in in_vcf:
            if line[0] == '#':
                out_vcf.write(line)
            else:
                last_line = current_line
                current_line = line.split('\t')
                if last_line and int(last_line[1])+min_dist >= int(current_line[1]):
                    add_filter(last_line)
                    add_filter(current_line)
                out_vcf.write("\t".join(last_line))
        out_vcf.write("\t".join(current_line))
