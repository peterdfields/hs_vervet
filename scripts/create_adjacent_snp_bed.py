#!/usr/bin/env python
import os, sys
import vcf
import argparse

parser = argparse.ArgumentParser(description="Creates a bed file with the intervals of all indels in the input vcf.")
parser.add_argument('in_vcf', type = argparse.FileType('r'), default = '-',help="INDEL vcf filename.")
parser.add_argument("-o","--out-fn",type = argparse.FileType('w'), default = '-',help="name of the output bed.")
parser.add_argument("-e","--extend_interval",default=0,type=int,help="extend the interval by n bases to left and right")
parser.add_argument("--min_dist",default=2,type=int,help="Minimum distance between SNPs not to be masked.\n"
                                                           "default=2, meaning only directly adjacent snps included"
                                                                                                        " in the mask")

args = parser.parse_args()

reader = vcf.Reader(args.in_vcf)


#attention bed has 0-based coordinates while VCF has 1-based coordinates
last_record = None
intv_start = None
intv_end = None
for record in reader:
    if last_record is not None:
        if record.CHROM != last_record.CHROM:
            if intv_start is not None:
                intv_end = last_record.POS
        elif record.POS - last_record.POS < args.min_dist:
            if intv_start is None:
                intv_start = last_record.POS
        elif intv_start is not None:
            intv_end = last_record.POS
        if intv_start is not None and intv_end is not None:
            #convert start/end to 0-indexed right open bed format
            bed_start = max(0, intv_start - 1 - args.extend_interval)
            bed_end = min(reader.contigs[last_record.CHROM].length, intv_end + args.extend_interval)
            intv_start = None
            intv_end = None
            args.out_fn.write("{}\t{}\t{}\n".format(last_record.CHROM,bed_start,bed_end))
    last_record = record
 


