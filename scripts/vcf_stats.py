#!/usr/bin/env python
import os, sys
import vcf
import argparse, json

parser = argparse.ArgumentParser(description="Print some statistics about the variants in the vcf.")
parser.add_argument('in_vcf', type = argparse.FileType('r'), default = '-',help="VCF filename to parse.")
parser.add_argument("-o","--out-fn",type = argparse.FileType('w'), default = '-',help="name of the outgroup")

args = parser.parse_args()

reader = vcf.Reader(args.in_vcf)


var_stats = {'total':0,'snps':0,'indels':0,'other_variants':0,'pass':0,'filters':{k:0 for k in reader.filters.keys()}}

if 'AA' in reader.infos.keys():
    var_stats.update({'ancestral_known':0})
    var_stats.update({'pass_ancestral_known':0})
    var_stats.update({'pass_ancestral_is_ref':0})
    var_stats.update({'pass_ancestral_is_alt':0})
    var_stats.update({'pass_ancestral_third_allele':0})

for record in reader:
    s = var_stats
    s['total']+=1
    if record.is_snp:
        s['snps'] += 1
    elif record.is_indel:
        s['indels'] += 1
    else:
        s['other_variants'] += 1
    if not record.FILTER:
        s['pass'] += 1
    else:
        for ft in record.FILTER:
            s['filters'][ft]+=1
            
    if 'ancestral_known' in s.keys():
        try:
            aa = record.INFO['AA']
        except KeyError:
            continue
        if aa in ['A','C','T','G']:
            s['ancestral_known'] += 1
            if not record.FILTER:
                s['pass_ancestral_known'] += 1
            if aa == record.REF:
                s['pass_ancestral_is_ref'] += 1
            elif aa == record.ALT[0]:
                s['pass_ancestral_is_alt'] += 1
            else:
                s['pass_ancestral_third_allele'] += 1
        else:
            raise ValueError(record.CHROM+' '+str(record.POS)+' alternative allele has unknown state {}'.format(aa))

json.dump(s,args.out_fn)

