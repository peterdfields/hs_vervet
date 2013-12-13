#!/usr/bin/env python
import os, sys
import pandas as pd
import vcf
import argparse

parser = argparse.ArgumentParser(description="Add ancestral allele to INFO column in vcf. Ancestral state is taken\
                                 from the outgroup in the input tsv. Input should be produced with the script \
                                 vervet_snp_other_species_state.py ")
parser.add_argument("outgroup_state_tsv",help="tsv with outgroup state produced by vervet_snp_other_species_state.py")
parser.add_argument("in_vcf",help="Input vcf filename.")
parser.add_argument("out_vcf",help="Output vcf filename.")
parser.add_argument("-o","--outgroup",default='macaque',help="name of the outgroup")

args = parser.parse_args()

def get_absolute_path(fn):
    if fn[0]=='/':
        return fn
    else:
        return os.path.join(os.getcwd(),fn)

anc_df_fn = get_absolute_path(args.outgroup_state_tsv)

anc_df = pd.read_csv(anc_df_fn,sep='\t',index_col=0)

vcf_reader = vcf.Reader(filename=get_absolute_path(args.in_vcf))
vcf_reader.infos.update({'AA':vcf.parser._Info('AA', '1', 'String', 'Ancestral Allele as derived from {}'.format(args.outgroup))})
vcf_writer = vcf.Writer(open(get_absolute_path(args.out_vcf), 'w'), vcf_reader)

for record in vcf_reader:
    try:
        aa = anc_df.loc[record.POS-1][args.outgroup]
        if type(aa)==str:
            record.INFO.update({'AA':aa})
    except KeyError:
        pass
    vcf_writer.write_record(record)
vcf_writer.flush()
vcf_writer.close()

