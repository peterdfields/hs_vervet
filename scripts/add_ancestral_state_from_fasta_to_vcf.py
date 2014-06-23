#!/usr/bin/env python
import os, sys
import vcf as vcf
from pyfasta import Fasta


def add_ancestral_state(vcf_reader,vcf_writer):
    for i,record in enumerate(vcf_reader):
        if record.var_type == "snp":
            aa = out_grp["Chr"+record.CHROM][record.POS-1]
            if aa == "-":
                pass
            else:
                record.INFO.update({'AA':aa})
            vcf_writer.write_record(record)

def update_vcf_header(vcf_reader):
    vcf_reader.infos.update({'AA':\
                            vcf.parser._Info('AA',
                                             '1', 
                                            'String',
                                             'Ancestral allele as derived from {}'\
                                                            .format(args.outgroup))})

if __name__ == "__main__":
    import argparse
    from hs_vervet.tools import bioparallel

    parser = argparse.ArgumentParser(description="Add ancestral allele to INFO column in vcf. Ancestral state is taken\
                                     from the outgroup in the input fasta.")
    parser.add_argument("outgroup_fasta",help="fasta with outgroup state")
    parser.add_argument("in_vcf",help="Input vcf filename.")
    parser.add_argument("out_vcf",help="Output vcf filename.")
    parser.add_argument("--chrom",nargs='+',help="Chromosome")
    parser.add_argument("--chrom_len",type=int,default=None,nargs='*',help="Chromosome length")
    parser.add_argument("-o","--outgroup",required=True,help="name of the outgroup")
    parser.add_argument("--ncpus",type=int,default=None,help="Number of processed to spawn.")
    parser.add_argument("--tmp_dir",default=".",help="Directory to write temporary files to.")

    args = parser.parse_args()

    out_grp = Fasta(args.outgroup_fasta)

    parser =  bioparallel.VCFParser(args.in_vcf,add_ancestral_state,
                                    chromosomes=args.chrom,chrom_len=args.chrom_len,mode="vcf_write",
                                    out_fn=args.out_vcf,
                                    update_vcf_header_fun=update_vcf_header,
                                    tmp_dir=args.tmp_dir)
    parser.run(ncpus=args.ncpus)
    

