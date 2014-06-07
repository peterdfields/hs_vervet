#!/usr/bin/env python
import os, sys, math
#import pandas as pd
import vcf
import pyfasta

def parse_delta_file(fn):
    """
    Read a nucmer delta file and 
    produce a dictionary that gives the position
    in outgroup reference (chom,pos) if the position (pos) in 
    species reference is given as key.
    """
    query_ref_dic = {}
    with open(fn,'r') as f:
        f.next()
        f.next()
        for line in f:
            if line[0] == ">":
                ref_chrom = line[1:].split()[0]
                new_chrom = True
                new_align = True
            elif new_align:
                pos_ref = int(line.split()[0])
                pos_query = int(line.split()[2])
                end_ref = int(line.split()[1])
                end_query = int(line.split()[3].strip())
                new_align = False
            elif line.strip() == '0':
                query_ref_dic.update({pq:(ref_chrom,pr) for pq,pr in zip(range(pos_query,end_query),
                                                                         range(pos_ref,end_ref))})
                new_align = True
            else:
                next_indel = abs(int(line.strip()))
                in_or_del = math.copysign(1,int(line.strip()))
                query_ref_dic.update({pq:(ref_chrom,pr) for pq,pr in \
                                      zip(range(pos_query,pos_query+next_indel-1),
                                               range(pos_ref,pos_ref+next_indel-1))})
                pos_ref = pos_ref + next_indel - 1  + (1 if in_or_del > 0 else 0)
                pos_query = pos_query + next_indel - 1 + (1 if in_or_del < 0 else 0)
    return query_ref_dic

def write_outgroup_state_to_vcf(in_vcf,out_vcf,query_ref_dic,outgroup_name,outgr_ref_fn):
    
    outgr_ref = pyfasta.Fasta(outgr_ref_fn)
    vcf_reader = vcf.Reader(in_vcf)
    vcf_reader.infos.update({'AA':vcf.parser._Info('AA', '1', 'String', 
                                                   'Ancestral Allele as '
                                                   'derived from {} state.'.format(outgroup_name))})

    vcf_writer = vcf.Writer(args.out_vcf, vcf_reader)

    for record in vcf_reader:
        try:
            outgr_chrom, outgr_pos = query_ref_dic[record.POS-1]
            outgroup_state = outgr_ref[outgr_chrom][outgr_pos]
            record.INFO.update({'AA':outgroup_state})
        except KeyError:
            pass
        vcf_writer.write_record(record)
    vcf_writer.flush()
    vcf_writer.close()




if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Add ancestral allele to INFO column in vcf.")
    parser.add_argument("nucmer_delta_file",help="Nucmer .delta file from species (query) "
                                                 "vs. outgroup (ref) alignment.")
    #parser.add_argument("species_reference",help="Reference fasta file of the target species.")
    parser.add_argument("outgroup_reference",help="Reference fasta file of the outgroup species.")
    parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-',help="Input vcf filename.")
    parser.add_argument("out_vcf",type = argparse.FileType('w'), default = '-',help="Output vcf filename.")
    parser.add_argument("--outgroup_name",required='True',help="name of the outgroup")
    args = parser.parse_args()

    #test that
    for path in [args.nucmer_delta_file,args.outgroup_reference]:
        assert os.path.isfile(path)

    query_ref_dic = parse_delta_file(args.nucmer_delta_file)
    write_outgroup_state_to_vcf(args.in_vcf,args.out_vcf,query_ref_dic,args.outgroup_name,args.outgroup_reference)
