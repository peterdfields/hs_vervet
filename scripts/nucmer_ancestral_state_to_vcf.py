#!/usr/bin/env python
import os, sys, math
#import pandas as pd
import vcf
import pyfasta

def complement(base):
    bases = "ACGTN"
    comp = "TGCAN"
    return comp[bases.index(base)]

def print_align(str0,str1,line_length):
    j=-1
    for j in range(len(str0)/line_length):
        print str0[j*line_length:(j+1)*line_length]
        print str1[j*line_length:(j+1)*line_length]
        print
    print str0[(j+1)*line_length:]
    print str1[(j+1)*line_length:]
    print

def parse_delta_file(fn,ref_fasta,output_aligns=False,query_fasta=None):
    """
    Read a nucmer delta file and.
    produce a dictionary that gives the position
    in outgroup reference (chom,pos) if the position (pos) in.
    species reference is given as key.
    """
    c = complement
    ix = 0
    query_ref_dic = {}
    if output_aligns:
        aligns = []
        start_ref_ls = []
        start_query_ls = []
        q = query_fasta
    with open(fn,'r') as f:
        f.next()
        f.next()
        for line in f:
            if line[0] == ">":
                #print "starting alignment",line
                ref_chrom = line[1:].split()[0]
                r = ref_fasta[ref_chrom]
                new_chrom = True
                new_align = True
            elif new_align:
                #print "new alinment",line
                if output_aligns:
                    macaque_seq = ""
                    vervet_seq = ""
                al_dic = {}          
                start_ref = int(line.split()[0])
                start_query = int(line.split()[2])
                end_ref = int(line.split()[1])
                end_query = int(line.split()[3].strip())
                assert start_ref < end_ref, 
                        "referse strand in reference not implemented"
                ref_indices = range(start_ref-1,end_ref)
                pos_ref = 0

                if start_query < end_query:
                    query_indices = range(start_query-1,end_query)
                    query_reverse = False 
                else:
                    query_indices = range(end_query-1,start_query)[::-1]
                    query_reverse = True
                pos_query = 0
                new_align = False
            elif line.strip() == '0':
                macaque_str = "".join([(r[pr] if not query_reverse \
                                                else c(r[pr]))
                                                    for pr in ref_indices[pos_ref:]])
                al_dic.update({pq:br for pq,br in zip(query_indices[pos_query:],
                                                                    macaque_str)})
                query_ref_dic.update(al_dic)
                if output_aligns:
                    macaque_seq += macaque_str
                    vervet_seq += "".join([q[pq] for pq in \
                                    query_indices[pos_query:]])
                    aligns.append((macaque_seq,vervet_seq))
                    start_ref_ls.append(start_ref)
                    start_query_ls.append(start_query)
                new_align = True
            else:
                next_indel = abs(int(line.strip()))-1
                #
                in_or_del = math.copysign(1,int(line.strip()))
                macaque_str = "".join([(r[pr] if not query_reverse \
                                                else c(r[pr]))
                            for pr in ref_indices[pos_ref:pos_ref+next_indel]])
                al_dic.update({pq:br for pq,br in \
                                      zip(query_indices[pos_query:pos_query+next_indel],
                                               macaque_str)})
                if output_aligns:
                    macaque_seq += macaque_str
                    vervet_seq += "".join([q[pq] for pq in \
                                    query_indices[pos_query:pos_query+next_indel]])
                    macaque_seq += (r[ref_indices[pos_ref+ next_indel]]\
                                    if not query_reverse else \
                                    c(r[ref_indices[pos_ref+ next_indel]]))\
                                                       if in_or_del > 0 else "."
                    vervet_seq += q[query_indices[pos_query+ next_indel]] \
                                                if in_or_del < 0 else "."
                pos_ref = pos_ref + next_indel + (1 if in_or_del > 0 else 0)
                pos_query = pos_query + next_indel + (1 if in_or_del < 0 else 0)
    if output_aligns:
        return (query_ref_dic,start_ref_ls,start_query_ls,aligns)
    else:
        return query_ref_dic



def write_outgroup_state_to_vcf(in_vcf,out_vcf,query_ref_dic,outgroup_name):
    vcf_reader = vcf.Reader(in_vcf)
    vcf_reader.infos.update({'AA':vcf.parser._Info('AA', '1', 'String', 
                                                   'Ancestral allele as '
                                                   'derived from {} state.'.format(outgroup_name))})

    vcf_writer = vcf.Writer(args.out_vcf, vcf_reader)

    for record in vcf_reader:
        try:
            outgroup_state = query_ref_dic[record.POS-1]
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
    parser.add_argument("--debug",action="store_true")
    args = parser.parse_args()

    #test that
    for path in [args.nucmer_delta_file,args.outgroup_reference]:
        assert os.path.isfile(path)

    outgr_ref = pyfasta.Fasta(args.outgroup_reference)
    query_ref_dic = parse_delta_file(args.nucmer_delta_file,outgr_ref)

    write_outgroup_state_to_vcf(args.in_vcf,args.out_vcf,
                                query_ref_dic,args.outgroup_name)
