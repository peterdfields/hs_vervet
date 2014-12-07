#!/usr/bin/env python
"""
Creates 0,1,2 genotype matrix from VCF
"""
#import sys, os
import warnings

#--------support functions------------
def get_012(gt_str):
    if gt_str[:3] in ["0/0","0|0"]:
        return "0"
    elif gt_str[:3] in ["1/1","1|1"]:
        return "2"
    elif gt_str[:3] in ["0/1","0|1","1|0"]:
        return "1"
    elif gt_str[:3] == "./.":
        return "N"
    else:
        raise ValueError("Unsupported genotype " + gt_str)


def get_AA(info_str):
    aa = info_str.split("AA=")[1][0]    
    return aa


def generic_parser(parse_fun,header_fun,vcf_fh, tsv_fh,*args,**kwa):
    """
    parse_fun... applied to each split data line of the vcf
    header_fun... applied to each line in the vcf header
    """
    prev_chrom = None
    prev_pos = -1
    h = None
    p = None
    for line in vcf_fh:
        if line[0] == "#":
            h = header_fun(line, tsv_fh, h)
            continue
        d = line.strip().split("\t")
        chrom = d[0]
        pos = int(d[1])
        if chrom == prev_chrom:
            assert pos >= prev_pos, "vcf positions not in "\
                                    "ascending order at: {}:{},{}".format(chrom,prev_pos,pos) 
            if pos == prev_pos:
                warnings.warn("Warning, multiple entries for pos {}:{}.\n"
                              "Skipping all but the first.".format(chrom,pos))
                continue
        p = parse_fun(d,tsv_fh, p, h, *args,**kwa)
        prev_chrom = chrom
        prev_pos = pos

#--------parsing functions-------------

def vcf_to_012_header(line, tsv_fh,*args):
    if line[1:6] == "CHROM":
        tsv_fh.write("chrom\tpos\t"+line.split("\t",9)[-1])

def vcf_to_012_parse_fun(d, tsv_fh,*args):
    gt = map(get_012,d[9:])
    tsv_fh.write("\t".join(d[:2])+"\t"+"\t".join(gt)+"\n")


#-----------

def ref_alt_anc_header(line, tsv_fh,*args):
    if line[1:6] == "CHROM":
        tsv_fh.write("chrom\tpos\tref\talt\tanc\n")

def ref_alt_anc_parse_fun(d, tsv_fh,*args):
    try:
        aa = get_AA(d[7])
    except IndexError:
        aa = "N"
    tsv_fh.write("\t".join(d[:2])+"\t"+"\t".join(d[3:5])+"\t"+aa+"\n")

#-----------
import itertools

def msmc_header(line, tsv_fh,*args):
    if line[1:6] == "CHROM":
        return line[1:].strip().split()

def get_alleles(allele_dic,gt_str):
    allele0 = allele_dic[gt_str[0]]
    allele1 = allele_dic[gt_str[2]]
    phased = (True if gt_str[1]=="|" else False)
    return (allele0, allele1, phased)

def msmc_input_parse_fun(d, tsv_fhs, p, h, ind_tuples):
    """
    *attention, input must be a VCF that includes all sites,
     segregating and non-segregating, filtered and non-filtered
    *attention, tsv_fh can be a list of multiple file handles
    here, one per output file to be written
    """
    try:
        len(tsv_fhs)
    except TypeError:
        tsv_fhs = [tsv_fhs]
    assert len(ind_tuples) == len(tsv_fhs)
    if p is None:
        p = [0]*len(ind_tuples)
    #check for filter
    #if d[6] or len(d[3])>1 or len(d[4])>1:
    #    print "filtered", d[6]
    #    return p
    #else:
    #    p += 1
    if d[4] in "ACGT": #if SNP
        if d[6] in ["PASS","","."]: #check whether not filtered
            allele_dic = {"0":d[3],"1":d[4],".":"?"}
            for j,(tsv_fh,ids) in enumerate(zip(tsv_fhs,ind_tuples)):
                genotypes = [get_alleles(allele_dic,d[h.index(id)]) for id in ids]
                gts = []
                for gt in genotypes:
                   if gt[2] or (gt[0]=='?' and gt[1]=='?') or (gt[0] == gt[1]):
                        gts.append([gt[0]+gt[1]])
                   else:
                        gts.append([gt[0]+gt[1],gt[1]+gt[0]])
                gt_str = ",".join(["".join(i) for i in itertools.product(*gts)])
                if "?" in gt_str: #don't print site with missing gt
                    pass
                elif gt_str == len(gt_str) * gt_str[0]: #don't print if all alleles equal
                    p[j] += 1
                else:
                    p[j] += 1
                    tsv_fh.write(d[0]+"\t"+d[1]+"\t"+str(p[j])+"\t"+gt_str+"\n")
                    p[j] = 0
            return p
        else: 
            return p
    else:
        if d[6] in ["PASS","",".","LowQual"] and len(d[3])==1 and (len(d[4])==1 or d[4]=="None"):
            p = [i+1 for i in p]
        return p


if __name__ == "__main__":
    import argparse
    import gzip

    parser = argparse.ArgumentParser(description="Extract data from vcf file.")
    parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-', help="Input vcf filename.")
    parser.add_argument("out_fn",type = argparse.FileType('w'), default = '-', help="Output tsv filename")
    #parser.add_argument("--ancestral_derived",action = 'store_true', help=  "0,1,2 corresponds to ancestral vs\n"
    #                                                                        "derived instead of ref vs alt.\n"
    #                                                                        "Using AA tag from VCF info.\n"
    #                                                                        "Sites for which this tag is "
    #                                                                        "unavailable or neither ref \n"
    #                                                                        "nor alt are ignored.")
    subparsers = parser.add_subparsers(dest='mode',
                                        help='E.g. "to012", or "ref_alt_anc"')
    parser0 = subparsers.add_parser("to012",help="Convert a vcf into a tsv with 0,1,2 genotypes.\n"
                                                   "0..reference,1..het,2..alt.")

    parser1 = subparsers.add_parser("ref_alt_anc",help="Extract tsv with chrom, pos, "
                                                       "ref_allele, alt_allele, anc_allele. "
                                                       "Requires annotation AA in info column")

    parser2 = subparsers.add_parser("msmc",help="Create msmc input files from VCF which includes "
                                                                        "fixed and filtered sites.")
    parser2.add_argument("--add_out_fns",nargs="+",type = argparse.FileType('w'),
                                                help="List of additional output files used with"
                                                                    " second to last ind_tuples")
    parser2.add_argument("--ind_tuples",required=True,nargs='+',help="List of individual ids that should be "
                                                            "outputted into the output files."
                                                        "Each n-tuple of ids for the same file should be "
                                                        "quoted in the command line, e.g., 'id1 id2 id3' 'id2 id4') "
                                                        "produces two output files.")

    args = parser.parse_args()
    #for k,v in vars(args).iteritems():
    #    print k,v,
    if args.mode == "to012":
        generic_parser(vcf_to_012_parse_fun,vcf_to_012_header,
                                        args.in_vcf, args.out_fn)
    elif args.mode == "ref_alt_anc":
        generic_parser(ref_alt_anc_parse_fun, ref_alt_anc_header,
                                           args.in_vcf, args.out_fn)
    elif args.mode == "msmc":
        out_fns = [args.out_fn] if args.add_out_fns is None else [args.out_fn]+args.add_out_fns
        assert len(out_fns) == len(args.ind_tuples), "There are {} out files but {} tuples.".format(len(out_fns),len(args.ind_tuples))
        generic_parser(msmc_input_parse_fun, msmc_header, args.in_vcf, out_fns,
                                                ind_tuples = [it.split() for it in args.ind_tuples])

    else:
        raise UserException("Unknown mode.")
