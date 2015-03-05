#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.
Todo:
-- parallel
-- pack into classes
"""
import sys
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)



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

class VCFParser(object):
    """
    Parse a VCF file.
    """
    def __init__(self,vcf_fh, parse_fun=None, header_fun=None, cleanup_fun=None,
                    parse_arg_dic=None, header_arg_dic=None, cleanup_arg_dic=None,
                                                    skip_multiple_entries=True,**kwa):
        self.vcf_fh = vcf_fh
        if parse_fun is None:
            parse_fun = lambda *x,**y: None
        self.parse_fun = parse_fun
        if header_fun is None:
            header_fun = lambda *x,**y: None
        self.header_fun = header_fun
        if cleanup_fun is None:
            cleanup_fun = lambda *x,**y: None
        self.cleanup_fun = cleanup_fun
        if parse_arg_dic is None:
            parse_arg_dic = {}
        self.parse_arg_dic = parse_arg_dic
        if header_arg_dic is None:
            header_arg_dic = {}
        self.header_arg_dic = header_arg_dic
        if cleanup_arg_dic is None:
            cleanup_arg_dic = {}
        self.cleanup_arg_dic = cleanup_arg_dic

        self.skip_multiple_entries = skip_multiple_entries

    def parse(self):
        """
        header_out and parser_out
        can be used to store, append and reuse information.
        If information from previous lines should be propagated,
        make sure that header_out and parser_out are mutable types
        that are changed by side-effects
        such as dicts or lists.
        """
        header_out = {}
        parser_out = {}
        prev_chrom = None
        prev_pos = -1
        logging.info("Starting parsing.")
        header_line_nr = 0
        for i,line in enumerate(self.vcf_fh):
            #logging.debug(line)
            if line[0] == "#":
                header_line_nr = i
                #logging.debug("Parsing header line {}".format(i))
                self.header_fun(line,header_out=header_out,**self.header_arg_dic)
                continue
            #logging.debug("Parsing main line {}".format(i-header_line_nr))
            d = line.strip().split("\t")
            chrom = d[0]
            pos = int(d[1])
            if chrom == prev_chrom:
                assert pos >= prev_pos, "vcf positions not in "\
                                        "ascending order at: {}:{},{}".format(chrom,prev_pos,pos)
                if pos == prev_pos:
                    if not skip_multiple_entries:
                        logging.warning("Multiple entries for pos {}:{}.\n"
                                  "Keeping all entries.".format(chrom,pos))
                    else:
                        logging.warning("Warning, multiple entries for pos {}:{}.\n"
                                  "Skipping all but the first.".format(chrom,pos))
                        continue
            self.parse_fun(d,header_out=header_out,parser_out=parser_out,**self.parse_arg_dic)
            prev_chrom = chrom
            prev_pos = pos
            #if i>10000:
            #    break
        self.parser_out = parser_out

    def cleanup(self):
        """
        """
        logging.info("Starting cleanup.")
        self.cleanup_fun(self.parser_out,**self.cleanup_arg_dic)


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

def msmc_input_parse_fun(d, tsv_fhs, p, h, ind_tuples,haplotypes):
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
    if d[4] in ['A','C','G','T']: #if SNP
        if d[6] in ["PASS","","."] and len(d[3])==1: #check whether not filtered or deletion
            allele_dic = {"0":d[3],"1":d[4],".":"?"}
            for j,(tsv_fh,ids) in enumerate(zip(tsv_fhs,ind_tuples)):
                genotypes = [get_alleles(allele_dic,d[h.index(id)]) for id in ids]
                phases = [g[2] for g in genotypes]
                if haplotypes==0:
                    genotype_strs = [g[0] for g in genotypes]
                elif haplotypes==1:
                    genotype_strs = [g[1] for g in genotypes]
                elif haplotypes==2:
                    genotype_strs = [g[0]+g[1] for g in genotypes]
                gts = []
                for gt, phase in zip(genotype_strs,phases):
                   if phase or ('?' in gt) or (gt == gt[0]*len(gt)):
                        gts.append([gt])
                   else:
                        gts.append([gt,gt[::-1]])
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

#--------------


def add_outgroup_header(line, out_vcf_fh,*args, **kwa):
    #if line[:6]=='#CHROM':
    #    line = line.strip() + "\t" + sample_name + '\n'
    out_vcf_fh.write(line)



def add_outgroup_parse_fun(line,out_vcf_fh, *args,**kwa):
    """
    """
    bases = ["A",'C','T','G']
    no_var = ['.','X']
    og_base = og_fasta[line[0]][int(line[1])-1]
    if line[4] in no_var and og_base.upper() in bases and og_base.upper() != line[3]:
        line[4] = og_base.upper()
        #print "added SNP at", line[0], line[1]
    #print out_vcf_fh
    #print "\t".join(line)
    out_vcf_fh.write("\t".join(line)+'\n')


#-----------------------------


def get_filter_stats_parse_fun(line,parser_out,**kwa):
    def add_to_countdic(dic,key):
        try:
            dic[key] += 1
        except KeyError:
            dic[key] = 1
    filters = line[6].split(';')
    filters.sort()
    filters = tuple(filters)
    add_to_countdic(parser_out,'n_sites')
    add_to_countdic(parser_out,filters)

def get_filter_stats_cleanup_fun(countdic,out_fn=sys.stdout):
    import pandas as pd
    filter_info = pd.Series(countdic.values(),index=countdic.keys())
    filter_info.sort(inplace=True,ascending=False)
    filter_info.to_csv(out_fn,sep='\t')
    #json.dump(countdic,open(out_fn,'w'))


#---MAIN----FUNCTION-------------


def run(vcf_fh, parse_fun=None, header_fun=None, cleanup_fun=None,
        parse_arg_dic=None, header_arg_dic=None, cleanup_arg_dic=None,
                                        no_skip_multiple_entries=False):


    parser = VCFParser(vcf_fh,parse_fun,header_fun,cleanup_fun,
                        parse_arg_dic,header_arg_dic,cleanup_arg_dic,
                    skip_multiple_entries=not no_skip_multiple_entries)
    parser.parse()
    parser.cleanup()

    return parser.parser_out

if __name__ == "__main__":
    #import gzip
    import argparse
    parser = argparse.ArgumentParser(description="Parse a Variant Call Format (VCF) file.")
    parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-', help="Input vcf filepath.")
    parser.add_argument("--parse_fun",help="Name of function to be used to parse the vcf body.")
    parser.add_argument("--header_fun",help="Name of function to be used to parse the vcf header.")
    parser.add_argument("--cleanup_fun",help="Name of function to be used to handle the parser output.")
    parser.add_argument("--parse_kwa",help='key=value pairs to pass to parse fun. E.g., "foo=bar;test=hallo du"')
    parser.add_argument("--header_kwa",help='key=value pairs to pass to header fun. E.g., "foo=bar;test=hallo du"')
    parser.add_argument("--cleanup_kwa",help='key=value pairs to pass to cleanup fun. E.g., "foo=bar;test=hallo du"')
    parser.add_argument("--no_skip_multiple_entries",action='store_true',
                         help='Do not skip all but the first entry for the same site in the VCF.')


    args = parser.parse_args()


    if args.parse_fun is not None:
        pf = eval(args.parse_fun)
    else:
        pf = None
    if args.header_fun is not None:
        hf = eval(args.header_fun)
    else:
        hf = None
    if args.cleanup_fun is not None:
        cf = eval(args.cleanup_fun)
    else:
        cf = None
    if args.parse_kwa is not None:
        parse_arg_dic = {k:v for (k,v) in [t.split('=') for t in args.parse_kwa.split(';')]}
    else:
        parse_arg_dic = None
    if args.header_kwa is not None:
        header_arg_dic = {k:v for (k,v) in [t.split('=') for t in args.header_kwa.split(';')]}
    else:
        header_arg_dic = None
    if args.cleanup_kwa is not None:
        cleanup_arg_dic = {k:v for (k,v) in [t.split('=') for t in args.cleanup_kwa.split(';')]}
    else:
        cleanup_arg_dic = None

    run(args.in_vcf,pf,hf,cf,
          parse_arg_dic,header_arg_dic,cleanup_arg_dic,args.no_skip_multiple_entries)
    
