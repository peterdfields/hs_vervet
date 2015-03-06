#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.
Todo:
-- parallel
-- pack into classes
"""
import sys, os
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
eu = os.path.expanduser

analyses = {}
def add_analysis(name,funs,info='',always_req_params=None,
                    command_line_req_params=None,
                    opt_params=None):
    """
    Add analysis to dict of analyses.
    Input:
        funs: dict of {name:function} with
    """
    entry = {'funs':funs,
             'info':info,
               'command_line_req_params':\
                 command_line_req_params if command_line_req_params is not None else {},
               'always_req_params':\
                 always_req_params if always_req_params is not None else {},
               'opt_params': opt_params if opt_params is not None else {}}
    analyses.update({name:entry})



def check_params(arg_dic,req_param_dic):
    """
    Check whether params noted in req_param_dic are 
    present in arg dic
    """
    for k in req_param_dic:
        if k not in arg_dic:
            raise TypeError('Required parameter {} missing from arg_dic.'.format(k))

def get_parser(vcf_fh,analysis_name,arg_dic,skip_multiple_entries=None):
    check_params(arg_dic,analyses[analysis_name]['always_req_params'])
    return VCFParser(vcf_fh,arg_dic=arg_dic,skip_multiple_entries=skip_multiple_entries,
                                                **analyses[analysis_name]['funs'])



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

#---------main object-------------------

class VCFParser(object):
    """
    Generic parser for VCF file.
    Input:
    vcf_fh ... file handle to input vcf file
    parse_fun ... Function applied to each line of the vcf body.
                  Takes as input a split line (list) and two
                    dictionaries with additional parameters and variables,
                    e.g, to count something.
                    The first dictionary is its 
    header_fun ... Function applied to each line of the vcf header.
                   Takes as input a line (string) abd a dictionary
                    with additional parameters and variables.
    setup_fun ... Function applied to all input dictionaries to modify 
                   them if needed, e.g., to initialise the variables used
                    by parse fun.
    parse_fun ... Function applied to each line of the vcf body.

    Examples:
    
    

    """
    def __init__(self,vcf_fh, parse_fun=None, header_fun=None,
                               setup_fun=None, cleanup_fun=None,
                                arg_dic=None,
                        #parse_arg_dic=None, header_arg_dic=None,
                        #setup_arg_dic=None, cleanup_arg_dic=None,
                                    skip_multiple_entries=True,**kwa):
        self.vcf_fh = vcf_fh
        if parse_fun is None:
            parse_fun = lambda *x,**y: None
        self.parse_fun = parse_fun
        if header_fun is None:
            header_fun = lambda *x,**y: None
        self.header_fun = header_fun
        if setup_fun is None:
            setup_fun = lambda *x,**y: None
        self.setup_fun = setup_fun
        if cleanup_fun is None:
            cleanup_fun = lambda *x,**y: None
        self.cleanup_fun = cleanup_fun
        if arg_dic is None:
            arg_dic = {}
        self.arg_dic = arg_dic

        self.skip_multiple_entries = skip_multiple_entries

    def setup(self):
        logging.info("Starting setup.")
        self.setup_fun(self.arg_dic)
        #self.arg_dic = arg_dic

    def parse(self):
        """
        """
        prev_chrom = None
        prev_pos = -1
        logging.info("Starting parsing:.")
        multiple_count = 0
        before_header = True
        before_body = True
        for i,line in enumerate(self.vcf_fh):
            #logging.debug(line)
            if line[0] == "#":
                if before_header:
                    logging.info("Parsing header.")
                    before_header = False
                self.header_fun(line,self.arg_dic)
                continue
            else:
                if before_body:
                    logging.info("Parsing main body.")
                    before_body = False
            #logging.debug("Parsing main line {}".format(i-header_line_nr))
            d = line.strip().split("\t")
            chrom = d[0]
            pos = int(d[1])
            if chrom == prev_chrom:
                assert pos >= prev_pos, "vcf positions not in "\
                                        "ascending order at: {}:{},{}".format(chrom,prev_pos,pos)
                if pos == prev_pos:
                    multiple_count += 1
                    if mutiple_count > 100:
                        logging.warning("Omitting further multiple entry warnings.")
                    else:
                        if not skip_multiple_entries:
                            logging.warning("Multiple entries for pos {}:{}.\n"
                                      "Keeping all entries.".format(chrom,pos))
                        else:
                            logging.warning("Warning, multiple entries for pos {}:{}.\n"
                                  "Skipping all but the first.".format(chrom,pos))
                            continue
            self.parse_fun(d,self.arg_dic)
            prev_chrom = chrom
            prev_pos = pos
            #if i>10000:
            #    break

    def cleanup(self):
        """
        """
        logging.info("Starting cleanup.")
        self.result = self.cleanup_fun(self.arg_dic)


    def run(self):
        self.setup()
        self.parse()
        self.cleanup()
        #return self.result



#------support functions----------------


def add_to_countdic(dic,key):
    try:
        dic[key] += 1
    except KeyError:
        dic[key] = 1

#-----------------------------------------
#--------parsing functions---------------
#-----------------------------------------


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
info = \
    """Count occurences of all combinations
       of filters in the filter column."""

def get_filter_stats_setup_fun(arg_dic):
    try:
        arg_dic['out_fh'] = open(arg_dic['out_fn'],'w')
    except KeyError:
        logging.warning('No out filepath specified. Returning output directly.')
    arg_dic['count_dic'] = {}

def get_filter_stats_parse_fun(line,arg_dic):
    filters = line[6].split(';')
    filters.sort()
    filters = tuple(filters)
    add_to_countdic(arg_dic['count_dic'],'n_sites')
    add_to_countdic(arg_dic['count_dic'],filters)

def get_filter_stats_cleanup_fun(arg_dic):
    import pandas as pd
    filter_info = pd.Series(arg_dic['count_dic'].values(),
                                index=arg_dic['count_dic'].keys())
    filter_info.sort(inplace=True,ascending=False)
    try:
        filter_info.to_csv(arg_dic['out_fh'],sep='\t')
    except KeyError:
        return filter_info
    #json.dump(countdic,open(out_fn,'w'))

add_analysis('get_filter_stats',
             {'setup_fun':get_filter_stats_setup_fun,
              'parse_fun':get_filter_stats_parse_fun,
               'cleanup_fun':get_filter_stats_cleanup_fun},
                info,
                    command_line_req_params={'out_fn':"Filename to write to."})

#--------add_ancestral_from_fasta-----------
info = \
    """Add ancestral state from fasta to vcf.
        E.g.,
       Macaque state is taken from vervet-macaque
       alignment and added in the info column with
       tag AA.
    """

def add_ancestral_from_fasta_setup_fun(arg_dic):
    from pyfasta import Fasta
    arg_dic["out_vcf_fh"] = open(eu(arg_dic["out_vcf"]),'w')
    arg_dic["fasta"] = Fasta(eu(arg_dic["ancestral_fasta"]))
    arg_dic['info_parsed'] = False

def add_ancestral_from_fasta_header_fun(line,arg_dic):
    if line[:7] == '##INFO=':
         arg_dic['info_parsed'] = True
    elif arg_dic['info_parsed']:
        arg_dic["out_vcf_fh"].write(\
                                '##INFO=<ID=AA,Number=1,Type=Character'
                                ',Description="Ancestral allele as'
                                ' derived from {}">\n'.format(arg_dic['ancestral_source']))
        arg_dic['info_parsed'] = False
    arg_dic["out_vcf_fh"].write(line)

def add_ancestral_from_fasta_parse_fun(line,arg_dic):
    aa = arg_dic['fasta'][line[0]][int(line[1])-1]
    if aa not in ['N','n']:
        line[7] = line[7] + ';AA=' + aa
    arg_dic["out_vcf_fh"].write("\t".join(line)+'\n')


add_analysis('add_ancestral_fasta',
             {'setup_fun':add_ancestral_from_fasta_setup_fun,
              'header_fun':add_ancestral_from_fasta_header_fun,
              'parse_fun':add_ancestral_from_fasta_parse_fun},
                info,
                always_req_params={'ancestral_source':\
                                        "Name of source of ancestral allele info.",
                                   'ancestral_fasta':\
                                        'Filepath of fasta with ancestral state.',
                                   'out_vcf':"Filepath to write output vcf to."})

#------------------------

#def filter_missing_stats_parse_fun(line,parser_out,):
#    
#    snp = ["A","G","C","T"]
#    novar = [None,"X",'.']
#    sites_dic = {"total":0,"pass_nosnp":0,"pass_snp":0,"filter_nosnp":0,"filter_snp":0}
#    N_df = pd.DataFrame(0,columns=sites_dic.keys(),index=reader.samples)
#    Nxy = np.zeros((len(reader.samples),len(reader.samples)))
#
#    parser_out['sites_dic']["total"] += 1
#    ns = np.array([1 if s.data.GT is None else 0 for s in rec.samples]) #vector of Ns
#    parser_out['N_df']["total"] +=  ns
#
#        Nxy += np.outer(ns,ns)
#        try:
#            alt = rec.ALT[0].sequence
#        except AttributeError:
#            alt = rec.ALT[0].
#        is_variant = not alt in novar
#        if is_variant:
#            if not rec.FILTER:
#                category = "pass_snp"
#
#            else:
#                category = "filter_snp"
#        else:
#            if not rec.FILTER or (len(rec.FILTER)==1 and ("LowQual" in rec.FILTER) and rec.QUAL>5):
#                category = "pass_nosnp"
#            else:
#                category = "filter_nosnp"
#        sites_dic[category] += 1
#        N_df[category] += ns
#        prev_pos = rec.POS
#....
#    Nxy = pd.DataFrame(Nxy,index=reader.samples,columns=reader.samples)
#    return (sites_dic, N_df, Nxy)


#---MAIN----FUNCTION-------------



if __name__ == "__main__":
    #import gzip
    import argparse
    import select
    parser = argparse.ArgumentParser(description="Parse a Variant Call Format (VCF) file.")
    parser.add_argument("--variant",'-V',type = argparse.FileType('r'), default = '-', help="Input vcf filepath.")
    parser.add_argument("--analysis_type","-T",help="Name of type of analysis, "
                                                      "that defines the functions to use. "
                                                        "Run --show_analyses to see available tools.")
    
    parser.add_argument("--show_analyses",action='store_true',help="List available analyses and exit.")
    parser.add_argument("--analysis_info",help="Get info for specified analysis and exit.")
    parser.add_argument("--arg_dic",help='key=value pairs to pass as dict to functions. E.g., "foo=bar;test=hallo du"')
    parser.add_argument("--no_skip_multiple_entries",action='store_true',
                         help='Do not skip all but the first entry for the same site in the VCF.')


    args = parser.parse_args()

    if args.show_analyses:
        for k in analyses:
            print k+':', analyses['k']['info']
    elif args.analysis_info is not None:
        print analyses[args.analysis_info]
    elif args.analysis_type is not None:
        assert args.analysis_type in analyses, "Analysis {} does not exist."\
                                                "Possible analyses are {}.".format(args.analysis_type,
                                                                                   analyses.keys())
        assert select.select([args.variant,],[],[],0.0)[0], "Input vcf has no data."
        if args.arg_dic is not None:
            arg_dic = {k:v for (k,v) in [t.split('=') for t in args.arg_dic.split(';')]}
        else:
            arg_dic = {}
        check_params(arg_dic,analyses[args.analysis_type]['command_line_req_params'])
        parser = get_parser(args.variant,args.analysis_type,arg_dic,not args.no_skip_multiple_entries)
        parser.run()

    else:
        logging.warning("No analysis specified. Run with flag -h for options.")




    
