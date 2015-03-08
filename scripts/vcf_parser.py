#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.
Todo:
-- parallel
"""
import sys, os, json
import logging
import numpy as np
import pandas as pd
logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s  %(message)s')
logger.setLevel(logging.DEBUG)
eu = os.path.expanduser

nucleotides = ['A','C','T','G']
no_var = ['N','.','X']

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
        self.parse_fun = parse_fun
        self.header_fun = header_fun
        self.setup_fun = setup_fun
        self.cleanup_fun = cleanup_fun
        if arg_dic is None:
            arg_dic = {}
        self.arg_dic = arg_dic
        if 'in_vcf_fh' not in self.arg_dic:
            self.arg_dic['in_vcf_fh'] = vcf_fh
        else:
            logging.warning("Key 'in_vcf_fh' of arg_dic already in use. "
                                                     " Won't overwrite.")
        self.skip_multiple_entries = skip_multiple_entries

    def setup(self):
        logging.info("Starting setup.")
        self.setup_fun(self.arg_dic)
        #self.arg_dic = arg_dic


    def parse_header(self):
        logging.info("Parsing header.")
        for line in self.vcf_fh:
            if line[0] == '#':
                if self.header_fun is not None:
                    self.header_fun(line,self.arg_dic)
            else:
                break

    def parse(self):
        """
        """
        logging.info("Parsing vcf body.")
        prev_chrom = None
        prev_pos = -1
        multiple_count = 0
        for line in self.vcf_fh:
            if line[0] == "#":
                logging.warning("Skipping comment line in vcf: {}".format(line))
                continue
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
        if self.setup_fun is not None:
            self.setup()
        self.parse_header()
        if self.parse_fun is not None:
            self.parse()
        if self.cleanup_fun is not None:
            self.cleanup()
        else:
            self.result = None
        logging.info("Run finished.")



#------support functions----------------

#=====support setup======

def try_add_out_fh(arg_dic,key):
    try:
        arg_dic[key+'_fh'] = open(eu(arg_dic[key]),'w')
    except KeyError:
        pass

#====support parsing=====

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

def add_to_countdic(dic,key):
    try:
        dic[key] += 1
    except KeyError:
        dic[key] = 1

def get_info_dic(line):
    return {k:v for (k,v) in [t.split('=') for t in line[7].split(';')]}

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

#-------add_filter_info---------

info = ("Adds info header entries for each pair "
        "in expressions/descriptions iterables to "
        "in_vcf and writes it to out_vcf ")

def add_filter_info_setup_fun(arg_dic):
    arg_dic["out_vcf_fh"] = open(eu(arg_dic["out_vcf"]),'w')
    arg_dic["line_counter"] = 0
    arg_dic["filter_found"] = False
    arg_dic["filters_added"] = False

def add_filter_info_header_fun(line,arg_dic):
    arg_dic['line_counter'] += 1
    b = len("##FILTER=<ID=")
    if line[:b] == "##FILTER=<ID=":
        logging.debug("Exisitng filter found: {}".format(line))
        arg_dic["filter_found"] = True
        for j,e in enumerate(arg_dic['expressions']):
            if line[b:b+len(e)] == e:
                line = line[:b+len(e)+len(',Description="')] + arg_dic['descriptions'][j] + '">\n'
                logging.info("Updating filter expression {}.".format(e))
                del arg_dic['expressions'][j]
                del arg_dic['descriptions'][j]
                break
    elif not arg_dic["filters_added"] and (arg_dic["filter_found"] or (line[:6] == '#CHROM')):
        for e,d in zip(arg_dic['expressions'],arg_dic['descriptions']):
            arg_dic['out_vcf_fh'].write('##FILTER=<ID=' + e + ',Description="' + d + '">\n')
            logging.info("Adding filter expression {}.".format(e))
        arg_dic["filters_added"] = True
    arg_dic['out_vcf_fh'].write(line)


def add_filter_info_cleanup_fun(arg_dic):
    import subprocess
    arg_dic['out_vcf_fh'].flush()
    logging.info("Starting to cat vcf body.")
    command = ["tail","-n +" + str(arg_dic['line_counter']+1),arg_dic['in_vcf_fh'].name]
    p = subprocess.Popen(" ".join(command),shell=True,stdout=arg_dic['out_vcf_fh'])
    #p.wait()
    p.communicate()

add_analysis('add_filter_info',
             {'setup_fun':add_filter_info_setup_fun,
              'header_fun':add_filter_info_header_fun,
              'cleanup_fun':add_filter_info_cleanup_fun},
                info,
                always_req_params={'expressions':\
                                        "List of filter expressions to add to the vcf header.",
                                   'descriptions':\
                                        'List of filter descriptions to add to the vcf header.',
                                   'out_vcf':"Filepath to write output vcf to."})


#----------------vcf_stats-----------------

info = "Print some statistics about the variants in the vcf."

def vcf_stats_setup_fun(arg_dic):
    try_add_out_fh(arg_dic)
    arg_dic['var_stats'] = {}
    arg_dic['filters'] = {}

def vcf_stats_parse_fun(line,arg_dic):
    add = lambda s: add_to_countdic(arg_dic['var_stats'],s)
    addf = lambda s: add_to_countdic(arg_dic['filters'],s)
    add('total')
    info = get_info_dic(line)
    ref = line[3]
    alt = line[4].split(',')
    pass0 = (line[6] in ['PASS','Pass'])
    if len(alt) > 1:
        add('multiallelic')
    elif len(ref) > 1 or len(alt[0]) > 1:
        add('indels')
    elif alt[0] in nucleotides:
        add('snps')
    if pass0:
        add('pass')
    else:
        for f in line[6].split(';'):
            addf(f)

    try:
        aa = info['AA']
    except KeyError:
        return
    if aa in nucleotides:
        add('ancestral_known')
        if pass0:
            add('pass_ancestral_known')
            if aa == ref:
                add('pass_ancestral_is_ref')
            elif aa == alt[0]:
                add('pass_ancestral_is_alt')
            else:
                add('pass_ancestral_is_third_allele')


def vcf_stats_cleanup_fun(arg_dic):
    arg_dic['var_stats']['filters'] = arg_dic['filters']
    try:
        json.dump(arg_dic['var_stats'],arg_dic['out_fh']) 
    except KeyError:
        return arg_dic['var_stats']

add_analysis('vcf_stats',
             {'setup_fun':vcf_stats_setup_fun,
              'parse_fun':vcf_stats_parse_fun,
              'cleanup_fun':vcf_stats_cleanup_fun},
                info,
                command_line_req_params={'out_fn':\
                                        "File path to write output json to."})

#-------filter_missing_stats---------

info = ("Parse a whole genome (all sites) VCF "
        "to get statistics about the number of called "
        "SNPs and nonsnps. This gives information on the number "
        "of accessible (i.e., non Filtered, non missing genotype) "
        " sites both global and per individual.")

out_fns = ["out_filter_count","out_N_count","out_N_corr"]


def filter_missing_stats_setup_fun(arg_dic):
    for fn in out_fns:
        try_add_out_fh(arg_dic,fn)
    arg_dic['sites_dic'] = {"total":0,"pass_nosnp":0,
                            "pass_snp":0,"filter_nosnp":0,"filter_snp":0}

def filter_missing_stats_header_fun(line,arg_dic):
    if line[:6] == '#CHROM':
        arg_dic['samples'] = line.strip().split('\t')[9:]
        arg_dic['N_df'] = pd.DataFrame(0,
                                        columns=arg_dic['sites_dic'].keys(),
                                        index=arg_dic['samples'])
        arg_dic['Nxy'] = np.zeros((len(arg_dic['samples']),len(arg_dic['samples'])))

def filter_missing_stats_parse_fun(line,arg_dic):
    add = lambda s: add_to_countdic(arg_dic['sites_dic'],s)

    #sites_dic = {"total":0,"pass_nosnp":0,"pass_snp":0,"filter_nosnp":0,"filter_snp":0}

    add("total")
    ns = np.array([1 if '.' in s.split(':')[0] else 0 for s in line[9:]]) #vector of Ns
    arg_dic['N_df']["total"] +=  ns
    arg_dic['Nxy'] += np.outer(ns,ns)

    ref = line[3]
    alt = line[4].split(',')
    pass0 = (line[6] in ['PASS','Pass'])

    if pass0:
        category = 'pass_'
    else:
        category = 'filter_'

    if len(alt) == 1 and len(ref) == 1 and len(alt[0]) == 1 and alt[0] in nucleotides:
        category += 'snp'
    else:
        category += 'nosnp'

    add(category)
    arg_dic['N_df'][category] += ns

def filter_missing_stats_cleanup_fun(arg_dic):
    N_df = arg_dic['N_df']
    sites_dic = arg_dic['sites_dic']
    Nxy = pd.DataFrame(arg_dic['Nxy'],index=arg_dic['samples'],columns=arg_dic['samples'])
    corr = (Nxy-1./sites_dic["total"]*np.outer(N_df["total"],N_df["total"]))/\
            np.sqrt(np.outer(N_df["total"]*(1-1./sites_dic["total"]*N_df["total"]),N_df["total"]*(1-1./sites_dic["total"]*N_df["total"])))
    try:
        json.dump(sites_dic,arg_dic['out_filter_count_fh'])
        N_df.to_csv(arg_dic['out_N_count'],sep='\t')
        corr.to_csv(arg_dic['out_N_corr'],sep='\t')
    except KeyError:
        logging.info("At least one output filename not found. Returning results.")
        return (sites_dic, N_df, corr)

add_analysis('filter_missing_stats',
             {'setup_fun':filter_missing_stats_setup_fun,
              'header_fun':filter_missing_stats_header_fun,
              'parse_fun':filter_missing_stats_parse_fun,
              'cleanup_fun':filter_missing_stats_cleanup_fun},
                info,
                command_line_req_params={'out_filter_count':\
                                        "File path to write count of filtered sites as json to.",
                                         'out_N_count':\
                                        "Tsv path to write per individual missing genotype count to.",
                                        'out_N_corr':\
                                        "Tsv path to write cross individual correlations "
                                        "of missing genotypes to."})

if __name__ == "__main__":
    #import gzip
    import argparse
    import select
    parser = argparse.ArgumentParser(description="Parse a Variant Call Format (VCF) file.")
    parser.add_argument("--variant",'-V',type = argparse.FileType('r'), default = '-', help="Input vcf filepath.")
    parser.add_argument("--analysis_type","-T", choices = analyses.keys(),
                                                 help="Name of type of analysis, "
                                                      "that defines the functions to use. "
                                                      "Run --show_analyses to see available tools.")
    
    parser.add_argument("--show_analyses",action='store_true',help="List available analyses and exit.")
    parser.add_argument("--analysis_info",help="Get info for specified analysis and exit.")
    #parser.add_argument("--arg_dic",help='key=value pairs to pass as dict to functions. E.g., "foo=bar;test=hallo du"')
    parser.add_argument("--no_skip_multiple_entries",action='store_true',
                         help='Do not skip all but the first entry for the same site in the VCF.')
    parser.add_argument('--logging_level','-l',
                        choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],default='INFO',
                        help='Minimun level of logging.')

    args, additional_args = parser.parse_known_args()


    logger.setLevel(getattr(logging,args.logging_level))

    if args.show_analyses:
        for k in analyses:
            print k+':', analyses[k]['info']
    elif args.analysis_info is not None:
        print analyses[args.analysis_info]
    elif args.analysis_type is not None:
        assert args.analysis_type in analyses, "Analysis {} does not exist."\
                                                "Possible analyses are {}.".format(args.analysis_type,
                                                                                   analyses.keys())
        assert select.select([args.variant,],[],[],0.0)[0], "Input vcf has no data."
        ana = analyses[args.analysis_type]
        if additional_args:
            assert additional_args[0][:2] == '--', "First additional argument does not "\
                                                   "start with flag '--': {}".format(additional_args[0])
            arg_dic = {}
            key = None
            for arg in additional_args:
                if arg[:2] == '--':
                    if key is not None:
                        if len(arg_dic[key]) == 1:
                            arg_dic[key] = arg_dic[key][0]
                        elif len(arg_dic[key]) == 0:
                            arg_dic[key] = True
                    key = arg[2:]
                    arg_dic[key] = []
                else:
                    arg_dic[key].append(arg)
            if key is not None:
                if len(arg_dic[key]) == 1:
                    arg_dic[key] = arg_dic[key][0]
                elif len(arg_dic[key]) == 0:
                    arg_dic[key] = True
            #additional_args = "&|?@".join(additional_args).split('--')
            #print additional_args
            #print [a.split(' ',1) for a in additional_args]
            #arg_dic = {k:v.strip() for (k,v) in [a.split('&|?@',1) for a in additional_args if a]}
        else:
            arg_dic = {}
        non_recognized_args = {}
        for arg in arg_dic:
            if arg not in ana['always_req_params'] and \
                arg not in ana['command_line_req_params'] and \
                arg not in ana['opt_params']:
                non_recognized_args.update({arg:arg_dic[arg]})
        logging.info("Standard arguments recognized: {}".format(vars(args)))
        logging.info("Additional arguments: {}".format(arg_dic))
        if non_recognized_args:
            logging.warning("Arguments that were not recognized "
                            "for this analysis_type: {}".format(non_recognized_args))
        check_params(arg_dic,ana['command_line_req_params'])
        
        parser = get_parser(args.variant,args.analysis_type,arg_dic,not args.no_skip_multiple_entries)
        parser.run()

    else:
        logging.warning("No analysis specified. Run with flag -h for options.")




