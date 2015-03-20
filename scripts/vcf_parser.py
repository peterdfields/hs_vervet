#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.
Todo:
-- call parser classer walkers
-- package parse functions in parser classes
-- handle multiple input within argdic 
-- the current line should be stored in arg_dic
-- and the next() method should be called withn the parse function
 (or _yield generator for each addidional input in the future parser class)
Make class to package parse funs:
arg_dic --> self
"""
import sys, os, json, uuid
import logging
import numpy as np
import pandas as pd
import subprocess
import multiprocessing as mp
logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
logger.setLevel(logging.DEBUG)
eu = os.path.expanduser
jn = os.path.join
try:
    import tabix
except ImportError:
    logging.warning('pytabix could not be imported. No support for interval use.')

# assure that all errors get logged
def excepthook(*args):
  logging.getLogger().error('Uncaught exception:', exc_info=args)
sys.excepthook = excepthook



nucleotides = ['A','C','T','G']
no_var = ['N','.','X']





#---------main object-------------------


class Parser(object):
    """
    Generic parser for text file.
    Input:
    in_fh ... file handle of input text file
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
    def __init__(self,in_fh, sep='\t',parse_fun=None, header_fun=None,
                               setup_fun=None, cleanup_fun=None,
                                output_fun = None,
                                arg_dic=None,id='',
                                line_write_vars = None,
                                    skip_multiple_entries=True,
                                    progress_report_interval=50000,**kwa):
        self.in_fh = in_fh
        self.sep = sep
        self.parse_fun = parse_fun
        self.header_fun = header_fun
        self.setup_fun = setup_fun
        self.cleanup_fun = cleanup_fun
        self.output_fun = output_fun
        if line_write_vars is None:
            line_write_vars = []
        self.line_write_vars = line_write_vars
        self.id = id
        if arg_dic is None:
            arg_dic = {}
        self.arg_dic = arg_dic
        #if 'in_fh' in self.arg_dic:
        #    logging.warning("Key 'in_fh' of arg_dic already in use. "
        #                                             "Overwrite!")
        self.arg_dic['in_fh'] = in_fh
        self.skip_multiple_entries = skip_multiple_entries
        self.progress_report_interval = progress_report_interval
        self.finished = False
        self.i = 0
        self.prev_chrom = None
        self.prev_pos = -1
        self.multiple_count = 0
        self.print_warning = True

 
    def _split_line(self,line):
        self.i += 1
        return line.strip().split(self.sep)

    def _yield_split_line(self,fh):
        while True:
            line = self._split_line(fh.next())
            self._report_progress()
            while self._skip_duplicate_line(line) or self._skip_comment(line):
                line = self._split_line(fh.next())
                self._report_progress()
            yield line



    def _header_line_parser(self,line):
        if self.header_fun is not None:
            self.header_fun(line,self.arg_dic)


    def _skip_duplicate_line(self,line):
        chrom = line[0]
        pos = int(line[1])
        if chrom == self.prev_chrom:
            assert pos >= self.prev_pos, "vcf positions not in "\
                                            "ascending order at: {}:{},{}".format(chrom,prev_pos,pos)
            if pos == self.prev_pos:
                self.multiple_count += 1
                if self.multiple_count > 100 and print_warning:
                    logging.warning("Omitting further multiple entry warnings.")
                    self.print_warning = False
                elif self.print_warning:
                    if not self.skip_multiple_entries:
                        logging.warning("Multiple entries for pos {}:{}.\n"
                                  "Keeping all entries.".format(chrom,pos))
                        return False
                    else:
                        logging.warning("Warning, multiple entries for pos {}:{}.\n"
                              "Skipping all but the first.".format(chrom,pos))
                        return True
        self.prev_chrom = chrom
        self.prev_pos = pos
        return False

    def _report_progress(self):
        if self.i % self.progress_report_interval == 0:
            logging.info("{} Parsed {} lines: {} - {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def _skip_comment(self,line):
        if line[0][0] == "#":
            logging.warning("Skipping comment line in vcf: {}".format(line))
            return True
        else:
            return False

    def setup(self):
        if self.setup_fun is not None:
            logging.info("Starting setup.")
            self.setup_fun(self.arg_dic)


    def parse_header(self):
        logging.info("Parsing header.")
        for line in self.in_fh:
            if line[0] == '#':
                self._header_line_parser(line)
            else:
                break
        for v in self.line_write_vars:
            self.arg_dic[v+'_fh'].flush()

    def parse(self,fh):
        """
        """
        if self.parse_fun is not None:
            logging.info("{} Parsing vcf body of {}.".format(self.id,fh))
            line_it = self._yield_split_line(fh)
            for d in line_it:
                self.parse_fun(d,self.arg_dic)
            #not perfect implementation, prev_pos is not necessarily updated 
            #in children if _skip_duplicate_line is overidden
            logging.info("{} Finished: {} lines at {} {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def cleanup(self):
        """
        """
        if self.cleanup_fun is not None:
            logging.info("Starting cleanup.")
            self.cleanup_fun(self.arg_dic)
        for v in self.line_write_vars:
            self.arg_dic[v+'_fh'].close()

    def output(self):
        if self.output_fun is not None:
            logging.info("Creating output.")
            self.result = self.output_fun(self.arg_dic)
        else:
            self.result = None
        self.finished = True

    def run(self):
        self.setup()
        self.parse_header()
        self.parse(self.in_fh)
        self.cleanup()
        self.output()
        logging.info("Run finished.")


class SerialParser(Parser):
    """
    Parse several regions serially.
    Only the parse function is applied to each
    (tabix) file handle.
    Input:
    vcf_fh ... file handle of vcf 
               (must be tabixed and opened with gvcf)
    intervals ... list of string intervals, in Samtools format,
                                            such as Chr1:1-5000
    kwa ... same keyword arguments as VCFParser.
    
    Note that SerialParser does not support seperator specification,
    since pytabix does the line splitting automatically.
    """
    def __init__(self,in_fh,intervals,auto_tabix=False, **kwa):
        super(SerialParser, self).__init__(in_fh,**kwa)
        self.intervals = intervals
        #if 'line_write_vars' not in kwa or kwa['line_write_vars'] is None:
        #    kwa['line_write_vars'] = []
        #self.line_wire_vars = kwa['line_write_vars']
        self.tabix_fh = tabix.open(in_fh.name)
        try:
            self._query_tabix(intervals[0])
        except tabix.TabixError, e:
            logging.error("Tabix raised error: {}".format(str(e)))
            logging.warning("Compress file with bgzip and produce tabix index.")
            raise e
            #if auto_tabix:
            #    logging.warning("Tabix raised error: {}".format(str(e)))
            #    p=subrocess.Popen([tabix])


    def _query_tabix(self,interval):
        try:
            return self.tabix_fh.querys(interval)
        except TypeError:
            try:
                return self.tabix_fh.query(*interval)
            except TypeError, e:
                logging.error("Chromosome must be string and position integer."
                                               " e.g. ('Chr1',1000000,1005000). "
                                "Alternatively use string 'Chr:1000000-1005000'")
                raise e
#    def _split_line(self,fh):
#        """
#        The line in the tabix iterator is already split.
#        """
#        while True:
#            yield fh.next()


    def _split_line(self,line):
        return line



    def run_no_output(self):
        self.setup()
        self.parse_header()
        if self.parse_fun is not None:
            for interval in self.intervals:
                fh = self._query_tabix(interval)
                self.parse(fh)
        else:
            logging.warning("No vcf body parse function supplied.")
        self.cleanup()

    def run(self):
        self.run_no_output()
        self.output()
        logging.info("Run finished.")

#class MultiInputSerialParser(SerialParser):
#    """
#    Parse multiple vcf_files with corresponding lines.
#    """
#    import itertools
#    def __init__(self,in_fhs,intervals, **kwa):
#        #this is a bit sloppy, the further in_fhs are not tested whether
#        #they can be queried with tabix
#        super(MultiInputSerialParser, self).__init__(in_fhs[0],intervals,**kwa)
#        self.in_fh = in_fhs
#
#
#    def _yield_split_line(self,fhs):
#        """
#        The line in the tabix iterator is already split.
#        """
#        nxt = lambda fh: super(MultiInputSerialParser, self)._yield_split_line(fh)
#        while True:
#            lines = [nxt(fh) for fh in fhs]
#            #lines = [fhi.next() for fhi in fh]
#            poss = [int(l[1]) for l in lines]
#            assert len(set([l[0] for l in lines])) <= 1, ("Unequal chromosomes in "
#                                                "MultiInSerialParser not implemented.")
#            while max(poss)>min(poss):
#                lines = [lines[i]  if poss[i]==max(poss) else nxt(fhs[i]) for i in range(len(lines))]
#                poss = [int(l[1]) for l in lines]
#            yield tuple(lines)

class MultiRegionParallelParser(SerialParser):
    def __init__(self, in_fh, intervals,
                            ncpus='auto', tmp_dir='.',**kwa):
        super(MultiRegionParallelParser, self).__init__(in_fh,intervals,**kwa)
        self.reduce_fun = kwa['reduce_fun']
        arg_dic = kwa['arg_dic']
        del kwa['arg_dic']
        kwa['header_fun'] = None
        self.parsers_kwa = kwa
        if ncpus == 'auto':
            ncpus = mp.cpu_count()
        self.ncpus = ncpus
        self.tmp_dir = tmp_dir

    def set_ncpus(self):
        self.ncpus = min(len(self.intervals),self.ncpus)
        logging.info("{} regions specified.  Parallelising by region. "
                     "Using {} processes.".format(len(self.intervals),self.ncpus))

    def setup_parsers(self):
        parsers = []
        temp_files = []
        for i,interval in enumerate(self.intervals):
            p_arg_dic = self.arg_dic.copy()
            for var in self.line_write_vars:
                tmp_fn = jn(self.tmp_dir,var+str(uuid.uuid4()) + \
                                        "_" + str(i) + ".tmp")
                p_arg_dic[var] = tmp_fn
                temp_files.append(tmp_fn)
            parsers.append(SerialParser(self.in_fh, [interval],
                                        arg_dic=p_arg_dic, id="Interval {}:".format(interval), **self.parsers_kwa))
        self.parsers = parsers
        self.temp_files = temp_files

    def reduce(self):
        if self.reduce_fun is not None:
            logging.info("Starting reduce step.")
            self.reduce_fun(self.arg_dic,[ad for ad in self.out_arg_dics])
        else:
            logging.warning("No reduce function given. Won't do anything.")

    def del_temp_files(self):
        logging.info("Removing temp files.")
        while self.temp_files:
            os.remove(self.temp_files.pop())

    def run_parser(self,i):
        self.parsers[i].run_no_output()
        return self.parsers[i].arg_dic

    def run(self):
        self.setup()
        self.parse_header()
        self.set_ncpus()
        self.setup_parsers()
        out_arg_dics = parmap(self.run_parser,range(len(self.parsers)))
        self.out_arg_dics = out_arg_dics
        logging.info("Reducing.")
        self.reduce()
        logging.info("Creating output.")
        self.output()
        self.del_temp_files()
        logging.info("Run finished.")



class SingleRegionParallelParser(MultiRegionParallelParser):
    def __init__(self, in_fh, intervals, **kwa):
        assert len(intervals) == 1, ("ParallelSingleRegionParser requires a "
                                                           "single interval.")
        #logging.info("One interval specified, starting single_region_parallel mode with {} cores.".format(ncpus))
        super(SingleRegionParallelParser, self).__init__(in_fh, intervals, **kwa)
        self.in_fh = in_fh
        region = intervals[0]
        try:
            #self.tabix_fh.querys(region)
            chrompos = region.split(":")
            chrom = chrompos[0]
            try:
                startend = chrompos[1].split('-')
                start = int(startend[0])
                try:
                    end = int(startend[1])
                except:
                    end = None
            except IndexError:
                start = None
                end = None
        except TypeError:
            chrom = region[0]
            start = region[1]
            end = region[2]
        self.chrom = chrom
        self.start = start
        self.end = end
        if self.start is None:
            self.start = 0
        if self.end is None:
            self._header_line_parser = self._header_line_parser_search_contig_len
            logging.info("No chromosome end specified, searching for 'contig'"
                         "tag in vcf header.")
        self.kwa = kwa


    def parse_header(self):
        super(SingleRegionParallelParser, self).parse_header()
        self.intervals = self.get_intervals()

#    def parse_header_contig_len(self):
#        logging.info("Parsing header.")
#        for line in self.in_fh:
#            if line[0] == '#':
#                if line[:9] == '##contig=':
#                    contig_dic = get_header_line_dic(line)
#                    if contig_dic['ID'] == self.chrom:
#                        length = int(contig_dic['length'])
#                        self.end = length
#                if self.header_fun is not None:
#                    self.header_fun(line,self.arg_dic)
#            else:
#                break

        self.intervals = self.get_intervals()

    def _header_line_parser_search_contig_len(self,line):
        if line[:9] == '##contig=':
                    contig_dic = get_header_line_dic(line)
                    if contig_dic['ID'] == self.chrom:
                        length = int(contig_dic['length'])
                        self.end = length
        if self.header_fun is not None:
            self.header_fun(line,self.arg_dic)

    def get_intervals(self):
        n_chunks = self.ncpus
        chunksize = int((self.end-self.start)/n_chunks) + 1
        starts = []
        for s in range(self.start,self.end,chunksize):
            starts.append(s)
        intervals = [(self.chrom,s,e) for s,e in zip(starts,[s-1 for s in starts[1:]]+[self.end])]
        return intervals



#--------------------------------------------------------

analyses = {}
def add_analysis(name,funs,info='',always_req_params=None,
                    command_line_req_params=None,
                    opt_params=None,line_write_vars=None):
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
               'opt_params': opt_params if opt_params is not None else {},
               'line_write_vars':line_write_vars}
    analyses.update({name:entry})



def check_params(arg_dic,req_param_dic):
    """
    Check whether params noted in req_param_dic are 
    present in arg dic
    """
    for k in req_param_dic:
        if k not in arg_dic:
            raise TypeError('Required parameter --{} <{}> '
                            'missing from arg_dic.'.format(k,k))

def get_parser(vcf_fh,analysis,**kwa):
    check_params(kwa['arg_dic'],analysis['always_req_params'])
    kwa.update(analysis['funs'])
    try:
        kwa['line_write_vars'] = analysis['line_write_vars']
    except KeyError:
        logging.warning('Analysis has no key "line_write_vars". Assuming '
                        'that there is no output written line by line.')
    if not 'intervals' in kwa or not kwa['intervals']:
        #if 'intervals' in kwa:
        #    del kwa['intervals']
        if 'ncpus' in kwa and kwa['ncpus'] > 1:
            logging.warning("For ncpus>1, specify at least one interval,  "
                            "e.g., -L Chr1. Falling back to single process...")
        #    del kwa['ncpus']
        logging.info("Initialising Parser.")
        parser = Parser(vcf_fh,**kwa)
    elif 'ncpus' in kwa and kwa['ncpus'] > 1:
        assert 'reduce_fun' in analysis['funs'],("This analysis does not support "
                                                 "parallel execution, remove "
                                                 "--ncpus option.")
        if len(kwa['intervals']) == 1:
            logging.info("Initialising SingleRegionParallelParser.")
            parser = SingleRegionParallelParser(vcf_fh,**kwa)
        else:
            logging.info("Initialising MultiRegionParallelParser.")
            parser = MultiRegionParallelParser(vcf_fh,**kwa)
    else:
        logging.info("Initialising SerialParser.")
        parser = SerialParser(vcf_fh,**kwa)
    return parser

#parallel support


def fun(f,q_in,q_out):
    while True:
        i,x = q_in.get()
        if i is None:
            break
        q_out.put((i,f(x)))

def parmap(f, X, nprocs = mp.cpu_count()):
    q_in   = mp.Queue(1)
    q_out  = mp.Queue()

    proc = [mp.Process(target=fun,args=(f,q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i,x in sorted(res)]


#------support functions----------------

#=====support setup======

def try_add_out_fh(arg_dic,key):
    try:
        arg_dic[key+'_fh'] = open(eu(arg_dic[key]),'w')
    except KeyError:
        logging.warning("Could add filehandle for {}, no file specified. ".format(key))
        #pass

#====support parsing=====

def get_header_line_dic(line):
    dic = {k:v for k,v in [t.split('=') for t  in line.strip().split('=<')[1][:-1].split(',')]}
    return dic


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

def sum_countdics(countdics):
    tot_countdic = {}
    for countdic in countdics:
        for k,v in countdic.iteritems():
            try:
                tot_countdic[k] += v
            except KeyError:
                tot_countdic[k] = v
    return tot_countdic


def get_info_dic(line):
    return {k:v for (k,v) in [t.split('=') for t in line[7].split(';')]}

#-----------------------------------------
#--------parsing functions---------------
#-----------------------------------------

#-----vcf_to_012--------

info = ("Extract genotype information into a tsv file. "
        "Coding 0 for homozygous reference, 1 for heterozygote "
        "and 2 for homozygous alternative allele. "
        "Missing genotypes are coded by 'N'")

def vcf_to_012_setup_fun(arg_dic):
    try_add_out_fh(arg_dic,'out_tsv')

def vcf_to_012_header_fun(line, arg_dic):
    if line[1:6] == "CHROM":
        arg_dic['out_tsv_fh'].write("chrom\tpos\t"+line.split("\t",9)[-1])

def vcf_to_012_parse_fun(d, arg_dic):
    gt = map(get_012,d[9:])
    arg_dic['out_tsv_fh'].write("\t".join(d[:2])+"\t"+"\t".join(gt)+"\n")



def vcf_to_012_reduce_fun(arg_dic,arg_dics):
    command = ["cat"]+[ad['out_tsv'] for ad in arg_dics]
    p = subprocess.Popen(command, stdout=arg_dic['out_tsv_fh'])
    p.communicate()

add_analysis('vcf_to_012',
             {'setup_fun':vcf_to_012_setup_fun,
              'parse_fun':vcf_to_012_parse_fun,
              'header_fun':vcf_to_012_header_fun,
              'reduce_fun':vcf_to_012_reduce_fun},
                info,
                always_req_params={'out_tsv':\
                                        "File path to write output tsv to."},
                line_write_vars=['out_tsv'])



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
    filter_info = pd.Series(arg_dic['count_dic'].values(),
                                index=arg_dic['count_dic'].keys())
    filter_info.sort(inplace=True,ascending=False)
    arg_dic['filter_info'] = filter_info

def get_filter_stats_reduce_fun(arg_dic,arg_dics):
    fi = arg_dics[0]['filter_info']
    for ad in arg_dics[1:]: 
        fi = fi.add(ad['filter_info'], fill_value=0) 
    arg_dics['filter_info'] = fi

def get_filter_stats_output_fun(arg_dic):
    try:
        arg_dic['filter_info'].to_csv(arg_dic['out_fh'],sep='\t')
    except KeyError:
        return filter_info


add_analysis('get_filter_stats',
             {'setup_fun':get_filter_stats_setup_fun,
              'parse_fun':get_filter_stats_parse_fun,
              'reduce_fun':get_filter_stats_parse_fun,
               'cleanup_fun':get_filter_stats_cleanup_fun,
               'output_fun':get_filter_stats_cleanup_fun},
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
    try_add_out_fh(arg_dic,'out_fn')
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
    else:
        logging.warning("Unknown variant type at {} - {}: {}".format(line[0],line[1],line[4]))

    if float(info['AF'])>0 and float(info['AF'])<1:
        if pass0:
            add('pass')
        else:
            for f in line[6].split(';'):
                addf(f)
    else:
        addf('non_segregating')

    try:
        aa = info['AA']
    except KeyError:
        return
    if aa in nucleotides:
        add('ancestral_known')
        if pass0 and float(info['AF'])>0 and float(info['AF'])<1:
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
        json.dump(arg_dic['var_stats'],arg_dic['out_fn_fh']) 
    except KeyError:
        logging.info("No output file supplied returning result within python.")
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
    for fn in ["out_filter_count","out_N_count","out_N_corr"]:
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
    try:
        af = float(get_info_dic(line)['AF'].split(',')[0])
        if len(alt) == 1 and len(ref) == 1 and len(alt[0]) == 1 and alt[0] in nucleotides and af>0 and af<1:
            category += 'snp'
        else:
            category += 'nosnp'
    except KeyError:
        category += 'nosnp'

    add(category)
    arg_dic['N_df'][category] += ns

def filter_missing_stats_reduce_fun(arg_dic,arg_dics):
    arg_dic['N_df'] = sum([d['N_df'] for d in arg_dics])
    arg_dic['sites_dic'] = sum_countdics([d['sites_dic'] for d in arg_dics])
    arg_dic['Nxy'] = sum([d['Nxy'] for d in arg_dics])

def filter_missing_stats_output_fun(arg_dic):
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
              'reduce_fun':filter_missing_stats_reduce_fun,
              'output_fun':filter_missing_stats_output_fun},
                info,
                command_line_req_params=\
                {'out_filter_count':\
                        "File path to write count of filtered sites as json to.",
                 'out_N_count':\
                        "Tsv path to write per individual missing genotype count to.",
                 'out_N_corr':\
                        "Tsv path to write cross individual correlations "
                                                "of missing genotypes to."})


#-----------create_indel_bed--------------
#attention bed has 0-based coordinates while VCF has 1-based coordinates

info = ("Creates a bed file with the intervals "
       "of all indels in the input vcf.")

def create_indel_bed_setup_fun(arg_dic):
    arg_dic["out_bed_fh"] = open(arg_dic["out_bed"],'w')
    if 'extend_interval' not in arg_dic:
        arg_dic['extend_interval'] = 0
    else:
        arg_dic['extend_interval'] = int(arg_dic['extend_interval'])
    arg_dic['contig_dic'] = {}

def create_indel_bed_header_fun(line, arg_dic):
    if line[:9] == '##contig=':
        c_dic = get_header_line_dic(line)
        arg_dic['contig_dic'].update({c_dic['ID']:int(c_dic['length'])})

def create_indel_bed_parse_fun(line,arg_dic):
    ref = line[3]
    alt = line[4].split(',')
    pos = int(line[1])
    chrom = line[0]
    if len(ref)>1 or [a for a in alt if len(a)>1]:
        start = max(0, pos - 1 - arg_dic['extend_interval'])
        end = min(arg_dic['contig_dic'][chrom],start + len(ref) + arg_dic['extend_interval'])
        arg_dic['out_bed_fh'].write("{}\t{}\t{}\n".format(chrom,start,end))


add_analysis("create_indel_bed",
            {'setup_fun': create_indel_bed_setup_fun,
            'header_fun': create_indel_bed_header_fun,
            'parse_fun': create_indel_bed_parse_fun,
            },
            info,
            always_req_params=\
            {'out_bed':'Filepath of the output bed.'},
            opt_params={'extend_interval':'Extend the interval by n bases '
                                           'to left and right. Default 0.'},
            line_write_vars=['out_bed'])

#----------create_adjacent_snp_bed--------------------
#attention bed has 0-based coordinates while VCF has 1-based coordinates

info = ("Creates a bed file with the intervals "
       "of all adjacent or close snps in the input vcf.")

def create_adjacent_snp_bed_setup_fun(arg_dic):
    arg_dic["out_bed_fh"] = open(arg_dic["out_bed"],'w')
    if 'extend_interval' not in arg_dic:
        arg_dic['extend_interval'] = 0
    else:
        arg_dic['extend_interval'] = int(arg_dic['extend_interval'])
    if 'min_dist' not in arg_dic:
        arg_dic['min_dist'] = 2
    else:
        arg_dic['min_dist'] = int(arg_dic['min_dist'])
    arg_dic['contig_dic'] = {}
    arg_dic['last_chrom'] = None
    arg_dic['last_pos'] = -inf
    arg_dic['intv_start'] = False
    arg_dic['intv_end'] = False

def create_adjacent_snp_bed_header_fun(line, arg_dic):
    if line[:9] == '##contig=':
        c_dic = get_header_line_dic(line)
        arg_dic['contig_dic'].update({c_dic['ID']:int(c_dic['length'])})

def create_adjacent_snp_bed_parse_fun(line,arg_dic):
    ref = line[3]
    alt = line[4].split(',')
    pos = int(line[1])
    chrom = line[0]
    
    if alt[0] in nucleotides: #SNP?
        if arg_dic['last_chrom'] != chrom:
            if arg_dic['intv_start']:
                arg_dic['intv_end'] = arg_dic['last_pos']
        elif pos - arg_dic['last_pos'] < arg_dic['min_dist']:
            if not arg_dic['intv_start']:
                arg_dic['intv_start'] = arg_dic['last_pos']
        elif arg_dic['intv_start']:
            arg_dic['intv_end'] = arg_dic['last_pos']

        if arg_dic['intv_start'] and arg_dic['intv_end']:
            start = max(0, arg_dic['intv_start'] - 1 - arg_dic['extend_interval'])
            end = min(arg_dic['contig_dic'][arg_dic['last_chrom']], arg_dic['intv_end'] + arg_dic['extend_interval'])
            arg_dic['out_bed_fh'].write("{}\t{}\t{}\n".format(arg_dic['last_chrom'],start,end))
            arg_dic['intv_start'] = False
            arg_dic['intv_end'] = False
        arg_dic['last_chrom'] = chrom
        arg_dic['last_pos'] = pos

def create_adjacent_snp_bed_cleanup_fun(arg_dic):
    if arg_dic['intv_end']:
        start = max(0, arg_dic['intv_start'] - 1 - arg_dic['extend_interval'])
        end = min(arg_dic['contig_dic'][arg_dic['last_chrom']], arg_dic['intv_end'] + arg_dic['extend_interval'])
        arg_dic['out_bed_fh'].write("{}\t{}\t{}\n".format(arg_dic['last_chrom'],start,end))


add_analysis("create_adjacent_snp_bed",
            {'setup_fun': create_adjacent_snp_bed_setup_fun,
            'header_fun': create_adjacent_snp_bed_header_fun,
            'parse_fun': create_adjacent_snp_bed_parse_fun,
            'cleanup_fun':create_adjacent_snp_bed_cleanup_fun
            },
            info,
            always_req_params=\
            {'out_bed':'Filepath of the output bed.'},
            opt_params={'extend_interval':'Extend the interval by n bases '
                                           'to left and right. Default 0.',
                        'min_dist':'Minimun distance between SNPs not to be masked.'
                                    'Default 2 meaning that only directly adjacent SNPs are included.'},
            line_write_vars=['out_bed'])


#---------filter_by_bed----------

info = ("Add filter tag to sites in intervals of bed file. "
        "This parser is experimental, only working for special cases.")

def filter_by_bed_setup_fun(arg_dic):
    if isinstance(arg_dic['in_beds'], basestring):
        arg_dic['in_beds'] = [arg_dic['in_beds']]
    if isinstance(arg_dic['filter_names'], basestring):
        arg_dic['filter_names'] = [arg_dic['filter_names']]
    assert len(arg_dic['filter_names']==arg_dic['in_beds']), "There must be as many filter names as beds."
    arg_dic['in_beds_fh'] = [open(b) for b in arg_dic['in_beds']] 
    arg_dic['out_vcf_fh'] = open(arg_dic['out_vcf'])
    arg_dic['last_rec'] = [None for b in arg_dic['in_beds']]


def filter_by_bed_header_fun(arg_dic):
    arg_dic['out_vcf_fh'].write(line)


def filter_by_bed_parse_fun(line,arg_dic):
    def get_rec(fh):
        rec = fh.next().strip().split()
        rec[1] = int(rec[1])
        rec[2] = int(rec[2])
    return rec

    ref = line[3]
    alt = line[4].split(',')
    pos = int(line[1])
    chrom = line[0]
    for i in range(len(arg_dic['last_rec'])):
        if arg_dic['last_rec'][i] is None or arg_dic['last_rec'][i][2]<pos:
            arg_dic['last_rec'][i] = get_rec(arg_dic['in_bed_fh'][i])
        if arg_dic['last_rec'][i][0] == chrom \
            and arg_dic['last_rec'][i][1] < pos \
            and arg_dic['last_rec'][i][2] > pos:
            if line[6] in ['.','PASS']:
                line[6] = arg_dic['filter_names'][i]
            else:
                line[6] = line[6] + ','  + arg_dic['filter_names'][i]
    arg_dic['out_vcf_fh'].write("\t".join(line)+'\n')


add_analysis("filter_by_bed",
            {'setup_fun': filter_by_bed_setup_fun,
            'header_fun': filter_by_bed_header_fun,
            'parse_fun': filter_by_bed_parse_fun,
            },
            info,
            always_req_params=\
            {'in_beds':'Filepath of the input beds.',
             'filter_names':'Filter names corresponding to beds.',
             #'overwrite':'Overwrite all other filters '
             #           'and mark sites as PASS if not '
             #          'covered by the filter. Default False.',
             'out_vcf':'Filepath to output vcf.'},
            line_write_vars=['out_vcf'])


#------------samtools_to_gatk-------------


info = ("Samtools to GATK format. !!!Experimental!!!")


def samtools_to_gatk_setup_fun(arg_dic):
    arg_dic['out_vcf_fh'] = open(arg_dic['out_vcf'],'w')

def samtools_to_gatk_parse_fun(line,arg_dic):
    alt = line[4].split(',')
    alt = [a for a in alt if a != 'X']
    if not alt:
        alt = ['.']
    line[4] = ",".join(alt)
    arg_dic['out_vcf_fh'].write("\t".join(line)+'\n')


add_analysis("samtools_to_gatk",
            {'setup_fun':  samtools_to_gatk_setup_fun,
            'parse_fun':  samtools_to_gatk_parse_fun,
            },
            info,
            always_req_params=\
            {'out_vcf':'Filepath to output vcf.'},
            line_write_vars=['out_vcf'])




if __name__ == "__main__":
    import gzip
    import argparse
    import select
    import time

    parser = argparse.ArgumentParser(description="Parse a Variant Call Format (VCF) file.")
    parser.add_argument("--variant",'-V',type = argparse.FileType('r'), default = '-', help="Input vcf filepath.")
    parser.add_argument("--analysis_type","-T", choices = analyses.keys(),
                                                 help="Name of type of analysis, "
                                                      "that defines the functions to use. "
                                                      "Run --show_analyses to see available tools.")
    
    parser.add_argument("--show_analyses",action='store_true',help="List available analyses and exit.")
    parser.add_argument("--analysis_info",help="Get info for specified analysis and exit.")
    parser.add_argument("--intervals",'-L', nargs='*', dest='intervals', action='append',
                            help='Specify intervals to consider e.g. Chr1:1-50000. '
                                 'Input vcf must be compressed with bgzip and indexed '
                                                                          'with tabix.')
    parser.add_argument("--ncpus", '-nct',
                        type=int, default=1,
                                  help='Number of processes for parallel parsing. '
                                       ' Requires at least one interval to be '
                                       'specified '
                                       'with -L. \n'
                                       '1) If a single interval is '
                                       'specified (e.g. -L Chr1:1-5000000), '
                                       'this interval will be split into equal parts. '
                                       'If start/end are not given, we try to infer '
                                       'them from the VCF header (tag contig=...). \n'
                                       '2) If multiple intervals are specified, '
                                       'we parallelise across intervals. \n'
                                       'Parallel parsing is not implemented for all '
                                       'analyses.')
    parser.add_argument("--no_skip_multiple_entries",action='store_true',
                         help='Do not skip all but the first entry for the same site in the VCF.')
    parser.add_argument('--logging_level','-l',
                        choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                        default='INFO',
                                                help='Minimun level of logging.')
    parser.add_argument('--progress_report_interval',
                        type=int, default=200000,
                        help='Number of lines after which to report progress. '
                                            'Only output if logging level >= INFO')
    args, additional_args = parser.parse_known_args()

    #print args.intervals
    #for line in args.variant:
    #    print line
    #    break
    #print args.variant.name
    #sys.exit()

    logger.setLevel(getattr(logging,args.logging_level))

    if args.show_analyses:
        for k in analyses:
            print k+':', analyses[k]['info']
    elif args.analysis_info is not None:
        print analyses[args.analysis_info]
    elif args.analysis_type is not None:
        assert args.analysis_type in analyses, "Analysis {} does not exist."\
                                               "Possible analyses are {}."\
                                               .format(args.analysis_type,
                                                            analyses.keys())
        try:
            assert select.select([args.variant,],[],[],0.0)[0], "Input vcf has no data."
        except TypeError: #the above check does not work for tabix
            pass

        args.intervals = [el for elements in args.intervals for el in elements] \
                                                if args.intervals is not None else []
        extension = os.path.splitext(args.variant.name)[-1]
        if extension == '.gz':
            args.variant = gzip.open(args.variant.name)
        else:
            assert not args.intervals,("Interval mode (-L) only supported on bgzipped "
                                 "files with tabix index. But input has not ending .gz")
            if  extension != '.vcf':
                logging.warning("Unrecognized file extension. "
                                "Assuming this is a vcf: {}".format(args.variant.name))


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

        parser = get_parser(args.variant,analyses[args.analysis_type],
                             arg_dic=arg_dic, intervals = args.intervals,
                             ncpus=args.ncpus,
                                skip_multiple_entries=not args.no_skip_multiple_entries,
                                progress_report_interval=args.progress_report_interval)
        start = time.time()
        parser.run()
        end = time.time()
        delta = end - start
        logging.info("This run took {} seconds = {} minutes = {} hours.".format(delta,delta/60.,delta/3600.))

    else:
        logging.warning("No analysis specified. Run with flag -h for options.")




