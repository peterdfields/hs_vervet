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
import logging, argparse
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


class Walker(object):
    """
    Generic walker through delimited text file.
    Applies methods of a parser object to each
    line in the file.

    Input:
    in_fh ... file handle of input text file
    parser... object instance that derives 
              from Parser class
    Examples:

    """
    def __init__(self,in_fh, parser, sep='\t',
                               id='',
                                    skip_multiple_entries=True,
                                    progress_report_interval=50000,**kwa):
        self.in_fh = in_fh
        self.sep = sep
        self.parser = parser
        self.id = id
        self.skip_multiple_entries = skip_multiple_entries
        self.progress_report_interval = progress_report_interval
        self.finished = False
        self.i = 0
        self.prev_chrom = None
        self.prev_pos = -1
        self.multiple_count = 0
        self.print_warning = True

 
    def _split_line(self,line):
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
        if self.parser.header_fun is not None:
            self.parser.header_fun(line)


    def _skip_duplicate_line(self,line):
        chrom = line[0]
        pos = int(line[1])
        if chrom == self.prev_chrom:
            assert pos >= self.prev_pos, "vcf positions not in "\
                                            "ascending order at: {}:{},{}".format(chrom,prev_pos,pos)
            if pos == self.prev_pos:
                self.multiple_count += 1
                if self.multiple_count > 100 and self.print_warning:
                    logging.warning("Omitting further multiple entry warnings.")
                    self.print_warning = False
                if not self.skip_multiple_entries:
                    if self.print_warning:
                        logging.warning("Multiple entries for pos {}:{}.\n"
                                  "Keeping all entries.".format(chrom,pos))
                    return False
                else:
                    if self.print_warning:
                        logging.warning("Warning, multiple entries for pos {}:{}.\n"
                              "Skipping all but the first.".format(chrom,pos))
                    return True
        self.prev_chrom = chrom
        self.prev_pos = pos
        return False


    def _report_progress(self):
        self.i += 1
        if self.i % self.progress_report_interval == 0:
            logging.info("{} Parsed {} lines: {} - {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def _skip_comment(self,line):
        if line[0][0] == "#":
            logging.warning("Skipping comment line in vcf: {}".format(line))
            return True
        else:
            return False



    def parse_header(self):
        logging.info("Parsing header.")
        for line in self.in_fh:
            if line[0] == '#':
                self._header_line_parser(line)
            else:
                break
        for v in vars(self.parser):
            #to be sure that header is written before body
            #personally I find duck-typing dangerous...
            try:
                getattr(self.parser,v).flush()
            except AttributeError:
                pass

    def parse(self,fh):
        """
        """
        if self.parser.parse_fun is not None:
            logging.info("{} Parsing vcf body of {}.".format(self.id,fh))
            line_it = self._yield_split_line(fh)
            for d in line_it:
                self.parser.parse_fun(d)
            #not perfect implementation, prev_pos is not necessarily updated 
            #in children if _skip_duplicate_line is overidden
            logging.info("{} Finished: {} lines at {} {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def cleanup(self):
        """
        """
        if self.parser.cleanup_fun is not None:
            logging.info("Starting cleanup.")
            self.parser.cleanup_fun()
        for v in vars(self.parser):
            #to be sure that header is written before body
            #personally I find duck-typing dangerous...
            try:
                getattr(self.parser,v).close()
            except AttributeError:
                pass

    def output(self):
        if self.parser.output_fun is not None:
            logging.info("Creating output.")
            self.result = self.parser.output_fun()
        else:
            self.result = None
        self.finished = True

    def run(self):
        self.parse_header()
        self.parse(self.in_fh)
        self.cleanup()
        self.output()
        logging.info("Run finished.")


class SerialWalker(Walker):
    """
    Parse several regions serially.
    Only the parse function of parser 
    is applied to each
    (tabix) file handle.
    Input:
    vcf_fh ... file handle of vcf 
               (must be tabixed and opened with gvcf)
    intervals ... list of string intervals, in Samtools format,
                                            such as Chr1:1-5000
    kwa ... same keyword arguments as VCFWalker.
    
    Note that SerialWalker does not support seperator specification,
    since pytabix does the line splitting automatically.
    """
    def __init__(self,in_fh,parser,intervals,auto_tabix=False, **kwa):
        super(SerialWalker, self).__init__(in_fh,parser,**kwa)
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

    def _split_line(self,line):
        """
        Tabix handles the split automatically.
        """
        return line

    def run_no_output(self):
        self.parse_header()
        if self.parser.parse_fun is not None:
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

#class MultiInputSerialWalker(SerialWalker):
#    """
#    Parse multiple vcf_files with corresponding lines.
#    """
#    import itertools
#    def __init__(self,in_fhs,intervals, **kwa):
#        #this is a bit sloppy, the further in_fhs are not tested whether
#        #they can be queried with tabix
#        super(MultiInputSerialWalker, self).__init__(in_fhs[0],intervals,**kwa)
#        self.in_fh = in_fhs
#
#
#    def _yield_split_line(self,fhs):
#        """
#        The line in the tabix iterator is already split.
#        """
#        nxt = lambda fh: super(MultiInputSerialWalker, self)._yield_split_line(fh)
#        while True:
#            lines = [nxt(fh) for fh in fhs]
#            #lines = [fhi.next() for fhi in fh]
#            poss = [int(l[1]) for l in lines]
#            assert len(set([l[0] for l in lines])) <= 1, ("Unequal chromosomes in "
#                                                "MultiInSerialWalker not implemented.")
#            while max(poss)>min(poss):
#                lines = [lines[i]  if poss[i]==max(poss) else nxt(fhs[i]) for i in range(len(lines))]
#                poss = [int(l[1]) for l in lines]
#            yield tuple(lines)

class MultiRegionParallelWalker(SerialWalker):
    def __init__(self, in_fh, intervals,
                            ncpus='auto', tmp_dir='.',**kwa):
        super(MultiRegionParallelWalker, self).__init__(in_fh,intervals,**kwa)
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
            parsers.append(SerialWalker(self.in_fh, [interval],
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



class SingleRegionParallelWalker(MultiRegionParallelWalker):
    def __init__(self, in_fh, intervals, **kwa):
        assert len(intervals) == 1, ("ParallelSingleRegionWalker requires a "
                                                           "single interval.")
        #logging.info("One interval specified, starting single_region_parallel mode with {} cores.".format(ncpus))
        super(SingleRegionParallelWalker, self).__init__(in_fh, intervals, **kwa)
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
        super(SingleRegionParallelWalker, self).parse_header()
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


##------Parsers---------

#analyses = {}


parser = argparse.ArgumentParser(description="Parse a Variant Call Format (VCF) file.")
parser.add_argument("--variant",'-V',required=True,type = argparse.FileType('r'), 
                            default = '-', help="Input vcf filepath.")

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
parser.add_argument("--skip_multiple_entries",action='store_true',
                     help='Skip all but the first entry for the same site in the VCF.')
parser.add_argument('--logging_level','-l',
                    choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                    default='INFO',
                                            help='Minimun level of logging.')
parser.add_argument('--progress_report_interval',
                    type=int, default=200000,
                    help='Number of lines after which to report progress. '
                                        'Only output if logging level >= INFO')

subparsers = parser.add_subparsers(dest='parser',
                                        help='Type of analysis/parser to run on input file.')

class MetaParser(type):
    """
    This meta-class handles the creation of subparsers
    for each new parser class on class definition.
    """
    def __new__(cls, clsname, bases, attrs):
        newclass = super(MetaParser, cls).__new__(cls, clsname, bases, attrs)
        if clsname != 'Parser':
            
            new_subparser = subparsers.add_parser(clsname,help=newclass.__doc__)
            args = getattr(newclass,'args')
            if args is not None:
                for arg, pars in args.iteritems():
                    new_subparser.add_argument("--"+arg,**pars)
            #Add original Parser init to Child init.
            def parserinit(init):
                def newinit(self,**kwa):
                    Parser.__init__(self,**kwa)
                    init(self,**kwa)
                return newinit
            setattr(newclass,'__init__',parserinit(newclass.__init__))
        return newclass




class Parser(object):
    """
    This is the basic parser object.
    Not to be used directly.
    All parsing tools should derive from this class.
    Parsing tools are supplied to a walker
    to be used to parse the lines of the file.
    """
    __metaclass__ = MetaParser
    args = None
    def __init__(self,**kwa):
        known_args = self.__class__.args if self.__class__.args is not None else {}
        for arg in kwa:
            assert arg in known_args, ("Unknown argument {},"
                                       "possible args are {}".format(arg,known_args.keys()))
        for arg in known_args:
            try:
                a = kwa[arg]
                if a is None:
                    try:
                        nargs = known_args[arg]['nargs']
                        if nargs in ['*','+']:
                            a = []
                    except KeyError:
                        pass
                setattr(self, arg, a)
            except KeyError:
                try:
                    d = known_args[arg]['default']
                except KeyError:
                    req = False
                    try:
                        if known_args[arg]['required']:
                            req = True
                    except KeyError:
                        pass
                    if not req:
                        d = None
                        try:
                            nargs = known_args[arg]['nargs']
                            if nargs in ['*','+']:
                                d = []
                        except KeyError:
                            pass
                    else:
                        raise TypeError("Argument {} not supplied but required.".format(arg))
                #print "Argument {} not supplied, using default {}.".format(arg,d)
                logging.info("Argument {} not supplied, using default {}.".format(arg,d))
                setattr(self,arg,d)
    header_fun = None
    parse_fun = None
    cleanup_fun = None
    reduce_fun = None
    output_fun = None


class VCFTo012(Parser):
    """
    Extract genotype information into a tsv file.
    Coding 0 for homozygous reference, 1 for heterozygote
    and 2 for homozygous alternative allele.
    Missing genotypes are coded by 'N'
    """
    args ={
        'out_tsv':
            {'required':True,
             'type':argparse.FileType('w'),
             'help':"File path to write output tsv to."}
        }

    def header_fun(self,line):
        if line[1:6] == "CHROM":
            self.out_tsv.write("chrom\tpos\t"+line.split("\t",9)[-1])

    def parse_fun(self,sline):
        gt = map(get_012,sline[9:])
        self.out_tsv.write("\t".join(sline[:2])+"\t"+"\t".join(gt)+"\n")

    def reduce_fun(self,selfs):
        command = ["cat"]+[s.out_tsv.name for s in selfs]
        p = subprocess.Popen(command, stdout=self.out_tsv)
        p.communicate()



class FilterByBed(Parser):
    """
    Add filter tag to sites in intervals of bed file.
    This parser is experimental, only working for special cases.
    """
    args = {'in_beds':
                      {'required':True,
                       'nargs':'+','type':argparse.FileType('r'),
                       'help':"List of filepathes to the input beds."},
            'filter_names':{'required':True,'nargs':'+',
                            'help':'Filter names corresponding to beds.'},
            'out_vcf':{'required':True,
                       'type':argparse.FileType('w'),
                       'help':'Filepath to output vcf.'}
            }
    def __init__(self,**kwa):
        #print self.in_beds
        #print self.filter_names
        assert len(self.filter_names)==len(self.in_beds), \
                       "There must be as many filter names as beds."
        self.last_rec = [None for _ in self.in_beds]
        self.seen_chroms = set()

    def header_fun(self,line):
        self.out_vcf.write(line)

    def parse_fun(self,sline):
        def get_rec(fh):
            rec = fh.next().strip().split()
            rec[1] = int(rec[1])
            rec[2] = int(rec[2])
            return rec
        ref = sline[3]
        alt = sline[4].split(',')
        pos = int(sline[1])
        chrom = sline[0]
        self.seen_chroms.update((chrom,))
        for i in range(len(self.last_rec)):
            while self.last_rec[i] is None \
                    or (self.last_rec[i][0]!=chrom and self.last_rec[i][0] in self.seen_chroms) \
                                             or (self.last_rec[i][0]==chrom and self.last_rec[i][2]<pos):
                try:
                    self.last_rec[i] = get_rec(self.in_beds_fh[i])
                except StopIteration:
                    break
            if self.last_rec[i][0] == chrom \
                and self.last_rec[i][1] < pos \
                and self.last_rec[i][2] + 1 > pos:
                if sline[6] in ['.','PASS']:
                    sline[6] = self.filter_names[i]
                else:
                    sline[6] = sline[6] + ','  + self.filter_names[i]
        self.out_vcf_fh.write("\t".join(sline)+'\n')




#--------------------------------------------------------


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
        logging.info("Initialising Walker.")
        parser = Walker(vcf_fh,**kwa)
    elif 'ncpus' in kwa and kwa['ncpus'] > 1:
        assert 'reduce_fun' in analysis['funs'],("This analysis does not support "
                                                 "parallel execution, remove "
                                                 "--ncpus option.")
        if len(kwa['intervals']) == 1:
            logging.info("Initialising SingleRegionParallelWalker.")
            parser = SingleRegionParallelWalker(vcf_fh,**kwa)
        else:
            logging.info("Initialising MultiRegionParallelWalker.")
            parser = MultiRegionParallelWalker(vcf_fh,**kwa)
    else:
        logging.info("Initialising SerialWalker.")
        parser = SerialWalker(vcf_fh,**kwa)
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





if __name__ == "__main__":
    import gzip
    import select
    import time


    args = parser.parse_args()


    logger.setLevel(getattr(logging,args.logging_level))

    args.intervals = [a for b in args.intervals for a in b]

    try:
        assert select.select([args.variant,],[],[],0.0)[0], "Input vcf has no data."
    except TypeError: #the above check does not work for tabix
        pass

    extension = os.path.splitext(args.variant.name)[-1]
    if extension == '.gz':
        args.variant = gzip.open(args.variant.name)
    else:
        assert not args.intervals,("Interval mode (-L) only supported on bgzipped "
                             "files with tabix index. But input has not ending .gz")
        if  extension != '.vcf':
            logging.warning("Unrecognized file extension. "
                            "Assuming this is a vcf: {}".format(args.variant.name))



    parser_class = globals()[args.parser]
    #get attributes specific to subparser
    parser_attr = [a.dest for a in subparsers.choices[args.parser]._actions if a.dest != 'help']
    parser = parser_class(**{arg:getattr(args,arg) for arg in vars(args) if arg in parser_attr})

    if args.intervals is None:

        walker = Walker(args.variant, parser, sep='\t',
                                                  id='',
                                        skip_multiple_entries=args.skip_multiple_entries,
                                        progress_report_interval=args.progress_report_interval)
    else:
        print args.intervals
        walker = SerialWalker(args.variant, parser, args.intervals, sep='\t',
                                                                       id='',
                                            skip_multiple_entries=args.skip_multiple_entries,
                                            progress_report_interval=args.progress_report_interval)


    walker.run()
    #print parser

#        walker = get_parser(args.variant,analyses[args.analysis_type],
#                             arg_dic=arg_dic, intervals = args.intervals,
#                             ncpus=args.ncpus,
#                                skip_multiple_entries=not args.no_skip_multiple_entries,
#                                progress_report_interval=args.progress_report_interval)
#        start = time.time()
#        parser.run()
#        end = time.time()
#        delta = end - start
#        logging.info("This run took {} seconds = {} minutes = {} hours.".format(delta,delta/60.,delta/3600.))





