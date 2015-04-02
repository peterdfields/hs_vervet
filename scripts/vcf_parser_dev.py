#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.

ATTENTION:
If ever adding a reference to the walker as an attribute of the parser,
then the deepcopy in MultiRegionParallelWalker might make problems.

"""
import sys, os, json, uuid, gzip
import logging, argparse, inspect, copy
import numpy as np
import subprocess
import multiprocessing as mp
logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
logger.setLevel(logging.DEBUG)
eu = os.path.expanduser
jn = os.path.join
try:
    import pandas as pd
except ImportError:
    logging.warning('Python pandas could not be imported. Several parsers will not work.')

try:
    import tabix
except ImportError:
    logging.warning('pytabix could not be imported. No support for interval or parallel use.')

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
    in_fh ... file handle or filename of input delimited
              text file. Can be compressed with gzip or bgzip.
    parser... object instance that derives 
              from Parser class
    Examples:

    """
    def __init__(self,in_fh, parser, sep='\t',
                               id='',
                                    skip_multiple_entries=True,
                                    progress_report_interval=50000,**kwa):
        self.id = id
        if not hasattr(in_fh, 'read'):
            logging.info("Input file has no .read method. Assuming it is a filepath string.")
            extension = os.path.splitext(in_fh)[-1]
            if extension in ['.gz','.bgz']:
                logging.info("{}Input file has extension {}. Opening with gzip.".format(self.id,extension))
                in_fh = gzip.open(eu(in_fh))
            else:
                in_fh = open(eu(in_fh))
        else:
            extension = os.path.splitext(in_fh.name)[-1]
            if extension in ['.gz','.bgz']:
                logging.info("{}Input file has extension {}. Opening with gzip.".format(self.id,extension))
                in_fh = gzip.open(eu(in_fh.name))
        self.in_fh = in_fh
        self.sep = sep
        self.parser = parser
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
                                            "ascending order at: {}:{},{}".format(chrom,self.prev_pos,pos)
            if pos == self.prev_pos:
                self.multiple_count += 1
                if self.multiple_count > 10 and self.print_warning:
                    logging.warning("Omitting further multiple entry warnings.")
                    self.print_warning = False
                if not self.skip_multiple_entries:
                    if self.print_warning:
                        logging.warning("Multiple entries for pos {}:{}. "
                                  "Keeping all entries.".format(chrom,pos))
                    return False
                else:
                    if self.print_warning:
                        logging.warning("{}Warning, multiple entries for pos {}:{}. "
                              "Skipping all but the first.".format(self.id,chrom,pos))
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
        for a in self.parser.line_write_attrs:
            try:
                getattr(self.parser,a).flush()
            except AttributeError:
                pass

    def parse(self,fh):
        """
        """
        if self.parser.parse_fun is not None:
            logging.info("{}Parsing vcf body of {}.".format(self.id,fh))
            line_it = self._yield_split_line(fh)
            for d in line_it:
                self.parser.parse_fun(d)
            #not perfect implementation, prev_pos is not necessarily updated 
            #in children if _skip_duplicate_line is overidden
            logging.info("{}Finished: {} lines at {} {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def cleanup(self):
        """
        """
        if self.parser.cleanup_fun is not None:
            logging.info("Starting cleanup.")
            self.parser.cleanup_fun()
        for a in self.parser.line_write_attrs:
            try:
                getattr(self.parser,a).close()
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
    Walk trough several regions serially.
    Same arguments as Walker, with the addition:
    Input:
    in_fh ... file handle of vcf 
               (must be tabixed and opened with gvcf)
    intervals ... list of string intervals, in Samtools format,
                                            such as Chr1:1-5000
    Note that SerialWalker does not support seperator specification,
    since pytabix does the line splitting automatically.
    """
    def __init__(self,in_fh,parser,intervals,auto_tabix=False, **kwa):
        super(SerialWalker, self).__init__(in_fh, parser, **kwa)
        self.intervals = intervals
        self.auto_tabix = auto_tabix
        self._tabix_init()

    def _tabix_init(self):
        self.tabix_fh = tabix.open(self.in_fh.name)
        self._check_tabix()

    def _check_tabix(self):
        try:
            self._query_tabix(self.intervals[0])
        except tabix.TabixError, e:
            logging.warning("Tabix raised error: {}".format(str(e)))
            if self.auto_tabix:
                logging.warning("Trying to (re-)index file. This can take a while.")
                p = subprocess.Popen(['tabix','-p','vcf','-f',self.in_fh.name],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                if p.returncode != 0:
                    if "was bgzip used to compress this file" in err:
                        base, extension  = os.path.splitext(self.in_fh.name)
                        if extension in ['.gz','.bgz']:
                            logging.warning("File seems not to be compressed with bgzip but ends in .gz or.bgz. "
                                             "Trying to decompress. This can take a while.")

                            p = subprocess.Popen(['gzip','-d','-f',self.in_fh.name],
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            out, err = p.communicate()
                            if p.returncode != 0:
                                logging.error(err)
                                raise e
                            name  = base
                            logging.warning("Trying to compress. This can take a while.")
                        else:
                            logging.warning("File seems not to be compressed with bgzip. "
                                             "Trying to compress. This can take a while.")
                            name = self.in_fh.name

                        p = subprocess.Popen(['bgzip','-f',name],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        out, err = p.communicate()
                        if p.returncode != 0:
                            logging.error(err)
                            raise
                        logging.warning("Trying to index file. This can take a while.")
                        self.in_fh = gzip.open(name+'.gz')
                        p = subprocess.Popen(['tabix','-p','vcf','-f',self.in_fh.name],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        out, err = p.communicate()
                        if p.returncode != 0:
                            logging.error(err)
                            raise
                    else:
                        logging.error("Reindexing failed with unhandeled error: {}".format(err))
                        raise e
                self.tabix_fh = tabix.open(self.in_fh.name)
                try:
                    self._query_tabix(intervals[0])
                except tabix.TabixError, e:
                    logging.error("Failed to auto-tabix input file.")
                    logging.error("Is the interval {} in the vcf? ".format(self.intervals[0]))
                    raise e
                logging.info("Auto-tabix successful.")
            else:
                logging.error("Is file compressed with bgzip and does it have tabix index? "
                                " If not, produce it or run with flag/argument auto_tabix.")
                logging.error("Is the interval {} in the vcf? ".format(self.intervals[0]))
                raise e


    def _query_tabix(self,interval):
        try:
            return self.tabix_fh.querys(interval)
        except TypeError:
            try:
                return self.tabix_fh.query(*interval)
            except TypeError, e:
                logging.error("Interval: {}: Chromosome must be string and position integer."
                                               " e.g. ('Chr1',1000000,1005000). "
                                "Alternatively use string 'Chr:1000000-1005000'".format(interval))
                raise e
            except tabix.TabixError, e:
                logging.error("Is the interval {} in the vcf? ".format(interval))
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


#class MultiRegionParallelWalker(SerialWalker):
#    def __init__(self, in_fh, parser, intervals,
#                            ncpus='auto', tmp_dir='.',**kwa):
#        super(MultiRegionParallelWalker, self).__init__(in_fh,parser,intervals,**kwa)
#        assert hasattr(self.parser,'reduce_fun'), ("Parser {} has no reduce_fun method. "
#                                                 "No support for parallel execution.")
#        if self.parser.reduce_fun is None:
#            logging.warning("Reduce function is None. Is this really intended?")
#        self.kwa = kwa
#        if ncpus == 'auto':
#            ncpus = mp.cpu_count()
#        self.ncpus = ncpus
#        self.tmp_dir = os.path.expanduser(tmp_dir)
#
#    def set_ncpus(self):
#        self.ncpus = min(len(self.intervals),self.ncpus)
#        logging.info("{} regions specified.  Parallelising by region. "
#                     "Using {} processes.".format(len(self.intervals),self.ncpus))
#
#    def setup_subwalkers(self):
#        subwalkers = []
#        temp_fns = []
#        for i,interval in enumerate(self.intervals):
#            subparser = copy.deepcopy(self.parser)
#            subparser.header_fun = None #don't parse header in subparser
#            parse_source = inspect.getsource(subparser.parse_fun)
#            for a in subparser.line_write_attrs:
#                tmp_fn = jn(self.tmp_dir,a+'_'+str(uuid.uuid4()) + \
#                                         "_" + str(i) + ".tmp")
#                temp_fns.append(tmp_fn)
#                setattr(subparser,a+'_fn',tmp_fn)
#            subwalkers.append(SerialWalker(self.in_fh, subparser, [interval],
#                                         id="Interval {}:".format(interval), **self.kwa))
#        self.subwalkers = subwalkers
#        self.temp_fns = temp_fns
#
#    def reduce(self):
#        if self.parser.reduce_fun is not None:
#            logging.info("Starting reduce step.")
#            print [p.out_tsv for p in self.subparsers]
#            self.parser.reduce_fun([p for p in self.subparsers])
#        else:
#            logging.warning("No reduce function given. Won't do anything.")
#
#    def del_temp_files(self):
#        logging.info("Removing temp files.")
#        while self.temp_fns:
#            os.remove(self.temp_fns.pop())
#
#    def run_parser(self,i):
#        s = self.subwalkers[i]
#        for a in s.parser.line_write_attrs:
#            setattr(s.parser,a,open(getattr(s.parser,a+'_fn'),'w'))
#        s.run_no_output()
#        return s.parser
#
#    def run(self):
#        self.parse_header()
#        self.set_ncpus()
#        self.setup_subwalkers()
#        subparsers = parmap(self.run_parser,range(len(self.subwalkers)),self.ncpus)
#        #this is a hacky solution acconting for the fact that open files cannot be
#        #pickled and sent to the child processes by multiprocessing
#        for p in subparsers:
#            for a in p.line_write_attrs:
#                setattr(p,a,open(getattr(p,a+'_fn'),'a'))
#        self.subparsers = subparsers
#        self.reduce()
#        logging.info("Creating output.")
#        self.output()
#        self.del_temp_files()
#        logging.info("Run finished.")


class ParallelWalker(SerialWalker):
    """
    chunk ...        if chunk is False, intervals are not divided into chunks, i.e.,
                     multiprocessing runs with min(ncpus,len(intervals))
    chunk_factor ... Multiply ncpus by this factor to get the number of chunks.
                     Larger value will imply overhead,
                     but helps to use all processors if some chromosomes are very short.
    """
    def __init__(self, in_fh, parser, intervals=None, ncpus='auto', auto_tabix=False,
                                                            tmp_dir='.',chunk=True,
                                                             chunk_factor=2, **kwa):
        if intervals is None:
            #setup without any interval specific things
            super(SerialWalker, self).__init__(in_fh, parser, **kwa)
            self.intervals = None
            self.auto_tabix = auto_tabix
        else:
            super(ParallelWalker, self).__init__(in_fh, parser,intervals, **kwa)
        self.chunk = chunk
        assert hasattr(self.parser,'reduce_fun'), \
                                            ("Parser {} has no reduce_fun method. "
                                             "No support for parallel execution."\
                                                .format(self.parser.__class__.__name__))
        if self.parser.reduce_fun is None:
            logging.warning("Reduce function is None. Is this really intended?")
        self.kwa = kwa
        if ncpus == 'auto':
            ncpus = mp.cpu_count()
        self.ncpus = ncpus
        self.tmp_dir = os.path.expanduser(tmp_dir)
        if chunk:
            self.n_chunks = chunk_factor * self.ncpus
            self.missing = True
            if self.intervals is None:
                self.replace_missing = self.replace_missing_no_interval
                logging.info("No intervals specificed for parallel parsing. "
                             "Trying to infer intervals by searching vcf header. "
                                "Hence, tabix checks will happen at runtime only.")
            else:
                self.intervals = [self._parse_interval(i) for i in self.intervals]
                if not any([i[2] is None for i in self.intervals]):
                    self.missing = False
            if self.missing:
                self._header_line_parser = self._header_line_parser_search_contig_len
        else:
            self.ncpus = min(len(self.intervals),self.ncpus)
            logging.info("Chunking is turned off. {} intervals specified. "
                                                "Parallelising by region. "
                     "Using {} processes.".format(len(self.intervals),self.ncpus))
        self.contic_dic = {}


    def _parse_interval(self,interval):
        try:
            chrompos = interval.split(":")
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
            chrom = interval[0]
            start = interval[1]
            end = interval[2]
        if start is None:
            start = 0
        return [chrom, start, end]


    def _header_line_parser_search_contig_len(self,line):
        if line[:9] == '##contig=':
                    dic = get_header_line_dic(line)
                    self.contic_dic.update({dic['ID']:int(dic['length'])})
        if self.parser.header_fun is not None:
            self.parser.header_fun(line)

    def replace_missing_no_interval(self):
        logging.info("No intervals given, considering all contigs given in VCF header.")
        self.intervals = [(k,1,v) for k,v in  self.contic_dic.iteritems()]
        super(ParallelWalker, self)._tabix_init()

    def replace_missing(self):
        for interval in self.intervals:
            if interval[2] is None:
                interval[2] = self.contic_dic[interval[0]]
            if interval[1] is None:
                interval[1] = 1

    def get_chunks(self):
        chunks = self.intervals
        while len(chunks) < self.n_chunks:
            lengths = [i[2]-i[1] for i in chunks]
            idx = np.argmax(lengths)
            longest = chunks[idx]
            midpoint = longest[1]+(longest[2]-longest[1])/2
            left_chunk = [longest[0],longest[1],midpoint]
            right_chunk = [longest[0],midpoint+1,longest[2]]
            chunks[idx] = left_chunk
            chunks.insert(idx+1,right_chunk)
        #chunks.sort()
        logging.debug("Chunks used: {}".format(chunks))
        return chunks

    def setup_subwalkers(self):
        subwalkers = []
        temp_fns = []
        for i, chunk in enumerate(self.chunks):
            subparser = copy.deepcopy(self.parser)
            subparser.header_fun = None #don't parse header in subparser
            parse_source = inspect.getsource(subparser.parse_fun)
            for a in subparser.line_write_attrs:
                tmp_fn = jn(self.tmp_dir,a+'_'+str(uuid.uuid4()) + \
                                         "_" + str(i) + ".tmp")
                temp_fns.append(tmp_fn)
                setattr(subparser,a+'_fn',tmp_fn)
            subwalkers.append(SerialWalker(self.in_fh, subparser, [chunk],
                                         id="Chunk {}: ".format(chunk), **self.kwa))
        self.subwalkers = subwalkers
        self.temp_fns = temp_fns

    def reduce(self):
        if self.parser.reduce_fun is not None:
            logging.info("Starting reduce step.")
            self.parser.reduce_fun([p for p in self.subparsers])
        else:
            logging.warning("No reduce function given. Won't do anything.")

    def del_temp_files(self):
        logging.info("Removing temp files.")
        while self.temp_fns:
            os.remove(self.temp_fns.pop())

    def run_parser(self,i):
        s = self.subwalkers[i]
        for a in s.parser.line_write_attrs:
            setattr(s.parser,a,open(getattr(s.parser,a+'_fn'),'w'))
        s.run_no_output()
        return s.parser

    def run(self):
        self.parse_header()
        if self.chunk:
            if self.missing:
                self.replace_missing()
            self.chunks = self.get_chunks()
        else:
            self.chunks = self.intervals
        self.setup_subwalkers()
        subparsers = parmap(self.run_parser,range(len(self.subwalkers)),self.ncpus)
        #this is a hacky solution acconting for the fact that open files cannot be
        #pickled and sent to the child processes by multiprocessing
        for p in subparsers:
            for a in p.line_write_attrs:
                #attention the files are opened in read mode now
                setattr(p,a,open(getattr(p,a+'_fn'),'r'))
        self.subparsers = subparsers
        self.reduce()
        logging.info("Creating output.")
        self.output()
        self.del_temp_files()
        logging.info("Run finished.")

#class SingleRegionParallelWalker(MultiRegionParallelWalker):
#    def __init__(self, in_fh, parser, intervals, **kwa):
#        assert len(intervals) == 1, ("ParallelSingleRegionWalker requires a "
#                                                           "single interval.")
#        #logging.info("One interval specified, starting single_region_parallel mode with {} cores.".format(ncpus))
#        super(SingleRegionParallelWalker, self).__init__(in_fh, parser, intervals, **kwa)
#        chrom, start, end = self._parse_interval(intervals[0])
#        self.contigs
#
#    def _parse_interval(self,interval):
#        try:
#            chrompos = interval.split(":")
#            chrom = chrompos[0]
#            try:
#                startend = chrompos[1].split('-')
#                start = int(startend[0])
#                try:
#                    end = int(startend[1])
#                except:
#                    end = None
#            except IndexError:
#                start = None
#                end = None
#        except TypeError:
#            chrom = interval[0]
#            start = interval[1]
#            end = interval[2]
#        if start is None:
#            start = 0
#        if end is None:
#            self._header_line_parser = self._header_line_parser_search_contig_len
#            logging.info("No chromosome end specified, searching for 'contig'"
#                         "tag in vcf header.")
#        return (chrom, start, end)
#
#
#    def parse_header(self):
#        super(SingleRegionParallelWalker, self).parse_header()
#        self.intervals = self.get_intervals()
#
#    def _header_line_parser_search_contig_len(self,line):
#        if line[:9] == '##contig=':
#                    contig_dic = get_header_line_dic(line)
#                    if contig_dic['ID'] == self.chrom:
#                        length = int(contig_dic['length'])
#                        self.end = length
#        if self.parser.header_fun is not None:
#            self.parser.header_fun(line)
#
#    def get_intervals(self):
#        n_chunks = self.ncpus
#        chunksize = int((self.end-self.start)/n_chunks) + 1
#        starts = []
#        for s in range(self.start,self.end,chunksize):
#            starts.append(s)
#        intervals = [(self.chrom,s,e) for s,e in zip(starts,[s-1 for s in starts[1:]]+[self.end])]
#        return intervals



#--------------SUPPORT FUNCTIONS-------------------------

def get_walker(in_fh,parser,intervals=None,ncpus=None,**kwa):
    if ncpus is None or ncpus <= 1:
        if intervals is None:
            return Walker(in_fh,parser,**kwa)
        else:
            if not intervals:
                logging.warning("Intervals given but empty, Walker won't do anything.")
            return SerialWalker(in_fh, parser, intervals,**kwa)
    else:
        return ParallelWalker(in_fh,parser,intervals=intervals,ncpus=ncpus,**kwa)

#parallel support

def fun(f,q_in,q_out):
    while True:
        i,x = q_in.get()
        if i is None:
            break
        q_out.put((i,f(x)))

def parmap(f, X, nprocs):
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


#---------support functions used by several parsers------------

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

##--------------------Parsers-----------------------


void_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = void_parser.add_subparsers(dest='parser')


parser_classes = {}
def register(cls):
    parser_classes[cls.__name__] = cls



class MetaParser(type):
    """
    This meta-class handles the creation of subparsers
    for each new parser class on class definition.
    """
    def __new__(cls, clsname, bases, attrs):
        newclass = super(MetaParser, cls).__new__(cls, clsname, bases, attrs)
        if 'register' not in attrs.keys():
            attrs['register'] = True
        if attrs['register']:
            register(newclass)
            new_subparser = subparsers.add_parser(clsname, description=newclass.__doc__)
            args = getattr(newclass,'args')
            newclass.line_write_attrs = []
            newclass.end_write_attrs = []
            newclass.line_read_attrs = []
            if args is not None:
                for arg, pars in args.iteritems():
                    new_subparser.add_argument("--"+arg,**pars)
                    #if argument is a file open in write/read mode 
                    #and it is used in parse_fun --> add it to line_write_attr/ line_read_attrs
                    try:
                        t = pars['type']
                        try:
                            mode = t._mode
                        except AttributeError:
                            try:
                                mode = t.mode
                            except AttributeError:
                                continue
                        if  newclass.parse_fun is not None:
                            parse_code = inspect.getsource(newclass.parse_fun)
                            if arg in parse_code:
                                if mode == 'w':
                                    newclass.line_write_attrs.append(arg)
                                elif mode == 'r':
                                    newclass.line_read_attrs.append(arg)
                                else:
                                    logging.warning("Arg {} of parser class {} has _mode or mode "
                                                    "attribute but mode is not 'w' or 'r'. "
                                                    "This will cause problems in parallel mode "
                                                    "if there is any read or write "
                                                    "happening within the parse method.".format(arg,clsname))
                            #check for files that are written in the end
                            elif mode == 'w' and newclass.output_fun is not None:
                                parse_code = inspect.getsource(newclass.output_fun)
                                if arg in parse_code:
                                    newclass.end_write_attrs.append(arg)
                    except KeyError:
                        pass

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
    register = False
    def __init__(self,**kwa):
        known_args = self.__class__.args if self.__class__.args is not None else {}
        for arg in kwa:
            assert arg in known_args, ("Unknown argument {},"
                                       "possible args are {}".format(arg,known_args.keys()))
        #this is not really necessary, but useful if line_write_attr are changed in instance
        self.line_write_attrs = copy.copy(self.__class__.line_write_attrs)
        self.end_write_attrs = copy.copy(self.__class__.end_write_attrs)
        self.line_read_attrs = copy.copy(self.__class__.line_read_attrs)
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
            except KeyError:
                try:
                    a = known_args[arg]['default']
                    logging.info("Argument {} not supplied, using default {}.".format(arg,a))
                except KeyError:
                    req = False
                    try:
                        if known_args[arg]['required']:
                            req = True
                    except KeyError:
                        pass
                    if not req:
                        a = None
                        try:
                            nargs = known_args[arg]['nargs']
                            if nargs in ['*','+']:
                                a = []
                        except KeyError:
                            pass
                    elif arg in self.end_write_args:
                        a = None
                        logging.info("Argument {} not supplied. It looks like a file that is used for final output. "
                                     "Setting it None and assuming that method output_fun is returning to variable. ")
                    else:
                        raise TypeError("Argument {} not supplied but required.".format(arg))
            setattr(self,arg,a)


        for attr in self.line_write_attrs:
            fh = getattr(self,attr)
            if not hasattr(fh, 'write'):
                logging.info("{} has no .write method. Assuming it is a filepath string: {}.".format(attr,fh))
                fh = open(eu(fh),'w')
            base, extension = os.path.splitext(fh.name)
            if extension in ['.gz','.bgz']:
                logging.info("{} ends in {}, compressing with tabix on the fly.".format(attr, extension))
                base1, ext1 = os.path.splitext(base)
                if ext1 and ext1[1:] in['gff', 'bed', 'sam', 'vcf', 'psltbl', 'gvcf']:
                    logging.info("{} has ending {}, supposing that it is reference ordered data."
                                                                    " Will create tabix index.".format(attr,ext1))
                    index = 'vcf' if ext1[1:] == 'gvcf' else ext1[1:]
                else:
                    index = False
                setattr(self,attr,TabixWrite(fh,index))
    header_fun = None
    parse_fun = None
    cleanup_fun = None
    output_fun = None

class TabixWrite(object):
    """
    A roughly file-like object,
    That pipes everything that is
    writen to it to bgzip in order to
    compress it.

    Input:
    fh ... file handle
    index ... (bool) Also create tabix index on closing the file.
    """

    def __init__(self,fh, index):
        #idea: inherit all attributes from fh?
        self.fh = fh
        self.name = self.fh.name
        bgzip = subprocess.Popen(['bgzip','-c'], stdin=subprocess.PIPE, stdout=fh,stderr=subprocess.PIPE)
        self.bgzip = bgzip
        self.index = index

    def write(self,string):
        self.bgzip.stdin.write(string)

    def flush(self):
        self.bgzip.stdin.flush()

    def close(self):
        out, err = self.bgzip.communicate()
        self.bgzip.stdin.close()
        if self.bgzip.returncode:
            logging.error("Failed to bgzip {}: {}".format(self.fh.name, err))
            raise IOError
        if self.index:
            logging.info("Creating tabix index for {}: {}".format(self.fh.name, err))
            tabix = subprocess.Popen(['tabix','-p',self.index,self.fh.name],stderr=subprocess.PIPE)
            out, err = tabix.communicate() 
            if self.bgzip.returncode:
                logging.warning("Failed to create tabix index for {}: {}".format(self.fh.name, err))

class ReduceError(Exception):
    pass

def line_write_reduce_cat(filenames, out_fh):
    command = ["cat"] + filenames
    p = subprocess.Popen(command, stdout=out_fh)
    out, err = p.communicate()
    if p.returncode:
        raise ReduceError("Cat reduce step had non-zero exit status: {}.".format(err))

def line_write_reduce_python(file_handles, out_fh):
    for fh in file_handles:
        for line in fh:
            out_fh.write(line)

class LineWriteParser(Parser):
    """
    Standard parser that reads from VCF
    and writes to one file line by line 
    (vcf or any other).

    Subclass it and add custom parse method!
    """
    register = False
    args ={
        'out_file':
            {'required':True,
             'type':argparse.FileType('w'),
             'help':"File path to write output to."}
        }
    def reduce_fun(self,selfs):
        try:
            line_write_reduce_cat([s.out_tsv.name for s in selfs], self.out_tsv)
        except ReduceError:
            logging.warning("Cat reduce step had non-zero exit status: {}."
                            "Trying pytonic reduce.".format(err))
            line_write_reduce_python([s.out_tsv for s in selfs], self.out_tsv)

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
                    self.last_rec[i] = get_rec(self.in_beds[i])
                except StopIteration:
                    break
            if self.last_rec[i][0] == chrom \
                and self.last_rec[i][1] < pos \
                and self.last_rec[i][2] + 1 > pos:
                if sline[6] in ['.','PASS']:
                    sline[6] = self.filter_names[i]
                else:
                    sline[6] = sline[6] + ','  + self.filter_names[i]
        self.out_vcf.write("\t".join(sline)+'\n')


class GetFilterStats(Parser):
    """
    Count occurences of all combinations
    of filters in the filter column.
    """
    args = {'out_fn':{'required':True,
                      'type':argparse.FileType('w'),
                      'help':"Filename to write to."}}

    def __init__(self,**kwa):
        try:
            pd
        except NameError, e:
            logging.error("{} required python pandas, but it seems not to be imported.")
            raise e
        self.count_dic = {}

    def parse_fun(self,sline):
        filters = sline[6].split(';')
        filters.sort()
        filters = tuple(filters)
        add_to_countdic(self.count_dic,'n_sites')
        add_to_countdic(self.count_dic,filters)

    def cleanup_fun(self):
        filter_info = pd.Series(self.count_dic.values(),
                                    index=self.count_dic.keys())
        filter_info.sort(inplace=True,ascending=False)
        self.filter_info = filter_info

    def reduce_fun(self,selfs):
        fi = selfs[0].filter_info
        for ad in selfs[1:]:
            fi = fi.add(ad.filter_info, fill_value=0)
        self.filter_info = fi

    def output_fun(self):
        if self.out_fn is not None:
            self.filter_info.to_csv(self.out_fn,sep='\t')
        else:
            return self.filter_info


class SNPEFFParser(Parser):
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


    CHR_IX =0
    POS_IX =1
    REF_IX =3
    ALT_IX =4
    INFO_IX =7
    AT_START_IX = 9
    AT_LENGTH = 1135
    AT_END_IX = AT_START_IX+ AT_LENGTH
    LYR_START_IX = 1171
    LYR_END_IX = LYR_START_IX +2
    REF = "0"
    ALT = "1"
    DOUBLE_REF = "00"
    DOUBLE_ALT = "11"

    def header_fun(self,line):
        if line[1:6] == "CHROM":
            self.out_tsv.write("chrom\tpos\tref\talt\tlyr\tgenotype\tinfo\n")

    def parse_fun(self,sline):
        output = self._filter_(sline)
        if output is not None:
            self.out_tsv.write("\t".join(output)+"\n")

    def reduce_fun(self,selfs):
        command = ["cat"]+[s.out_tsv.name for s in selfs]
        p = subprocess.Popen(command, stdout=self.out_tsv)
        p.communicate()


    def _filter_(self,fields):
        if fields[self.ALT_IX] == ".":
            return None
        if len(fields[self.REF_IX]) > 1 :
            return None
        lookupMap = self._create_set_(fields[self.AT_START_IX:self.AT_END_IX])
        if len(lookupMap) != 2:
            return None
        alt = fields[self.ALT_IX]
        genotype = "1"
        if len(alt) > 1:
            alt_split = alt.split(",")
            if len(alt_split)  == 1:
                return None
            keys = set(lookupMap.keys())
            keys = sorted(keys)
            genotype = keys[-1][0]
            ix = int(genotype)
            if len(alt_split) >= ix:
                alt = alt_split[ix-1]
                if len(alt) != 1:
                    return None
            else:
                return None
        lookupMap = self._create_set_(fields[self.LYR_START_IX:self.LYR_END_IX])
        lyr = self._get_lyr_allele_(lookupMap)
        return [fields[self.CHR_IX],fields[self.POS_IX],fields[self.REF_IX],alt,lyr,genotype,fields[self.INFO_IX]]



    def _create_set_(self,fields):
        lookupMap = {}
        for item in fields:
            sep = item[1]
            isValid = False
            if sep == '|':
                isValid = True
            elif sep == '/':
                if item[0] != '.' and item[2] != '.':
                    isValid = True
            if not isValid:
                continue
            genotype = item[0]+item[2]
            lookupMap[genotype] = True
        return lookupMap

    def _get_lyr_allele_(self,lookupMap):
        size = len(lookupMap)
        if size == 2:
            return 'S'
        elif size == 1:
            if self.DOUBLE_REF in lookupMap:
                return '0'
            if self.DOUBLE_ALT in lookupMap:
                return '1'
            if self.REF+self.ALT in lookupMap:
                return 'S'
            if self.ALT+self.REF in lookupMap:
                return 'S'
        return 'NA'







def get_argparser():
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description="Parse a Variant Call Format (VCF) file.",
                                                                                 add_help=False)


    argparser.add_argument("--variant",'-V', type = argparse.FileType('r'),
                                #default = '-',
                                 help="Input vcf filepath.")

    argparser.add_argument("--parser",'-P', choices = subparsers.choices.keys(),
                                                                help="Subparser to be used.")

    argparser.add_argument("--intervals",'-L', nargs='*', dest='intervals', action='append',
                            help='Specify intervals to consider e.g. Chr1:1-50000. '
                                 'Input vcf must be compressed with bgzip and indexed '
                                                                         'with tabix.')
    argparser.add_argument("--auto_tabix",action='store_true',
                                help="Automatically compress and/or index on the fly.")
    argparser.add_argument("--ncpus", '-nct',
                            type=int, default=1,
                                  help='Number of processes for parallel parsing. '
                                       'If no intervals are specified with '
                                       'with -L (e.g. -L Chr1:1-5000000), '
                                       'or if intervals lack the end position (e.g. -L Chr1), '
                                       'then we try to infer all regions present in the VCF '
                                       'from the VCF header (tag contig=...). '
                                       'Make the header contig tags are consisten with the '
                                       'body in this case. '
                                       'Parallel parsing is not implemented for all '
                                       'parsers.')
    argparser.add_argument("--no_chunk", action='store_true',
                            help="If no_chunk is set, intervals are not divided into chunks, i.e., "
                                 "multiprocessing runs with min(ncpus,len(intervals))")
    argparser.add_argument("--chunk_factor", default=6,
                            help="Multiply ncpus by this factor to get the number of chunks. "
                                  "Larger value will imply overhead, "
                                  "but helps to use all processors for the whole runtime "
                                  "when chromosomes or intervals are of unequal length.")
    argparser.add_argument("--temp_dir",
                                 default='.',
                                  help= "Path to folder to write temporary files to. "
                                        "Only used when parallel processing.")
    argparser.add_argument("--skip_multiple_entries",action='store_true',
                         help='Skip all but the first entry for the same site in the VCF.')
    argparser.add_argument('--logging_level','-l',
                        choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                        default='INFO',
                                                help='Minimun level of logging.')
    argparser.add_argument('--progress_report_interval',
                        type=int, default=200000,
                        help='Number of lines after which to report progress. '
                                            'Only output if logging level >= INFO')
    argparser.add_argument("--help",'-h', action='store_true',
                                                     help="Print help message and exit.")


    return argparser

def parse_args(argparser):
    args, unknown = argparser.parse_known_args()

    if args.help:
        help_str = "\n"+argparser.format_help()
        logger.setLevel(logging.INFO)
    try:
        subargparser = subparsers.choices[args.parser]
        if args.help:
            help_str += "\n\n-----------Help for parser {} ------------\n\n".format(args.parser) \
                                                                         + subargparser.format_help()
    except KeyError, e:
        if args.help:
            pass
        elif args.parser is None:
            argparser.error("argument --parser/-P is required")
        else:
            raise e

    if args.help:
       logging.info(help_str)
       sys.exit(0)

    if args.variant is None:
        argparser.error("argument --variant/-V is required")

    #if args.variant is sys.stdin:
    #    logging.warning("argument --variant/-V not supplied, trying to read from stdin")

    sub_args = subargparser.parse_args(unknown)
    return args, sub_args

def parse(args,sub_args):
    import select
    import time
    logger.setLevel(getattr(logging,args.logging_level))
    if args.intervals is not None:
        args.intervals = [a for b in args.intervals for a in b]

    try:
        assert select.select([args.variant,],[],[],0.0)[0], "Input vcf has no data."
    except TypeError: #the above check does not work for tabix
        pass

    parser_class = parser_classes[args.parser]
    parser = parser_class(**{arg:getattr(sub_args,arg) for arg in vars(sub_args)})


    walker = get_walker(args.variant, parser, intervals=args.intervals,
                                            auto_tabix=args.auto_tabix,
                            chunk=not args.no_chunk, chunk_factor=args.chunk_factor,
                                 tmp_dir=args.temp_dir, ncpus=args.ncpus,
                           skip_multiple_entries=args.skip_multiple_entries,
                            progress_report_interval=args.progress_report_interval)

#    if args.intervals is None:
#        if args.ncpus > 1:
#            logging.warning("Ignoring --ncpus! Multiprocessing is only supported "
#                                            "if at least one interval is specified, "
#                                                                    "e.g., -L Chr1.")
#        if args.auto_tabix:
#            logging.warning("Flag --auto_tabix given, but not in "
#                            "interval (-L) mode, won't auto_tabix.")
#        walker = Walker(args.variant, parser, sep='\t',
#                                        skip_multiple_entries=args.skip_multiple_entries,
#                                        progress_report_interval=args.progress_report_interval)
#    else:
#        if args.ncpus <= 1:
#            walker = SerialWalker(args.variant, parser, args.intervals, sep='\t',auto_tabix=args.auto_tabix,
#                                                skip_multiple_entries=args.skip_multiple_entries,
#                                                progress_report_interval=args.progress_report_interval)
#        elif len(args.intervals) > 1:
#            walker = MultiRegionParallelWalker(args.variant, parser, args.intervals, sep='\t', auto_tabix=args.auto_tabix,
#                                                                                   tmp_dir=args.temp_dir,ncpus=args.ncpus,
#                                                    skip_multiple_entries=args.skip_multiple_entries,
#                                                    progress_report_interval=args.progress_report_interval)
#        else:
#            walker = SingleRegionParallelWalker(args.variant, parser, args.intervals, sep='\t', auto_tabix=args.auto_tabix,
#                                                                                     tmp_dir=args.temp_dir,ncpus=args.ncpus,
#                                                    skip_multiple_entries=args.skip_multiple_entries,
#                                                    progress_report_interval=args.progress_report_interval)
    logging.info("Using {} to traverse the file.".format(walker.__class__.__name__))

    start = time.time()
    walker.run()
    end = time.time()
    delta = end - start
    logging.info("This run took {} seconds = {} minutes = {} hours.".format(delta,delta/60.,delta/3600.))



def main():
    argparser = get_argparser()
    args, sub_args = parse_args(argparser)
    try:
        parse(args,sub_args)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        logging.exception(e)
        return 2


if __name__ == "__main__":
    sys.exit(main())








