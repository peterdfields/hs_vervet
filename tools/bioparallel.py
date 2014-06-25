#!/usr/bin/env python
import os
import vcf as pyvcf
from hs_vervet.tools import hs_vervet_basics as hvb
import multiprocessing as mp
import subprocess
import uuid
import collections, itertools, csv
eu = os.path.expanduser
jn = os.path.join


class HeaderlessWriter(pyvcf.Writer):
    """
    This allows to use the pyvcf writer
    without writing the header.
    This class could brake with newer
    versions of pyvcf
    """
    def __init__(self,stream, template, lineterminator="\r\n"):
        self.writer = csv.writer(stream, 
                                 delimiter="\t", 
                                 lineterminator=lineterminator)
        self.template = template
        self.stream = stream
        # Order keys for INFO fields defined in the header (undefined fields
        # get a maximum key).
        self.info_order = collections.defaultdict(
            lambda: len(template.infos),
            dict(zip(template.infos.iterkeys(), itertools.count())))


class VCFParser(object):
    """
    Parse a vcf file using pyvcf
    allowing for parallelisation.
    Chops the chromosome in equally
    sized chunks and processed them in
    parallel.
    If output is a file, then the intermediate
    results are written to tmp files.
    You can supply any parse and
    reduce function.
    ---------------------------------------------
    INPUT:
    vcf_fn: vcf filename
        (should bgzipped and tabix indexed)
    parse_fun: function that walks through the records of the vcf
            it takes the arguments:
                in mode "count": parse_fun(reader)
                in mode "write" of "vcf_write": parse_fun(reader,out_file)
    chrom: chromosome
    chrom_len: length of chromosome, if None, we check whether length
                is defined in the VCF
    reduce_fun: function to reduce the parallel processed output
                there are default reduce_funs, see below
            it takes the arguments:
                in mode "count": reduce_fun(result_list)
                                default is to add up the results
                in mode "write" of "vcf_write": reduce_fun(fn_list,out_file)
                                default is to cat the files together
    mode: 'count': iterate over the vcf without writing to tempfiles
                    e.g., for counting a certain type of entries
          'write': iterate over records and write feature to a file
          'vcf_write': iterate over records and write new vcf
    out_fn: (mode 'write' or 'vcf_write') output filename to write to
    update_vcf_header_fun: (mode 'vcf_write') function to apply to reader to
                                                update vcf header
                it takes the arguments:
                    update_vcf_header_fun(reader)
    tmp_dir: directory to write temporal files to
    ----------------------------------------------
    Examples:
 
    1) extract tsv with columns chrom,pos,ref,alt
        mode = "write"
        
    #define parse function
    def extract_ref_alt(reader, out_file):
        for i,record in enumerate(reader):
            out_file.write(record.CHROM+\
                            '\t'+str(record.POS)+\
                            '\t'+record.REF+'\t'+\
                            str(record.ALT[0])+"\n")
    #intitialise parallel parser
    parser = VCFParser(fn,extract_ref_alt,
                    "CAE29",mode="write",out_fn="parallel_parse.test")
    #parallel run parser with 8 cores
    parser.run(ncpus=8)
    #output can be found in "parallel_parse.test"
    
    2) write a (modified vcf)
        mode = "vcf_write"
        
    #define parse function
    def write_vcf(reader, writer):
        for i,record in enumerate(reader):
            #do something on reader (stupid example)
            record.ALT = "X"
            if record.REF == "A":
                writer.write_record(record)
                
    #define function to update the vcf header
    def update_vcf_info(reader):
        reader.infos.update({'AA':\
                            pyvcf.parser._Info('AA', '1', 'String',
                            'Ancestral Allele")})
    #initialise parallel parser
    parser = VCFParser(fn,write_vcf,
                    chromosome,mode="vcf_write",
                    out_fn="parallel_parse.test",
                    update_vcf_header_fun=update_vcf_info)
    #parallel run parser
    parser.run()
    
        
    3) "Count the number of sites with missing calls for all individuals"
        mode = "count"
        
    
    #define parse function
    def count_accessible_genome(reader):
        genome_dict = {s:{"accessible":0,"non_accessible":0} \
                                        for s in reader.samples}
        for record in reader:
            for call in record.samples:
                if call["GT"] is None:
                    genome_dict[call.sample]["non_accessible"] +=1
                else:
                    genome_dict[call.sample]["accessible"] +=1
        return genome_dict
        
    #define custom reduce function    
    def reduce_accessible_genome(dicts):
        # make copy of sub-dicts
        total_dict = {k:v.copy() for k,v in dicts[0].iteritems()}
        for id in total_dict.keys():
            for cg in total_dict[id].keys():
                total_dict[id][cg] = sum([dic[id][cg] for dic in dicts])
        return total_dict
    
    #intitialise parallel parser
    parser = VCFParser(fn,count_accessible_genome,
                    chromosome,reduce_accessible_genome)
    #parallel run parser
    parser.run()
    #retrieve results
    print parser.out
    
    TODO:
    -check for bgzip/tabix index
    -bgzip/tabix index on the fly
    -define region of the chromosome to be used
    -allow to analyse multiple chromosomes
    -use more standard multiprocessing tools
    -raise on child error
    
    """
    def __init__(self,vcf_fn,parse_fun,chromosomes,chrom_len=None,reduce_fun=None,
                    mode="count",out_fn=None,update_vcf_header_fun=
                    lambda x: None,tmp_dir="."):
        #implement a check whether file is bgzip
        #and has index, ask to create if not
        modes = ["count","write","vcf_write"]
        if mode in ["write","vcf_write"]:
            assert out_fn is not None, "pass out_fn in mode write"
            self.out_fn = out_fn
            if reduce_fun is None:
                reduce_fun = self._default_write_reduce
        elif mode == "count":
            if reduce_fun is None:
                reduce_fun = self._default_count_reduce
        else:
            raise ValueError("mode must be in" + str(modes))
        self.in_fn = vcf_fn
        self.parse_fun = parse_fun
        if type(chromosomes) == str:
            self.chromosomes = [chromosomes]
        else:
            self.chromosomes = chromosomes[:]
        if chrom_len is None:
            self.chrom_len = []
            try:
                reader = pyvcf.Reader(filename=self.in_fn)
                for chrom in self.chromosomes:
                    self.chrom_len.append(reader.\
                                    contigs[chrom].length)
            except:
                print "VCF has no contig length info, provide chrom_len"
                raise
        else:
            if not hasattr(chrom_len, '__iter__'):
                chrom_len = [chrom_len]
            self.chrom_len = chrom_len[:]
        
        
        self.reduce_fun = reduce_fun
        
        self.mode = mode
        self.update_vcf_header = update_vcf_header_fun
        self.tmp_fns = []
        self.tmp_dir = tmp_dir
        
    def __exit__(self, type, value, traceback):
        self.del_tmp_files()

    def check_file(self):
        pass
    
    def prepare_file(self,keep_original=False):
        """
        Prepare vcf for random access:
        bgzip vcf
        tabix index vcf
        """
        #check whether file is not gz yet
        #check whether tabix index exists....
        pass
    #----------default reduce functions--------
    @staticmethod
    def _default_count_reduce(results):
        return reduce(lambda x,y: x+y,results)
    @staticmethod        
    def _default_write_reduce(fns,out_file):
        p = subprocess.Popen(["cat"]+fns,stdout=out_file)
        p.communicate()
    @staticmethod
    def _pythonic_write_reduce(fns,out_file):
        """
        should be equivalent to 
        _default_write_reduce
        but nativly in python
        probably slower
        """
        for fn in fns:
            with open(fn,"r") as f:
                for line in f:
                    out_file.write(line)
    #----------------------------------------
    
    def fetch(self,chrom,start,end):
        chunk = pyvcf.Reader(filename=self.in_fn).fetch(chrom,start,end)
        return chunk
    
    def get_intervals(self,ncpus):
        min_chunks = ncpus
        chunksize = int(sum(self.chrom_len)/min_chunks) + 1
        intervals = []
        for chrom,chrom_len in zip (self.chromosomes,self.chrom_len):
            intervals += [[chrom,start,start+chunksize] \
                            for start in range(0,chrom_len+1,
                                                        chunksize)]
        print intervals
        return intervals
    
    def parse(self,chrom,start,end,fn=None):
        reader = self.fetch(chrom,start,end)
        if self.mode == "count":
            return self.parse_fun(reader)
        elif self.mode == "write":
            with open(fn,"w") as tmp_file:     
                return self.parse_fun(reader,tmp_file)
        elif self.mode == "vcf_write":
            with open(fn,"w") as tmp_file:
                writer =  HeaderlessWriter(tmp_file,reader)
                return self.parse_fun(reader,writer)
        else:
            raise
                
    def del_tmp_files(self):
        while self.tmp_fns:
            os.remove(self.tmp_fns.pop())
        
    def reduce(self):
        if self.mode == "count":
            self.out = self.reduce_fun(self.parallel_out)
        elif self.mode == "write":
            with open(self.out_fn,"w") as out_handle:
                self.out = self.reduce_fun(self.tmp_fns,out_handle)
            #self.del_tmp_files()
        elif self.mode == "vcf_write":
            with open(self.out_fn,"w") as out_file:
                reader = pyvcf.Reader(filename=self.in_fn)
                self.update_vcf_header(reader)
                header_maker = pyvcf.Writer(out_file,reader)                      
            with open(self.out_fn,"a") as out_handle:
                self.out = self.reduce_fun(self.tmp_fns,out_handle)
        else:
            raise
        self.del_tmp_files()
        
    def run(self,ncpus=None):
        try:
            if ncpus is None:
                ncpus = mp.cpu_count()
            intervals = self.get_intervals(ncpus)
            #add temp filename to parse parameters
            if self.mode in ["write","vcf_write"]:
                params = []
                for intv in intervals:
                    tmp_fn = jn(self.tmp_dir,str(uuid.uuid4()) + \
                                        "_" + str(intv[1]) + ".tmp")
                    self.tmp_fns.append(tmp_fn)
                    params.append(intv+[tmp_fn])
            else:
                params = intervals
            print "starting parse"
            self.parallel_out = hvb.parmap(lambda par: 
                                            self.parse(*par),
                                            params,ncpus)
            print "starting reduce"
            self.reduce()
            print "reduce finished"
        finally:
            self.del_tmp_files()


if __name__ == "__main__":
    pass
