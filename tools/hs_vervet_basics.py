#!/usr/bin/env python
"""
This file is a collection of basic functionality
used by several of hs' vervet tools.
I purged all content on 20130903 to get rid of redundent functions.
If something is missing here, copy it from hs_vervet_basics_old.py.
"""
from __future__ import print_function
import os,gzip,re
import numpy as np

from ast import literal_eval

db_dir=os.path.expanduser("~/vervet_data/data/db/")

#------------
# general
#------------

class HsError( Exception ): 
    """
    a general error class to show that this error is raised by myself
    """
    pass

def try_make_dirs(direc):
    """
    make directory only if it does not exist
    """
    import errno
    direc = os.path.expanduser(direc)
    try:
        os.makedirs(direc)
        #print 'made', direc
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(direc):
            pass
            #print 'exists:', direc
        else: raise

def v_print(text,min_verbosity=0,verbosity=0,file=None,append=True):
    """
    verbose printing with different verbosity levels
    """
    file = os.path.expanduser(file)   
    if verbosity >= min_verbosity:
        if file is None:
            print(text)
        else:
            try:
                print(text,file=file)
            except AttributeError:
                if append:
                    print(text,file=open(file,'a'))
                else:
                    print(text,file=open(file,'w'))
                    
    

#-----------
# tsv I/0
#-----------

def save_genotype(df,fname,*args,**kwargs):
    df.to_csv(fname,sep='\t',index_label='snp_pos',na_rep='N',float_format='%i',*args,**kwargs)

def read_genotype(fname,*args,**kwargs):
    import pandas as pd
    df1=pd.read_csv(fname,sep='\t',index_col=0,na_values='N',*args,**kwargs)
    df1.columns=pd.MultiIndex.from_tuples(df1.columns.map(literal_eval), 
                                              names=['ucla_id','haplotype'])
    return df1.astype(np.float16)
    
def read_metadata(fname,*args,**kwargs):
    import pandas as pd
    df1=pd.read_csv(fname,sep='\t',index_col=0,*args,**kwargs)
    for allele in ['ALT','REF']:
        parsed=df1[allele].apply(lambda s: list(s.replace('[','').replace(']','')))
        df1[allele] = parsed
    df1['POS']=df1['POS']
    return df1

#-------------
# parallel mapping of class methods
#------------

import multiprocessing

def spawn(f):
    def fun(q_in,q_out):
        while True:
            i,x = q_in.get()
            if i == None:
                break
            q_out.put((i,f(x)))
    return fun

def parmap(f, X, nprocs = multiprocessing.cpu_count()):
    q_in   = multiprocessing.Queue(1)
    q_out  = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=spawn(f),args=(q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i,x in sorted(res)]


#----------
# misc
#----------

def get_contig_length(contig_ls,lookup_fname):
    """
    hacky solution to get the contig-length from the header in a VCF file
    """
    vcffname=os.path.join(db_dir,lookup_fname)
    vcfdatafile = gzip.open(vcffname)
    contiglen=[]
    contigname=[]
    occurrence=False
    for line in vcfdatafile.readlines():   
        if line[0:13]=='##contig=<ID=':
            occurrence=True
            splitline=re.split('=|,|>',line)
            contigname.append(splitline[2])
            contiglen.append(splitline[4])
        if occurrence and line[0:13]!='##contig=<ID=':
            break
    return [contiglen[contigname.index(contig)] for contig in contig_ls]

def country_dict():
    countryDict={1:"South Africa",\
     129:"Zambia",\
     130:"Ghana",\
     131:"Ethiopia",\
     132:"Botswana",\
     133:"Central African Republic",\
     134:"Tanzania",\
     135:"United States of America",\
     136:"Barbados",\
     144:"Saint Kitts",\
     148:"Nevis",\
     151:"Gambia",\
     152:"Kenya"}
    return countryDict

def taxonomic_short_dict():
    taxDict={60711:"sab",\
              60712:"tan",\
              101841:"aet",\
              460674:"pyg",\
              460675:"cyn"}
    return taxDict
 
def taxShortToLong(tax_short): 
    return {"sab":"Chlorocebus sabaeus",\
              "tan":"Chlorocebus tantalus",\
              "aet":"Chlorocebus aethiops aethiops",\
              "pyg":"Chlorocebus pygerythrus pygerythrus",\
              "cyn":"Chlorocebus pygerythrus cynosurus"}[tax_short] 



              
                     


