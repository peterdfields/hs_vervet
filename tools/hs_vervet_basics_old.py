#!/usr/bin/python
"""
This file is a collection of basic functionality
used by several of hs' vervet tools.
"""
import socket, os, csv, string, re , gzip, random
import numpy as np
import stats_on_data as sod
import pandas as pd

"""
general low level support functions used in several places
"""

def flatten(list1):
    """
    flattens out all levels of list1
    """
    if isinstance(list1,(np.ndarray,tuple)):
        list1=list(list1)
    if isinstance(list1,list):
        return [a for i in list1 for a in flatten(i)]
    else:
        return [list1]    

def runCommand(cmd_ls,cwd=None):
    import subprocess
    if cwd==None:
        cwd=os.getcwd()
    p = subprocess.Popen(cmd_ls, stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=cwd)
    out, err = p.communicate()
    return (out,err)

#---------------------
#parallel mapping of class methods:
"""
from multiprocessing import Process, Pipe
from itertools import izip

def spawn(f):
    def fun(ppipe, cpipe,x):
        ppipe.close()
        cpipe.send(f(x))
        cpipe.close()
    return fun

def parmap(f,X):
    pipe=[Pipe() for x in X]
    proc=[Process(target=spawn(f),args=(p,c,x)) for x,(p,c) in izip(X,pipe)]
    [p.start() for p in proc]
    ret = [p.recv() for (p,c) in pipe]
    [p.join() for p in proc]
    return ret
"""
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

#----------------------


def try_make_dirs(direc):
    """
    make directory only if it does not exist
    """
    import errno
    try:
        os.makedirs(direc)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(direc):
            pass
        else: raise
        

def copyScriptToOutdir(filepath,outDir):
    import datetime,shutil
    timestamp=datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    scriptname, extension = os.path.splitext(os.path.basename(filepath))
    scriptname=scriptname+'_'+timestamp
    shutil.copyfile(filepath, os.path.join(outDir,scriptname+extension)) 
    
    
def get_data_dir():
    """
    Get the location of the data dir depending on the machine 
    I am working on.
    """    
    if socket.gethostname()=='hsT420s':
        dataDir='/media/Data/Akademisches_data/vervetpopgen'
    else:
        dataDir='/net/gmi.oeaw.ac.at/nordborg/lab/vervetpopgen'
    
    return dataDir


def get_metadata_fname(genotypeMethodID):
    dataDir=get_data_dir()
    anaPath='analyses/_general_input_files/method_{}'.format(genotypeMethodID)
    metaDir=os.path.join(dataDir,anaPath)
    try_make_dirs(metaDir)
    metaFname=os.path.join(metaDir,'metadata_m{}.tsv'.format(genotypeMethodID))
    return metaFname   

def getMetaDF(genotypeMethodID):
    return pd.read_csv(getMetadataFname(genotypeMethodID),sep='\t',index_col=0)

def get_VCFfnames_fname(genotypeMethodID):
    """
    get the full path + filename where the metadata-file with the VCF flilenames and
    chromosome info for a given genotype method is stored or will be stored 
    """
    dataDir=get_data_dir()
    anaPath='analyses/_general_input_files/method_{}'.format(genotypeMethodID)
    metaDir=os.path.join(dataDir,anaPath)
    try_make_dirs(metaDir)
    Fname=os.path.join(metaDir,'VCFfnames_m{}.tsv'.format(genotypeMethodID))
    return Fname

def removeSexContigs(contig_ls):
    sex_contigs=['Contig83','Contig149','Contig193']
    return [contig for contig in contig_ls if contig not in sex_contigs]

def getSexContigs(genotypeMethodID):
    if genotypeMethodID<=40:
        return ['Contig83','Contig149','Contig193']




"""
file i/o
"""

def loadTSVraw(Fname):
    """
    loads a tsv without seperating the row headers
    """
    data=[]
    dataF=open(Fname).readlines()
    for line in dataF:
        data.append(map(string.rstrip,line.split('\t')))
    return data

def dictFromTSV(Fname):
    """
    transforms tsv into a dict of the form  {first row element: rest of row}
    """
    rawdata=loadTSVraw(Fname)
    return {line[0]:line[1:] for line in rawdata}

def writeTSV_direct(Fname,mat):
    """
    writes a tsv file where the header of each row is given in headers
    columns don't have headers
    """
    with open(Fname, 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        for i,row in enumerate(mat):
            writer.writerow(row)
        del writer
    

def writeTSV(Fname,mat,headers):
    """
    writes a tsv file where the header of each row is given in headers
    columns don't have headers
    """
    with open(Fname, 'w') as f:
        writer=csv.writer(f, delimiter='\t')
        for i,row in enumerate(mat):
            writer.writerow([headers[i]]+row)
        del writer

"""
file i/o pandas
"""


def df_to_tsv(Fname,dataFrame):
    with open(Fname, 'w') as f:
        dataFrame.to_csv(f,sep='\t')


def df_from_tsv(Fname):
    with open(Fname, 'r') as f:
        dataFrame=pd.read_csv(f,sep='\t',index_col=0,header=0)
    return dataFrame

def dipl(ucla_id):
    if not hasattr(ucla_id,'__iter__'):
        ucla_id=[ucla_id]
    return np.reshape([(i+'.0',i+'.1') for i in ucla_id],-1)

def ucla_id(dipl):
    if not hasattr(dipl,'__iter__'):
        dipl=[dipl]
    return np.array([i.split(".")[0] for i in dipl])[::2]

"""
specialised functions only used in few places
"""


def countryDict():
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

def taxonomicShortDict():
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

    

def get_contig_length(contig_ls,lookup_fname):
    """
    hacky solution to get the contig-length from the header in a VCF file
    """
    vcffname=os.path.join(get_data_dir(),'db/',lookup_fname)
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



"""
parsing data and metadata dicts
"""
def reduceDict(dic,key,sub_val):
    """
    takes a dictionary where all values are lists of equal length and
    removes the entries that are not in sub_val
    """
    new_dic=dic.copy()
    retain_idx=[i for i in range(len(dic[key])) if dic[key][i] in sub_val]
    for key1,val in dic.items():
        new_dic[key1]=[val[i] for i in range(len(val)) if i in retain_idx]
    return new_dic

def reduceDF(df,remove_idx):
    checkfun=lambda x:x not in remove_idx
    return df[df.index.map(checkfun)] 


def sProp(propname):
    """
    returns the appropriate entry in the propertyList for some defined short-hands
    """
    shorthandDict={'Africa':['country',[el for el in countryDict().values() if el not in ["Barbados","Saint Kitts","Nevis","United States of America"]]],\
                   'Caribbean':['country',["Barbados","Saint Kitts","Nevis"]],\
                   'pygcyn':['species',['pyg','cyn']]}
    return shorthandDict[propname]

def idLSfromPropertyLs(metadataDict,propertyList):
    """
    translates a list of properties, e.g.
    [['country',['South Africa','Gambia']],...,['species',['pyg','sab']]
    into
    [[ucla_id11, ucla_id12,...],...,[ucla_idn1,ucla_idn2,...]]
    """
    idLS=[]
    for prop in propertyList:
        idLS.append([metadataDict['ucla_id'][i] for i,k in enumerate(metadataDict[prop[0]]) if k in prop[1]])
    return idLS
    

def getSubPopMeta(metadataDict,propertyList):
    """
    returns a metadataDict for a sub-population that meets the criteria in property list
    propertyList is of the form [[property1,[possible_val1,possible_val2]],...,[]]
    where each property has to be a key of metadataDict and each value should be a possible value to this key
    e.g.
    [['country',['South Africa','Gambia']],...,['species',['pyg','sab']]
    the inner-most list can be understood as logical OR, the outer layer as logical AND
    """
    
    #by default: remove the contaminated individual:
    noncontaminated_ids=metadataDict['ucla_id'][:]
    noncontaminated_ids.remove('VGHC1006')
    propertyList.append(['ucla_id',noncontaminated_ids])
    
    satisfied=True
    newMetadataDict={key:[] for key in metadataDict.keys()}
    for ind in range(len(metadataDict.values()[0])):
            for prop in propertyList:
                if not metadataDict[prop[0]][ind] in prop[1]:
                    satisfied=False
                    break
            if satisfied:
                for meta in newMetadataDict.keys():
                    newMetadataDict[meta].append(metadataDict[meta][ind])
            satisfied=True
    return newMetadataDict        

def afrSP(metaDF,species):
    #this allows for species identifiers such as ['aet','pygcyn']
    if isinstance(species, list):
        species=",".join(species)
    africa=["South Africa","Zambia","Ghana",\
     "Ethiopia","Botswana","Central African Republic",\
     "Tanzania","Gambia","Kenya"]
    contaminated_id='VGHC1006'
    checkfun=lambda x,y:x in y
    return metaDF[metaDF['species'].apply(checkfun,args=(species,))*\
           metaDF['country'].apply(checkfun,args=(africa,))*(metaDF.index!=contaminated_id)]        

def speciesSubPop(metaDF,species):
    #this allows for species identifiers such as ['aet','pygcyn']
    contaminated_id='VGHC1006'
    if isinstance(species, list):
        species=",".join(species)
    checkfun=lambda x,y:x in y
    return metaDF[metaDF['species'].apply(checkfun,args=(species,))*\
           (metaDF.index!=contaminated_id)]    

#def metaSubPop():

def subSample(metaDF,sample_size):
    sample=random.sample(metaDF.index,sample_size)
    return metaDF.ix[sample]

def getSubPopData(df,ucla_id):
    #try:
    #print ucla_id
    ucla_id=np.array(ucla_id)
    #except:
    #    print ucla_id
    return df[dipl(ucla_id)]

def randomSubSampleMeta(metadataDict,sample_size):
    newMetadataDict={key:[] for key in metadataDict.keys()}
    sample_ind=random.sample(range(len(metadataDict['ucla_id'])),sample_size)
    for key, value in metadataDict.items():
        newMetadataDict[key]=[value[i] for i in sample_ind]
    return newMetadataDict    


def freqFilterDF(df,cutoff):
    aaf=sod.calcAF(df.values)    
    return df[(aaf>cutoff)*(aaf<1-cutoff)]

                             
                             
"""
def getSubPopData(dataDict,ucla_id):

    newDataDict=dataDict.copy()
    idx_ls=np.array([k for k in range(len(dataDict['ucla_id'])) if dataDict['ucla_id'][k] in ucla_id])
    newDataDict['ucla_id']=dataDict['ucla_id'][idx_ls]
    newDataDict['data']=dataDict['data'][:,np.reshape(np.array([idx_ls*2,idx_ls*2+1]).T,-1)]
    return newDataDict
"""

"""
def freqFilter(dataDict,cutoff):
    aaf=sod.calcAF(dataDict['data'])    
    newDataDict=dataDict.copy()
    abovecutoffpos=np.nonzero((aaf>cutoff)*(aaf<1-cutoff))
    abovecutoffpos=abovecutoffpos[0]
    for key in dataDict.keys():
        if isinstance(dataDict[key],(list,tuple,np.ndarray)) and  len(dataDict[key])==len(dataDict['snp_pos']):
            newDataDict[key]=dataDict[key][abovecutoffpos]

    return newDataDict
"""



              
                     


