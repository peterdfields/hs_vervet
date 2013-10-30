#!/usr/bin/env python

import sys, os, errno
import pandas as pd

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper
from vervet.src import VervetDB
from pymodule.yhio.VCFFile import VCFFile 
import hs_vervet_basics as hvb

class DBFetcher(AbstractVervetMapper):       
    __doc__ = __doc__
    """
    get meta-information from Yu's data base and write it to files
    """                        
    
    def __init__(self,genotype_method):
        """
        """
        self.genotype_method=genotype_method
        self.out_dir=os.path.expanduser("~/vervet_data/metadata/genotype_method_{}/".format(genotype_method))
        self.db_dir=os.path.expanduser("~/vervet_data/db/")
        
    
    def get_db_object(self, drivername='postgresql', hostname='crocea.mednet.ucla.edu', \
                db_user='hannes', db_passwd='1o9p2a4', schema='public', dbname='vervetdb', port=None):
        """
        returns a data-base object to connect to Yu's vervet database
        """
        db_vervet = VervetDB.VervetDB(drivername=drivername, username=db_user, password=db_passwd, \
                                    hostname=hostname, database=dbname, schema=schema, port=port)
        db_vervet.setup(create_tables=False)        
        return db_vervet
            
    
        
    def get_vcf_ind(self, uclaidlist, chromosome, format1='VCF'):
        """
        2012.9.19
            get entries of VCF-file that correspond to a sub-population with ucla_id in uclaidlist
            and return genotype matrix
        """
        db_vervet = self.get_db_object()
        session = db_vervet.session
        
        session.begin()
 
        
        genotypeFile = db_vervet.getGenotypeFile(genotype_method_id=self.genotype_method, chromosome=chromosome, format=format1)
        
        if not genotypeFile:
            sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n" % (self.genotype_method, chromosome))
            sys.exit(2)
        filename = os.path.join(self.db_dir, genotypeFile.path)
        if os.path.isfile(filename):
            from pymodule.yhio.VCFFile import VCFFile
            
            vcfFile = VCFFile(inputFname=filename, minDepth=0)
            # this is a list with the read-group names
            readgroupIDList = vcfFile.getSampleIDList()
            
            new_ucla_id_ls = []
            columnIndexList = []
    
            for i in xrange(len(readgroupIDList)):
                readgroupID = readgroupIDList[i]
                # this is the first part of the read group
                individualAlignment = db_vervet.parseAlignmentReadGroup(readgroupID).individualAlignment
                uclaid = individualAlignment.individual_sequence.individual.ucla_id
                if uclaid in uclaidlist:            
                    # header.append(readgroupID)
                    columnIndexList.append(i)
                    new_ucla_id_ls.append(str(uclaid))
            session.close()        
            return (columnIndexList, new_ucla_id_ls)        
                      
    
    def create_chromosome_metadata(self):
        """
        create a list of VCF Filenames for each contig
        """
        db_vervet = self.get_db_object()
        session = db_vervet.session
        query = VervetDB.GenotypeFile.query.filter_by(genotype_method_id=self.genotype_method)
        contig_ls=[]
        VCFfilename_ls=[]
        for entry in query:
            contig_ls.append(str(entry.chromosome))
            VCFfilename_ls.append(str(entry.path))   
        session.close()
        # get list of contig-length from Contig917's (choose any small file) VCF-file header 
        contiglen=hvb.get_contig_length(contig_ls,VCFfilename_ls[contig_ls.index('CAE19')])
        df=pd.DataFrame([[VCFfilename_ls[i],contiglen[i]] for i in range(len(contig_ls))],index=contig_ls,columns=['VCFfname','length']) 
        df.index.name='chromosome'
        return df
    
    def save_chromosome_metadata(self):
        """
        create and save a list of VCF Filenames for each contig
        """
        df=self.create_chromosome_metadata()
        hvb.try_make_dirs(self.out_dir)
        fname = os.path.join(self.out_dir,"chromosome_metadata_m{}.tsv".format(self.genotype_method))
        df.to_csv(fname,sep='\t')   
  
    
    def create_individual_metadata_df(self,chromosome='CAE19'):
        """
        creates a data-frame containing some useful metadata for each individual (see header below)
        """  
        db_vervet = self.get_db_object()
        session = db_vervet.session      
        session.begin()
        try:      
            genotypeFile = db_vervet.getGenotypeFile(genotype_method_id=self.genotype_method, chromosome=chromosome, format='VCF')
                
            if not genotypeFile:
                sys.stderr.write("Error: genotype_method_id %s, chromosome %s does not exist.\n" % (self.genotype_method, chromosome))
                sys.exit(2)
               
            filename = os.path.join(self.db_dir, genotypeFile.path)
                    
            if os.path.isfile(filename):
                           
                #allow 0 depth-> no missing data
                vcfFile = VCFFile(inputFname=filename,minDepth=0)
                sampleIDList = vcfFile.getSampleIDList()
                
                dataMat=[]
                uclaIDList=[]
                
                taxDict=hvb.taxonomic_short_dict()
                countryDict=hvb.country_dict()
                
                header=['VCF_idx','species','country','site_name','longitude','latitude','readgroup','sex','coverage','mean_depth','perc_mapped']
                
                for i in xrange(len(sampleIDList)):
                    sampleID = sampleIDList[i]
                    individualAlignment = db_vervet.parseAlignmentReadGroup(sampleID).individualAlignment
                    if not individualAlignment.individual_sequence.is_contaminated:
                        dataRow=[]
                        #'VCF_idx'
                        dataRow.append(i)
                        
                        
                        species=taxDict[int(individualAlignment.individual_sequence.individual.tax_id)]
                        dataRow.append(species)
                        
                        country=countryDict[int(individualAlignment.individual_sequence.individual.site.country_id)]
                        dataRow.append(country)
                        
                        dataRow.append(individualAlignment.individual_sequence.individual.site.short_name)
                        
                        dataRow.append(individualAlignment.individual_sequence.individual.site.longitude)
                        dataRow.append(individualAlignment.individual_sequence.individual.site.latitude)
                                      
                        dataRow.append(sampleID)
                        
                        dataRow.append(individualAlignment.individual_sequence.individual.sex)
                        
                        dataRow.append(individualAlignment.individual_sequence.coverage)
                        
                        dataRow.append(individualAlignment.mean_depth)
                        
                        dataRow.append(individualAlignment.perc_reads_mapped)
                        
                        uclaIDList.append(individualAlignment.individual_sequence.individual.ucla_id)
                        
                        dataMat.append(dataRow)
                    
                
                metadata=pd.DataFrame(dataMat,index=uclaIDList,columns=header)
                metadata.index.name='ucla_id'
                #[uclaIDList,columnIndexList,species,country,site_row,longitudeList,latitudeList,sampleIDlist]
                return metadata
            else: 
                raise IOError("{} does not exist".format(filename))
        finally:
            session.close()
                

    
    def save_individual_metadata(self):
        """
        saves the metadata matrix to the folder <data_dir>/metadata/method_{}/
        this metadata matrix should be produced once and used for all analyses
        """
        hvb.try_make_dirs(self.out_dir)
        out_fname=os.path.join(self.out_dir,"individual_metadata_m{}.tsv".format(self.genotype_method))
        metadata=self.create_individual_metadata_df()
        metadata.to_csv(out_fname,sep='\t')
        
    def run(self):
        """
        2012.7.13
        """
        pass



if __name__ == '__main__':
    """
    main_class = CalculateStatsForSubPop
    po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
    print po.arguments
    print po.long_option2value
    instance = main_class(po.arguments, **po.long_option2value)
    instance.run()
    """
    #x=hsDBtools()
    #x.saveMetadataMat(40)

