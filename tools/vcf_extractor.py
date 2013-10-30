#!/usr/bin/env python
"""
"""

import sys, os
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

#from vervet.src.mapper.AbstractVervetMapper import AbstractVervetMapper

import numpy as np
import pandas as pd
import hs_vervet_basics as hvb

#AbstractVervetMapper

class VCFExtractor():       
    __doc__ = __doc__
    """
    object that extracts genotypes or other information from VCF files
    """                        
    
    def __init__(self,genotype_method=233):
        self.genotype_method=genotype_method
        self.db_dir=os.path.expanduser("~/vervet_data/db/")
        self.data_dir=os.path.expanduser("~/vervet_data/")
        self.meta_dir=os.path.expanduser("~/vervet_data/metadata/genotype_method_{}/".format(genotype_method))
        ind_meta_fname=os.path.join(self.meta_dir,"individual_metadata_m{}.tsv".format(genotype_method))
        self.individual_metadata=pd.read_csv(ind_meta_fname,sep='\t',index_col=0,header=0)
        
        chrom_meta_fname=os.path.join(self.meta_dir,"chromosome_metadata_m{}.tsv".format(genotype_method))
        self.chromosome_metadata=pd.read_csv(chrom_meta_fname,sep='\t',index_col=0,header=0)
        
    def extract_genotype_data(self,chromosome,qualityfilter=0,minDepth=0):
        """
        2012.9.19
            get entries of VCF-file that correspond to a sub-population with vcf file index in vcf_indx_ls
            and return genotype dictionary
        """
        #import pdb
        vcf_indx_ls=self.individual_metadata['VCF_idx'].values
        ucla_id_ls=self.individual_metadata.index.values
        genotype_method=self.genotype_method
        vcf_fname=self.chromosome_metadata['VCFfname'][chromosome]
        
        filename = os.path.join(self.db_dir,vcf_fname)
        if os.path.isfile(filename):
            counter= 0
            from pymodule.yhio.VCFFile import VCFFile
            
            vcfFile = VCFFile(inputFname=filename,minDepth=minDepth)
            
                
            #this is a list with the read-group names
            #readgroupIDList = vcfFile.getSampleIDList()
            #writer = csv.writer(open(self.outputFname, 'w'), delimiter='\t')
            #header = ['Chromosome', 'position', 'ref','alt']
            chrom_ls=[]; ref_ls=[]; snp_pos_ls=[]; alt_ls=[]; quality_ls=[];  info_ls=[]
            columnIndexList = map(int,vcf_indx_ls)
            datalist=[]
            for vcfRecord in vcfFile:
                #objects in VCFRecord: [vcfRecord.chr, vcfRecord.pos, vcfRecord.vcf_locus_id, vcfRecord.refBase, vcfRecord.altBase, \
#                    vcfRecord.quality, vcfRecord.filter,\
#                    vcfRecord.info, vcfRecord.format] + vcfRecord.row[self.sampleStartingColumn:]                
                if float(vcfRecord.quality)>=qualityfilter:
                    data_row=[]
                    chrom_ls.append(vcfRecord.chr)
                    snp_pos_ls.append(vcfRecord.pos)
                    refBase = vcfRecord.refBase
                    nonRefBase = vcfRecord.altBase
                    ref_ls.append(refBase)
                    alt_ls.append(nonRefBase)
                    quality_ls.append(vcfRecord.quality)
                    info_ls.append(vcfRecord.info)
                    for columnIndex in columnIndexList:
                        #for vcfCall in vcfRecord.data_row[1:]: #data_row is a list of dictionary {'GT': base-call, 'DP': depth}, or None if missing.
                        #it includes one extra sample in index 0, which is the reference individual (from the ref column of VCF).
                        vcfCall = vcfRecord.data_row[columnIndex+1]
                        if vcfCall:
                            if vcfCall['GT'][0]==refBase:
                                gt=0.
                            elif vcfCall['GT'][0]==nonRefBase:
                                gt=1.    
                            else:
                                #gt=float('nan')
                                gt=np.NaN
                            data_row.append(gt)
                            if vcfCall['GT'][1]==refBase:
                                gt=0.
                            elif vcfCall['GT'][1]==nonRefBase:
                                gt=1.    
                            else:
                                #gt=float('nan')
                                gt=np.NaN
                            data_row.append(gt)
                        else:
                            #data_row.append(float('nan'))
                            #data_row.append(float('nan'))
                            data_row.append(np.NaN)
                            data_row.append(np.NaN)
                    counter += 1
                    datalist.append(data_row)
            sys.stderr.write("%s: %s loci in %i individuals outputted.\n"%(chromosome,counter,len(columnIndexList)))
            #pdb.set_trace()
            data=np.array(datalist)
            
            #test wether all entries have the same chromosome
            #if chrom_ls.count(chrom_ls[0]) != len(chrom_ls):
            #    raise hvb.hsError("All SNPs should be on the same chromosome. Otherwise the output has to be modified.")
            
            datadict={'ucla_id':np.array(ucla_id_ls),'chromosome':np.array(chrom_ls),\
                                        'ref':np.array(ref_ls),'snp_pos':np.array(snp_pos_ls),\
                                        'alt':np.array(alt_ls),\
                                        'quality':np.array(map(float,quality_ls)),'info':np.array(info_ls),\
                                        'data':data,'genotype_method':genotype_method}
            return datadict        
    
    
    def create_genotype_df(self,datadict):
        """
        creates a genotype data frame
        with snp_pos as Index
        and 2 column levels, ucla_id and strand
        entries are 0/1 for ref/alt alleles
        """
        dd=datadict
        ucla_ids=[i for k in dd['ucla_id'] for i in [k]*2]
        strands=[0,1]*len(dd['ucla_id'])
        tuples=list(zip(ucla_ids,strands))
        cols=pd.MultiIndex.from_tuples(tuples)
        cols.names=['ucla_id','haplotype']
        df=pd.DataFrame(dd['data'],index=dd['snp_pos'],columns=cols)
        return df
    
    def get_genotype_df(self,chromosome,qualityfilter=0,minDepth=0):
        return self.create_genotype_df(self.extract_genotype_data(chromosome,qualityfilter=0,minDepth=0))
    
    def create_genotype_metadata(self,datadict):
        dd=datadict.copy()
        del dd['data']
        del dd['ucla_id']
        del dd['genotype_method']
        df=pd.DataFrame(dd)
        df.set_index('snp_pos',inplace=True)
        return df
    
    def add_ancestral_to_metadata(self,chromosome,genotype_metadata):
        fname=os.path.join(self.meta_dir,'align_outgroup/vervet_snps_macaque_{}.tsv'.format(chromosome))
        ancestral_df=pd.read_csv(fname,sep='\t',index_col=0,header=0)
        ancestral=ancestral_df['mac']
        ancestral[(ancestral_df['mac']!=ancestral_df['ref'])*( ancestral_df['mac']!=ancestral_df['alt'])]=np.NaN
        # correct for the fact that the snp_pos in the VCF_file is 1-indexed and not 0-indexed
        if np.all(ancestral.index+1==genotype_metadata.index):
            ancestral.index=ancestral.index+1
        elif np.any(ancestral.index!=genotype_metadata.index):
            raise hvb.hsError("Not the same SNP sets for ancestral state and metadata. Modify function to allow non-identical sets.")       
        genotype_metadata['ancestral']=ancestral
        return genotype_metadata
    
    def polarise_ancestral_genotype(self,genotype_df,genotype_metadata):
        #invert 0. and 1. for snps where the alternative state is the ancestral state
        genotype_df_pol=genotype_df.copy()
        genotype_df_pol[genotype_metadata['alt']==genotype_metadata['ancestral']]=(genotype_df[genotype_metadata['alt']==genotype_metadata['ancestral']]-1).abs()
        return genotype_df_pol
        
    

    def save_genotype_data(self,chromosome,qualityfilter=0,minDepth=0):
        
        meta_fname=os.path.join(self.meta_dir,"snp_metadata_m{}_{}.tsv".format(self.genotype_method,chromosome))
        hvb.try_make_dirs(os.path.dirname(meta_fname))
        geno_fname=os.path.join(self.data_dir,"hs_data/genotype_method_{}/".format(self.genotype_method),\
                                "genotype_mat_ancestral_m{}_{}.tsv".format(self.genotype_method,chromosome))
        hvb.try_make_dirs(os.path.dirname(geno_fname))
        
        dd=self.extract_genotype_data(chromosome,qualityfilter=qualityfilter,minDepth=minDepth)
        genotype_df=self.create_genotype_df(dd)
        metadata_df=self.add_ancestral_to_metadata(chromosome,self.create_genotype_metadata(dd))
        genotype_df=self.polarise_ancestral_genotype(genotype_df,metadata_df)
        
        metadata_df.to_csv(meta_fname,sep='\t')
        hvb.save_genotype(genotype_df,geno_fname)
        
        
    

    
    
    def parseVariantInfo(self,info_ls):
        """
        this function takes the np.array of the VCF info column and returns a dictionary of 
        VCF info lists, where each entry corresponds to a variant
        """
        dict_ls=[]
        for variant_entry in info_ls:
            variant_entry2=variant_entry.split(';')
            dict_creator=[]
            for entry in variant_entry2:
                entry2=entry.split("=")
                try:
                    entry2[1]=int(entry2[1])
                except:
                    try:
                        entry2[1]=float(entry2[1])
                    except:
                        pass
                dict_creator.append(tuple(entry2))
            dict_ls.append(dict(dict_creator))
#            dict_ls=map(lambda s: s.split("="),dict_ls)
        
        #convert list of dicts into dict of lists:
        #numpy array has problems with mixed data types...
        key_ls=dict_ls[0].keys()
        empt_array_ls=[np.empty(len(dict_ls),dtype=type(dict_ls[0][key])) for key in key_ls]
        for el in empt_array_ls:
            el[:]=np.NAN
#        empt_array[:]=np.NAN
        info_dic=dict([(key,empt_array) for key,empt_array in zip(key_ls,empt_array_ls)])
        for (indx,dic) in enumerate(dict_ls):
            for key in key_ls:
#                print key,indx,dic[key],info_dic[key][indx]
                try:
                    info_dic[key][indx]=dic[key]
                except:
                    pass        
        return info_dic
    
    
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
    print "main"
