#!/usr/bin/env python
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))





def vcf_to_tsv(in_vcf,out_fn,interval=None,dont_randomise_unphased=False,include_triallelic=False,keep_original_id=False,include_filtered=False,encode_ancestral_derived=False,verbose=False):
    """
    Reads a vcf and outputs a tsv where each row is one SNP and each column is one allele of an individual.
    The first two rows and columns are labels, chrom pos, and ucla_id first/second allele, respectively.
    First argument should be a file-stream or file.
    """ 
    import vcf
    import random
    
    def try_fetch(reader,*args):
        try:
            return reader.fetch(*args)
        except IOError:
            print "For interval functionality, the vcf needs to be tabix indexed. Use bgzip+tabix."
            raise

    
    #open if in_vcf is filename 
    try:
        open(in_vcf,'r')
    except TypeError, e:
        if 'file found' in str(e):
            pass
        else:
            raise

    reader = vcf.Reader(in_vcf)
    
    if encode_ancestral_derived:
        if 'AA' not in reader.infos.keys():
            raise Exception('For encoding ancestral/derived, there must be an info tag "AA" for the alternative allele.')

    if interval is not None:
        interval1  = interval.split(':')
        chrom = interval1[0]
        if len(interval1)==1:
            reader = try_fetch(reader,chrom)
        elif len(interval1)==2:
            start_end = interval1[1].split('-')
            start = start_end[0]
            if len(start_end) == 1:
                reader = try_fetch(reader,chrom,start)
            elif len(start_end) == 2:
                end = start_end[1]
                reader = try_fetch(reader,chrom,start,end)
            else:
                raise Exception("Interval should be of the form chrom[:start[-end]] but is {}".format(interval))
        else:
            raise Exception("Interval should be of the form chrom[:start[-end]] but is {}".format(interval))
            
    

    ids_original = reader.samples
    if not keep_original_id:
        id_dict={'pygerythrus':'AG5417','cynosurus':'VZA3008','aethiops':'VEC1016','tantalus':'AGM141'}
        ids = [(h.split('_')[2] if len(h.split('_'))>2 else h) for h in ids_original]
        for i,id0 in enumerate(ids):
            try:
                ids[i] = id_dict[id0] 
            except KeyError:
                pass

    #duplicate ids; first two columns reserved for index
    header0 =['',''] + [i for k in ids for i in [k]*2]
    #indicates first and second allele for each locus
    header1 =['',''] + ['0','1']*len(ids)
    header2 = ['chrom','position'] + ['']*len(ids)
    with open(out_fn,'w') as f:
        f.write('\t'.join(header0)+'\n')
        f.write('\t'.join(header1)+'\n')
        f.write('\t'.join(header2)+'\n')
        for record in reader:
            if dont_randomise_unphased:
                het0 = 0
                het1 = 1
            else:
                het0 = random.randint(0,1)
                het1 = abs(het0-1)
            genotype_dict={None:['N','N'],'0/0':[0,0],'0/1':[het0,het1],'1/1':[1,1],'0|0':[0,0],'0|1':[0,1],'1|0':[1,0],'1|1':[1,1]}
            if not include_triallelic:
                if len(record.ALT)>1:
                    if verbose:
                        print record.CHROM,record.POS,"skipped because more than 2 alleles."
                    continue
            if not include_filtered:
                if record.FILTER:
                    if verbose:
                        print record.CHROM,record.POS,"skipped because of filters {}.".format(record.FILTER)
                    continue
            if encode_ancestral_derived:
                try:
                    aa = record.INFO['AA']
                except KeyError:
                    if verbose:
                        print record.CHROM,record.POS,"skipped because no ancestral state info."
                    continue
                if (aa != record.REF) and (aa != record.ALT[0].sequence):
                    if verbose:
                        print record.CHROM,record.POS,"skipped because of ancestral state={0} neither ref={1} or alt={2}".format(aa,record.REF,record.ALT[0])
                    continue
                if  aa == record.ALT[0].sequence:                    
                    genotype_dict={None:['N','N'],'0/0':[1,1],'0/1':[het0,het1],'1/1':[1,1],'0|0':[1,1],'0|1':[1,0],'1|0':[0,1],'1|1':[0,0]}
                elif aa == record.REF:
                    pass
                else:
                    raise Exception('Alternative allele not recognised.')
            line = []
            line.append(record.CHROM)
            line.append(str(record.POS))
            genotypes_row = [str(k) for id1 in ids_original for k in genotype_dict[record.genotype(id1)['GT']]]
            line = line + genotypes_row
            f.write('\t'.join(line)+'\n')
            if verbose:
                print record.CHROM,record.POS, 'written to output'
        


def vcf_to_df(in_vcf,include_triallelic=False,keep_original_id=False):
    pass


def load(fn,**kwargs):
    import pandas as pd
    df = pd.read_csv(fn,sep='\t',index_col=[0,1],header=[0,1],tupleize_cols=False,na_values='N',**kwargs)
    return df    

def save(df,fn,**kwargs):
    df.to_csv(fn,sep='\t',index_label=['chrom','position'],na_rep='N',float_format='%i',tupleize_cols=False,**kwargs)





if __name__ == '__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Parse VCF into a 0/1 (pseudo-)haplotype tsv file. (If not phased, the alleles are randomly assigned to haplotypes in heterozygotes.)")
    #the following allows piping the vcf as input stream
    parser.add_argument('in_vcf', type = argparse.FileType('r'), default = '-',help="VCF filename to parse.")
    parser.add_argument("out_fn",help="Filename of tsv to output.")
    parser.add_argument("-L","--interval",default=None,help="Target intervals of the form: 'chrom:[start-end]'")
    parser.add_argument("--dont-randomise-unphased",action='store_true',help='By default, unphased hets are are randomly assigned 0 1 or 1 0.')
    parser.add_argument("--include-triallelic",action="store_true")
    parser.add_argument("--keep-original-id",action="store_true",help='If not specified, the script tries to convert read-group id into ucla_id')
    parser.add_argument("--include-filtered",action="store_true",help='By default SNPs which do not have PASS in the filter column are excluded.')
    parser.add_argument("--encode-ancestral-derived",action="store_true",help='0/1 encodes ancestral/derived instead of ref/alt. Option only uses SNPs that have an AA tag in the INFO column. SNPs where AA is a third allele are excluded.')
    parser.add_argument("-v","--verbose",action="store_true",help="Verbose mode.")
    
    arg_dict = vars(parser.parse_args())
    vcf_to_tsv(**arg_dict)
