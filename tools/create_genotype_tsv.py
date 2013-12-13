#!/usr/bin/env python
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


genotype_dict={None:['N','N'],'0/0':[0,0],'0/1':[0,1],'1/0':[1,0],'1/1':[1,1],'0|0':[0,0],'0|1':[0,1],'1|0':[1,0],'1|1':[1,1]}


def vcf_to_tsv(in_vcf,out_fn,interval=None,include_triallelic=False,keep_original_id=False,include_filtered=False):
    """
    First argument should be a file-stream or file
    """ 
    import vcf
    
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
    
    if interval is not None:
        interval1  = interval.split(':')
        chrom = interval1[0]
        if len(interval1)==1:
            reader = try_fetch(reader,chrom)
        elif len(interval1)==2:
            start_end = interval1[1].split('-')
            start = start_end[0]
            if len(start_end)==1:
                reader = try_fetch(reader,chrom,start)
            elif len(start_end)==2:
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
            if not include_triallelic:
                if len(record.ALT)>1:
                    continue
            if not include_filtered:
                if record.FILTER:
                    continue
            line = []
            line.append(record.CHROM)
            line.append(str(record.POS))
            genotypes_row = [str(k) for id1 in ids_original for k in genotype_dict[record.genotype(id1)['GT']]]
            line = line + genotypes_row
            f.write('\t'.join(line)+'\n')
        


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

    parser=argparse.ArgumentParser(description="Parse VCF into a 0/1 genotype tsv file.")
    #the following allows piping the vcf as input stream
    parser.add_argument('in_vcf', type = argparse.FileType('r'), default = '-',help="VCF filename to parse.")
    parser.add_argument("out_fn",help="Filename of tsv to output.")
    parser.add_argument("-L","--interval",default=None,help="Target intervals of the form: 'chrom:[start-end]'")
    parser.add_argument("--include-triallelic",action="store_true")
    parser.add_argument("--keep-original-id",action="store_true",help='If not specified, the script tries to convert read-group id into ucla_id')
    parser.add_argument("--include-filtered",action="store_true",help='By default SNPs which do not have PASS in the filter column are excluded.')
    #parser.add_argument("input_vcf",help="VCF filename to parse.")
    arg_dict = vars(parser.parse_args())
    vcf_to_tsv(**arg_dict)
