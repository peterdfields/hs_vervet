#!/usr/bin/env python
"""
This script infers ancestral states of vervet snps from alignments 
of vervet reference against macaque reference. 
"""

import os,sys
import numpy as np
#needs pandas version >=0.12 to support regex replace
import pandas as pd



class AlignmentLine():
    def __init__(self,name,start,length,strand,sequence):
        self.name=name
        self.start=start
        self.length=length
        self.strand=strand
        self.sequence=sequence
    def __len__(self):
        return self.length

class Alignment(list):
    """
    """
    def __init__(self,score,alignment_lines):
        self.score=score
        for line in alignment_lines:
            self.append(line)
    def append(self,alignment_line):
        list.append(self,alignment_line)


class MultiAlign(file):
    def __init__(self,name,mode):
        file.__init__(self,name,mode)
    
    def _read_maf_entry(self):
        """
        read one entry (one line starting with "a" and several starting with "s") from a .maf file
        error if the first non-empty non-comment line does not start with an "a" 
        returns a tuple (entry score), where entry is [align0,align1,...]
        where align is a "s" line from a maf file without the initial s
        align0=[name,start_pos,bases,strand,total_sequence_len,aligned_sequence]
        see maf-file docu http://last.cbrc.jp/doc/last-tutorial.html
        """
        def try_to_int(s):
            try: 
                return int(s)
            except ValueError:
                return s
        line=file.next(self)
        while line==[] or line[0]=='#':
            file.next(self)    
        if line[0]=='a':
            score=line.split('score=')[1].strip()
        else:
            raise IOError('first non-empty line should start with "a", but is: {}'.format(line[0]))
        entries=[]
        line=file.next(self)
        while line[0]=='s':
            entries.append(map(lambda st: try_to_int(st.strip()),line.split()[1:]))
            line=file.next(self)
        return Alignment(int(score),[AlignmentLine(name=el[0],start=el[1],length=el[2],strand=el[3],sequence=el[5]) for el in entries])  
        
    def next(self):
        return self._read_maf_entry()

class HsError( Exception ): pass    


def create_align_df(align_ls,window_start,window_sz):
    """
     Creates a DataFrame which contains all the macaque sequences from the alignment list
     index: postition (integers) from window_start to window_start+window_sz
     columns: one column 'ref' for the reference sequence, then for each macaque alignment
      in this region one column labeled with (i,score) where i=0,1,...
    """
    df=pd.DataFrame(index=range(window_start,window_start+window_sz))
    for i,al in enumerate(align_ls):
        df0=pd.DataFrame({'ver':list(al[0].sequence),'mac':list(al[1].sequence)})
        df0=df0[df0['ver']!='-']
        df1=pd.DataFrame({(i,al.score):df0['mac'].values},index=range(al[0].start,al[0].start+al[0].length))
        df=df.merge(df1,how='left',left_index=True,right_index=True)
    return df


def polarise_snp(align_df,snp_df):
    pol_df=align_df.merge(snp_df,how='left',left_index=True,right_index=True)
    pol_df=pol_df[pol_df['ref'].notnull()]
    align_df2=align_df.ix[pol_df.index]
    ancestral=decide_macaque_base(align_df2)
    return ancestral
        


def decide_macaque_base(align_df):
    #test if there are no snps in the region
    if not align_df:
        return pd.Series()
    bases=np.array([np.nan,'A','C','G','T'],dtype=object)
    #arrays of len (pol_df)=n_snp giving for each SNP the scores of A,C,G,T
    score_nan=np.zeros(align_df.shape[0])
    score_A=(align_df.replace(['[Aa]','[^A^a]'],[1,0],regex=True).fillna(0).values*np.array(map(lambda x: int(x[1]),align_df.columns.values))).sum(axis=1)
    score_C=(align_df.replace(['[Cc]','[^C^c]'],[1,0],regex=True).fillna(0).values*np.array(map(lambda x: int(x[1]),align_df.columns.values))).sum(axis=1)
    score_G=(align_df.replace(['[Gg]','[^G^g]'],[1,0],regex=True).fillna(0).values*np.array(map(lambda x: int(x[1]),align_df.columns.values))).sum(axis=1)
    score_T=(align_df.replace(['[Tt]','[^T^t]'],[1,0],regex=True).fillna(0).values*np.array(map(lambda x: int(x[1]),align_df.columns.values))).sum(axis=1)
    score_mat=np.array([score_nan,score_A,score_C,score_G,score_T]).T
    #the base with the maximum sum of scores is chosen as true base for each position:
    return pd.Series(bases[score_mat.argmax(axis=1)],index=align_df.index)


def get_ancestral(align_fname,snp_df,window_sz=10000,chromosome_length=False,verbose=False):
    """
    Gives an ancestral state for all the snps on a given chromosome.
    input...
    top aligned sequence must be the one from which the snp list originates (vervet).
    second aligned sequence must be the one used to infer the ancestral state (e.g. macaque).
    the '.maf' sequence file must be sorted with respect to positions in the top aligned sequence
    
    """
    window_last=None
    with MultiAlign(align_fname,'r') as aligns:
        align_ls=[]
        ancestral=pd.Series()
        for align in aligns:
            if align[0].strand!='+':
                raise HsError('Strand is {}. Alignment of {} not on forward strand, not implemented.'.format(align[0].strand,align[0].name))
            if align_ls and align[0].start<align_ls[-1][0].start:
                raise HsError('{} alignment must be sorted. Alignment starting at {} after one starting at {}.'.\
                                    format(align[0].name,align[0].start,align_ls[-1][0].start))
            #integer division, gives the current window
            window_current=align[0].start/window_sz
            #check if alignments are in the same window
            if window_last!=None:
                while window_current>window_last:
                    if chromosome_length and verbose:
                        n_window=chromosome_length/window_sz
                        if round(window_last*10./n_window)/10. != round((window_last+1)*10./n_window)/10.:
                            print "window", window_last,',',round(window_last*10./n_window)/10.*100, "% polarised"
                    if snp_df.ix[window_last*window_sz:(window_last+1)*window_sz]:        
                        align_df=create_align_df(align_ls,window_last*window_sz,window_sz)
                        ancestral=ancestral.append(polarise_snp(align_df,snp_df))
                    #keep alignments that overlap the next window
                    align_ls=[al for al in align_ls if (al[0].start+al[0].length)>= (window_last+1)*window_sz]
                    window_last+=1
            align_ls.append(align)
            window_last=window_current
        if align_ls:
            #add the final alignments
            align_df=create_align_df(align_ls,window_current*window_sz,window_sz)
            ancestral=ancestral.append(polarise_snp(align_df,snp_df))
    return ancestral

if __name__ == '__main__':
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
    #from hs_vervet.tools import hs_vervet_basics as hvb
    #from hs_vervet.tools import stats_on_data as sod
    import vcf
    import argparse
    
    parser=argparse.ArgumentParser(description="Return macaque state for each SNP. Output  \
                                                written to outfile!")
    parser.add_argument("alignment_file",help="Alignment file in .maf format,\
                                         sorted by vervet genome position with vervet as top sequence.")
    parser.add_argument("vcf_file",help="File with the SNP calls for the same chromosome as the\
                                        alignment file")
    parser.add_argument("out_file",help="filename to which the output is written")
    parser.add_argument("--outgroup_name",default="outgroup",help="name of the other species that vervet is aligned to")
    parser.add_argument("-v","--verbose",action="store_true")
    parser.add_argument("-w", "--window", default=10000, type=int,
                    help="change the size of the window for which ancestral state inference is done \
                    simultaneously (expert use)")
    args=parser.parse_args()
    if args.verbose:
        import time
        start=time.clock()

    out_file=open(args.out_file,'w')
    vcf = vcf.Reader(open(args.vcf_file,'r'))

    if args.verbose:
        print args.vcf_file, "loaded"
    snp_pos=[]
    ref=[]
    alt=[]
    for record in vcf:
        #if record.FILTER:
        #    continue
        if not record.is_snp:
            continue
        snp_pos.append(record.POS)
        ref.append(record.REF)
        alt.append(record.ALT[0])

    snp_df=pd.DataFrame({'alt':alt,'ref':ref},index=np.array(snp_pos)-1)
    
    if args.verbose:
        print len(snp_df), " SNPs read, time elapsed:", time.clock()-start
    ancestral=get_ancestral(args.alignment_file,snp_df,window_sz=args.window,\
                            chromosome_length=snp_df.index.values[-1],verbose=args.verbose)
    ancestral.name=args.outgroup_name
    
    if args.verbose:
        print "polarising finished. Saving to file. Time elapsed:", time.clock()-start  
    out_df=snp_df.join(ancestral,how='outer')
    out_df.to_csv(out_file,sep='\t',index_label='position')
    out_file.close()
    if args.verbose:
        print time.clock()-start, "seconds elapsed"
    
