#!/usr/bin/env python
"""
This script infers macaque states for vervet reference from alignments 
of vervet reference against macaque reference. 
"""



if __name__ == '__main__':
    import os,sys
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
    import numpy as np
    import pandas as pd
    from hs_vervet.scripts import vervet_snp_other_species_state as vsoss
    from hs_vervet.tools import hs_vervet_basics as hvb
    from pyfasta import Fasta
    import argparse
    
    parser=argparse.ArgumentParser(description="Return differences between one reference (e.g. vervet)\
                                                and another species (macaque). Output \
                                                saved to file!")
    parser.add_argument("alignment_file",help="Alignment file in .maf format,\
                                         sorted by vervet genome position with vervet as top sequence.")
    parser.add_argument("ref_file",help="Reference file. This file must contain 1 Output coordinates will be in the \
                                        coordinates of this species")
    parser.add_argument("out_file",help="filename to which the output is written")
    parser.add_argument("chrom_name",help="name of the chromosome on which the analysis runs\
                                    as given in header starting with '>' the reference file, e.g. CAE1")
    parser.add_argument("--outgroup_name",default="outgroup",help="name of the other species that vervet is aligned to")
    parser.add_argument("-v","--verbose",action="store_true")
    parser.add_argument("-w", "--window", default=10000, type=int,
                    help="change the size of the window for which ancestral state inference is done \
                    simultaneously (expert use)")
    args=parser.parse_args()
    
    if args.verbose:
        import time
        start=time.clock()

    hvb.try_make_dirs(os.path.dirname(args.out_file))
    out_file=open(args.out_file,'w')
    
    ref=Fasta(args.ref_file)[args.chrom_name]
    
    if args.verbose:
        print args.ref_file, "loaded"
    
    #these are not just snps, but the whole chromosome
    snp_df=pd.DataFrame({'ref':[c for c in ref]})
   
    if args.verbose:
        print "loading reference finished, time elapsed:", time.clock()-start
    
    other_species=vsoss.get_ancestral(args.alignment_file,snp_df,window_sz=args.window,\
                            chromosome_length=snp_df.index.values[-1],verbose=args.verbose)
    other_species.name=args.outgroup_name
    
    if args.verbose:
        print "creating consensus alignment for "+args.chrom_name+" finished, time elapsed:", time.clock()-start
        
    align_df=snp_df.join(other_species,how='outer')
    #replace 'N' in the vervet reference with NaNs 
    align_df['ref'][align_df['ref']=='N']=np.nan
    align_df_no_nan=align_df.dropna()    
    ver_mac_diff=align_df_no_nan[align_df_no_nan['ref']!=align_df_no_nan[args.outgroup_name]]    
    
    ver_mac_diff.rename(columns=lambda x: 'vervet' if x=='ref' else x,inplace=True)
    
    ver_mac_diff.to_csv(out_file,sep='\t',index_label='position')
        
    out_file.close()
    if args.verbose:
        print time.clock()-start, "seconds elapsed"