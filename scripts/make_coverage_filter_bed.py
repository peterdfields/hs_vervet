#!/usr/bin/env python
info = """
Uses the output of GATK DepthOfCoverage as input.
Writes bed files that mark low and high coverage percentiles.
Output files are of the form
coverage_<top/low>_<cutoff_percent>pct_<chrom>.bed
Assumes that there is only one chromosome per file!!
Not true:Assumes that coordinates are consecutive (increment of 1).
"""
import numpy as np

def get_start_end_pos(zero_one_vec,positions):
    """
    find intervals of consecutive ones in zero_one_vec
    return the the corresponding coordinates at the start end end positions
    of these intervals
    the vector <positions> gives the coordinate corresponding to the indices 
    """
    zero_one_vec = zero_one_vec*1
    first_pos = (1 if zero_one_vec[0]==1 else 0)
    last_pos = (-1 if zero_one_vec[1]==1 else 0)
    interval_start_end = np.append(np.insert(np.diff(zero_one_vec),0,first_pos),
                                                                        last_pos)
    start = np.where(interval_start_end==1)[0]
    end = np.where(interval_start_end==-1)[0]
    start_pos = position[start]
    end_pos = position[end]
    return start_pos, end_pos

def make_bed(fname,chrom,start,end,info=None):
    """
    write a bedfile
    fname ... filename to write to
    start, end ... are arrays of equal length containing 
            the start and end coordinates to write
    chrom ... string or list of strings of chromosome names
            if single string it is assumed that all entries 
            have the same chromosome
    """
    #if type(chrom) == str:
    #    chrom = [chrom]*len(start)
    #if type(info) == str:
    #    info = [info]*len(start)
    with open(fname,"w") as f:
        if info is None:
            for s,e in zip(start,end):
                line = [chrom,str(s),str(e)]
                f.write("\t".join(line)+"\n")
        else:
            for s,e in zip(start,end):
                line = [chrom,str(s),str(e),info]
                f.write("\t".join(line)+"\n")



if __name__ == '__main__':
    import os, sys
    import argparse

    parser = argparse.ArgumentParser(description=info)
    subparsers = parser.add_subparsers(dest='mode',help='Coverage can be filtered by \n'
                                                    'percentile or by fold higher/lower '
                                                                   'than median.')
    parser_a = subparsers.add_parser('percentile', help='Cutoff by percentile.')
    parser_a.add_argument("--excess_cutoff_percent",type=float,default=None)
    parser_a.add_argument("--low_cutoff_percent",type=float,default=None)
   
    parser_b = subparsers.add_parser('fold', help='Cutoff by fold higher/lower '
                                                  'than median.')
    parser_b.add_argument("--low_fold_cutoff",type=float,default=None)
    parser_b.add_argument("--excess_fold_cutoff",type=float,default=None)



    parser.add_argument("in_fn", type = argparse.FileType('r'),
                                        help="Input depth_of_coverage filename.")
    parser.add_argument("out_folder", help="Output_folder.")
    parser.add_argument("--chrom",default="",
                        help="Chromosome name to add to filename\n"
                             "and inside bed file.")

    args = parser.parse_args()

    assert os.path.isdir(args.out_folder)

    #print args
    #sys.exit(0)
    
    if args.mode == "percentile":
        excess_cutoff_percent = args.excess_cutoff_percent
        low_cutoff_percent = args.low_cutoff_percent
        assert (excess_cutoff_percent is not None) or (low_cutoff_percent is not None)
        excess_out_fn = os.path.join(args.out_folder,"coverage_top_"+\
                                   str(excess_cutoff_percent)+"pct_"+\
                                        args.chrom+".bed")

        low_out_fn = os.path.join(args.out_folder,"coverage_low_"+\
                              str(low_cutoff_percent)+"pct_"+\
                                    args.chrom+".bed")
    elif args.mode == "fold":
        assert (args.low_fold_cutoff is not None) or (args.excess_fold_cutoff is not None)
        excess_out_fn = os.path.join(args.out_folder,"coverage_above_"+\
                                   str(args.excess_fold_cutoff)+"fold_"+\
                                        args.chrom+".bed")

        low_out_fn = os.path.join(args.out_folder,"coverage_below_"+\
                              str(args.low_fold_cutoff)+"fold_"+\
                                    args.chrom+".bed")
        


    #with open(coverage_fn,'r') as cf:
    depth = []
    position = []
    chrom = args.chrom
    args.in_fn.readline()
    for line in args.in_fn:
        depth.append(int(line.split()[1]))
        coord = line.split()[0].split(":")
        #chrom.append(coord[0])
        position.append(int(coord[1]))
        #if position[-1]>1010000:
        #    break
    position = np.array(position)
    position = np.append(position,position[-1]+1)
    depth = np.array(depth)

    if args.mode == "percentile":

        if excess_cutoff_percent is not None:
            excess_cutoff = np.percentile(depth,100-excess_cutoff_percent)
            #print excess_cutoff
            excess_start, excess_end = get_start_end_pos(depth>excess_cutoff,position)
            make_bed(excess_out_fn,chrom,excess_start,excess_end,"COVgt"+str(round(excess_cutoff*10)/10))

        del excess_start
        del excess_end

        if low_cutoff_percent is not None:
            low_cutoff = np.percentile(depth,low_cutoff_percent)
            #print low_cutoff
            low_start, low_end = get_start_end_pos(depth<low_cutoff,position)
            make_bed(low_out_fn,chrom,low_start,low_end,"COVlt"+str(round(low_cutoff*10)/10))

    elif args.mode == "fold":
        median = np.median(depth)
        if args.excess_fold_cutoff is not None:
            excess_cutoff = median * args.excess_fold_cutoff
            excess_start, excess_end = get_start_end_pos(depth>excess_cutoff,position)
            make_bed(excess_out_fn,chrom,excess_start,excess_end,"COVgt"+str(round(excess_cutoff*10)/10))

        del excess_start
        del excess_end

        if args.low_fold_cutoff is not None:
            low_cutoff = median * 1. / args.low_fold_cutoff
            #print low_cutoff
            low_start, low_end = get_start_end_pos(depth<low_cutoff,position)
            make_bed(low_out_fn,chrom,low_start,low_end,"COVlt"+str(round(low_cutoff*10)/10))
