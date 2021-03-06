#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.
Todo:
-- parallel
-- pack into classes
"""
#import sys, os
import warnings

#--------support functions------------
def get_012(gt_str):
    if gt_str[:3] in ["0/0","0|0"]:
        return "0"
    elif gt_str[:3] in ["1/1","1|1"]:
        return "2"
    elif gt_str[:3] in ["0/1","0|1","1|0"]:
        return "1"
    elif gt_str[:3] == "./.":
        return "N"
    else:
        raise ValueError("Unsupported genotype " + gt_str)


def get_AA(info_str):
    aa = info_str.split("AA=")[1][0]    
    return aa


def generic_parser(parse_fun,header_fun,vcf_fh, tsv_fh,no_skip_multiple_entries,*args,**kwa):
    """
    parse_fun... applied to each split data line of the vcf
    header_fun... applied to each line in the vcf header
    """
    iwarn = 0
    max_warn = 50
    prev_chrom = None
    prev_pos = -1
    h = None
    p = None
    for line in vcf_fh:
        if line[0] == "#":
            #try:
            h = header_fun(line, tsv_fh, h,**kwa)
            #except Exception:
                
            continue
        d = line.strip().split("\t")
        chrom = d[0]
        pos = int(d[1])
        if chrom == prev_chrom:
            assert pos >= prev_pos, "vcf positions not in "\
                                    "ascending order at: {}:{},{}".format(chrom,prev_pos,pos) 
            if pos == prev_pos:
                if no_skip_multiple_entries:
                    warnings.warn("Warning, multiple entries for pos {}:{}.\n"
                              "Keeping all entries.".format(chrom,pos))
                else:
                    if iwarn < max_warn:
                        warnings.warn("Warning, multiple entries for pos {}:{}.\n"
                                  "Skipping all but the first.".format(chrom,pos))
                        iwarn += 1
                    continue
        p = parse_fun(d,tsv_fh, p, h, *args,**kwa)
        prev_chrom = chrom
        prev_pos = pos

#--------parsing functions-------------

def vcf_to_012_header(line, tsv_fh,*args):
    if line[1:6] == "CHROM":
        tsv_fh.write("chrom\tpos\t"+line.split("\t",9)[-1])

def vcf_to_012_parse_fun(d, tsv_fh,*args):
    gt = map(get_012,d[9:])
    tsv_fh.write("\t".join(d[:2])+"\t"+"\t".join(gt)+"\n")


#-----------

def ref_alt_anc_header(line, tsv_fh,*args):
    if line[1:6] == "CHROM":
        tsv_fh.write("chrom\tpos\tref\talt\tanc\n")

def ref_alt_anc_parse_fun(d, tsv_fh,*args):
    try:
        aa = get_AA(d[7])
    except IndexError:
        aa = "N"
    tsv_fh.write("\t".join(d[:2])+"\t"+"\t".join(d[3:5])+"\t"+aa+"\n")

#-----------
import itertools

def msmc_header(line, tsv_fh,*args,**kw):
    if line[1:6] == "CHROM":
        return line[1:].strip().split()

#def get_alleles(allele_dic,gt_str):
#    allele0 = allele_dic[gt_str[0]]
#    allele1 = allele_dic[gt_str[2]]
#    phased = (True if gt_str[1]=="|" else False)
#    return (allele0, allele1, phased)


def get_alleles(gt_str,ref,alts):
    """
    alts is a list obtained by alt.split(',')
    """
    def get_allele(num_gt):
        if num_gt == '0':
            return ref
        elif num_gt == '.':
            return '?'
        else:
            try:
                allele = alts[int(num_gt)-1]
            except IndexError, e:
                print ref, alts, gt_str
                raise e
            return allele
    gts = gt_str.split(':')[0]
    if '|' in gts:
        phased = True
        gts_ls = gts.split('|')
    elif '/' in gts:
        phased = False
        gts_ls = gts.split('/')
    else:
        raise TypeError('Unsupported genotype {}'.format(gt_str))

    #phased = (True if gt_str[1]=="|" else False)
    try:
        a0 = get_allele(gts_ls[0])
    except ValueError,e:
        print("Problem with gts_ls[0]",gts_ls[0])
        raise e
    try:
        a1 = get_allele(gts_ls[1])
    except ValueError, e:
        print("Problem with gts_ls[1]",gts_ls[1])
        raise e
    return(a0, a1, phased)


def msmc_input_parse_fun(d, tsv_fhs, p, h, ind_tuples,haplotypes):
    """
    *attention, input must be a VCF that includes all sites,
     segregating and non-segregating, filtered and non-filtered
    *attention, tsv_fh can be a list of multiple file handles
    here, one per output file to be written

    d ... split vcf line
    tsv_fhs ... list of output filehandles
    p ... return value of parse fun, counts the number of informative sites for each 
        comparison
    h ... header line
    ind_tuples ...  list of tuples of individual ids, same length as tsv_fhs
    haplotypes ... which haplotypes to take
    """
    alts = d[4].split(',')


    #allele_dic = {"0":d[3],"1":alt[0],".":"?",'2':alt[1],'3':alt[]}

    try:
        len(tsv_fhs)
    except TypeError:
        tsv_fhs = [tsv_fhs]
    assert len(ind_tuples) == len(tsv_fhs)
    if p is None:
        p = [0]*len(ind_tuples)
    #check for filter
    #if d[6] or len(d[3])>1 or len(d[4])>1:
    #    print "filtered", d[6]
    #    return p
    #else:
    #    p += 1

    #if pass and not indel or multiallelic
    if d[6] in ["PASS","",".","LowQual"]:
        for j,(tsv_fh,ids) in enumerate(zip(tsv_fhs,ind_tuples)):
            try:
                genotypes = [get_alleles(d[h.index(id)],d[3],alts) for id in ids]
            except ValueError, e:
                print(d)
                raise(e)
            phases = [g[2] for g in genotypes]
            if haplotypes==0:
                genotype_strs = [g[0] for g in genotypes]
            elif haplotypes==1:
                genotype_strs = [g[1] for g in genotypes]
            elif haplotypes==2:
                genotype_strs = [g[0]+g[1] for g in genotypes]
            gts = []
            for gt, phase in zip(genotype_strs,phases):
               if phase or ('?' in gt) or (gt == gt[0]*len(gt)):
                    gts.append([gt])
               else:
                    gts.append([gt,gt[::-1]])

            gt_len = 2 if haplotypes==2 else 1
            skip = False
            for gt in gts:
                #skip if missing data or if more than biallelic
                #for some reason, gt is a list (with one (?) string element)
                for g in gt:
                    if len(g) > gt_len or \
                                   '?' in g:
                        skip = True
                        break
                   # else:
                    #    print('noskip:',gt)
            if not skip:
                try:
                    gt_str = ",".join(["".join(i) for i in itertools.product(*gts)])
                except TypeError, e:
                    print(d[1],d[3],d[4], genotypes, gts)
                    raise e
                    
                #if ids == ['14319', '18694', '15560', '18696']:
                #    print(d[1],ids,gts,gt_str)
                #if "?" in gt_str: #don't print site with missing gt
                #    pass
                #elif gt_str == len(gt_str) * gt_str[0]: #don't print if all alleles equal
                p[j] += 1 #count site as informative
                if gt_str != len(gt_str) * gt_str[0]: #print if not all alleles equal
                    #print(d[1],genotypes,gts,gt_str)
                    tsv_fh.write(d[0]+"\t"+d[1]+"\t"+str(p[j])+"\t"+gt_str+"\n")
                    #why is the following check needed?
                    #if len(gt_str) > 2:
                    #    print(d)
                    #    print gt_str
                    #    raise Exception('bla')
                    p[j] = 0
    return p


def msmc_input_parse_fun_old(d, tsv_fhs, p, h, ind_tuples,haplotypes):
    """
    *attention, input must be a VCF that includes all sites,
     segregating and non-segregating, filtered and non-filtered
    *attention, tsv_fh can be a list of multiple file handles
    here, one per output file to be written

    d ... split vcf line
    tsv_fhs ... list of output filehandles
    p ... return value of parse fun, counts the number of informative sites for each.
        comparison
    h ... header line
    ind_tuples ...  list of tuples of individual ids, same length as tsv_fhs
    haplotypes ... which haplotypes to take
    """
    try:
        len(tsv_fhs)
    except TypeError:
        tsv_fhs = [tsv_fhs]
    assert len(ind_tuples) == len(tsv_fhs)
    if p is None:
        p = [0]*len(ind_tuples)
    #check for filter
    #if d[6] or len(d[3])>1 or len(d[4])>1:
    #    print "filtered", d[6]
    #    return p
    #else:
    #    p += 1
    if d[4] in ['A','C','G','T']: #if SNP
        if d[6] in ["PASS","","."] and len(d[3])==1: #check whether not filtered or deletion
            allele_dic = {"0":d[3],"1":d[4],".":"?"}
            for j,(tsv_fh,ids) in enumerate(zip(tsv_fhs,ind_tuples)):
                genotypes = [get_alleles(allele_dic,d[h.index(id)]) for id in ids]
                phases = [g[2] for g in genotypes]
                if haplotypes==0:
                    genotype_strs = [g[0] for g in genotypes]
                elif haplotypes==1:
                    genotype_strs = [g[1] for g in genotypes]
                elif haplotypes==2:
                    genotype_strs = [g[0]+g[1] for g in genotypes]
                gts = []
                for gt, phase in zip(genotype_strs,phases):
                   if phase or ('?' in gt) or (gt == gt[0]*len(gt)):
                        gts.append([gt])
                   else:
                        gts.append([gt,gt[::-1]])
                gt_str = ",".join(["".join(i) for i in itertools.product(*gts)])
                if "?" in gt_str: #don't print site with missing gt
                    pass
                elif gt_str == len(gt_str) * gt_str[0]: #don't print if all alleles equal
                    p[j] += 1
                else:
                    p[j] += 1
                    tsv_fh.write(d[0]+"\t"+d[1]+"\t"+str(p[j])+"\t"+gt_str+"\n")
                    p[j] = 0
            return p
        else:
            return p
    else:
        if d[6] in ["PASS","",".","LowQual"] and len(d[3])==1 and (len(d[4])==1 or d[4]=="None"):
            p = [i+1 for i in p]
        return p



#--------------


def add_outgroup_header(line, out_vcf_fh,*args, **kwa):
    #if line[:6]=='#CHROM':
    #    line = line.strip() + "\t" + sample_name + '\n'
    out_vcf_fh.write(line)



def add_outgroup_parse_fun(line,out_vcf_fh, *args,**kwa):
    bases = ["A",'C','T','G']
    no_var = ['.','X']
    og_base = og_fasta[line[0]][int(line[1])-1]
    if line[4] in no_var and og_base.upper() in bases and og_base.upper() != line[3]:
        line[4] = og_base.upper()
        #print "added SNP at", line[0], line[1]
    #print out_vcf_fh
    #print "\t".join(line)
    out_vcf_fh.write("\t".join(line)+'\n')

if __name__ == "__main__":
    import argparse
    import gzip

    parser = argparse.ArgumentParser(description="Extract data from vcf file.")
    parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-', help="Input vcf filename.")
    parser.add_argument("out_fn",type = argparse.FileType('w'), default = '-', help="Output filename")
    parser.add_argument("--no_skip_multiple_entries",action="store_true",help="Do not skip multiple entries for the same "
                                                                              "position can be the case for indels")
    #parser.add_argument("--ancestral_derived",action = 'store_true', help=  "0,1,2 corresponds to ancestral vs\n"
    #                                                                        "derived instead of ref vs alt.\n"
    #                                                                        "Using AA tag from VCF info.\n"
    #                                                                        "Sites for which this tag is "
    #                                                                        "unavailable or neither ref \n"
    #                                                                        "nor alt are ignored.")
    subparsers = parser.add_subparsers(dest='mode',
                                        help='E.g. "to012", or "ref_alt_anc"')
    parser0 = subparsers.add_parser("to012",help="Convert a vcf into a tsv with 0,1,2 genotypes.\n"
                                                   "0..reference,1..het,2..alt.")

    parser1 = subparsers.add_parser("ref_alt_anc",help="Extract tsv with chrom, pos, "
                                                       "ref_allele, alt_allele, anc_allele. "
                                                       "Requires annotation AA in info column")

    parser2 = subparsers.add_parser("msmc",help="Create msmc input files from VCF which includes "
                                                                        "fixed and filtered sites.")
    parser2.add_argument("--add_out_fns",nargs="+",type = argparse.FileType('w'),
                                                help="List of additional output files used with"
                                                                    " second to last ind_tuples")
    parser2.add_argument("--ind_tuples",required=True,nargs='+',help="List of individual ids that should be "
                                                            "outputted into the output files."
                                                        "Each n-tuple of ids for the same file should be "
                                                        "quoted in the command line, e.g., 'id1 id2 id3' 'id2 id4') "
                                                        "produces two output files.")
    parser2.add_argument("--haplotypes",type=int,default=2,help="Which haplotypes to use:0..first,1..second,2..both")

    parser3 = subparsers.add_parser("add_outgroup_from_fasta", help="Add outgroup difference from "
                                                                    "outgroup fasta (in vcf ref coord) as SNPs to VCF"
                                                                    ". The genotype itself is not added "
                                                                    " Attention: might violate "
                                                                    "VCF specs. "
                                                                    "Read function docu for "
                                                                    "more info.")

    parser3.add_argument("--in_fasta",help="Differences to be added as SNPs into VCF.")

    #parser3.add_argument("--sample_name",help="Name of the sample to be added.")

    args = parser.parse_args()
    #for k,v in vars(args).iteritems():
    #    print k,v,
    if args.mode == "to012":
        generic_parser(vcf_to_012_parse_fun,vcf_to_012_header,
                                        args.in_vcf, args.out_fn,args.no_skip_multiple_entries)
    elif args.mode == "ref_alt_anc":
        generic_parser(ref_alt_anc_parse_fun, ref_alt_anc_header,
                                           args.in_vcf, args.out_fn,args.no_skip_multiple_entries)
    elif args.mode == "msmc":
        out_fns = [args.out_fn] if args.add_out_fns is None else [args.out_fn]+args.add_out_fns
        assert len(out_fns) == len(args.ind_tuples), "There are {} out files but {} tuples.".format(len(out_fns),len(args.ind_tuples))
        generic_parser(msmc_input_parse_fun, msmc_header, args.in_vcf, out_fns,args.no_skip_multiple_entries,
                                                ind_tuples = [it.split() for it in args.ind_tuples],
                                                                               haplotypes=args.haplotypes)
    elif args.mode == "add_outgroup_from_fasta":
        import pyfasta
        og_fasta = pyfasta.Fasta(args.in_fasta)
        generic_parser(add_outgroup_parse_fun, add_outgroup_header,  args.in_vcf, args.out_fn,
                                                args.no_skip_multiple_entries, og_fasta=og_fasta)
    else:
        raise UserException("Unknown mode.")
