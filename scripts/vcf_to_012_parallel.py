#!/usr/bin/env python
"""
Creates 0,1,2 genotype matrix from VCF
"""
import sys, os

def write_header(reader,tf):
    print "running write header in vcf_to_012", reader, tf
    tf.write("chrom"+"\t"+"pos"+"\t"+"\t".join(reader.samples)+"\n")

def write_012_tsv(reader,tf,ancestral_derived):
    allele0 = ["0/0","0|0"]
    allele1 = ["1/1","1|1"]

    for record in reader:

        if ancestral_derived:
            try:
                aa = record.INFO['AA']
            except KeyError:
                continue
            if aa == record.REF:
                allele0 = ["0/0","0|0"]
                allele1 = ["1/1","1|1"]
            elif aa == record.ALT[0]:
                allele0 = ["1/1","1|1"]
                allele1 = ["0/0","0|0"]
            else:
                continue

        tf.write(record.CHROM+"\t"+str(record.POS))

        for sample in record.samples:
            tf.write("\t")
            if sample["GT"] is None:
                tf.write("N")
            elif sample["GT"] in allele0:
                tf.write("0")
            elif sample["GT"] in ["0/1","0|1","1|0"]:
                tf.write("1")
            elif sample["GT"] in allele1:
                tf.write("2")
            else:
                raise ValueError("Unsupported genotype " + sample["GT"])
        tf.write("\n")

if __name__ == "__main__":
    import argparse
    from hs_vervet.tools import bioparallel

    parser = argparse.ArgumentParser(description="Convert a vcf into a tsv "
                                                 "with 0,1,2 genotypes.\n"
                                                 "0..reference,1..het,2..alt.")
    parser.add_argument("in_vcf",help="Input vcf filename.")
    parser.add_argument("out_fn",help="Output tsv filename.")
    parser.add_argument("--ancestral_derived",action = 'store_true',
                         help=  "0,1,2 corresponds to ancestral vs\n"
                                "derived instead of ref vs alt.\n"
                                "Using AA tag from VCF info.\n"
                                "Sites for which this tag is "
                                "unavailable or neither ref \n"
                                               "nor alt are ignored.")

    #parallel specific arguments
    parser.add_argument("--chrom",nargs='+',help="Chromosome")
    parser.add_argument("--chrom_len",type=int,default=None,nargs='*',
                                                    help="Chromosome length")
    parser.add_argument("--ncpus",type=int,default=None,
                            help="Number of processed to spawn.")
    parser.add_argument("--tmp_dir",default=".",
                        help="Directory to write temporary files to.")

    args = parser.parse_args()

    write_fun = lambda r,f: write_012_tsv(r,f,args.ancestral_derived)

    parser =  bioparallel.VCFParser(args.in_vcf,write_fun,
                                    chromosomes=args.chrom,
                                    chrom_len=args.chrom_len,
                                    mode="write",
                                    out_fn=args.out_fn,
                                    write_header_fun=write_header,
                                    tmp_dir=args.tmp_dir)
    parser.run(ncpus=args.ncpus)

