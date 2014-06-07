#!/usr/bin/env python
import os, sys
import pandas as pd
import vcf
import argparse

##todo: extend to multiple info tags added at the same time

parser = argparse.ArgumentParser(description=   
                                "Add one or more tags to INFO column in vcf.\n"
                                "Tag values taken from columns 3+ form input tsv.")
parser.add_argument("in_tsv",type = argparse.FileType('r'), default = '-',
                             help="tsv with chrom,pos,value_tag1,value_tag2,...")
parser.add_argument("in_vcf",type = argparse.FileType('r'), default = '-',
                                                    help="Input vcf filename.")
parser.add_argument("out_vcf",type = argparse.FileType('w'), default = '-',
                                                    help="Output vcf filename.")
parser.add_argument("--tag_name",nargs='+',
                            help="Name for info tags. Number of arguments "
                                                "must match cols in tsv-2.")
parser.add_argument("--tag_description",nargs="+",
                                help="Description for tag to be but in VCF header.")
parser.add_argument("--tag_type",nargs="+",
                                help="Tag value type: Float, String, Integer")
#parser.add_argument("--tsv_col",type=int,nargs='+',help="Column of the tsv to use for annotation")

args = parser.parse_args()

tag_df = pd.read_csv(args.in_tsv,sep='\t',index_col=(0,1))
vcf_reader = vcf.Reader(args.in_vcf)

assert len(args.tag_name) == len(args.tag_description) == len(args.tag_type),\
                                        "--tag_name, --tag_type and "\
                                        "--tag_description must have "\
                                        "same numer of arguments."\

assert len(args.tag_name) == len(tag_df.columns),\
                            "number of --tag_name arguments must match "\
                            "number of tags in tsv, but they are "+\
                            str(args.tag_name) + " and "\
                            + str(tag_df.columns.values)


vcf_reader.infos.update(
        {tn: vcf.parser._Info(tn, '1', tt, td) for tn,td,tt in 
                            zip(args.tag_name,args.tag_description,args.tag_type)}
            )

vcf_writer = vcf.Writer(args.out_vcf, vcf_reader)

for record in vcf_reader:
    try:
        tag_col = tag_df.loc[record.CHROM,int(record.POS)]
        for name, tag, tag_type in zip(args.tag_name,tag_col,args.tag_type):
            if tag_type == "Float":
                tag = round(tag*10000.)/10000
            #if type(aa)==str:
            record.INFO.update({name:tag})
    except KeyError:
        print "key", record.POS, "not in tsv"
    vcf_writer.write_record(record)
vcf_writer.flush()
vcf_writer.close()


