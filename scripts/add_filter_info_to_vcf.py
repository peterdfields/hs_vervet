#!/usr/bin/env python
import vcf


def add_filter_info_to_vcf(in_vcf,out_vcf,expressions,descriptions):
    """
    adds info header entries for each pair in expressions/descriptions iterables to in_vcf and writes it to out_vcf
    """


    try:
        open(in_vcf,'r')
    except TypeError, e:
        if 'file found' in str(e):
            pass
        else:
            raise

    reader = vcf.Reader(in_vcf)

    for expr, desc in zip(expressions,descriptions):
        reader.filters.update({expr:vcf.parser._Filter(expr,desc)})
            
    vcf_writer = vcf.Writer(open(out_vcf, 'w'), reader) 
    for record in reader:
        vcf_writer.write_record(record)
    vcf_writer.flush()
    vcf_writer.close()

    
    




if __name__ == '__main__':
    import os, sys
    import argparse

    parser = argparse.ArgumentParser(description="Add or update filter information in vcf header.")
    parser.add_argument("in_vcf", type = argparse.FileType('r'), default = '-',help="Input vcf filename.")
    parser.add_argument("out_vcf",help="Output vcf filename.")
    parser.add_argument("-f","--filter-expression",default=None, help="name of the outgroup",nargs="+")
    parser.add_argument("-d","--filter-description",default=None,help="name of the outgroup",nargs="+")

    args = parser.parse_args()

    def get_absolute_path(fn):
        if fn[0]=='/':
            return fn
        else:
            return os.path.join(os.getcwd(),fn)

    if args.filter_expression is None or args.filter_description is None:
        parser.print_help()
        raise Exception('filter expressions and descriptions must be provided.')
    elif len(args.filter_expression) != len(args.filter_description):
        raise Exception('filter expressions and descriptions must be of equal length.')

    add_filter_info_to_vcf(args.in_vcf,args.out_vcf,args.filter_expression,args.filter_description)

