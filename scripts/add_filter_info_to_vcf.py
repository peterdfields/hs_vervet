#!/usr/bin/env python
import subprocess


def add_filter_info_to_vcf(in_vcf,out_vcf,expressions,descriptions):
    """
    adds info header entries for each pair in expressions/descriptions iterables to in_vcf and writes it to out_vcf
    """

    for i,line in enumerate(in_vcf):
        if line[:2] == '##':
            b = len("##FILTER=<ID=")
            if line[:b] == "##FILTER=<ID=":
                for j,e in enumerate(expressions):
                    match = False
                    if line[b:b+len(e)] == e:
                        newline = line[:b+len(e)+len(',Description="')] + descriptions[j] + '">\n'
                        del expressions[j]
                        del descriptions[j]
                        match = True
                        break
                if not match:
                    newline = line
            else:
                newline = line
            out_vcf.write(newline)
        else:
            for e,d in zip(expressions,descriptions):
                out_vcf.write('##FILTER=<ID=' + e + ',Description="' + d + '">\n')
            out_vcf.flush()
            command = ["tail","-n +" + str(i+1),in_vcf.name]
            print command
            p = subprocess.Popen(" ".join(command),shell=True,stdout=out_vcf)
            #p.wait()
            p.communicate()
            break


if __name__ == '__main__':
    import os, sys
    import argparse

    parser = argparse.ArgumentParser(description="Add or update filter information in vcf header.")
    parser.add_argument("in_vcf", type = argparse.FileType('r'), default = '-',help="Input vcf filename.")
    parser.add_argument("-o",'--out_vcf',type = argparse.FileType('w'),default = '-',help="Output vcf filename.")
    parser.add_argument("-f","--filter-expression",default=None, help="name of the outgroup",nargs="+")
    parser.add_argument("-d","--filter-description",default=None,help="name of the outgroup",nargs="+")

    args = parser.parse_args()


    if args.filter_expression is None or args.filter_description is None:
        parser.print_help()
        raise Exception('filter expressions and descriptions must be provided.')
    elif len(args.filter_expression) != len(args.filter_description):
        raise Exception('filter expressions and descriptions must be of equal length.')

    add_filter_info_to_vcf(args.in_vcf,args.out_vcf,args.filter_expression,args.filter_description)

