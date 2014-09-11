#!/usr/bin/env python
from hs_vervet.tools import bioparallel
import vcf as pyvcf
import os



def write_vcf(reader, writer):
    for i,record in enumerate(reader):
        #do something on reader
        record.ALT = "B"
        writer.write_record(record)
        if i>1000:
            break
    return fn
def update_vcf_header_1(reader):
    reader.infos.update({'AA':pyvcf.parser._Info('AA', '1', 'String', 'Ancestral Allele as derived from {}'.format("a cow"))})


fn = os.path.expanduser("~/1001genomes_project/final/VCFpopIntersectionGMIandMPI/eventOnly/"
                        "1001genomes_snp-short-indel.vcf.bgzip.gz")

parser =  bioparallel.VCFParser(fn,write_vcf,
                    chrom="1",chrom_len=30427671,mode="vcf_write",
                    out_fn="parallel_parse.test",
                            update_vcf_header_fun=update_vcf_header_1)

print "running parser"
parser.run(ncpus=)
print "finished"
