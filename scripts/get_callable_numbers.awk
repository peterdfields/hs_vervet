#!/bin/awk -f

BEGIN {OFS="\t";callable_nonvar=0;pass_snp=0;noncallable_nonvar=0;filter_snp=0;total=0}
NR>1{if ($2==last) next;last=$2; # skip over indels
              total++;
            n=split($4,filters,","); # split filter column
            if ($6!="NA"){if (filters[1]=="PASS") pass_snp++ #check whether SNP
                          else filter_snp++} 
            else {if (filters[1]=="PASS") callable_nonvar++
                  else {if (n<2 && filters[1]=="LowQual") callable_nonvar++
                        else noncallable_nonvar++}}}
END {print "total","pass_nonvar","filter_nonvar","pass_snp","filter_snp";
     print total,callable_nonvar,noncallable_nonvar,pass_snp,filter_snp}
