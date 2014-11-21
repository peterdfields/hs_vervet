#!/bin/awk -f
#get the correlation coefficient of Ns between the first samlple and all other samples
#(using the first is arbitary, but it would be too much to do all pairwise correlations)
#awk -v countfile="${countfile}" '
BEGIN {nsites=0;
      for(i=1;i<=163;++i) {nx[i]=0;
                           for(j=1;j<=i;j++) n1x[i,j]=0
                           }
       }
/^#[C]/ {printf "\t";
         printf "\t" > countfile;
         for(i=1;i<=163;i++) {id[i]=$(i+9);printf $(i+9)"\t";printf $(i+9)"\t" > countfile};
         printf "\n" > countfile
         printf "\n"}
!/^#/   {nsites++;
         for(i=1;i<=163;++i){
            if ($(i+9) == "./.") nx[i]++;
            for(j=1;j<=i;j++) {if ($(i+9) == "./." && $(j+9) == "./.") n1x[i,j]++}}}
#END{print nsites;for(i=1;i<=1;i++) printf nx[i]";"n1x[i]-1/nsites*nx[1]*nx[i]";"sqrt(nx[1](1-1/nsites*nx[1])*nx[i]*(1-1/nsites*nx[i]))"\t"}
END     {printf "missing_num" > countfile;
         for(i=1;i<=163;i++) {printf id[i];
                              printf "\t"nx[i] > countfile;
                              for(j=1;j<=i;j++) {printf "\t"(n1x[i,j]-1/nsites*nx[j]*nx[i])/(sqrt(nx[j]*(1-1/nsites*nx[j])*nx[i]*(1-1/nsites*nx[i])))};
                               printf "\n"
                             };
         printf "\n""missing_percent" > countfile;
         for(i=1;i<=163;i++) printf "\t"nx[i]/nsites > countfile;
         }
#'
