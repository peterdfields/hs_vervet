#!/bin/awk -f
#get the correlation coefficient of Ns between the first samlple and all other samples
#(using the first is arbitary, but it would be too much to do all pairwise correlations)
BEGIN{nsites=0;for(i=1;i<=163;++i) {nx[i]=0;n1x[i]=0}}
/^#[C]/{for(i=1;i<163;i++) printf $(i+9)"\t";printf "\n"}
!/^#/{nsites++;for(i=1;i<=163;++i) {if ($(i+9) == "./.") {nx[i]++;if ($10 == "./.") n1x[i]++}}}
#END{print nsites;for(i=1;i<=1;i++) printf nx[i]";"n1x[i]-1/nsites*nx[1]*nx[i]";"sqrt(nx[1](1-1/nsites*nx[1])*nx[i]*(1-1/nsites*nx[i]))"\t"}
END{for(i=1;i<=163;i++) printf (n1x[i]-1/nsites*nx[1]*nx[i])/(sqrt(nx[1]*(1-1/nsites*nx[1])*nx[i]*(1-1/nsites*nx[i])))"\t"}
