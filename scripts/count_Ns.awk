#!/bin/awk -f
#count the number of Ns for each sample in a VCF
BEGIN{for(i=1;i<=163;++i) x[i]=0} 
/^#[C]/{for(i=1;i<163;i++) printf $(i+9)"\t";printf "\n"} 
!/^#/{for(i=1;i<=163;++i) {if ($(i+9) == "./.") x[i]++}} 
END{for(i=1;i<=163;i++) printf x[i]"\t"}
