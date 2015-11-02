library("rjson")
library(phangorn)
library(RSvgDevice)

figure_dir <- "/home/GMI/hannes.svardal/Akademisches/Projects/VervetPopgen/Figures"

fn <-"/home/GMI/hannes.svardal/vervet_project/analyses/20150504_163_UnifiedGenotyper/_data/subspecies_Africa_pi_minus_within.tsv"

pops3 = c("aet","cyn","pyn","pys","sab","sac","sar","tan")
colors3 = fromJSON(file="/home/GMI/hannes.svardal/vervet_project/metadata/colors3.json")


#split mat in kilo years(2* because difference is 2 times the coalescent time)
diffmat <- as.dist(as.matrix(read.table(fn, head=T, row.names=1, sep='\t')))
#make Unweighted Pair Group Method with Arithmetic Mean tree (hirarchical clustering, assumes equal rates of evolution)
treeUPGMA2 <- upgma(diffmat)


tree.order <- c("sab","tan","pyn",'pys','cyn','aet')
treeUPGMA2 =rotateConstr(treeUPGMA2, tree.order)

#write tree to file
write.tree(treeUPGMA2,file="~/vervet_project/analyses/20150520_popgen_analysis_UG/_data/UPGMA_tree_subspecies_Africa_pi_minus_within.newick")
#devSVG(paste(figure_dir,"163_subspecies_Africa_UPGMA.svg",sep='/'))
plot(treeUPGMA2,cex=1.8,
     label.offset=0.00005)
#dev.off()