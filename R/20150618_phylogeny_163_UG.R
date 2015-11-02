library(phangorn)


fn <- "/home/GMI/hannes.svardal/vervet_project/analyses/20150504_163_UnifiedGenotyper/_data/pi_per_site.tsv"


#split mat in kilo years(2* because difference is 2 times the coalescent time)
diffmat <- as.dist(as.matrix(read.table(fn, head=T, row.names=1)))
#make Unweighted Pair Group Method with Arithmetic Mean tree (hirarchical clustering, assumes equal rates of evolution)
treeUPGMA2 <- upgma(diffmat)
#write tree to file
#write.tree(treeUPGMA2,file="~/vervet_project/analyses/20140607_vervetpopgen_manuscript/_data/UPGMA_tree_163_pi_per_site.newick")

plot(treeUPGMA2,)

