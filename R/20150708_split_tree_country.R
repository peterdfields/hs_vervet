library(phangorn)
library(RSvgDevice)

figure_dir <- "/home/GMI/hannes.svardal/Akademisches/Projects/VervetPopgen/Figures"

fn <- "/home/GMI/hannes.svardal/vervet_project/analyses/20150504_163_UnifiedGenotyper/_data/country_pi_minus_within.tsv"




#split mat in kilo years(2* because difference is 2 times the coalescent time)
diffmat <- as.dist(as.matrix(read.table(fn, head=T, row.names=1, sep='\t')))
#make Unweighted Pair Group Method with Arithmetic Mean tree (hirarchical clustering, assumes equal rates of evolution)
treeUPGMA2 <- upgma(diffmat)


tree.order <- c("Nevis","Saint Kitts","Gambia","Barbados","Ghana","Ethiopia","Central African Republic","Kenya","Tanzania","Zambia","Boswana","South Africa")
treeUPGMA2 =rotateConstr(treeUPGMA2, tree.order)

#write tree to file
write.tree(treeUPGMA2,file="~/vervet_project/analyses/20150520_popgen_analysis_UG/_data/UPGMA_tree_country_pi_minus_within.newick")
devSVG(paste(figure_dir,"163_country_UPGMA.svg",sep='/'))
plot(treeUPGMA2)
dev.off()
