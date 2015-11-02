library(phangorn)

#geographic distance:

fn <- "~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/geographic_dist_10x.tsv"
#split mat in kilo years 
distmat <- as.dist(as.matrix(read.table(fn, head=T, row.names=1)))
#make Unweighted Pair Group Method with Arithmetic Mean tree (hirarchical clustering, assumes equal rates of evolution)
treeUPGMA = upgma(distmat)
#write tree to file
write.tree(treeUPGMA,file="~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/UPGMA_tree_10x_geographic_distance.newick")
plot(treeUPGMA)


# fn <- "~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/split_time_10x.tsv"
# 
#   
# #split mat in kilo years 
# split_mat <- as.dist(as.matrix(read.table(fn, head=T, row.names=1)))/1000
# #distance is 2*splittimes
# distmat = split_mat*2
# #make Unweighted Pair Group Method with Arithmetic Mean tree (hirarchical clustering, assumes equal rates of evolution)
# treeUPGMA = upgma(distmat)
# #write tree to file
# write.tree(treeUPGMA,file="~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/UPGMA_tree_10x.newick")
# 
# #optional:
# #plot(treeUPGMA)
# 
# #for comparison: tree on the pairwise distance matrix (not correcting for within-species differences)
# 
# fn <- "~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/coal_time_10x.tsv"
# 
# 
# #split mat in kilo years(2* because difference is 2 times the coalescent time)
# diffmat <- 2*as.dist(as.matrix(read.table(fn, head=T, row.names=1)))/1000
# #make Unweighted Pair Group Method with Arithmetic Mean tree (hirarchical clustering, assumes equal rates of evolution)
# treeUPGMA2 <- upgma(diffmat)
# #write tree to file
# write.tree(treeUPGMA2,file="~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/UPGMA_tree_10x_coal_time.newick")
