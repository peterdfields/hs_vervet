library(ape)

fn_no_mac <- "~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/split_time_10x.tsv"
fn_mac <- "~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/split_time_10x_macaque.tsv"

#par(mfrow=c(2,1))

make_tree <- function(fn, type='nj'){
  split_mat <- as.dist(as.matrix(read.table(fn, head=T, row.names=1)))/1000
  if (type == 'nj') {
    split_tree <- nj(split_mat)
    }
  else if (type == 'hc'){
    split_tree <- hclust(split_mat)
    }
  return(split_tree)
}

save_tree <- function(infn, type='nj',outfn){
  tree <- make_tree(fn,type)
  tree.write(out_fn)
}

t <- make_tree(fn_mac,type='nj')

nj.out <- root(t,"mac",resolve.root=TRUE)

plot(nj.out)

t.write("~/vervet_project/analyses/20140524_SNP_calling_ref3500_10X_individuals/output/nj_tree_10x_root_macaque.tree")


#root tree
#t$root.edge <- 2 #rpde = NULL)#,resolve.root=TRUE)
#make tree ultrametric
#t <- chronopl(t, lambda=0)
#drop.tip(t,3)
#show(t)
#t <- as.hclust.phylo(t)

#plot(t)

#t2 <- make_tree(fn_no_mac,type='hc')
#plot(t2,hang=-1) #hang=-1 plots labels down to present
#png('hclust_tree_dendrogram.png',#width=2200,height=2000,units="px"
#    width=1200,height=1000,units="px",pointsize = 50)
#png('hclust_tree_dendrogram.png',
#    width=3.25,height=3,units="in",res=800)
#plot(as.dendrogram(t2),ylab = "split time (kya)")
#dev.off()
#write.tree(split_tree,file="split_tree.newick")

#plot(split_tree,type='unrooted')

#why doesn't rhis work this tree is still unrooted
#rooted_split_tree <- root(split_tree, "mac",resolve.root=TRUE)

#write.tree(rooted_split_tree,file="rooted_split_tree.newick")

#split_dendro <- as.dendrogram(as.hclust(multi2di(split_tree)))

#plot(split_dendro)

#plot(rooted_split_tree)

#split_hclust <- hclust(as.dist(split_mat)) 

#plot(split_hclust)

#chronopl(split_tree, lambda=0.2)
#split_tree$root.edge <- 0
#print(is.ultrametric(split_tree))
#print(is.binary.tree(split_tree))
#print(is.rooted(split_tree))
#plot(split_tree)
#hc <- as.hclust.phylo(split_tree)
#plot(hc)