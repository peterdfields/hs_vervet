tree = treeUPGMA

man_ana_dir="/home/GMI/hannes.svardal/vervet_project/analyses/20150520_popgen_analysis_UG/_data"

plot_ids = fromJSON(file=paste(man_ana_dir,"matrix_plot_ids.json",sep='/'))

tree =rotateConstr(tree, plot_ids)


edge.colors = tree$edge[,2]

#color for external nodes corresponds to pop3
edge.colors[edge.colors%in%1:length(tree$tip)]<-mt[tree$tip[edge.colors[edge.colors%in%1:length(tree$tip)]],]$colors3
#set internal nodes black
edge.colors[edge.colors%in%length(tree$tip):length(edge.colors)]<- rgb(0,0,0)
#devSVG(paste(figure_dir,"163_rooted.svg",sep='/'))
png(paste(figure_dir,"163_rooted.png",sep='/'))
plot.phylo(tree,
           edge.color=edge.colors,
           show.tip.label=FALSE)
#nodelabels(1:324,1:324)

tree_order <- tree$tip.label[tree$edge[1:324,2]]
tree_order = toJSON(tree_order[!is.na(tree_order)])
fileConn<-file(paste(man_ana_dir,"matrix_plot_ids_tree.json",sep='/'))
write(tree_order, fileConn)
close(fileConn)