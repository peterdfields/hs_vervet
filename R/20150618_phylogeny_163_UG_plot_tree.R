library("rjson")
library(RSvgDevice)



figure_dir <- "/home/GMI/hannes.svardal/Akademisches/Projects/VervetPopgen/Figures"

mt = read.csv("/home/GMI/hannes.svardal/vervet_project/metadata/163_master_table.csv",header=TRUE)
rownames(mt) <- mt$ucla_id

pops3 = c("aet","cyn","pyn","pys","sab","sac","sar","tan")
colors3 = fromJSON(file="/home/GMI/hannes.svardal/vervet_project/metadata/colors3.json")

mt$colors3='N'

for (i in 1:length(pops3)) {
  mt[mt$pop3==pops3[i],]$colors3 = rgb(colors3[[i]][1],colors3[[i]][2],colors3[[i]][3])
}



# plot.phylo(tree,type='unrooted',
#            edge.color=edge.colors,
#            cex=0.5,lab4ut="axial",
#            label.offset=0.00005,
#            tip.color=mt[tree$tip.label,]$colors3)
plot_rooted <- function(tree,save=FALSE){
  edge.colors = tree$edge[,2]
  
  #color for external nodes corresponds to pop3
  edge.colors[edge.colors%in%1:length(tree$tip)]<-mt[tree$tip[edge.colors[edge.colors%in%1:length(tree$tip)]],]$colors3
  #set internal nodes black
  edge.colors[edge.colors%in%length(tree$tip):length(edge.colors)]<- rgb(0,0,0)
  if (save) {
    devSVG(paste(figure_dir,"163_rooted.svg",sep='/'))}
  plot.phylo(tree,
           edge.color=edge.colors,
           cex=0.5,
           show.node.label=TRUE,
           label.offset=0.00005,
           tip.color=mt[tree$tip.label,]$colors3)
  nodelabels(1:324, 1:324)
  if (save) {
    dev.off()}
}


plot_unrooted <- function(tree2,save=FALSE){
  
  #Unrooted plot from neighbor joining tree
  #make tree with sampling locations as labels (for temporary use)
  sampling.locations <- as.vector(mt[tree2$tip.label,]$site)
  edge.colors2 = tree2$edge[,2]
  #color for external nodes corresponds to pop3
  edge.colors2[edge.colors2%in%1:length(tree2$tip)]<-mt[tree2$tip[edge.colors2[edge.colors2%in%1:length(tree2$tip)]],]$colors3
  #set internal nodes black
  edge.colors2[edge.colors2%in%length(tree2$tip):length(edge.colors2)]<- rgb(0,0,0)
  tree2$tip.label <- sampling.locations
  if (save) {
    devSVG(paste(figure_dir,"163_unrooted_site.svg",sep='/'))}
  plot.phylo(tree2,type='unrooted',
            edge.color=edge.colors2,
            cex=0.5,lab4ut="axial",
            label.offset=0.00005)
  if (save) {
    dev.off()}
}


