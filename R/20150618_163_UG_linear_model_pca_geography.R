library("rjson")

man_ana_dir <- "~/vervet_project/analyses/20150520_popgen_analysis_UG/_data"
pca = read.csv(paste(man_ana_dir,"pcadapt_linear_model_africa.csv",sep='/'),header=TRUE)

#library(lme4)
#model1 <- lmer(PCA~long+lat)
#model2<-lmer(PCA~long+lat+spec)
#summary(model1)$loglik

for (i in 1:6) {
  f <- paste(paste('X',i,sep=''), "~", "longitude + latitude")
  f2 <- paste(f,'+ population')
  model1 <- lm(f, data=pca)
  model2 <- lm(f2, data=pca)
  anov <- anova(model1, model2)
  print(anov)
  capture.output(anov, file=paste(man_ana_dir,'/',"pcadapt_ibd_anova_K",i,".txt",sep='')) 
}