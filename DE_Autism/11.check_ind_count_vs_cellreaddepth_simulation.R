#this code is modified from 7.2_fit_zinb.R
#this code is used for checking the relationship between cell level read depth and counts. 
# we use nb to fit  cell-level-count ~ log-read depth
# if the coefficient ~1, then we say the relathionship between cell-level-count ~ log-read depth are around 1.

#---------------------------log of 7.2-----------------------
#this code fit raw count data with ZINB/NB for each individuals per gene.

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with the original cells of certain type,certain individual
# cluster_tag=1
# file_tag="3k10"
# pre_tag="dca" #c("dca","scvi")
# 
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/10.Result/sim_v6/sim_data/")


source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/7.0_ZINB_fit_functions.R")

file_tag="1.2_1.2_0.6_0.2_1"

sim_matrix=readRDS(paste0("sim_matrix_",file_tag,".rds"))
sim_meta=readRDS(paste0("sim_meta_",file_tag,".rds"))
sim_cell_readdepth=readRDS(paste0("sim_cell_readdepth_",file_tag,".rds"))

de.mean=readRDS(paste0("../de_label/sim_de.mean_",file_tag,".rds"))
de.var=readRDS(paste0("../de_label/sim_de.var_",file_tag,".rds"))
de.mult=readRDS(paste0("../de_label/sim_de.mult_",file_tag,".rds"))
de.dp=readRDS(paste0("../de_label/sim_de.dp_",file_tag,".rds"))


hist(sim_cell_readdepth,breaks=100)

#adjust1, devide directly
sim_matrix_adj=sim_matrix
sim_matrix_adj[]=c(sim_matrix_adj)/rep(sim_cell_readdepth,each=nrow(sim_matrix))*10000

#adjust2, lm regression
log_sim_matrix=log(sim_matrix+1)
log_sim_cell_readdepth=log(sim_cell_readdepth)
log_sim_resid=matrix(ncol=ncol(sim_matrix),nrow=nrow(sim_matrix))
rownames(log_sim_resid)=rownames(sim_matrix)
colnames(log_sim_resid)=colnames(sim_matrix)

for(ig in 1:nrow(log_sim_matrix)){
  log_sim_resid[ig,]=lm(log_sim_matrix[ig,]~log_sim_cell_readdepth)$residuals
  print(ig)
}

library(ggplot2)

phenotype=sim_meta$phenotype
phenotype[phenotype==1]="case"
phenotype[phenotype==0]="ctrl"
phenotype=as.factor(phenotype)
ind=as.factor(sim_meta$individual)

pdf(paste0("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/11.check_ind_count_vs_cellreaddepth_simulation_",file_tag,".pdf"),width = 5,height=3)
#for mean-DE
sub_ind=1:nrow(sim_meta)
sub_ind_num=10
ncell=800
sub_ind=((nrow(sim_meta)/2)-ncell*sub_ind_num+1):((nrow(sim_meta)/2)+ncell*sub_ind_num)

for(ig in which(de.mean==1)[1:5]){
  #count=sim_matrix[ig,]
  #count_adj=sim_matrix_adj[ig,]
  log_count=log(sim_matrix[ig,]+1)
  log_count_adj=log(sim_matrix_adj[ig,]+0.1)
  log_count_resid=log_sim_resid[ig,]
  cur_df=data.frame(log_count,log_count_adj,log_count_resid,ind,phenotype)
  cur_df=cur_df[sub_ind,]
  cur_df$sim_meta.phenotype=as.factor(cur_df$phenotype)
  p1=ggplot(cur_df, aes(x=log_count, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("raw log_count: mean-DE gene", ig))
  p2=ggplot(cur_df, aes(x=log_count_adj, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("adj log_count: mean-DE gene", ig))
  p3=ggplot(cur_df, aes(x=log_count_resid, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("resid log_count: mean-DE gene", ig))
  print(p1)
  print(p2)
  print(p3)
}

#for var-DE
for(ig in which(de.var==1)[1:5]){
  #count=sim_matrix[ig,]
  #count_adj=sim_matrix_adj[ig,]
  log_count=log(sim_matrix[ig,]+1)
  log_count_adj=log(sim_matrix_adj[ig,]+0.1)
  log_count_resid=log_sim_resid[ig,]
  cur_df=data.frame(log_count,log_count_adj,log_count_resid,ind,phenotype)
  cur_df=cur_df[sub_ind,]
  cur_df$sim_meta.phenotype=as.factor(cur_df$phenotype)
  p1=ggplot(cur_df, aes(x=log_count, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("raw log_count: var-DE gene", ig))
  p2=ggplot(cur_df, aes(x=log_count_adj, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("adj log_count: var-DE gene", ig))
  p3=ggplot(cur_df, aes(x=log_count_resid, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("resid log_count: var-DE gene", ig))
  print(p1)
  print(p2)
  print(p3)
}

#for dp-DE
for(ig in which(de.dp==1)[1:5]){
  #count=sim_matrix[ig,]
  #count_adj=sim_matrix_adj[ig,]
  log_count=log(sim_matrix[ig,]+1)
  log_count_adj=log(sim_matrix_adj[ig,]+0.1)
  log_count_resid=log_sim_resid[ig,]
  cur_df=data.frame(log_count,log_count_adj,log_count_resid,ind,phenotype)
  cur_df=cur_df[sub_ind,]
  cur_df$sim_meta.phenotype=as.factor(cur_df$phenotype)
  p1=ggplot(cur_df, aes(x=log_count, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("raw log_count: dp-DE gene", ig))
  p2=ggplot(cur_df, aes(x=log_count_adj, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("adj log_count: dp-DE gene", ig))
  p3=ggplot(cur_df, aes(x=log_count_resid, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("resid log_count: dp-DE gene", ig))
  print(p1)
  print(p2)
  print(p3)
}

#for mult-DE
for(ig in which(de.mult==1)[1:5]){
  #count=sim_matrix[ig,]
  #count_adj=sim_matrix_adj[ig,]
  log_count=log(sim_matrix[ig,]+1)
  log_count_adj=log(sim_matrix_adj[ig,]+0.1)
  log_count_resid=log_sim_resid[ig,]
  cur_df=data.frame(log_count,log_count_adj,log_count_resid,ind,phenotype)
  cur_df=cur_df[sub_ind,]
  cur_df$sim_meta.phenotype=as.factor(cur_df$phenotype)
  p1=ggplot(cur_df, aes(x=log_count, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("raw log_count: mult-DE gene", ig))
  p2=ggplot(cur_df, aes(x=log_count_adj, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("adj log_count: mult-DE gene", ig))
  p3=ggplot(cur_df, aes(x=log_count_resid, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("resid log_count: mult-DE gene", ig))
  print(p1)
  print(p2)
  print(p3)
}


#for none-DE
for(ig in which(de.mean+de.var+de.dp+de.mult==0)[1:5]){
  #count=sim_matrix[ig,]
  #count_adj=sim_matrix_adj[ig,]
  log_count=log(sim_matrix[ig,]+1)
  log_count_adj=log(sim_matrix_adj[ig,]+0.1)
  log_count_resid=log_sim_resid[ig,]
  cur_df=data.frame(log_count,log_count_adj,log_count_resid,ind,phenotype)
  cur_df=cur_df[sub_ind,]
  cur_df$sim_meta.phenotype=as.factor(cur_df$phenotype)
  p1=ggplot(cur_df, aes(x=log_count, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("raw log_count: none-DE gene", ig))
  p2=ggplot(cur_df, aes(x=log_count_adj, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("adj log_count: none-DE gene", ig))
  p3=ggplot(cur_df, aes(x=log_count_resid, color=phenotype, group=ind)) +geom_density(bw=.2) + ggtitle(paste0("resid log_count: none-DE gene", ig))
  print(p1)
  print(p2)
  print(p3)
}

dev.off()
