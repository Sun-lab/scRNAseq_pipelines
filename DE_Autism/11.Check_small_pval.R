dist_method="JSD"        #c("mean","JSD")
fit_method="empirical"   #c("empirical","zinb")
file_tag=1
r_mean=1.5
r_var=1.5
r_disp=1.5
r_change_prop=0.4
n=10      #c(20,15,10,5)
ncell=20 #c(100,80,60,40,20)



cluster_tag=11
file_tag="3k10"
pre_tag="dca"
dist_method="jsd"
fit_method="nbzinb"
F_method="p"
perm_label=1

ind_covariate_flag="ind" 
perm_num=500
tol=0.2


if(pre_tag=="dca"){
  t_mean=read.table(paste0("res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
  t_dispersion=read.table(paste0("res_dca_rawM",file_tag,"/dispersion.tsv"),stringsAsFactors = FALSE,row.names = 1)
  t_dropout=read.table(paste0("res_dca_rawM",file_tag,"/dropout.tsv"),stringsAsFactors = FALSE,row.names = 1)
}

library("ggplot2")
library("emdbook")

###########functions#############
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/7.0_ZINB_fit_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/9.0_Fstat_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/8.0_kl_divergence_functions.R")
setwd("/Volumes/SpecialPass/fh_data/Data_PRJNA434002/")

dist_array=readRDS(paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

dist_pval=readRDS(paste0("../Data_PRJNA434002/8.Result/kl_pval/p",perm_label,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

pval_order=order(dist_pval)

fit_data=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
sim_data=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

#input phenotype
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}
total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)
phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0



case_names=cur_individual[phenotype==1]
control_names=cur_individual[phenotype==0]

pdf(paste0("../Data_PRJNA434002/8.Result/check/check_p",perm_label,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),width = 30,height = 20)
for(iplot in pval_order[c(1:5)]){
  cur_gene=dimnames(sim_data)[[1]][iplot]
  print(c(cur_gene,dist_pval[iplot]))
  print(" case fit ")
  print(fit_data[dimnames(sim_data)[[1]]==cur_gene,phenotype==1,])
  print(" ctrl fit")
  print(fit_data[dimnames(sim_data)[[1]]==cur_gene,phenotype==0,])
  print(" case matrix")
  print(dist_array[dimnames(sim_data)[[1]]==cur_gene,phenotype==1,phenotype==1])
  print(" ctrl matrix")
  print(dist_array[dimnames(sim_data)[[1]]==cur_gene,phenotype==0,phenotype==0])
  
  op=par(mfrow = c(4, 6))
  for(i in case_names){
    hist(sim_data[dimnames(sim_data)[[1]]==cur_gene,meta$individual==i,],col=rgb(1,0,0,0.8),main=paste0(cur_gene," case ",i))}
  for(i in control_names){
    hist(sim_data[dimnames(sim_data)[[1]]==cur_gene,meta$individual==i,],col=rgb(0,0,1,0.8),main=paste0(cur_gene," ctrl ",i))}
  par(op)
  op=par(mfrow = c(2, 1))
  hist(sim_data[dimnames(sim_data)[[1]]==cur_gene,meta$diagnosis!="Control",],main=paste0(cur_gene," case ",i))
  hist(sim_data[dimnames(sim_data)[[1]]==cur_gene,meta$diagnosis=="Control",],main=paste0(cur_gene," ctrl ",i))
  par(op)
}

dev.off()




