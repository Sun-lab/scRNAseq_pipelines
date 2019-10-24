#this code calculate the KL divergence between each individual.

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with each other using the KL divergence

library("ggplot2")
library("emdbook")

#cluster_tag=1
#file_tag="3k10"
#pre_tag="dca" #c("dca","scvi")

sim_n=10
n_bin=20
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

###########functions#############
source("./Command/8.0_kl_divergence_functions.R")
###########input###############
if(pre_tag=="dca"){
  t_mean=read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
  t_dispersion=read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_tag,"/dispersion.tsv"),stringsAsFactors = FALSE,row.names = 1)
  t_dropout=read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_tag,"/dropout.tsv"),stringsAsFactors = FALSE,row.names = 1)
}
if(pre_tag=="scvi"){
  t_mean=t(read.table(paste0("../Data_PRJNA434002/scvi_",file_tag,"/scvi_imputed_values_",file_tag,".txt"),stringsAsFactors = FALSE,sep=",")) 
  t_dispersion=t(read.table(paste0("../Data_PRJNA434002/scvi_",file_tag,"/scvi_lnDispersion_",file_tag,".txt"),stringsAsFactors = FALSE,sep=",")) #originaly in log scale
  t_dropout=t(read.table(paste0("../Data_PRJNA434002/scvi_",file_tag,"/scvi_dropout_rate_",file_tag,".txt"),stringsAsFactors = FALSE,sep=",")) #originaly in logit scale
  t_dispersion=exp(t_dispersion)
  t_dispersion=matrix(rep(t_dispersion,times=ncol(t_mean)),nrow=length(t_dispersion),ncol=ncol(t_mean))
  t_dropout=exp(t_dropout)/(1+exp(t_dropout))
}
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}

##########data processing############
#1. restrict the cluster
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
sub_mean=as.matrix(t_mean[,tmeta$cluster==cur_cluster])

sub_dispersion=as.matrix(t_dispersion[,tmeta$cluster==cur_cluster])
sub_dropout=as.matrix(t_dropout[,tmeta$cluster==cur_cluster])

meta=tmeta[tmeta$cluster==cur_cluster,]

#2.generate idividual label
cur_individual=unique(meta$individual)

#3. calculate KL divergence.

kl_array=array(dim=c(nrow(sub_mean),length(cur_individual),length(cur_individual),2),
                      dimnames = list(rownames(sub_mean),cur_individual,cur_individual,c("equal","quantile")))

for(i_g in 1:nrow(sub_mean)){
  cur_sim=matrix(ncol=sim_n,nrow=ncol(sub_mean))
  for(i_s in 1:ncol(sub_mean)){
    cur_sim[i_s,]=emdbook::rzinbinom(sim_n,sub_mean[i_g,i_s], sub_dispersion[i_g,i_s], sub_dropout[i_g,i_s])
  }
  for(i_ind_a in 1:length(cur_individual)){
    for(i_ind_b in 1:length(cur_individual)){
      cur_ind_a=cur_individual[i_ind_a]
      cur_ind_b=cur_individual[i_ind_b]
      #fit sim
      cur_sim_ind_a=as.numeric(cur_sim[meta$individual==cur_ind_a,])
      cur_sim_ind_b=as.numeric(cur_sim[meta$individual==cur_ind_b,])
      kl_array[i_g,i_ind_a,i_ind_b,1]=mean_KL(cur_sim_ind_a,cur_sim_ind_b,bmeth="equal",nbins=n_bin)
      kl_array[i_g,i_ind_a,i_ind_b,2]=mean_KL(cur_sim_ind_a,cur_sim_ind_b,bmeth="quantile",nbins=n_bin)
    }
  }
  print(i_g)
}


saveRDS(kl_array,paste0("../Data_PRJNA434002/kl_array_impirical_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))


sessionInfo()
q(save="no")
