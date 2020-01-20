#this code calculate the KL divergence between each individual. #this version calculate directly from zinb distribution, output of DCA

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with each other using the KL divergence

library("ggplot2")
library("emdbook")

#cluster_tag=1
#file_tag="3k10"
#pre_tag="dca" #c("dca","scvi")



#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########functions#############
source("./Command/8.0_kl_divergence_functions.R")
###########input###############
if(pre_tag=="dca"){
  t_mean=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
  t_dispersion=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/dispersion.tsv"),stringsAsFactors = FALSE,row.names = 1)
  t_dropout=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/dropout.tsv"),stringsAsFactors = FALSE,row.names = 1)
}
if(pre_tag=="scvi"){
  t_mean=t(read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/scvi_",file_tag,"/scvi_imputed_values_",file_tag,".txt"),stringsAsFactors = FALSE,sep=",")) 
  t_dispersion=t(read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/scvi_",file_tag,"/scvi_lnDispersion_",file_tag,".txt"),stringsAsFactors = FALSE,sep=",")) #originaly in log scale
  t_dropout=t(read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/scvi_",file_tag,"/scvi_dropout_rate_",file_tag,".txt"),stringsAsFactors = FALSE,sep=",")) #originaly in logit scale
  t_dispersion=exp(t_dispersion)
  t_dispersion=matrix(rep(t_dispersion,times=ncol(t_mean)),nrow=length(t_dispersion),ncol=ncol(t_mean))
  t_dropout=exp(t_dropout)/(1+exp(t_dropout))
}
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=read.table("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}




#adjust t_mean with total read depth per 1000 counts.
read_depth=readRDS(paste0("../Data_PRJNA434002/rawM_read_depth_per_1Kcell_ind.rds"))

#Ajust option: adjust means with read depth
#read_depth_median=median(read_depth)
#individual_index=match(tmeta$individual,rownames(read_depth))

#read_depth_weight=read_depth_median/read_depth[individual_index]
#t_mean=t(t(t_mean)*read_depth_weight)

##########data processing############
#1. restrict the cluster
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
sub_mean=as.matrix(t_mean[,tmeta$cluster==cur_cluster])

sub_dispersion=as.matrix(t_dispersion[,tmeta$cluster==cur_cluster])
sub_dropout=as.matrix(t_dropout[,tmeta$cluster==cur_cluster])

meta=tmeta[tmeta$cluster==cur_cluster,]

#2. generate idividual label
cur_individual=unique(meta$individual)



###################calculation t#################################
print("start calculation: Part II: direct KLmean and JSD")

klmean_direct_array=array(dim=c(nrow(sub_mean),length(cur_individual),length(cur_individual)),
               dimnames = list(rownames(sub_mean),cur_individual,cur_individual))
jsd_direct_array=array(dim=c(nrow(sub_mean),length(cur_individual),length(cur_individual)),
                dimnames = list(rownames(sub_mean),cur_individual,cur_individual))

for(i_g in 1:nrow(sub_mean)){
  
  cell_param=cbind(sub_mean[i_g,], sub_dispersion[i_g,], sub_dropout[i_g,])
  klmean_direct_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter="mean"), error = function(e) {NA} )
  jsd_direct_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter="JSD"), error = function(e) {NA} )
  print(i_g)
}



###################calculation end, output#################################
print("calculation end")


saveRDS(klmean_direct_array,paste0("../Data_PRJNA434002/8.Result/klmean_direct_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(jsd_direct_array,paste0("../Data_PRJNA434002/8.Result/jsd_direct_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))


sessionInfo()
q(save="no")
