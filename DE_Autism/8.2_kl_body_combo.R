#this code calculate the KL divergence between each individual.

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with each other using the KL divergence

library("ggplot2")
library("emdbook")

#cluster_tag=1
#file_tag="3k10"
#pre_tag="dca" #c("dca","scvi")
fit_tag="" # ""(zinb) or "nb" 
fit_tag2="raw" # "" or "raw"(directely use dca res)
covariate_flag=NA #c(NA, "quantile99")
dataset_folder="Data_PRJNA434002"  #Data_PRJNA434002   MS

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########functions#############
source("./Command/8.0_kl_divergence_functions.R")
###########input###############
#input phenotype
if(length(grep("PFC",file_tag))>0){
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../",dataset_folder,"/meta_PFC.rds"))
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],"_PFC.rds"))
  }
}
if(length(grep("PFC",file_tag))==0){
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=read.table(paste0("../",dataset_folder,"/meta.tsv"),header = TRUE, sep = "\t")
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
  }
}

#name match for MS samples
colnames(tmeta)[grep("cell_type",names(tmeta))]="cluster"
colnames(tmeta)[grep("sample",names(tmeta))]="individual"

cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]

cur_individual=unique(meta$individual)

#input counts
if(!is.na(covariate_flag)){
  sim_fit=readRDS(paste0("../",dataset_folder,"/7.Result/fit_ind_",covariate_flag,"_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  # if(fit_tag=="raw"){
  #   sim_data=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/",dataset_folder,"/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
  #   dim(sim_data)=c(dim(sim_data),1)
  # }
  if(fit_tag2!="raw"){
    sim_data=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",covariate_flag,"_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  }
  
}
if(is.na(covariate_flag)){
  sim_fit=readRDS(paste0("../",dataset_folder,"/7.Result/fit_ind_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  if(fit_tag2=="raw"){
    sim_data=as.matrix(read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/",dataset_folder,"/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE))
    sim_data=array(sim_data,dim=c(dim(sim_data),1))
  }
  if(fit_tag2!="raw"){
    sim_data=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    
  }
}


###################calculation t#################################
print("start calculation: Part I: Empirical KLmean and JSD")

klmean_empirical_array=array(dim=c(nrow(sim_data),length(cur_individual),length(cur_individual)),
               dimnames = list(rownames(sim_data),cur_individual,cur_individual))
jsd_empirical_array=array(dim=c(nrow(sim_data),length(cur_individual),length(cur_individual)),
                dimnames = list(rownames(sim_data),cur_individual,cur_individual))

for(i_g in 1:nrow(sim_data)){
  cur_sim=sim_data[i_g,,,drop=FALSE]
  for(i_ind_a in 1:length(cur_individual)){
    for(i_ind_b in 1:length(cur_individual)){
      cur_ind_a=cur_individual[i_ind_a]
      cur_ind_b=cur_individual[i_ind_b]
      #fit sim
      cur_sim_ind_a=as.numeric(cur_sim[,meta$individual==cur_ind_a,])
      cur_sim_ind_b=as.numeric(cur_sim[,meta$individual==cur_ind_b,])
      klmean_empirical_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_sim_ind_a,cur_sim_ind_b,alter="mean",fit_model="empirical"), error = function(e) {NA} )
      jsd_empirical_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_sim_ind_a,cur_sim_ind_b,alter="JSD",fit_model="empirical"), error = function(e) {NA} )
    }
  }
  print(i_g)
}

print("start calculation: Part II: model fitted KLmean and JSD")

klmean_nbzinb_array=array(dim=c(nrow(sim_fit),length(cur_individual),length(cur_individual)),
                             dimnames = list(rownames(sim_fit),cur_individual,cur_individual))
jsd_nbzinb_array=array(dim=c(nrow(sim_fit),length(cur_individual),length(cur_individual)),
                          dimnames = list(rownames(sim_fit),cur_individual,cur_individual))

sim_fit[,,1]=exp(sim_fit[,,1]) #change the log mean to mean!!!

for(i_g in 1:nrow(sim_fit)){
 cur_fit=sim_fit[i_g,,]
  for(i_ind_a in 1:length(cur_individual)){
    for(i_ind_b in 1:length(cur_individual)){
      cur_a=cur_fit[i_ind_a,]
      cur_b=cur_fit[i_ind_b,]
      #kl and jsd
      klmean_nbzinb_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_a,cur_b,alter="mean",zinb.quantile=0.975,fit_model="zinb"), error = function(e) {NA} )
      jsd_nbzinb_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_a,cur_b,alter="JSD",zinb.quantile=0.975,fit_model="zinb"), error = function(e) {NA} )
    }
  }
  print(i_g)
}


###################calculation end, output#################################
print("calculation end")

if(fit_tag2=="raw"){
  fit_tag=fit_tag2
}

if(!is.na(covariate_flag)){
  saveRDS(klmean_empirical_array,paste0("../",dataset_folder,"/8.Result/klmean_empirical_array_",covariate_flag,"_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(jsd_empirical_array,paste0("../",dataset_folder,"/8.Result/jsd_empirical_array_",covariate_flag,"_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(klmean_nbzinb_array,paste0("../",dataset_folder,"/8.Result/klmean_nbzinb_array_",covariate_flag,"_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(jsd_nbzinb_array,paste0("../",dataset_folder,"/8.Result/jsd_nbzinb_array_",covariate_flag,"_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(covariate_flag)){
  saveRDS(klmean_empirical_array,paste0("../",dataset_folder,"/8.Result/klmean_empirical_array_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(jsd_empirical_array,paste0("../",dataset_folder,"/8.Result/jsd_empirical_array_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(klmean_nbzinb_array,paste0("../",dataset_folder,"/8.Result/klmean_nbzinb_array_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(jsd_nbzinb_array,paste0("../",dataset_folder,"/8.Result/jsd_nbzinb_array_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}

sessionInfo()
q(save="no")
