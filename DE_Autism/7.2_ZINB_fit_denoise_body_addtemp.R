
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

file_tag_seq=c("PFC3k","PFC5k")
covariate_flag_seq=c("readdepth","quantile99")
dataset_folder="Data_PRJNA434002" 
pre_tag="dca"
fit_tag=""

library("ggplot2")
library("emdbook")

###########functions#############
source("./Command/7.0_ZINB_fit_functions.R")
###########input###############

for(file_tag in file_tag_seq){
  for(covariate_flag in covariate_flag_seq){
    if(pre_tag=="dca"){
      t_mean=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/",dataset_folder,"/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
    }
    
    if(is.na(unlist(strsplit(file_tag,"k"))[2])){
      tmeta=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/",dataset_folder,"/meta.tsv"),header = TRUE, sep = "\t")
    }
    if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
      tmeta=readRDS(paste0("/fh/scratch/delete90/sun_w/mengqi/",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
    }
    if(!is.na(grep("PFC",file_tag))){
      tmeta=tmeta[tmeta$region=="PFC",]
    }
    
    ##########data processing############
    #1. restrict the cluster
    for(cluster_tag in 1:length((unique(tmeta$cluster)))){
      cur_cluster=as.character(unique(tmeta$cluster))[cluster_tag]
      sub_mean=as.matrix(t_mean[,tmeta$cluster==cur_cluster])
      meta=tmeta[tmeta$cluster==cur_cluster,]
      
      #2. generate idividual label
      cur_individual=unique(meta$individual)
      
      if("quantile99" %in% covariate_flag){ #only 1 covaraite allowed
        quantile99=apply(sub_mean,2,function(x)return(quantile(x,0.99)+1))
        covariate=as.matrix(quantile99)
        log_covariate=log(covariate)
        covariate_ratio=covariate/median(covariate)
      }
      if("readdepth" %in% covariate_flag){  #only 1 covaraite allowed
        covariate=apply(sub_mean,2,function(x)return(sum(x,na.rm = TRUE)))
        log_covariate=log(covariate)
        covariate_ratio=covariate/median(covariate)
      }
      
      saveRDS(covariate,paste0("../",dataset_folder,"/7.Result/covariate_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
      saveRDS(log_covariate,paste0("../",dataset_folder,"/7.Result/log_covariate_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
      saveRDS(covariate_ratio,paste0("../",dataset_folder,"/7.Result/covariate_ratio_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  }
  }
}





sessionInfo()
#q(save="no")
