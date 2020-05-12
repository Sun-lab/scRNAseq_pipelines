#this code fit raw count data with ZINB/NB for each individuals per gene.

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with the original cells of certain type,certain individual
# cluster_tag=1
file_tag_seq=c("PFC3k","PFC5k")
pre_tag="dca" #c("dca","scvi")
# 
# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

fit_tag="" # ""(zinb) or "nb"
sim_n=10
covariate_flag="readdepth" #c(NA, "quantile99","quantile99_readdepth","readdepth")
dataset_folder="Data_PRJNA434002"  #Data_PRJNA434002   MS

library("ggplot2")
library("emdbook")

###########functions#############
source("./Command/7.0_ZINB_fit_functions.R")
###########input###############

for(file_tag in file_tag_seq){
  for(cluster_tag in 1:17){
    covariate=readRDS(paste0("../",dataset_folder,"/7.Result/covariate_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    covariate_ratio=readRDS(paste0("../",dataset_folder,"/7.Result/covariate_ratio_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    log_covariate=log(covariate)
    
    sim_ind=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    sim_ind_adj=sim_ind #adj for covariate
    sim_ind_logresid=sim_ind #adj for covariate
    
    for(i_g in 1:dim(sim_ind)[1]){
      cur_sim=sim_ind[i_g,,]
      if(!is.na(covariate_flag)){
        cur_sim_adj=cur_sim
        cur_sim_logresid=cur_sim
        #resid
        for(i_sim in 1:sim_n){
          cur_sim_logresid[,i_sim]=lm(log(cur_sim_logresid[,i_sim]+1)~log_covariate)$residuals
        }
        #adj
        for(i_s in 1:dim(sim_ind)[2]){
          cur_sim_adj[i_s,]=cur_sim_adj[i_s,]/covariate_ratio[i_s] 
        }
        sim_ind_logresid[i_g,,]=cur_sim_logresid
        sim_ind_adj[i_g,,]=cur_sim_adj
      }
      if(i_g%%1000==0){
        print(i_g)
      }
    }
    
    saveRDS(sim_ind_adj,paste0("../",dataset_folder,"/7.Result/sim_ind_adj",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds")) #if covariate flag works, then adj can be defined from here
    saveRDS(sim_ind_logresid,paste0("../",dataset_folder,"/7.Result/sim_ind_logresid",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    print(i_g)
  }
}

sessionInfo()
q(save="no")
