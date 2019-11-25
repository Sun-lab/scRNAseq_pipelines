#this code fit raw count data with ZINB/NB for each individuals per gene.

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with the original cells of certain type,certain individual
# cluster_tag=1
# file_tag="3k10"
# pre_tag="dca" #c("dca","scvi")
# 
# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

sim_n=10
covariate_flag=NA #c(NA, "quantile99","quantile99_readdepth","readdepth")


library("ggplot2")
library("emdbook")

###########functions#############
source("./Command/7.0_ZINB_fit_functions.R")
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

#3. fit simulated samples

fit_ind_sub_sim=array(dim=c(nrow(sub_mean),length(cur_individual),3),
                      dimnames = list(rownames(sub_mean),cur_individual,c("logmean","dispersion","dropout_rate")))

sim_ind=array(dim=c(nrow(sub_mean),ncol(sub_mean),sim_n),
                      dimnames = list(rownames(sub_mean),colnames(sub_mean),1:sim_n))

if(!is.na(covariate_flag)){
  if("quantile99" %in% covariate_flag){
    quantile99=log(apply(sub_mean,2,function(x)return(quantile(x,0.99)+1)))
    covariate=as.matrix(quantile99)
  }
  if("readdepth" %in% covariate_flag){
    quantile99=log(apply(sub_mean,2,function(x)return(quantile(x,0.99)+1)))
    covariate=log(apply(sub_mean,2,function(x)return(sum(x,na.rm = TRUE))))
  }
  
  pdf(paste0("../Data_PRJNA434002/7.Result/hist_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 4,width = 6)
}
if(is.na(covariate_flag)){
  pdf(paste0("../Data_PRJNA434002/7.Result/hist_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 4,width = 6)
}

for(i_g in 1:nrow(sub_mean)){
  cur_sim=matrix(ncol=sim_n,nrow=ncol(sub_mean))
  for(i_s in 1:ncol(sub_mean)){
    cur_sim[i_s,]=emdbook::rzinbinom(sim_n,sub_mean[i_g,i_s], sub_dispersion[i_g,i_s], sub_dropout[i_g,i_s])
  }
  sim_ind[i_g,,]=cur_sim
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]

    #fit sim
    cur_sim_ind=as.numeric(cur_sim[meta$individual==cur_ind,])
    
    if(!is.na(covariate_flag)){
      cur_covariate=rep(covariate[meta$individual==cur_ind,],sim_n)
      fit_ind_sub_sim[i_g,i_ind,]=fit_nbzinb(cur_sim_ind,cur_covariate)
    }
    if(is.na(covariate_flag)){
      fit_ind_sub_sim[i_g,i_ind,]=fit_nbzinb(cur_sim_ind)
    }
    
    #plot
    if(i_g<=10 & i_ind<=5){
      cur_sim_ind=data.frame(cur_sim_ind)
      colnames(cur_sim_ind)="simulated_count"
      ggplot(cur_sim_ind, aes(x=simulated_count),stat="count") + geom_histogram(fill="lightblue")+
        labs(title=paste0("Histogram of ",pre_tag,", ",rownames(sub_mean)[i_g]," of ",cur_ind,", ",cur_cluster),x="Simulated Count", y = "Frequency")
      #+theme_classic()
    }
  }
  print(i_g)
}
dev.off()

if(!is.na(covariate_flag)){
  saveRDS(fit_ind_sub_sim,paste0("../Data_PRJNA434002/7.Result/fit_ind_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(sim_ind,paste0("../Data_PRJNA434002/7.Result/sim_ind_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(covariate_flag)){
  saveRDS(fit_ind_sub_sim,paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(sim_ind,paste0("../Data_PRJNA434002/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}



sessionInfo()
q(save="no")
