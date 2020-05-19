#this code fit raw count data with ZINB/NB for each individuals per gene.

library("ggplot2")
#cluster_tag=1
#file_tag="3k10"
covariate_flag="readdepth" #c(NA, "quantile99","quantile99_readdepth","readdepth")
dataset_folder="Data_PRJNA434002"  #Data_PRJNA434002   MS

###########functions#############
source("./Command/7.0_ZINB_fit_functions.R")
###########input###############
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


trawM=readRDS(paste0("../",dataset_folder,"/rawM",file_tag,".rds"))

if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=read.table(paste0("../",dataset_folder,"/meta.tsv"),header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}

if(!is.na(grep("PFC",file_tag))){
  tmeta=tmeta[tmeta$region=="PFC",]
}

#name match for MS samples
if(dataset_folder=="MS"){
  colnames(tmeta)[grep("cell_type",names(tmeta))]="cluster"
  colnames(tmeta)[grep("sample",names(tmeta))]="individual"
}

##########data processing############
for(i_c in 1:17){
  cur_cluster=as.character(unique(tmeta$cluster)[i_c])
  print(c(i_c,sum(tmeta$cluster==cur_cluster)))
}

#1. restrict the cluster
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
rawM=trawM[,tmeta$cluster==cur_cluster]
meta=tmeta[tmeta$cluster==cur_cluster,]

#2.generate idividual label
cur_individual=unique(meta$individual)

#3.fit raw count data
#logsum_count=log(apply(rawM,2,sum))
if("quantile99" %in% covariate_flag){ #only 1 covaraite allowed
  quantile99=apply(rawM,2,function(x)return(quantile(x,0.99)+1))
  covariate=as.matrix(quantile99)
  log_covariate=log(covariate)
  covariate_ratio=covariate/median(covariate)
}
if("readdepth" %in% covariate_flag){  #only 1 covaraite allowed
  covariate=apply(rawM,2,function(x)return(sum(x,na.rm = TRUE)))
  log_covariate=log(covariate)
  covariate_ratio=covariate/median(covariate)
}



if(covariate_flag!=""){
  #for real data adj of empirical
  covariate_ratio=covariate/median(covariate,na.rm = TRUE)
  rawM_adj=round(matrix(as.numeric(rawM)*rep(covariate_ratio,each=nrow(rawM)),nrow=nrow(rawM),ncol=ncol(rawM)))
  rownames(rawM_adj)=rownames(rawM)
  colnames(rawM_adj)=colnames(rawM)
  
  #for real data resid of empirical
  rawM_resid=rawM
  rawM_resid[]=NA
  for(i_g in 1:nrow(rawM)){
    rawM_resid[i_g,]=lm(log(rawM[i_g,]+0.1)~log(covariate))$residuals
  }
}

saveRDS(rawM,paste0("../",dataset_folder,"/7.Result/rawM_",covariate_flag,cluster_tag,"_",file_tag,".rds"))
if(covariate_flag!=""){
  saveRDS(rawM_adj,paste0("../",dataset_folder,"/7.Result/rawM_adj",covariate_flag,cluster_tag,"_",file_tag,".rds"))
  saveRDS(rawM_resid,paste0("../",dataset_folder,"/7.Result/rawM_logresid",covariate_flag,cluster_tag,"_",file_tag,".rds"))
}

pdf(paste0("../",dataset_folder,"/7.Result/hist_raw_",covariate_flag,cluster_tag,"_",file_tag,".pdf"),height = 4,width = 6)

fit_ind_org=array(dim=c(nrow(rawM),length(cur_individual),3),
                      dimnames = list(rownames(rawM),cur_individual,c("logmean","dispersion","dropout_rate")))

for(i_g in 1:nrow(rawM)){
  #for individual level model fit
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]
    #fit org
    cur_org_ind=rawM[i_g,meta$individual==cur_ind]
    if(covariate_flag!=""){
      #for git
      cur_covariate=covariate[meta$individual==cur_ind]
      fit_ind_org[i_g,i_ind,]=fit_nbzinb(cur_org_ind,cur_covariate)
      

    }
    if(covariate_flag==""){
      fit_ind_org[i_g,i_ind,]=fit_nbzinb(cur_org_ind)
    }
    
    if(i_g<=10 & i_ind<=5){
      cur_org_ind=data.frame(cur_org_ind)
      colnames(cur_org_ind)="raw_count"
      ggplot(cur_org_ind, aes(x=raw_count),stat="count") + geom_histogram(fill="lightblue")+
        labs(title=paste0("Histogram of rawcount, ",rownames(rawM)[i_g]," of ",cur_ind,", ",cur_cluster),x="Count", y = "Frequency")
      #+theme_classic()
    }
    
  }
  print(i_g)
}
dev.off()


saveRDS(fit_ind_org,paste0("../",dataset_folder,"/7.Result/fit_ind_rawM_",covariate_flag,cluster_tag,"_",file_tag,".rds"))


sessionInfo()
q(save="no")