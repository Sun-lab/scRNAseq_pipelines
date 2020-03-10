#this code do some stat analysis for the output of the simulation from 7.2. which is:

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.(7.2 did)
#step2: compare the distribution with the original cells of certain type,certain individual

library("emdbook")

#cluster_tag=1
file_tag="5k"
pre_tag="dca" #c("dca","scvi")

sim_n=10
covariate_flag=NA #c(NA, "quantile99")
dataset_folder="MS"  #Data_PRJNA434002   MS

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

###########functions#############
source("./Command/7.0_ZINB_fit_functions.R")
###########input###############


for(cluster_tag in 1:17){
  print(cluster_tag)
  sim_data=NA
  #input phenotype
  if(!is.na(covariate_flag) & file.exists(paste0("../",dataset_folder,"/7.Result/sim_ind_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))){
    sim_data=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  }
  if(is.na(covariate_flag) & file.exists(paste0("../",dataset_folder,"/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))){
    sim_data=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  }
  if(!is.na(sim_data)){
    
    if(is.na(unlist(strsplit(file_tag,"k"))[2])){
      tmeta=read.table("../",dataset_folder,"/meta.tsv",header = TRUE, sep = "\t")
      
    }
    if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
      tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
    }
    
    cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
    meta=tmeta[tmeta$cluster==cur_cluster,]
    
    cur_individual=unique(meta$individual)
    
    #input counts
    rawM=readRDS(paste0("../",dataset_folder,"/rawM",file_tag,".rds"))
    rawM=rawM[,tmeta$cluster==cur_cluster]
    
    
    #Use the simulated counts, calculate the percentage of 0 per gene across all the cells of one individual. 
    #Then compare the percentage of 0 in the observed data. Get a scatter plot across all the genes considered. 
    
    ob_zero_rate_ind=matrix(nrow=nrow(rawM),ncol=length(cur_individual))
    rownames(ob_zero_rate_ind)=rownames(rawM)
    colnames(ob_zero_rate_ind)=cur_individual
    
    sim_zero_rate_ind=matrix(nrow=nrow(rawM),ncol=length(cur_individual))
    rownames(sim_zero_rate_ind)=rownames(rawM)
    colnames(sim_zero_rate_ind)=cur_individual
    
    if(!is.na(covariate_flag)){
      pdf(paste0("../",dataset_folder,"/7.Result/zero_rate_scatter_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 4,width = 4)
    }
    if(is.na(covariate_flag)){
      pdf(paste0("../",dataset_folder,"/7.Result/zero_rate_scatter_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 4,width = 4)
    }
    
    for(i_ind in 1:length(cur_individual)){
      #cal raw count
      cur_ind=cur_individual[i_ind]
      cur_ind_m=rawM[,meta$individual==cur_ind,drop=FALSE]
      ob_zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(x,na.rm = TRUE))})/sum(meta$individual==cur_ind)
      
      #cal sim count
      cur_ind_sim_m=sim_data[,which(meta$individual==cur_ind),, drop=FALSE]
      sim_zero_rate_ind[,i_ind]=apply(cur_ind_sim_m==0,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})/(dim(sim_data)[3]*sum(meta$individual==cur_ind))
      plot(sort(ob_zero_rate_ind[,i_ind]),sort(sim_zero_rate_ind[,i_ind]),cex=.5,main=paste0("QQplot: zero_rate of individuals ",cur_ind,", ",file_tag))
      lines(c(0,1),c(0,1),col="red")
    }
    
    plot(sort(ob_zero_rate_ind),sort(sim_zero_rate_ind),cex=.5,main=paste0("QQplot:zero_rate of all individuals, ",file_tag))
    lines(c(0,1),c(0,1),col="red")
    
    dev.off()
  }
}





#sessionInfo()
#q(save="no")
