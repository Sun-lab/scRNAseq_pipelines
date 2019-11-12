#this code calculate the permutation based manova or permanovaS
#The input are the dist_array from step 8.

# cluster_tag=1
# file_tag="3k10"
# pre_tag="dca"
# dist_method="jsd"
# fit_method="nbzinb"
# F_method="p"


covariate_flag=NA #c(NA, "quantile99")


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########functions#############
source("./Command/9.0_Fstat_functions.R")
###########input###############

#input phenotype
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}

cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)
phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0

###################calculation t#################################
print("start calculation: Part I: Empirical KLmean and JSD")

if(!is.na(covariate_flag)){
  dist_array=readRDS(paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_array_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(covariate_flag)){
  dist_array=readRDS(paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}

dist_pval=cal_permanova_pval2(dist_array,phenotype,F_method=F_method)

second_index=which(is.na(dist_pval))
for(i2 in second_index){
  print(i2)
  x=dist_array[i2,,]
  flag=(!is.na(x[,1]))
  if(sum(flag)>0&& max(as.numeric(x),na.rm=TRUE)>0){
    dist_pval[i2]=cal_permanova_pval(x[flag,flag],phenotype[flag],F_method=F_method)
    }
}



if(!is.na(covariate_flag)){
  saveRDS(dist_pval,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(covariate_flag)){
  saveRDS(dist_pval,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}

sessionInfo()
q(save="no")
