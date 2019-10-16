
#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with the original cells of certain type,certain individual

library("pscl")
library("emdbook")
library("MASS")
cluster_tag=1
sim_n=10

###########functions#############
#fit_zinb fits the zinb model(if 0 included in data) and nb model (if 0 is not included)

#input:
#input_numeric: the numerical vectors as observation
#input_x: the covariates.If only need to estimate the observation distribution, input_x=1(default)

#output:
#fit_zinb returns 3 numbers, the fitted log transformed mean, fitted dispersion(theta)

#package requirement:
#fit_zinb needs package pscl::zeroinfl

input_numeric=cur_org_ind
input_x=cur_logsum_count

fit_nbzinb=function(input_numeric,input_x=NA){
  if((max(input_numeric,na.rm = TRUE)-min(input_numeric,na.rm = TRUE)==0)){
    warning("all observation shows the same, returns the log-transformed values themselves or 0 as mean")
    if(max(input_numeric,na.rm = TRUE)==0){
      fit_total=c(0,NA,1)
    }
    if(max(input_numeric,na.rm = TRUE)>0){
      fit_total=c(log(max(input_numeric,na.rm = TRUE)),NA,0)
    }
  }
  else{
    if(min(input_numeric)==0){
      fm_zinb=NA
      if(sum(!is.na(input_x))==0){
        fm_zinb=tryCatch(pscl::zeroinfl(as.numeric(input_numeric) ~  1, dist = "negbin"), error = function(e) {NA} )
        #print("zinb0")
      }
      if(sum(!is.na(input_x))>0){
        fm_zinb=tryCatch(pscl::zeroinfl(as.numeric(input_numeric) ~  input_x, dist = "negbin"), error = function(e) {print(e);NA} )
        #print("zinb1")
      }
      if(sum(!is.na(fm_zinb))>0){
        fit_logitdropout=as.numeric(fm_zinb$coefficients$zero[1])
        fit_dropout=exp(fit_logitdropout)
        fit_dropout=fit_dropout/(1+fit_dropout)
        fit_logmean=as.numeric(fm_zinb$coefficients$count[1])
        fit_dispersion=as.numeric(fm_zinb$theta)
        fit_total=c(fit_logmean,fit_dispersion,fit_dropout)
      }
      else{
        fit_total=c(NA,NA,NA)
      }
    }
    else{
      fm_nb=NA
      if(sum(!is.na(input_x))==0){
        fm_nb=tryCatch(MASS::glm.nb(as.numeric(input_numeric) ~  1), error = function(e) {print(e);NA} )
        #print("nb0")
      }
      if(sum(!is.na(input_x))>0){
        fm_nb=tryCatch(MASS::glm.nb(as.numeric(input_numeric) ~  input_x), error = function(e) {print(e);NA} )
        #print("nb1")
      }
      if(sum(!is.na(fm_nb))>0){
        fit_dropout=0
        fit_logmean=as.numeric(fm_nb$coefficients[1])
        fit_dispersion=as.numeric(fm_nb$theta)
        fit_total=c(fit_logmean,fit_dispersion,fit_dropout)
      }
      else{
        fit_total=c(NA,NA,NA)
      }
      
    }
  }
  names(fit_total)=c("logmean","dispersion","dropout_rate")
  return(fit_total)
}

###########input###############


setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


trawM=readRDS("../Data_PRJNA434002/rawM3k10.rds")
# tdca_mean=read.table("../Data_PRJNA434002/res_dca_rawM3k10/mean.tsv",stringsAsFactors = FALSE)
# tdca_dispersion=read.table("../Data_PRJNA434002/res_dca_rawM3k10/dispersion.tsv",stringsAsFactors = FALSE,row.names = 1)
# tdca_dropout=read.table("../Data_PRJNA434002/res_dca_rawM3k10/dropout.tsv",stringsAsFactors = FALSE,row.names = 1)

tmeta=readRDS("../Data_PRJNA434002/meta10.rds")


##########data processing############
for(cluster_tag in 1:17){
  cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
  print(c(cluster_tag,sum(tmeta$cluster==cur_cluster)))
}


#1. restrict the cluster
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
# dca_mean=as.matrix(tdca_mean[,tmeta$cluster==cur_cluster])
# dca_dispersion=as.matrix(tdca_dispersion[,tmeta$cluster==cur_cluster])
# dca_dropout=as.matrix(tdca_dropout[,tmeta$cluster==cur_cluster])
rawM=trawM[,tmeta$cluster==cur_cluster]
meta=tmeta[tmeta$cluster==cur_cluster,]

#2.generate idividual label
cur_individual=unique(meta$individual)


#3. fit dca simulated samples
# fit_ind_dca_sim=array(dim=c(nrow(dca_mean),length(cur_individual),3),
#                       dimnames = list(rownames(dca_mean),cur_individual,c("logmean","dispersion","dropout_rate")))
# 
# for(i_g in 1:nrow(dca_mean)){
#   cur_sim=matrix(ncol=sim_n,nrow=ncol(dca_mean))
#   for(i_s in 1:ncol(dca_mean)){
#     cur_sim[i_s,]=emdbook::rzinbinom(sim_n,dca_mean[i_g,i_s], dca_dispersion[i_g,i_s], dca_dropout[i_g,i_s])
#   }
#   for(i_ind in 1:length(cur_individual)){
#     cur_ind=cur_individual[i_ind]
#    
#     #fit sim
#     cur_sim_ind=as.numeric(cur_sim[meta$individual==cur_ind,])
#     fit_ind_dca_sim[i_g,i_ind,]=fit_nbzinb(cur_sim_ind)
#   }
#   print(i_g)
# }
# 
# saveRDS(fit_ind_dca_sim,paste0("../Data_PRJNA434002/fit_ind_dca_sim_",cluster_tag,"_3k10.rds"))


#4.fit raw count data
logsum_count=log(apply(rawM,2,sum))


fit_ind_lognorm_org=array(dim=c(nrow(rawM),length(cur_individual),3),
                      dimnames = list(rownames(rawM),cur_individual,c("logmean","dispersion","dropout_rate")))

for(i_g in 1:nrow(rawM)){
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]
    #fit org
    cur_logsum_count=logsum_count[meta$individual==cur_ind]
    cur_org_ind=rawM[i_g,meta$individual==cur_ind]
    fit_ind_lognorm_org[i_g,i_ind,]=fit_nbzinb(cur_org_ind,cur_logsum_count)
  }
  print(i_g)
}

saveRDS(fit_ind_lognorm_org,paste0("../Data_PRJNA434002/fit_ind_lognorm_",cluster_tag,"_3k10.rds"))









