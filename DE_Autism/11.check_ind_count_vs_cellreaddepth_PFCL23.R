#this code is modified from 7.2_fit_zinb.R
#this code is used for checking the relationship between cell level read depth and counts. 
# we use nb to fit  cell-level-count ~ log-read depth
# if the coefficient ~1, then we say the relathionship between cell-level-count ~ log-read depth are around 1.

#---------------------------log of 7.2-----------------------
#this code fit raw count data with ZINB/NB for each individuals per gene.

#step1: simulate 10 samples based on the parameters of distributions on each cell each gene.
#step2: compare the distribution with the original cells of certain type,certain individual
# cluster_tag=1
# file_tag="3k10"
# pre_tag="dca" #c("dca","scvi")
# 
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
seedtag=9995665
set.seed(seedtag)
file_tag="PFC3k"
cluster_tag=4 #L23
pre_tag="dca"
fit_tag=""

source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/7.0_ZINB_fit_functions.R")

t_mean=read.table(paste0("res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
t_dispersion=read.table(paste0("res_dca_rawM",file_tag,"/dispersion.tsv"),stringsAsFactors = FALSE,row.names = 1)
t_dropout=read.table(paste0("res_dca_rawM",file_tag,"/dropout.tsv"),stringsAsFactors = FALSE,row.names = 1)
tmeta=read.table(paste0("meta.tsv"),header = TRUE, sep = "\t")
tmeta=tmeta[tmeta$region=="PFC",]

sim_n=10
covariate_flag="readdepth" #c(NA, "quantile99","quantile99_readdepth","readdepth")
dataset_folder="Data_PRJNA434002"  #Data_PRJNA434002   MS

library("ggplot2")
library("emdbook")

##########data processing############
#1. restrict the cluster and top 200 genes
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
sub_mean=as.matrix(t_mean[1:200,tmeta$cluster==cur_cluster])

sub_dispersion=as.matrix(t_dispersion[1:200,tmeta$cluster==cur_cluster])
sub_dropout=as.matrix(t_dropout[1:200,tmeta$cluster==cur_cluster])

meta=tmeta[tmeta$cluster==cur_cluster,]

#2. generate idividual label
cur_individual=unique(meta$individual)

#3. fit simulated samples

fit_coefficient=array(dim=c(nrow(sub_mean),length(cur_individual),3),
                      dimnames = list(rownames(sub_mean),cur_individual,c("Intercept","covariate","dispersion")))

sim_ind=array(dim=c(nrow(sub_mean),ncol(sub_mean),sim_n),
                      dimnames = list(rownames(sub_mean),colnames(sub_mean),1:sim_n))
fit_coefficient_adj=fit_coefficient
sim_ind_adj=sim_ind #adj for covariate
sim_ind_logresid=sim_ind #adj for covariate
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

covariate_ratio=covariate/median(covariate)

#for(i_c in 1:dim(sim_ind)[2]){
#  sim_ind_adj[,i_c,]=sim_ind[,i_c,]/covariate_ratio[i_c]
#}

for(i_g in 1:nrow(sub_mean)){
  cur_sim=matrix(ncol=sim_n,nrow=ncol(sub_mean))
  cur_sim_adj=cur_sim
  cur_sim_logresid=cur_sim
  for(i_s in 1:ncol(sub_mean)){
    if(fit_tag==""){
      cur_sim[i_s,]=emdbook::rzinbinom(sim_n,sub_mean[i_g,i_s], sub_dispersion[i_g,i_s], sub_dropout[i_g,i_s])
    }
    if(fit_tag=="nb"){
      cur_sim[i_s,]=rnbinom(n=sim_n,mu=sub_mean[i_g,i_s], size=sub_dispersion[i_g,i_s])
    }
    cur_sim_adj[i_s,]=cur_sim[i_s,]/covariate_ratio[i_s]
  }
  #resid
  for(i_sim in 1:sim_n){
    cur_sim_logresid[,i_sim]=lm(log(cur_sim[,i_sim]+1)~log_covariate)$residuals
  }
  
  sim_ind[i_g,,]=cur_sim
  sim_ind_adj[i_g,,]=cur_sim_adj
  sim_ind_logresid[i_g,,]=cur_sim_logresid
  
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]
    
    #fit sim
    cur_sim_ind=as.numeric(cur_sim[meta$individual==cur_ind,])
    cur_covariate=rep(log_covariate[meta$individual==cur_ind],sim_n)
    cur_model=MASS::glm.nb(cur_sim_ind~cur_covariate)
    fit_coefficient[i_g,i_ind,]=c(cur_model$coefficients,cur_model$theta)
    #fit sim adj
    cur_sim_ind=as.numeric(cur_sim_adj[meta$individual==cur_ind,])
    cur_covariate=rep(log_covariate[meta$individual==cur_ind],sim_n)
    cur_model=MASS::glm.nb(cur_sim_ind~cur_covariate)
    fit_coefficient_adj[i_g,i_ind,]=c(cur_model$coefficients,cur_model$theta)
    
    #fit sim resid
    
    
  }
  print(i_g)
}

pdf("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/11.check_ind_count_vs_cellreaddepth_fit_coefficient_PFCL23.pdf",width=12,height=8)
op=par(mfrow = c(2, 3))
hist(fit_coefficient[,,2],main="fit NB intercept of sim_count") #if ratio=1
hist(fit_coefficient[,,1],main="fit NB beta log_covariate of sim_count") #intercept, showing gene level differences.
hist(fit_coefficient[,,3],breaks=100,main="fit NB dispersion(fitting quality) of sim_count") #fit NB dispersion, showing fittingquality, better not too big, not too small

hist(fit_coefficient_adj[,,2],main="fit NB intercept of sim_count_cellreaddepth_adj") #if ratio=1
hist(fit_coefficient_adj[,,1],main="fit NB beta log_covariate of sim_count_cellreaddepth_adj") #intercept, showing gene level differences.
hist(fit_coefficient_adj[,,3],breaks=100,main="fit NB dispersion(fitting quality) of sim_count_cellreaddepth_adj") #fit NB dispersion, showing fitting
par(op)
dev.off()

png("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/11.check_ind_count_vs_cellreaddepth_cor_scatter_PFCL23.png",width=16,height=24,units="in",res=300)
op=par(mfrow = c(6, 4))
for(i_g in 1:6){
  plot(log(sim_ind[i_g,,]+1),rep(log(covariate),sim_n),cex=.1,col=rgb(0,0,0,0.3),xlab="log sim_count",ylab="log covariate",main=paste0("count vs readdepth, log-transformed, gene ",i_g))
  #abline(0,1,col="red")
  plot(log(sim_ind_adj[i_g,,]+1),rep(log(covariate),sim_n),cex=.1,col=rgb(0,0,0,0.3),xlab="log sim_count_adj",ylab="log covariate",main=paste0("adj count vs readdepth, log-transformed, gene ",i_g))
  #abline(0,1,col="red")
  plot(log(sim_ind[i_g,,]+1),log(sim_ind_adj[i_g,,]+1),cex=.1,col=rgb(0,0,0,0.3),xlab="log sim_count",ylab="log sim_count_adj",main=paste0("count vs count_adj, log-transformed, gene ",i_g))
  abline(0,1,col="red")
  plot(log(sim_ind[i_g,,]+1),sim_ind_logresid[i_g,,],cex=.1,col=rgb(0,0,0,0.3),xlab="log sim_count",ylab="log sim_count_resid",main=paste0("count vs count_resid, log-transformed, gene ",i_g))
  abline(0,1,col="red")
}
par(op)
dev.off()



#conclusion: 
#from (fit_coefficient[,,3]) the model fitting quality is acceptable.
#from (fit_coefficient[,,1]) there does shows the differences between genes expression.
#from (fit_coefficient[,,2]) generally, the overall parameters is around 1 and it looks like perfect normal. Thus to adjust the cell level reading counts. simplely devide the counts with its cell level read reatio is acceptable. 

sessionInfo()

