#this code checks the permanova method

perm_label_seq=0:5
perm_method="" # "" or "s"(set cell threshold, any individuals with less cells would be removed)
cell_threshold=4

ind_covariate_flag="ind" 
perm_num=500
tol=0.2

source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/9.0_Fstat_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/8.0_kl_divergence_functions.R")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

cluster_tag_seq=1:17
file_tag_seq=c("3k","rdpadj5k")
pre_tag_seq="dca"
dist_method_seq=c("klmean","jsd")
fit_method_seq=c("empirical","nbzinb")
F_method_seq=c("p","ps")
perm_label_seq=1:5
perm_method=""

i_file=1
i_F=1
i_pre=1
i_perm_label=1
i_cluster=4


file_tag=file_tag_seq[i_file]
F_method=F_method_seq[i_F]
pre_tag=pre_tag_seq[i_pre]
perm_label=perm_label_seq[i_perm_label]
cluster_tag=cluster_tag_seq[i_cluster]


i_d=2
i_fit=1
dist_method=dist_method_seq[i_d]
fit_method=fit_method_seq[i_fit]


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

sim_fit=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
rawcount_data=readRDS(paste0("../Data_PRJNA434002/rawM",file_tag,".rds"))
rawcount_data=rawcount_data[,tmeta$cluster==cur_cluster,drop=FALSE]

dist_array=readRDS(paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)

if(perm_method=="s"){
  cell_num=table(meta$individual)
  cell_index=which(cell_num>=cell_threshold)
  cur_individual=cur_individual[cell_index]
  dist_array=dist_array[,cell_index,cell_index,drop=FALSE]
}

phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0
ob_phenotype=phenotype


if(perm_label>0){
  phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
  phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0
  
  perm_order=readRDS(paste0("../Data_PRJNA434002/7.Result/ind_perm_order.rds"))
  perm_order=as.numeric(perm_order[,perm_label])
  total_individual_ref=total_individual[perm_order]
  perm_order=match(total_individual_ref,cur_individual)
  perm_order=perm_order[!is.na(perm_order)]
  phenotype=phenotype[perm_order,drop=FALSE]
}

phenotype
sum(phenotype==ob_phenotype)


dist_array=dist_array[order(cur_pval),,]
sim_fit=sim_fit[order(cur_pval),,]
rawcount_data=rawcount_data[order(cur_pval),]
cur_pval2=cur_pval[order(cur_pval)]



calG_3a=function(m){#m is array dim1=n_genes dim2=dim3=n_individual
  m=m*m
  n=dim(m)[2]
  res=-1/2*apply(m,1,function(x){return(x+matrix(sum(x)/(n*n),n,n)-1/n *(matrix((rep(rowSums(x),times=n)+rep(colSums(x),each=n)),n,n)))})
  dim(res)=c(n,n,dim(m)[1])
  res
}
calH=function(x){
  x=as.matrix(x)
  x%*% solve(crossprod(x))%*%t(x)
}

G=calG_3a(dist_array[1:100,,])
H=calH(phenotype)
IH=diag(nrow(H))-H

GIG_dist=array(dimnames=dimnames(G),dim=dim(G))
for(i_g in 1:100){
  GIG_dist[,,i_g]=H%*%G[,,i_g]%*%H + IH%*%G[,,i_g]%*%IH -G[,,i_g]
}
hist(GIG_dist)
hist(G)
quantile(GIG_dist)
quantile(G)
quantile(abs(apply(G,3,max)))
quantile(abs(apply(GIG_dist,3,max)))
















#redo all the method

#levene test from Anderson 2005
library("ggplot2")
library("emdbook")
library("Rcpp")
library("Barycenter")

source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/7.0_ZINB_fit_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/8.0_kl_divergence_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/9.0_Fstat_functions.R")

library("vegan")


#dist_matrix=nxn symmetric matrix with diag =0
#label nx1 vector about categories.


cal_betadisper_pval=function(dist_matrix,label,disper_type=c("median","centroid"),div_method=c("anova","TukeyHSD")){
  flag=which(apply(dist_matrix,1,function(x){sum(is.na(x))==length(x)}))
  if(length(flag)>0){
    label=as.factor(label[-flag])
    dist_matrix=dist_matrix[-flag,-flag]
    dist_matrix[is.na(dist_matrix)]=median(dist_matrix)
  }
  
  dis=as.dist(dist_matrix)
  mod=betadisper(dis, label,type=disper_type)
  
  if(div_method=="anova"){
    pval=anova(mod)$"Pr(>F)"[1]
  }
  if(div_method=="TukeyHSD"){
    pval=TukeyHSD(mod)$group[4]
  }
  return(pval)
}

dist_pval_m=matrix(ncol=20,nrow=dim(dist_array)[1])
tukey_centroid_pval_m=matrix(ncol=20,nrow=dim(dist_array)[1])
anova_centroid_pval_m=matrix(ncol=20,nrow=dim(dist_array)[1])
tukey_median_pval_m=matrix(ncol=20,nrow=dim(dist_array)[1])
anova_median_pval_m=matrix(ncol=20,nrow=dim(dist_array)[1])

pdf(paste0("~/Desktop/11.check_simulation_",dist_method,"_",fit_method,"_",F_method,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 15,width = 25)


res=matrix(ncol=3,nrow=0)
for(ig in 1:3000){
  dist_matrix=dist_array[ig,,]
  for(i in 1:(nrow(dist_matrix)-1)){
    for(j in (i+1):nrow(dist_matrix)){
      dist_range=dist_matrix[i,-c(i,j)]-dist_matrix[j,-c(i,j)]
      dist_diff=abs(max(dist_range,na.rm = TRUE)-min(dist_range,na.rm = TRUE))
      if(dist_diff<0.00001){
        res=rbind(res,c(ig,i,j))
      }
    }
  }
}
tie_table=table(res[,1])


for(i_s in 1:20){
  dist_pval=NA
  tukey_centroid_pval=NA
  anova_centroid_pval=NA
  tukey_median_pval=NA
  anova_median_pval=NA
  
  set.seed(i_s)
  # part I homogenies
  cur_phenotype=phenotype[sample.int(length(phenotype),length(phenotype))]
  dist_res=cal_permanova_pval2(dist_array,cur_phenotype,perm_num.min =1000)
  dist_pval=dist_res$pval
  
  #part II ties

  
  if(sum(abs(0.5-dist_pval)<0.1)>500){
  
  print(i_s)
  
  # pval=matrix(ncol=1,nrow=3000)
  # for(i in 1:3000){
  #   print(i)
  #   pval[i]=cal_betadisper_pval(dist_array[i,,],label = cur_phenotype,disper_type = "centroid",div_method="TukeyHSD")
  # }
  tukey_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
  anova_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )})
  tukey_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
  anova_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype,disper_type = "median",div_method="anova"), error = function(e) {NA} )}) 
  
  
  op=par(mfrow = c(3, 5))
  tryCatch({hist(dist_pval,main="pval of dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_centroid_pval,main="pval of tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_centroid_pval,main="pval of anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_median_pval,main="pval of tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_median_pval,main="pval of anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({plot(dist_pval,tukey_centroid_pval,cex=.2, main=paste0("dist_pval vs tukey_centroid, cor |0.5-p| ",round(cor(abs(0.5-dist_pval),tukey_centroid_pval,use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval,anova_centroid_pval,cex=.2, main=paste0("dist_pval vs anova_centroid, cor |0.5-p| ",round(cor(abs(0.5-dist_pval),anova_centroid_pval,use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval,tukey_median_pval,cex=.2, main=paste0("dist_pval vs tukey_median, cor |0.5-p| ",round(cor(abs(0.5-dist_pval),tukey_median_pval,use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval,anova_median_pval,cex=.2, main=paste0("dist_pval vs anova_median, cor |0.5-p| ",round(cor(abs(0.5-dist_pval),anova_median_pval,use="complete.obs"),3)))}, error = function(e) {NA} )
 

  tryCatch({hist(dist_pval,main="pval of dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(dist_pval[anova_median_pval<=0.05],main="pval of anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(dist_pval[anova_median_pval>0.05],main="pval of anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(dist_pval[as.numeric(names(tie_table[tie_table>2]))],main="pval of tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

  tryCatch({hist(dist_pval[-as.numeric(names(tie_table[tie_table>2]))],main="pval of tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  par(op)
  
  dist_pval_m[,i_s]=dist_pval
  tukey_centroid_pval_m[,i_s]=tukey_centroid_pval
  anova_centroid_pval_m[,i_s]=anova_centroid_pval
  tukey_median_pval_m[,i_s]=tukey_median_pval
  anova_median_pval_m[,i_s]=anova_median_pval
  }
}
dev.off()



