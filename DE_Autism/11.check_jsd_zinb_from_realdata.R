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
pre_tag_seq=c("dca","scvi")
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

jsd_zinb_pval=NA
jsd_empirical_pval=NA
jsd_direct_pval=NA
klmean_zinb_pval=NA
klmean_empirical_pval=NA
klmean_direct_pval=NA
deseq2_pval=NA
MAST_pval=NA

jsd_zinb_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_nbzinb_pval/p",perm_label,perm_method,"_jsd_nbzinb_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
jsd_empirical_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_empirical_pval/p",perm_label,perm_method,"_jsd_empirical_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
jsd_direct_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_direct_pval/p",perm_label,perm_method,"_jsd_direct_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
klmean_direct_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_direct_pval/p",perm_label,perm_method,"_klmean_direct_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
klmean_zinb_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_nbzinb_pval/p",perm_label,perm_method,"_klmean_nbzinb_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
klmean_empirical_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_empirical_pval/p",perm_label,perm_method,"_klmean_empirical_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
deseq2_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_ob_pval_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )

MAST_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/MAST_pval/p",perm_label,perm_method,"_MAST_pval1_rawcount_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )

#histogram

op=par(mfrow = c(4, 2))
tryCatch({hist(deseq2_pval,main="pval of deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(MAST_pval,main="pval of MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_empirical_pval,main="pval of jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_empirical_pval,main="pval of klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_zinb_pval,main="pval of jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_zinb_pval,main="pval of klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_direct_pval,main="pval of jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_direct_pval,main="pval of klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

par(op)













i_d=2
i_fit=1
dist_method=dist_method_seq[i_d]
fit_method=fit_method_seq[i_fit]
cur_pval=jsd_empirical_pval


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

op=par(mfrow = c(2, 2))
for(i in 1:10){
  j=i+2000
  plot(sim_fit[i,,1],sim_fit[j,,1],main=paste0("param mean",i))
  lines(c(0,10),c(0,10),col="red")
  plot(log(sim_fit[i,,2]),log(sim_fit[j,,2]),main=paste0("param disp",i))
  lines(c(0,10),c(0,10),col="red")
  plot(sim_fit[i,,3],sim_fit[j,,3],main=paste0("param drop",i))
  lines(c(0,10),c(0,10),col="red")
  plot(rawcount_data[i,],rawcount_data[j,],main=paste0("count",i))
  lines(c(0,10),c(0,10),col="red")
}
par(op)



param12_quantile=t(apply(sim_fit[,,2]/sim_fit[,,1],1,function(x){quantile(x,na.rm = TRUE)}))
param12_var=t(apply(sim_fit[,,2]/sim_fit[,,1],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param12_mean=t(apply(sim_fit[,,2]/sim_fit[,,1],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

param1_quantile=t(apply(sim_fit[,,1],1,function(x){quantile(x,na.rm = TRUE)}))
param1_var=t(apply(sim_fit[,,1],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param1_mean=t(apply(sim_fit[,,1],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

param2_quantile=t(apply(sim_fit[,,2],1,function(x){quantile(x,na.rm = TRUE)}))
param2_var=t(apply(sim_fit[,,2],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param2_mean=t(apply(sim_fit[,,2],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

param3_quantile=t(apply(sim_fit[,,3],1,function(x){quantile(x,na.rm = TRUE)}))
param3_var=t(apply(sim_fit[,,3],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param3_mean=t(apply(sim_fit[,,3],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

op=par(mfrow = c(4, 2))
for(i in c(1,3,5)){
  plot(1:100,param1_quantile[1:100,i],ylim=c(0,5))
  plot(2001:2100,param1_quantile[2001:2100,i],ylim=c(0,5))
}
plot(1:100,param1_var[1:100],ylim=c(0,5))
plot(2001:2100,param1_var[2001:2100],ylim=c(0,5))

for(i in c(1,3,5)){
  plot(1:100,param2_quantile[1:100,i],ylim=c(0,40))
  plot(2001:2100,param2_quantile[2001:2100,i],ylim=c(0,40))
}
plot(1:100,param2_var[1:100],ylim=c(0,80))
plot(2001:2100,param2_var[2001:2100],ylim=c(0,80))

#plot(1:100,param2_quantile[1:100,5]-param2_quantile[1:100,1],ylim=c(0,60))
#plot(2001:2100,param2_quantile[2001:2100,5]-param2_quantile[2001:2100,1],ylim=c(0,60))

for(i in c(1,3,5)){
  plot(1:100,param3_quantile[1:100,i],ylim=c(0,1))
  plot(2001:2100,param3_quantile[2001:2100,i],ylim=c(0,1))
}
plot(1:100,param3_var[1:100],ylim=c(0,1))
plot(2001:2100,param3_var[2001:2100],ylim=c(0,1))

for(i in c(1,3,5)){
  plot(1:100,param12_quantile[1:100,i],ylim=c(0,300))
  plot(2001:2100,param12_quantile[2001:2100,i],ylim=c(0,300))
}
plot(1:100,param12_var[1:100],ylim=c(0,300))
plot(2001:2100,param12_var[2001:2100],ylim=c(0,300))

for(i in c(1,3,5)){
  plot(1:100,param2_quantile[1:100,5]-param2_quantile[1:100,1],ylim=c(0,10))
  plot(2001:2100,param2_quantile[2001:2100,5]-param2_quantile[2001:2100,1],ylim=c(0,10))
}
plot(1:100,param2_var[1:100],ylim=c(0,50))
plot(2001:2100,param2_var[2001:2100],ylim=c(0,50))
par(op)


plot(log(sort(param2_quantile[1:100,5])),log(sort(param2_quantile[2001:2100,5])))
lines(c(0,30),c(0,30))


op=par(mfrow = c(4, 2))
for(i in (order(param2_quantile[1:100,5],decreasing = TRUE)[1:8])){
  j=i+2000
  hist(sim_fit[i,,2])
  hist(sim_fit[j,,2])
}
par(op)

op=par(mfrow = c(2, 2))
for(i in 1:10){
  j=i+1500
  heatmap(dist_array[i,,])
  print(dist_array[i,,])
  print(quantile(dist_array[i,,]))
  hist(dist_array[i,,])
  hist(dist_array[j,,])
  plot(dist_array[i,,],dist_array[j,,])
  lines(c(0,10),c(0,10),col="red")
}
par(op)

max_count=apply(sim_matrix,1,max)
head(dist_array)
dist_quantile=t(apply(dist_array,1,function(x){quantile(x,na.rm = TRUE)}))
dist_var=t(apply(dist_array,1,function(x){var(as.numeric(x),na.rm = TRUE)}))
nlog10_dist_var=-log10(dist_var)
View(dist_quantile)

op=par(mfrow = c(4, 2))
for(i in c(1,3,5)){
  plot(1:100,dist_quantile[1:100,i],ylim=c(0,0.5))
  plot(2001:2100,dist_quantile[2001:2100,i],ylim=c(0,0.5))
}
plot(1:100,nlog10_dist_var[1:100],ylim=c(0,6))
plot(2001:2100,nlog10_dist_var[2001:2100],ylim=c(0,6))
par(op)


plot(sort(dist_quantile[1:100,5]),sort(dist_quantile[2001:2100,5]))
lines(c(0,30),c(0,30))


op=par(mfrow = c(2, 2))
#param1 max
plot(dist_quantile[1:100,3],param1_quantile[1:100,3],sub=paste0("cor=", round(cor(dist_quantile[1:100,3],param1_quantile[1:100,3]),3)))
plot(dist_quantile[1:1000,3],param1_quantile[1:1000,3],sub=paste0("cor=", round(cor(dist_quantile[1:1000,3],param1_quantile[1:1000,3]),3)))
plot(dist_quantile[1501:1600,3],param1_quantile[1501:1600,3],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,3],param1_quantile[1501:1600,3]),3)))
plot(dist_quantile[,3],param1_quantile[,3],sub=paste0("cor=", round(cor(dist_quantile[,3],param1_quantile[,3]),3)))

#param 1 median
plot(dist_quantile[1:100,5],param1_quantile[1:100,5],sub=paste0("cor=", round(cor(dist_quantile[1:100,5],param1_quantile[1:100,5]),3)))
plot(dist_quantile[1:1000,5],param1_quantile[1:1000,5],sub=paste0("cor=", round(cor(dist_quantile[1:1000,5],param1_quantile[1:1000,5]),3)))
plot(dist_quantile[1501:1600,5],param1_quantile[1501:1600,5],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,5],param1_quantile[1501:1600,5]),3)))
plot(dist_quantile[,5],param1_quantile[,5],sub=paste0("cor=", round(cor(dist_quantile[,5],param1_quantile[,5]),3)))


#param1 max
plot(dist_quantile[1:100,3],param2_quantile[1:100,3],sub=paste0("cor=", round(cor(dist_quantile[1:100,3],param2_quantile[1:100,3]),3)))
plot(dist_quantile[1:1000,3],param2_quantile[1:1000,3],sub=paste0("cor=", round(cor(dist_quantile[1:1000,3],param2_quantile[1:1000,3]),3)))
plot(dist_quantile[1501:1600,3],param2_quantile[1501:1600,3],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,3],param2_quantile[1501:1600,3]),3)))
plot(dist_quantile[2001:3000,3],param2_quantile[2001:3000,3],sub=paste0("cor=", round(cor(dist_quantile[,3],param2_quantile[,3]),3)))

#param 1 median
plot(dist_quantile[1:100,5],param2_quantile[1:100,5],sub=paste0("cor=", round(cor(dist_quantile[1:100,5],param2_quantile[1:100,5]),3)))
plot(dist_quantile[1:1000,5],param2_quantile[1:1000,5],sub=paste0("cor=", round(cor(dist_quantile[1:1000,5],param2_quantile[1:1000,5]),3)))
plot(dist_quantile[1501:1600,5],param2_quantile[1501:1600,5],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,5],param2_quantile[1501:1600,5]),3)))
plot(dist_quantile[2001:3000,5],param2_quantile[2001:3000,5],sub=paste0("cor=", round(cor(dist_quantile[,5],param2_quantile[,5]),3)))

#param1 max
plot(dist_quantile[1:100,3],param12_quantile[1:100,3],sub=paste0("cor=", round(cor(dist_quantile[1:100,3],param12_quantile[1:100,3]),3)))
plot(dist_quantile[1:1000,3],param12_quantile[1:1000,3],sub=paste0("cor=", round(cor(dist_quantile[1:1000,3],param12_quantile[1:1000,3]),3)))
plot(dist_quantile[1501:1600,3],param12_quantile[1501:1600,3],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,3],param12_quantile[1501:1600,3]),3)))
plot(dist_quantile[2001:3000,3],param12_quantile[2001:3000,3],sub=paste0("cor=", round(cor(dist_quantile[,3],param12_quantile[,3]),3)))

#param 1 median
plot(dist_quantile[1:100,5],param12_quantile[1:100,5],sub=paste0("cor=", round(cor(dist_quantile[1:100,5],param12_quantile[1:100,5]),3)))
plot(dist_quantile[1:1000,5],param12_quantile[1:1000,5],sub=paste0("cor=", round(cor(dist_quantile[1:1000,5],param12_quantile[1:1000,5]),3)))
plot(dist_quantile[1501:1600,5],param12_quantile[1501:1600,5],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,5],param12_quantile[1501:1600,5]),3)))
plot(dist_quantile[2001:3000,5],param12_quantile[2001:3000,5],sub=paste0("cor=", round(cor(dist_quantile[,5],param12_quantile[,5]),3)))
par(op)

res=matrix(ncol=3,nrow=0)
res_case=matrix(ncol=3,nrow=0)
res_ctrl=matrix(ncol=3,nrow=0)
res_out=matrix(ncol=3,nrow=0)
zero_num=matrix(ncol=1,nrow=3000)
dist_diff_vector=zero_num
for(ig in 1:3000){
  dist_matrix=dist_array[ig,,]
  for(i in 1:(nrow(dist_matrix)-1)){
    for(j in (i+1):nrow(dist_matrix)){
      dist_range=dist_matrix[i,-c(i,j)]-dist_matrix[j,-c(i,j)]
      dist_diff=abs(max(dist_range,na.rm = TRUE)-min(dist_range,na.rm = TRUE))
      if(dist_diff<0.003){
        res=rbind(res,c(ig,i,j))
        if(phenotype[i]==0 && phenotype[j]==0){
          res_ctrl=rbind(res_ctrl,c(ig,i,j))
        }
        if(phenotype[i]==0 && phenotype[j]==0){
          res_case=rbind(res_case,c(ig,i,j))
        }
        if(phenotype[i]!=phenotype[j]){
          res_out=rbind(res_out,c(ig,i,j))
        }
        
        #   print(c(ig,i,j))
        #   print(quantile(dist_range))
      }
    }
  }
  zero_num[ig]=sum(is.na(dist_matrix))
  dist_diff_vector[ig]=dist_diff
}

tie_table=table(res[,1])
tie_table_case=table(res_case[,1])
tie_table_ctrl=table(res_ctrl[,1])
tie_table_out=table(res_out[,1])

op=par(mfrow = c(2, 3))

hist(tie_table_case)
hist(tie_table_ctrl)
hist(tie_table_out)

plot(as.numeric(names(tie_table_case)),tie_table_case,cex=.1)
plot(as.numeric(names(tie_table_ctrl)),tie_table_ctrl,cex=.1)
plot(as.numeric(names(tie_table_out)),tie_table_out,cex=.1)

hist(cur_pval2)
hist(cur_pval2[as.numeric(names(tie_table))])
hist(cur_pval2[-as.numeric(names(tie_table))])

hist(cur_pval2)
hist(cur_pval2[as.numeric(names(tie_table_case))])
hist(cur_pval2[-as.numeric(names(tie_table_case))])

hist(cur_pval2)
hist(cur_pval2[as.numeric(names(tie_table_ctrl))])
hist(cur_pval2[-as.numeric(names(tie_table_ctrl))])

hist(cur_pval2)
hist(cur_pval2[as.numeric(names(tie_table_out))])
hist(cur_pval2[-as.numeric(names(tie_table_out))])

hist(cur_pval2)
hist(cur_pval2[which(zero_num>0)])
hist(cur_pval2[-which(zero_num>0)])

hist(cur_pval2)
hist(cur_pval2[unique(as.numeric(names(tie_table)),which(zero_num>0))])
hist(cur_pval2[-unique(as.numeric(names(tie_table)),which(zero_num>0))])
par(op)















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

pdf(paste0("~/Desktop/11.check_simulation_",dist_method,"_",fit_method,"_",F_method,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 25,width = 25)

for(i_s in 1:20){
  dist_pval=NA
  tukey_centroid_pval=NA
  anova_centroid_pval=NA
  tukey_median_pval=NA
  anova_median_pval=NA
  
  set.seed(i_s)
  # part I homogenies
  cur_phenotype_ind=phenotype_ind[sample.int(length(phenotype_ind),length(phenotype_ind))]
  dist_res=cal_permanova_pval2(dist_array,cur_phenotype_ind,perm_num.min =1000)
  dist_pval=dist_res$pval
  

  
  
  
  #part II ties
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
  
  
  if(sum(dist_pval<0.1)>500){
  
  print(i)
  
  # pval=matrix(ncol=1,nrow=3000)
  # for(i in 1:3000){
  #   print(i)
  #   pval[i]=cal_betadisper_pval(dist_array[i,,],label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD")
  # }
  tukey_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
  anova_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )})
  tukey_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
  anova_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="anova"), error = function(e) {NA} )}) 
  
  
  op=par(mfrow = c(6, 5))
  tryCatch({hist(dist_pval[mean_index==1],main="pval of mean-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[var_index==1],main="pval of var-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[disp_index==1],main="pval of disp-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[mult_index==1],main="pval of mult-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(tukey_centroid_pval[mean_index==1],main="pval of mean-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_centroid_pval[var_index==1],main="pval of var-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_centroid_pval[disp_index==1],main="pval of disp-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_centroid_pval[mult_index==1],main="pval of mult-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(anova_centroid_pval[mean_index==1],main="pval of mean-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_centroid_pval[var_index==1],main="pval of var-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_centroid_pval[disp_index==1],main="pval of disp-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_centroid_pval[mult_index==1],main="pval of mult-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(tukey_median_pval[mean_index==1],main="pval of mean-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_median_pval[var_index==1],main="pval of var-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_median_pval[disp_index==1],main="pval of disp-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_median_pval[mult_index==1],main="pval of mult-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(tukey_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(anova_median_pval[mean_index==1],main="pval of mean-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_median_pval[var_index==1],main="pval of var-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_median_pval[disp_index==1],main="pval of disp-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_median_pval[mult_index==1],main="pval of mult-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(anova_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  
  
  
  
  tryCatch({plot(dist_pval[mean_index==1],tukey_centroid_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],tukey_centroid_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[var_index==1],tukey_centroid_pval[var_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,var-DE genes, cor ",round(cor(dist_pval[var_index==1],tukey_centroid_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[disp_index==1],tukey_centroid_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],tukey_centroid_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mult_index==1],tukey_centroid_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],tukey_centroid_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs tukey_centroid,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
  
  tryCatch({plot(dist_pval[mean_index==1],anova_centroid_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],anova_centroid_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[var_index==1],anova_centroid_pval[var_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,var-DE genes, cor ",round(cor(dist_pval[var_index==1],anova_centroid_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[disp_index==1],anova_centroid_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],anova_centroid_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mult_index==1],anova_centroid_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],anova_centroid_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs anova_centroid,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
  
  
  tryCatch({plot(dist_pval[mean_index==1],tukey_median_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],tukey_median_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[var_index==1],tukey_median_pval[var_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,var-DE genes, cor ",round(cor(dist_pval[var_index==1],tukey_median_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[disp_index==1],tukey_median_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],tukey_median_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mult_index==1],tukey_median_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],tukey_median_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs tukey_median,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
  
  tryCatch({plot(dist_pval[mean_index==1],anova_median_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs anova_median,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],anova_median_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[var_index==1],anova_median_pval[var_index==1],cex=.2, main=paste0("dist_pval vs anova_median,var-DE genes, cor ",round(cor(dist_pval[var_index==1],anova_median_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[disp_index==1],anova_median_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs anova_median,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],anova_median_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mult_index==1],anova_median_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs anova_median,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],anova_median_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
  tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs anova_median,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
  par(op)
  
  #Part excusion
  op=par(mfrow = c(5, 5))
  tryCatch({hist(dist_pval[mean_index==1],main="pval of mean-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[var_index==1],main="pval of var-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[disp_index==1],main="pval of disp-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[mult_index==1],main="pval of mult-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(dist_pval[anova_median_pval<=0.05 & mean_index==1],main="pval of mean-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval<=0.05 & var_index==1],main="pval of var-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval<=0.05 & disp_index==1],main="pval of disp-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval<=0.05 & mult_index==1],main="pval of mult-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval<=0.05 & mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  
  tryCatch({hist(dist_pval[anova_median_pval>0.05 & mean_index==1],main="pval of mean-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval>0.05 & var_index==1],main="pval of var-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval>0.05 & disp_index==1],main="pval of disp-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval>0.05 & mult_index==1],main="pval of mult-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[anova_median_pval>0.05 & mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  
  
  
  tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(mean_index==1)))],main="pval of mean-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(var_index==1)))],main="pval of var-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(disp_index==1)))],main="pval of disp-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(mult_index==1)))],main="pval of mult-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which( mean_index==0 & var_index==0 & disp_index==0 & mult_index==0)))],main="pval of non-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  
  tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(mean_index==0)))],main="pval of mean-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(var_index==0)))],main="pval of var-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(disp_index==0)))],main="pval of disp-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(mult_index==0)))],main="pval of mult-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(mean_index==1 | var_index==1 | disp_index==1 | mult_index==1)))],main="pval of non-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  par(op)
  
  dist_pval_m[,i_s]=dist_pval
  tukey_centroid_pval_m[,i_s]=tukey_centroid_pval
  anova_centroid_pval_m[,i_s]=anova_centroid_pval
  tukey_median_pval_m[,i_s]=tukey_median_pval
  anova_median_pval_m[,i_s]=anova_median_pval
  }
}
dev.off()


op=par(mfrow = c(2, 3))
hist(dist_pval)
hist(dist_pval[as.numeric(names(tie_table[tie_table>2]))])
hist(dist_pval[-as.numeric(names(tie_table[tie_table>2]))])
hist(dist_pval[-which(anova_median_pval<0.05) ])
hist(dist_pval[which(anova_median_pval<0.05) ])


