#this code analysis the fit result based on the observation and permutation

#something from the header file
cluster_tag=7 #this tag indicate the clusters it can be choose in 1 to 17

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
perm_num=100

fit_ob=readRDS(paste0("../Data_PRJNA434002/zlms_3k10_",cluster_tag,"_0.rds"))$datatable
fit_ob=fit_ob[fit_ob$contrast=="diagnosisControl",c("component","primerid","Pr(>Chisq)")]


fit_perm=matrix(ncol=perm_num,nrow=nrow(fit_ob))
for(ip in 1:perm_num){
  if(file.exists(paste0("../Data_PRJNA434002/zlms_3k10_",cluster_tag,"_",ip,".rds"))){
    cur_fit=readRDS(paste0("../Data_PRJNA434002/zlms_3k10_",cluster_tag,"_",ip,".rds"))$datatable
    cur_fit=c(cur_fit[cur_fit$contrast=="diagnosisControl","Pr(>Chisq)"])
    fit_perm[,ip]=cur_fit
  }
}

fit_ob_h=fit_ob[fit_ob$component=="H",]
fit_perm_h=fit_perm[fit_ob$component=="H",]
fit_ob_h_pval=c(fit_ob_h[,"Pr(>Chisq)"])
zero_rate=sum(fit_ob_h_pval==1)/length(fit_ob_h_pval) #the proportion of non-expression data
zero_rate

pdf(paste0("histogram_observed_and_permutated_gene_pval_celltype_",cluster_tag,".pdf"),width=6,height=6)
for(i in sample.int(3000,50,replace=FALSE)){
  if(fit_ob_h_pval[i]!=1){
    hist(fit_perm_h[i,],breaks=20,main=paste0("permutated p-value distribution of gene#",i))
    abline(v=fit_ob_h_pval[i],col="red")
  }
}
rawp_h=rowSums(fit_perm_h-fit_ob_h_pval<=0)/perm_num
hist(fit_ob_h_pval[fit_ob_h_pval<1],breaks=20,main=paste0("absolute observed pval(pval<1),zero rate ",round(zero_rate,4)))
hist(rawp_h[rawp_h<1],breaks=20,main=paste0("permutated observed pval(pval<1),zero rate ",round(sum(rawp_h==1)/length(rawp_h),4)))
hist(fit_perm_h[fit_perm_h<1],breaks=20,main=paste0("absoluted permutated pval(pval<1),zero rate ",round(sum(fit_perm_h==1)/length(fit_perm_h),4)))
dev.off()
