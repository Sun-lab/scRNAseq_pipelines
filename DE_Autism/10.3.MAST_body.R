#this code care with the results of simulation data, and calculate the MAST

#note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.

#perm_tag=0
#file_tag=1
#sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
#splat.org--method 4, change the mean.shape and mean.rate originally
#zinb.naive--method 5, using naive zinb models to do so.

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#r_mean=1.5  #r_mean/r_var should < 1+mean.shape
#r_var=4

sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

######################Other method comparison: MAST #################################
print("start MAST calculation: Part II: ZINB KLmean and JSD")

library("moments")
sim_matrix_log = log2(1 + sim_matrix)


dim(sim_matrix_log)
sim_matrix_log[1:10,1:10]
cell_id=colnames(sim_matrix_log)
gene_id=rownames(sim_matrix_log)

diagnosis=as.character(meta$phenotype)
diagnosis[diagnosis==1]="Case"
diagnosis[diagnosis==0]="Control"

library("MAST")
library("lme4")

fData=data.frame(primerid=gene_id)
cData=data.frame(wellKey=cell_id)
colnames(meta)
length(fData)
length(cData)

sca=FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta$CDR) #from Chong and Paul
colData(sca)$diagnosis =as.factor(diagnosis)
colData(sca)$ind = as.factor(meta$individual)

colData(sca)


#do a permutation
if(perm_tag>0){
  #count cases and controls
  diag_info=paste0(colData(sca)$ind,":",colData(sca)$diagnosis)
  diag_kind=unique(diag_info)
  diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
  
  #permute
  diag_kind[,2]=diag_kind[sample.int(nrow(diag_kind),nrow(diag_kind),replace=F),2]
  
  #match back to each individuals
  ind_index=match(colData(sca)$ind,diag_kind[,1])
  colData(sca)$diagnosis=as.factor(diag_kind[ind_index,2])
}

b=zlm(formula=~diagnosis + ind + cngeneson, sca=sca)
bs=summary(b,logFC=TRUE,doLRT = paste0("diagnosis","Control"), level = 0.95)


bs$datatable


MAST_cur=bs$datatable
MAST_cur=MAST_cur[MAST_cur$contrast=="diagnosisControl",c("component","primerid","Pr(>Chisq)")]


#restriction for our analysis on model H
MAST_cur_pval=c(MAST_cur[MAST_cur$component=="H","Pr(>Chisq)"])
MAST_cur_pval=unlist(MAST_cur_pval)
saveRDS(MAST_cur_pval,paste0("../Data_PRJNA434002/10.Result/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",perm_tag,".rds"))


sessionInfo()
q(save="no")
