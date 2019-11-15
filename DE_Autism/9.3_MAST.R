#this code calculate the MAST on the simulated count data from DCA/scvi output

#cluster_tag=1
#file_tag="3k10"
#pre_tag="dca" #c("dca","scvi")

library("DESeq2")
perm_num=500
sim_n=10
covariate_flag=NA #c(NA, "quantile99")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########input###############
#input diagnosis
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}

cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)


#input counts
if(!is.na(covariate_flag)){
  sim_data=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(covariate_flag)){
  sim_data=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}

sim_matrix=sim_data
dim(sim_matrix)=c(dim(sim_data)[1],dim(sim_data)[2]*dim(sim_data)[3])
meta_org=meta
meta=matrix(ncol=ncol(meta_org),nrow=0)
for(im in 1:sim_n){
  cur_meta=meta_org
  cur_meta$cell=paste0(cur_meta$cell,"_",im)
  meta=rbind(meta,cur_meta)
}
dim(meta)
###################calculation t#################################

#individual level info
CDR=matrix(ncol=1,nrow=nrow(meta))
colnames(CDR)="CDR"
CDR_ind=matrix(ncol=1,nrow=length(cur_individual))
colnames(CDR_ind)="CDR_ind"
  
  
for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=sim_matrix[,meta$individual==cur_ind]
  cur_CDR=sum(cur_ind_m>0,na.rm = TRUE)/sum(cur_ind_m>-1,na.rm = TRUE)
  CDR[meta$individual==cur_ind]=cur_CDR
  CDR_ind[i_ind]=cur_CDR
}

#individual level info
cur_individual=unique(meta$individual)
cell_num=matrix(ncol=1,nrow=length(cur_individual))
rownames(cell_num)=cur_individual
colnames(cell_num)="cell_num"
read_depth=matrix(ncol=1,nrow=length(cur_individual))
rownames(read_depth)=cur_individual
colnames(read_depth)="read_depth"


zero_rate_ind=matrix(nrow=nrow(sim_matrix),ncol=length(cur_individual))
rownames(zero_rate_ind)=rownames(sim_matrix)
colnames(zero_rate_ind)=cur_individual
sim_matrix_bulk=matrix(nrow=nrow(sim_matrix),ncol=length(cur_individual))
rownames(sim_matrix_bulk)=rownames(sim_matrix)
colnames(sim_matrix_bulk)=cur_individual

for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=sim_matrix[,meta$individual==cur_ind]
  cell_num[i_ind]=ncol(cur_ind_m)
  read_depth[i_ind]=sum(cur_ind_m,na.rm = TRUE)/cell_num[i_ind]*1000
  
  zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(x,na.rm = TRUE))})/cell_num[i_ind]
  sim_matrix_bulk[,i_ind]=apply(cur_ind_m,1,function(x){return(sum(x,na.rm = TRUE))})
}

hist(read_depth)

plot(read_depth,CDR_ind)
cor(read_depth,CDR_ind)




######################Other method comparison: MAST #################################
print("start MAST calculation: Part II: ZINB KLmean and JSD")

library("moments")
sim_matrix_log = log2(1 + sim_matrix)


dim(sim_matrix_log)
sim_matrix_log[1:10,1:10]
cell_id=meta$cell
gene_id=dimnames(sim_data)[[1]]
rownames(sim_matrix_log)=gene_id
colnames(sim_matrix_log)=cell_id


library("MAST")
library("lme4")

fData=data.frame(primerid=gene_id)
cData=data.frame(wellKey=cell_id)
colnames(meta)
length(fData)
length(cData)

sca=FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(CDR) #from Chong and Paul
colData(sca)$diagnosis =as.factor(meta$diagnosis)
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

b=zlm(formula=~diagnosis + (1|ind) + cngeneson, sca=sca, method='glmer', ebayes=FALSE,strictConvergence = FALSE)
bs=summary(b,logFC=TRUE,doLRT = paste0("diagnosis","Control"), level = 0.95)


bs$datatable


MAST_cur=bs$datatable
MAST_cur=MAST_cur[MAST_cur$contrast=="diagnosisControl",c("component","primerid","Pr(>Chisq)")]


#restriction for our analysis on model H
MAST_cur_pval=c(MAST_cur[MAST_cur$component=="H","Pr(>Chisq)"])
MAST_cur_pval=unlist(MAST_cur_pval)
saveRDS(MAST_cur_pval,paste0("../Data_PRJNA434002/8.Result/MAST_org_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,"_",perm_tag,".rds"))


sessionInfo()
q(save="no")
