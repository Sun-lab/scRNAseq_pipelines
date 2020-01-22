#this code calculate the MAST on the simulated count data from DCA/scvi output

#cluster_tag=1
#file_tag="3k10"


perm_num=500
sim_n=10
covariate_flag=NA #c(NA, "quantile99")
perm_label=1 #perm_label =0 means calculate the observed data other wise, permutated data

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Volumes/SpecialPass/fh_data/Data_PRJNA434002/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########input###############
#input diagnosis
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}
total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)


#input counts
rawcount_data=readRDS(paste0("../Data_PRJNA434002/rawM",file_tag,".rds"))
rawcount_data=rawcount_data[,tmeta$cluster==cur_cluster,drop=FALSE]

dim(rawcount_data)
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
  cur_ind_m=rawcount_data[,meta$individual==cur_ind]
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

for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=rawcount_data[,meta$individual==cur_ind,drop=FALSE]
  cell_num[i_ind]=ncol(cur_ind_m)
  read_depth[i_ind]=sum(cur_ind_m,na.rm = TRUE)/cell_num[i_ind]*1000
}

hist(read_depth)

plot(read_depth,CDR_ind)
cor(read_depth,CDR_ind)

######################Other method comparison: MAST #################################
print("start MAST calculation: Part II: ZINB KLmean and JSD")

library("moments")


if(length(grep("norm",file_tag))==0){
  rawcount_data_log = log2(1 + rawcount_data)
}
if(length(grep("norm",file_tag))>0){
  rawcount_data_log = rawcount_data
}

dim(rawcount_data_log)
rawcount_data_log[1:10,1:10]
cell_id=meta$cell
gene_id=dimnames(rawcount_data)[[1]]
rownames(rawcount_data_log)=gene_id
colnames(rawcount_data_log)=cell_id


library("MAST")
library("lme4")

fData=data.frame(primerid=gene_id)
cData=data.frame(wellKey=cell_id)
colnames(meta)
length(fData)
length(cData)

sca=FromMatrix(rawcount_data_log, cData, fData)
colData(sca)$cngeneson = as.numeric(CDR) 
colData(sca)$diagnosis =as.factor(meta$diagnosis)
colData(sca)$ind = as.factor(meta$individual)

if(perm_label>0){
  #count cases and controls
  diag_info=paste0(colData(sca)$ind,":",colData(sca)$diagnosis)
  diag_kind=unique(diag_info)
  diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
  
  #permute
  perm_order=readRDS(paste0("../Data_PRJNA434002/7.Result/ind_perm_order.rds"))
  perm_order=as.numeric(perm_order[,perm_label])
  total_individual_ref=total_individual[perm_order]
  perm_order=match(total_individual_ref,cur_individual)
  perm_order=perm_order[!is.na(perm_order)]
  diag_kind[,2]=diag_kind[perm_order,2]
  
  #match back to each individuals
  ind_index=match(colData(sca)$ind,diag_kind[,1])
  colData(sca)$diagnosis=as.factor(diag_kind[ind_index,2])
}

colData(sca)



print(paste0("print system details, before b0"))
date()
gc()
b0 = zlm(formula = ~ diagnosis, sca = sca, parallel = TRUE)
print(paste0("print system details, before b1"))
date()
gc()
b1 = zlm(formula = ~ diagnosis + ( 1 | ind ), sca = sca, method = 'glmer', 
         ebayes = FALSE, parallel = TRUE)
print(paste0("print system details, after b1"))
date()
gc()

b0
b1
print(paste0("print system details, before lrt0"))
date()
gc()
lrt0 = lrTest(b0, "diagnosis")
print(paste0("print system details, before lrt1"))
date()
gc()
lrt1 = lrTest(b1, "diagnosis")
print(paste0("print system details, after lrt1"))
date()
gc()
dim(lrt1)
lrt1[1,,]

MAST_pval0 = apply(lrt0, 1, function(x){x[3,3]})
length(MAST_pval0)
MAST_pval0[1:4]

MAST_pval1 = apply(lrt1, 1, function(x){x[3,3]})
length(MAST_pval1)
MAST_pval1[1:4]

saveRDS(MAST_pval0,paste0("../Data_PRJNA434002/8.Result/MAST_pval/p",perm_label,"_MAST_pval0_rawcount_",cluster_tag,"_",file_tag,".rds"))
saveRDS(MAST_pval1,paste0("../Data_PRJNA434002/8.Result/MAST_pval/p",perm_label,"_MAST_pval1_rawcount_",cluster_tag,"_",file_tag,".rds"))

sessionInfo()
q(save="no")
