#this code calculate the MAST on the simulated count data from DCA/scvi output

#cluster_tag=1
#file_tag="3k10"
#pre_tag="dca" #c("dca","scvi")
#perm_label=1 #perm_label =0 means calculate the observed data other wise, permutated data


perm_num=500
sim_n=10
covariate_flag=NA #c(NA, "quantile99")
dataset_folder="MS"  #Data_PRJNA434002   MS

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

library("pryr")

########### Functions ###############
check_cal=function(){
  print(date())
  print(gc())
  print(mem_used())
}

###########input###############
#input diagnosis
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=read.table(paste0("../",dataset_folder,"/meta.tsv"),header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}
#name match for MS samples
colnames(tmeta)[grep("cell_type",names(tmeta))]="cluster"
colnames(tmeta)[grep("sample",names(tmeta))]="individual"

total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)


#input counts
if(!is.na(covariate_flag)){
  sim_data=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(covariate_flag)){
  sim_data=readRDS(paste0("../",dataset_folder,"/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
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
colData(sca)$cngeneson = as.numeric(CDR) 
colData(sca)$diagnosis =as.factor(meta$diagnosis)
colData(sca)$ind = as.factor(meta$individual)
colData(sca)$riboPercent = as.numeric(meta$RNA.ribosomal.percent)

if(perm_label>0){
  #count cases and controls
  diag_info=paste0(colData(sca)$ind,":",colData(sca)$diagnosis)
  diag_kind=unique(diag_info)
  diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
  
  #permute
  perm_order=readRDS(paste0("../",dataset_folder,"/7.Result/ind_perm_order.rds"))
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
check_cal()
b0 = zlm(formula = ~ diagnosis + cngeneson + age + sex + RIN + PMI + Capbatch + riboPercent, sca = sca, parallel = TRUE)
print(paste0("print system details, before b1"))
check_cal()
b1 = zlm(formula = ~ diagnosis + ( 1 | ind )+ cngeneson + age + sex + RIN + PMI + Capbatch + riboPercent, sca = sca, method = 'glmer',  ebayes = FALSE, parallel = TRUE)
print(paste0("print system details, after b1"))
check_cal()

b0
b1
print(paste0("print system details, before lrt0"))
check_cal()
lrt0 = lrTest(b0, "diagnosis")
print(paste0("print system details, before lrt1"))
check_cal()
lrt1 = lrTest(b1, "diagnosis")
print(paste0("print system details, after lrt1"))
check_cal()
dim(lrt1)
lrt1[1,,]

MAST_pval0 = apply(lrt0, 1, function(x){x[3,3]})
length(MAST_pval0)
MAST_pval0[1:4]

MAST_pval1 = apply(lrt1, 1, function(x){x[3,3]})
length(MAST_pval1)
MAST_pval1[1:4]

saveRDS(MAST_pval0,paste0("../",dataset_folder,"/8.Result/MAST_pval/p",perm_label,"_MAST_pval0_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(MAST_pval1,paste0("../",dataset_folder,"/8.Result/MAST_pval/p",perm_label,"_MAST_pval1_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

sessionInfo()
q(save="no")
