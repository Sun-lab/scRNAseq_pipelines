#this code is from 9.2_DESeq2_rawcount.R
#we do it with only changes the input from wei's file, and label the name as "org".

cluster_tag=4
file_tag="org"

library("DESeq2")
#perm_num=500
covariate_flag=NA #c(NA, "quantile99")
perm_label=0 #perm_label =0 means calculate the observed data other wise, permutated data

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

###########input###############
#input diagnosis
#input phenotype

tmeta=read.table(paste0("/fh/fast/sun_w/mengqi/Data_PRJNA434002/meta.tsv"),header = TRUE, sep = "\t")
tmeta=tmeta[tmeta$region=="PFC",]

total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
cur_cluster

meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)

dim(meta)
meta[1:2,]

#input counts
rawcount_data=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data/ct_mtx/PFC_L2_3.rds"))
rawcount_data=rawcount_data[,meta$cluster==cur_cluster,drop=FALSE]
dim(rawcount_data)
rawcount_data[1:2,1:5]
###################calculation t#################################
print("start calculation: Part I: Bulk data preparation")

#individual level info
cell_num=matrix(ncol=1,nrow=length(cur_individual))
rownames(cell_num)=cur_individual
colnames(cell_num)="cell_num"
read_depth=matrix(ncol=1,nrow=length(cur_individual))
rownames(read_depth)=cur_individual
colnames(read_depth)="read_depth"

zero_rate_ind=matrix(nrow=dim(rawcount_data)[1],ncol=length(cur_individual))
rownames(zero_rate_ind)=dimnames(rawcount_data)[[1]]
colnames(zero_rate_ind)=cur_individual
rawcount_matrix_bulk=matrix(nrow=dim(rawcount_data)[1],ncol=length(cur_individual))
rownames(rawcount_matrix_bulk)=dimnames(rawcount_data)[[1]]
colnames(rawcount_matrix_bulk)=cur_individual


for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=rawcount_data[,meta$individual==cur_ind,drop=FALSE]
  cell_num[i_ind]=dim(cur_ind_m)[2]
  read_depth[i_ind]=sum(as.numeric(cur_ind_m),na.rm = TRUE)/cell_num[i_ind]*1000
  
  zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})/(cell_num[i_ind]*10)
  rawcount_matrix_bulk[,i_ind]=apply(cur_ind_m,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})
}

saveRDS(cell_num,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/rawcount_ind_cellnum_",cluster_tag,"_",file_tag,".rds"))
saveRDS(read_depth,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/rawcount_ind_readdepth_",cluster_tag,"_",file_tag,".rds"))
saveRDS(zero_rate_ind,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/rawcount_ind_zero_rate_",cluster_tag,"_",file_tag,".rds"))
saveRDS(rawcount_matrix_bulk,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/rawcount_matrix_bulk_",cluster_tag,"_",file_tag,".rds"))

zero_rate=apply(rawcount_data==0,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})/(length(rawcount_data[1,]))
read_count_total=apply(rawcount_matrix_bulk,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})
saveRDS(zero_rate,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/rawcount_gene_zero_rate_",cluster_tag,"_",file_tag,".rds"))
saveRDS(read_count_total,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/rawcount_gene_read_count_total_",cluster_tag,"_",file_tag,".rds"))

cur_info=meta[,c("individual","diagnosis","age","sex","Seqbatch","RNA.Integrity.Number")]
cur_info=unique(cur_info)
rownames(cur_info)=cur_info$individual

###################calculation t#################################
print("start DESeq2 calculation")

#version1 with wei's covariates
dds=NA
dds=DESeqDataSetFromMatrix(countData = rawcount_matrix_bulk,
                           colData = cur_info,
                           design = ~ age + sex + Seqbatch + RNA.Integrity.Number + diagnosis)

dds=DESeq(dds)
de_pval=results(dds)$pvalue

quantile(de_pval)
saveRDS(de_pval,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/p",perm_label,"_DESeq2_pval_cov_",cluster_tag,"_",file_tag,".rds"))

#version2: without any covariates
dds=NA
dds=DESeqDataSetFromMatrix(countData = rawcount_matrix_bulk,
                           colData = cur_info,
                           design = ~ diagnosis)

dds=DESeq(dds)
de_pval0=results(dds)$pvalue

quantile(de_pval0)
saveRDS(de_pval0,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/p",perm_label,"_DESeq2_pval_no_cov_",cluster_tag,"_",file_tag,".rds"))


##############permutated version############
set.seed(20)
perm_label=1

if(perm_label>0){
  perm_order=readRDS(paste0("/fh/fast/sun_w/mengqi/Data_PRJNA434002/7.Result/ind_perm_order.rds"))
  perm_order=as.numeric(perm_order[,perm_label])
  total_individual_ref=total_individual[perm_order]
  perm_order=match(total_individual_ref,cur_info$individual)
  perm_order=perm_order[!is.na(perm_order)]
  cur_info[,"diagnosis"]=cur_info[perm_order,"diagnosis"]
}

print("start DESeq2 calculation")

#version1 with wei's covariates
dds=NA
dds=DESeqDataSetFromMatrix(countData = rawcount_matrix_bulk,
                           colData = cur_info,
                           design = ~ age + sex + Seqbatch + RNA.Integrity.Number + diagnosis)

dds=DESeq(dds)
de_pval=results(dds)$pvalue

quantile(de_pval)
names(de_pval)=dimnames(rawcount_data)[[1]]
saveRDS(de_pval,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/p",perm_label,"_DESeq2_pval_cov_",cluster_tag,"_",file_tag,".rds"))

#version2: without any covariates
dds=NA
dds=DESeqDataSetFromMatrix(countData = rawcount_matrix_bulk,
                           colData = cur_info,
                           design = ~ diagnosis)

dds=DESeq(dds)

quantile(de_pval0)
de_pval0=results(dds)$pvalue
names(de_pval0)=dimnames(rawcount_data)[[1]]
saveRDS(de_pval0,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/p",perm_label,"_DESeq2_pval_no_cov_",cluster_tag,"_",file_tag,".rds"))

sessionInfo()
#q(save="no")
