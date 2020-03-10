#this code calculate the DESeq2 on the simulated count data from DCA/scvi output

#cluster_tag=1
#file_tag="3k10"


library("DESeq2")
#perm_num=500
covariate_flag=NA #c(NA, "quantile99")
perm_label=1 #perm_label =0 means calculate the observed data other wise, permutated data
dataset_folder="MS"  #Data_PRJNA434002   MS

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

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
rawcount_data=readRDS(paste0("../",dataset_folder,"/rawM",file_tag,".rds"))
rawcount_data=rawcount_data[,tmeta$cluster==cur_cluster,drop=FALSE]

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

saveRDS(cell_num,paste0("../",dataset_folder,"/7.Result/rawcount_ind_cellnum_",cluster_tag,"_",file_tag,".rds"))
saveRDS(read_depth,paste0("../",dataset_folder,"/7.Result/rawcount_ind_readdepth_",cluster_tag,"_",file_tag,".rds"))
saveRDS(zero_rate_ind,paste0("../",dataset_folder,"/7.Result/rawcount_ind_zero_rate_",cluster_tag,"_",file_tag,".rds"))
saveRDS(rawcount_matrix_bulk,paste0("../",dataset_folder,"/7.Result/rawcount_matrix_bulk_",cluster_tag,"_",file_tag,".rds"))

zero_rate=apply(rawcount_data==0,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})/(length(rawcount_data[1,]))
read_count_total=apply(rawcount_matrix_bulk,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})
saveRDS(zero_rate,paste0("../",dataset_folder,"/7.Result/rawcount_gene_zero_rate_",cluster_tag,"_",file_tag,".rds"))
saveRDS(read_count_total,paste0("../",dataset_folder,"/7.Result/rawcount_gene_read_count_total_",cluster_tag,"_",file_tag,".rds"))

cur_info=meta[,c("individual","diagnosis")]
cur_info=unique(cur_info)
rownames(cur_info)=cur_info$individual

###################calculation t#################################
print("start DESeq2 calculation")

if(perm_label>0){
  perm_order=readRDS(paste0("../",dataset_folder,"/7.Result/ind_perm_order.rds"))
  perm_order=as.numeric(perm_order[,perm_label])
  total_individual_ref=total_individual[perm_order]
  perm_order=match(total_individual_ref,cur_info$individual)
  perm_order=perm_order[!is.na(perm_order)]
  cur_info[,"diagnosis"]=cur_info[perm_order,"diagnosis"]
}

dds=DESeqDataSetFromMatrix(countData = rawcount_matrix_bulk,
                           colData = cur_info,
                           design = ~ diagnosis)

dds=DESeq(dds)
de_ob_pval=results(dds)$pvalue

saveRDS(de_ob_pval,paste0("../",dataset_folder,"/8.Result/DESeq2_pval/p",perm_label,"_DESeq2_ob_pval_",cluster_tag,"_",file_tag,".rds"))


sessionInfo()
q(save="no")
