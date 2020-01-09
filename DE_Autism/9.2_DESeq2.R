#this code calculate the DESeq2 on the simulated count data from DCA/scvi output

#cluster_tag=1
#file_tag="3k10"
#pre_tag="dca" #c("dca","scvi")
#perm_label=1 #perm_label =0 means calculate the observed data other wise, permutated data

library("DESeq2")
#perm_num=500
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


###################calculation t#################################
print("start calculation: Part I: Bulk data preparation")

#individual level info
cell_num=matrix(ncol=1,nrow=length(cur_individual))
rownames(cell_num)=cur_individual
colnames(cell_num)="cell_num"
read_depth=matrix(ncol=1,nrow=length(cur_individual))
rownames(read_depth)=cur_individual
colnames(read_depth)="read_depth"

zero_rate_ind=matrix(nrow=dim(sim_data)[1],ncol=length(cur_individual))
rownames(zero_rate_ind)=dimnames(sim_data)[[1]]
colnames(zero_rate_ind)=cur_individual
sim_matrix_bulk=matrix(nrow=dim(sim_data)[1],ncol=length(cur_individual))
rownames(sim_matrix_bulk)=dimnames(sim_data)[[1]]
colnames(sim_matrix_bulk)=cur_individual


for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=sim_data[,meta$individual==cur_ind,,drop=FALSE]
  cell_num[i_ind]=dim(cur_ind_m)[2]
  read_depth[i_ind]=sum(as.numeric(cur_ind_m),na.rm = TRUE)/cell_num[i_ind]*1000
  
  zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})/(cell_num[i_ind]*10)
  sim_matrix_bulk[,i_ind]=apply(cur_ind_m,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})
}

saveRDS(cell_num,paste0("../Data_PRJNA434002/7.Result/sim_ind_cellnum_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(read_depth,paste0("../Data_PRJNA434002/7.Result/sim_ind_readdepth_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(zero_rate_ind,paste0("../Data_PRJNA434002/7.Result/sim_ind_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(sim_matrix_bulk,paste0("../Data_PRJNA434002/7.Result/sim_matrix_bulk_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

zero_rate=apply(sim_data==0,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})/(length(sim_data[1,,]))
read_count_total=apply(sim_matrix_bulk,1,function(x){return(sum(as.numeric(x),na.rm = TRUE))})
saveRDS(zero_rate,paste0("../Data_PRJNA434002/7.Result/sim_gene_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(read_count_total,paste0("../Data_PRJNA434002/7.Result/sim_gene_read_count_total_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

cur_info=meta[,c("individual","diagnosis")]
cur_info=unique(cur_info)
rownames(cur_info)=cur_info$individual

###################calculation t#################################
print("start DESeq2 calculation")

if(perm_label>0){
  perm_order=readRDS(paste0("../Data_PRJNA434002/7.Result/ind_perm_order.rds"))
  perm_order=as.numeric(perm_order[,perm_label])
  cur_info[,"diagnosis"]=cur_info[perm_order,"diagnosis"]
}

dds=DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                           colData = cur_info,
                           design = ~ diagnosis)

dds=DESeq(dds)
de_ob_pval=results(dds)$pvalue

saveRDS(de_ob_pval,paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_DESeq2_ob_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))


sessionInfo()
q(save="no")
