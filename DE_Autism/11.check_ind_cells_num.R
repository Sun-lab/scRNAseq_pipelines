
cluster_tag_seq=1:17
file_tag="3k10"


###########functions#############

setwd("/Volumes/SpecialPass/fh_data/Data_PRJNA434002/")

#input phenotype
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}
total_individual=unique(tmeta$individual)

phenotype=tmeta$diagnosis[match(total_individual,tmeta$individual)]

res=matrix(0,nrow=length(total_individual),ncol=length(cluster_tag_seq))
rownames(res)=paste0(phenotype,"_",total_individual)
colnames(res)=paste0(cluster_tag_seq,"_",unique(tmeta$cluster))
for(cluster_tag in cluster_tag_seq){
  cur_cell=readRDS(paste0("../Data_PRJNA434002/7.Result/rawcount_ind_cellnum_",cluster_tag,"_",file_tag,".rds"))
  res[match(rownames(cur_cell),total_individual),cluster_tag]=cur_cell
}

saveRDS(res,paste0("11.check_clusters_individual_cellnum_",file_tag,".rds"))
saveRDS(res,paste0("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/11.check_clusters_individual_cellnum_",file_tag,".rds"))
write.table(res,file=paste0("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/11.check_clusters_individual_cellnum_",file_tag,".csv"),sep="\t")
quantile(res)
