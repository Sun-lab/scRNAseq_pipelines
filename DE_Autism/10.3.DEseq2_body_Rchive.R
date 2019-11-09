#this code care with the results of simulation data, and calculate the DESeq2


#file_tag=1
#sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
#splat.org--method 4, change the mean.shape and mean.rate originally
#zinb.naive--method 5, using naive zinb models to do so.


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#r_mean=1.5  #r_mean/r_var should < 1+mean.shape
#r_var=4



perm_num=500





sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

cur_info=meta[,c("individual","phenotype")]
cur_info=unique(cur_info)
rownames(cur_info)=cur_info$individual
######################Other method comparison: DESeq2 #################################
print("start DESeq2 calculation")

perm_num=500

library("DESeq2")
dds=DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                           colData = cur_info,
                           design = ~ phenotype)

dds=DESeq(dds)
de_ob_pval=results(dds)$pvalue

de_perm_pval=matrix(ncol=perm_num,nrow=nrow(sim_matrix_bulk))
for(ip in 1:perm_num){
  respval=1
  cur_info$phenotype=cur_info$phenotype[sample.int(nrow(cur_info))]
  dds=DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                             colData = cur_info,
                             design = ~ phenotype)
  dds=DESeq(dds)
  respval=results(dds)$pvalue
  de_perm_pval[,ip]=respval
  print(ip)
}

de_pval=rowSums(de_perm_pval-de_ob_pval<=0)/perm_num



saveRDS(de_ob_pval,paste0("../Data_PRJNA434002/10.Result/DESeq2_ob_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
saveRDS(de_perm_pval,paste0("../Data_PRJNA434002/10.Result/DESeq2_perm_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
saveRDS(de_pval,paste0("../Data_PRJNA434002/10.Result/DESeq2_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

