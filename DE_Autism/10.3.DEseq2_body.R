#this code care with the results of simulation data, and calculate the DESeq2

# file_tag=1
# r_mean=1.5
# r_var=1.5
# r_disp=1.2
# r_change_prop=0.6
perm_label=1

n_seq=c(20,15,10,5)
ncell_seq=c(100,80,60,40,20)
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

library("DESeq2")

sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))


# We calculate bulk information by summing up raw counts of
# all cells(of certain cluster) of an individual within a genes

#inputs:
# meta file
# sim_matrix_bulk

cur_info = meta[, c("individual", "phenotype")]
cur_info = unique(cur_info)
rownames(cur_info) = cur_info$individual
cur_info$phenotype = as.factor(cur_info$phenotype)



for(n in n_seq){
  selected_index=sample.int(20,n)
  sub_sim_matrix_bulk=sim_matrix_bulk[,c(selected_index,(20+selected_index))]
  sub_cur_info=cur_info[c(selected_index,(20+selected_index)),]
  
  if(perm_label>0){
    sub_cur_info$phenotype=sub_cur_info$phenotype[sample.int(2*n)]
  }
  
  # object construction
  dds = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk,
                               colData = sub_cur_info,
                               design = ~ phenotype)
  
  # observed pvalue calculation
  dds = DESeq(dds)
  deseq_pval = results(dds)$pvalue
  saveRDS(deseq_pval,paste0("../Data_PRJNA434002/10.Result/p",perm_label,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),".rds"))
}




