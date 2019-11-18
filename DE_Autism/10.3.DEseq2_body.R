#this code care with the results of simulation data, and calculate the DESeq2

#file_tag=1
#sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
#splat.org--method 4, change the mean.shape and mean.rate originally
#zinb.naive--method 5, using naive zinb models to do so.
#r_mean=1.5  #r_mean/r_var should < 1+mean.shape
#r_var=4

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

library("DESeq2")

sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))


# We calculate bulk information by summing up raw counts of
# all cells(of certain cluster) of an individual within a genes

#inputs:
# meta file
# sim_matrix_bulk

cur_info = meta[, c("individual", "phenotype")]
cur_info = unique(cur_info)
rownames(cur_info) = cur_info$individual
cur_info$phenotype = as.factor(cur_info$phenotype)

# object construction
dds = DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                             colData = cur_info,
                             design = ~ phenotype)

# observed pvalue calculation
dds = DESeq(dds)
deseq_pval = results(dds)$pvalue

saveRDS(deseq_pval,paste0("../Data_PRJNA434002/10.Result/DESeq2_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

