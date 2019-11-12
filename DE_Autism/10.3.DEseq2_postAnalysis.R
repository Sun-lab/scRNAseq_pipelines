#this code care with the results of simulation data and DESeq2
library("ggplot2")

file_tag=1
sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
#splat.org--method 4, change the mean.shape and mean.rate originally
#zinb.naive--method 5, using naive zinb models to do so.


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


r_mean=1.2  #r_mean/r_var should < 1+mean.shape
r_var=1.5

perm_num=500

sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

cur_info=meta[,c("individual","phenotype")]
cur_info=unique(cur_info)
rownames(cur_info)=cur_info$individual
######################Other method comparison: DESeq2 #################################
print("start DESeq2 preparation")
#individual level info
cur_individual=unique(meta$individual)

sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
de_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/DESeq2_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

mean_index=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_de.mean_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
var_index=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_de.var_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

s_deMean_pval=de_pval[de_pval==0 & mean_index==1]
s_deMean_bulk=sim_matrix_bulk[de_pval==0 & mean_index==1,]

s_deVar_pval=de_pval[de_pval<0.02 & var_index==1]
s_deVar_bulk=sim_matrix_bulk[de_pval<0.02 & var_index==1,]

s_blank_pval=de_pval[de_pval==0 & mean_index==0 & var_index==0]
s_blank_bulk=sim_matrix_bulk[de_pval==0 & mean_index==0 & var_index==0,]


pdf(paste0("../Data_PRJNA434002/10.Result/counts_of_sig_DESeq2genes_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".pdf"),height = 4,width = 4)

for(ig in 1:5){
  meanDiff_sig_gene_count=s_deMean_bulk[ig,]
  cur_table=data.frame(cbind(meanDiff_sig_gene_count,cur_info))
  p=ggplot(cur_table, aes(x=phenotype, y=meanDiff_sig_gene_count,fill=phenotype)) + geom_boxplot()+ theme_classic()
  print(p)
  
  varDiff_sig_gene_count=s_deVar_bulk[ig,]
  cur_table=data.frame(cbind(varDiff_sig_gene_count,cur_info))
  p=ggplot(cur_table, aes(x=phenotype, y=varDiff_sig_gene_count,fill=phenotype)) + geom_boxplot()+ theme_classic()
  print(p)
  
  nonDE_sig_gene_count=s_blank_bulk[ig,]
  cur_table=data.frame(cbind(nonDE_sig_gene_count,cur_info))
  p=ggplot(cur_table, aes(x=phenotype, y=nonDE_sig_gene_count,fill=phenotype)) + geom_boxplot()+ theme_classic()
  print(p)
}

dev.off()
