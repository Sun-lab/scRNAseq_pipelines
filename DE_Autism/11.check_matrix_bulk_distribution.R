setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002")

seedtag=3673163
set.seed(seedtag)
file_tag=1
r_mean=1.2
r_var=1.2
r_change_prop=0.6
dp_minor_prop=0.2

perm_label=0
sim_folder="sim_v6"

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

library("DESeq2")

#Input and Construction ################
sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))


cur_info = meta[, c("individual", "phenotype")]
cur_info = unique(cur_info)
rownames(cur_info) = cur_info$individual
cur_info$phenotype = as.factor(cur_info$phenotype)


#DESeq2 ########################
# object construction
dds = DESeqDataSetFromMatrix(countData = sim_matrix_bulk,
                             colData = cur_info,
                             design = ~ phenotype)
# observed pvalue calculation
dds = DESeq(dds)
deseq_pval = results(dds)$pvalue

#Plot ##########################

pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/fig_bulk_count_hist/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".pdf"),height =15 ,width=10)

op=par(mfrow = c(3, 2))
for(ig in 1:30){
  cur_case=sim_matrix_bulk[ig,cur_info$phenotype==1]
  cur_ctrl=sim_matrix_bulk[ig,cur_info$phenotype==0]
  hist(cur_case,col=rgb(1,0,0,0.5),main=paste0("bulk of gene ",ig),sub=paste0(
    "mean-var case: ",round(mean(cur_case),2),"-",round(var(cur_case),2),
    " mean-var ctrl: ",round(mean(cur_ctrl),2),"-",round(var(cur_ctrl),2)))
  legend("topright",c("case","ctrl"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),pch=15)
  hist(cur_ctrl,add=T,col=rgb(0,0,1,0.5))
  plot(density(cur_case),col=rgb(1,0,0,0.5),main=paste0("bulk of gene ",ig),sub=paste0(
    "mean-var case: ",round(mean(cur_case),2),"-",round(var(cur_case),2),
    " mean-var ctrl: ",round(mean(cur_ctrl),2),"-",round(var(cur_ctrl),2)))
  legend("topright",c("case","ctrl"),col=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)),pch=15)
  lines(density(cur_ctrl),col=rgb(0,0,1,0.5))
}
par(op)
dev.off()
