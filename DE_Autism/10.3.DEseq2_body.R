#this code care with the results of simulation data, and calculate the DESeq2

# file_tag=1
# r_mean=1.5
# r_var=1.5
# dp_minor_prop=0.3
# r_change_prop=0.8
perm_label=0
perm_method="" #c("","b") 
sim_folder="sim_v6"

n_seq=c(50,30,20,10,5)
ncell_seq=c(800,400,200,100,50,20)

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

library("DESeq2")

sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))


# We calculate bulk information by summing up raw counts of
# all cells(of certain cluster) of an individual within a genes

#inputs:
# meta file
# sim_matrix_bulk

cur_info = meta[, c("individual", "phenotype")]
cur_info = unique(cur_info)
rownames(cur_info) = cur_info$individual
cur_info$phenotype = as.factor(cur_info$phenotype)


#additional analysis
mean_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
var_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.var_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
dp_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.dp_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
mult_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))


for(n in n_seq){
  selected_index=sample.int(max(n_seq),n)
  sub_sim_matrix_bulk=sim_matrix_bulk[,c(selected_index,(max(n_seq)+selected_index))]
  sub_cur_info=cur_info[c(selected_index,(max(n_seq)+selected_index)),]
  
  if(perm_label>0){
    if(perm_method=="b"){
      
      n_exchange=n/2
      n_exchange=floor(n_exchange)+(n_exchange-floor(n_exchange))*2*rbinom(1,1,0.5) #the number changed to other side
      
      i_exchange=sample.int(n,n_exchange)

      sub_cur_info$phenotype[(i_exchange)]=1
      sub_cur_info$phenotype[(i_exchange+n)]=0
      
    }
    if(perm_method!="b"){
      sub_cur_info$phenotype=sub_cur_info$phenotype[sample.int(2*n)]
    }
  }
  
  # object construction
  dds = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk,
                               colData = sub_cur_info,
                               design = ~ phenotype)
  
  # observed pvalue calculation
  dds = DESeq(dds)
  deseq_pval = results(dds)$pvalue
  saveRDS(deseq_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),".rds"))
  
  
  #additional analysis
  dds_mean = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk[which(mean_index==1),],colData = sub_cur_info,design = ~ phenotype)
  dds_var  = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk[which(var_index==1),], colData = sub_cur_info,design = ~ phenotype)
  dds_dp   = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk[which(dp_index==1),],  colData = sub_cur_info,design = ~ phenotype)
  dds_mult = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk[which(mult_index==1),],colData = sub_cur_info,design = ~ phenotype)
  dds_none = DESeqDataSetFromMatrix(countData = sub_sim_matrix_bulk[which(mean_index+var_index+dp_index+mult_index==0),],colData = sub_cur_info,design = ~ phenotype)
  dds_mean = DESeq(dds_mean)
  dds_var = DESeq(dds_var)
  dds_dp = DESeq(dds_dp)
  dds_mult = DESeq(dds_mult)
  dds_none = DESeq(dds_none)
  
  deseq_mean_pval = results(dds_mean)$pvalue
  deseq_var_pval  = results(dds_var)$pvalue
  deseq_dp_pval   = results(dds_dp)$pvalue
  deseq_mult_pval = results(dds_mult)$pvalue
  deseq_none_pval = results(dds_none)$pvalue
  
  saveRDS(deseq_mean_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_mean.rds"))
  saveRDS(deseq_var_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_var.rds"))
  saveRDS(deseq_dp_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_dp.rds"))
  saveRDS(deseq_mult_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_mult.rds"))
  saveRDS(deseq_none_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_none.rds"))
  
  deseq_sep_pval=matrix(1,nrow=length(deseq_pval),ncol=1)
  names(deseq_sep_pval)=names(deseq_pval)
  deseq_sep_pval[which(mean_index==1)]=results(dds_mean)$pvalue
  deseq_sep_pval[which(var_index==1)]=results(dds_var)$pvalue
  deseq_sep_pval[which(dp_index==1)]=results(dds_dp)$pvalue
  deseq_sep_pval[which(mult_index==1)]=results(dds_mult)$pvalue
  deseq_sep_pval[which(mean_index+var_index+dp_index+mult_index==0)]=results(dds_none)$pvalue
  saveRDS(deseq_sep_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_sep.rds"))
}




