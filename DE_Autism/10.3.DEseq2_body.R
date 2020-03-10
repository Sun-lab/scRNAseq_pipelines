#this code care with the results of simulation data, and calculate the DESeq2

# file_tag=1
# r_mean=1.5
# r_var=1.5
# r_disp=1.2
# r_change_prop=0.6
perm_label=0
perm_method="" #c("","b") 
sim_folder="sim_v5"

n_seq=c(50,30,20,10,5)
ncell_seq=c(200,100,50,20)

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
}




