
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

cluster_tag_seq=1:17
file_tag_seq=c("3k10","5k")
pre_tag_seq=c("dca","scvi")
dist_method_seq=c("klmean","jsd")
fit_method_seq=c("empirical","nbzinb")
F_method_seq=c("p","ps")


file_tag_seq="5k"
pre_tag_seq="scvi"

perm_label_seq=c(0,1)
ind_covariate_flag=NA

perm_method=""
###################functions###################

#power calculation
cal_power=function(x,threshold){
  return(sum(x<=threshold,na.rm = TRUE)/length(x))
}
#ks test: compared with uniform(0,1) distribution
cal_ks=function(x,method="two.sided"){
  return(ks.test(1:length(x)/length(x),x,alternative = method)$p.value)
}
###############################################
power_array=array(dim=c(
  length(file_tag_seq),
  length(F_method_seq),
  length(pre_tag_seq),
  length(cluster_tag_seq),
  length(perm_label_seq),
  8),
  dimnames = list(
    file_tag_seq,
    F_method_seq,
    pre_tag_seq,
    cluster_tag_seq,
    perm_label_seq,
    c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")))

ks_array=power_array
cor_nonexpres_ind_array=power_array
cor_zerorate_ind_mean_array=power_array
cor_zero_rate_array=power_array
cor_expression_level_array=power_array
cor_overdisp_median_array=power_array
cor_overdisp_max_array=power_array
cor_dropout_median_array=power_array
cor_dropout_max_array=power_array
count=1
zeros=matrix(ncol=8,nrow=length(file_tag_seq)*length(F_method_seq)*length(pre_tag_seq)*length(cluster_tag_seq)*length(perm_label_seq))
rownames_zeros=matrix(ncol=1,nrow=length(file_tag_seq)*length(F_method_seq)*length(pre_tag_seq)*length(cluster_tag_seq)*length(perm_label_seq))
colnames(zeros)=c("deseq2_pval","MAST_pval","jsd_empirical_pval","klmean_empirical_pval","jsd_zinb_pval","klmean_zinb_pval","jsd_direct_pval","klmean_direct_pval")

pval_list=list()

for(i_file in 1:length(file_tag_seq)){
  file_tag=file_tag_seq[i_file]
  pval_length=as.numeric(unlist(strsplit(file_tag,"k"))[1])*1000
  pval_list[[file_tag]]=array(dim=c(
    length(F_method_seq),
    length(pre_tag_seq),
    length(cluster_tag_seq),
    length(perm_label_seq),
    8,pval_length),
    dimnames = list(
      F_method_seq,
      pre_tag_seq,
      cluster_tag_seq,
      perm_label_seq,
      c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct"),1:pval_length))
  
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      for(i_cluster in 1:length(cluster_tag_seq)){
        for(i_perm_label in 1:length(perm_label_seq)){
          file_tag=file_tag_seq[i_file]
          F_method=F_method_seq[i_F]
          pre_tag=pre_tag_seq[i_pre]
          perm_label=perm_label_seq[i_perm_label]
          cluster_tag=cluster_tag_seq[i_cluster]
          
          jsd_zinb_pval=NA
          jsd_empirical_pval=NA
          jsd_direct_pval=NA
          klmean_zinb_pval=NA
          klmean_empirical_pval=NA
          klmean_direct_pval=NA
          deseq2_pval=NA
          MAST_pval=NA
          
          jsd_zinb_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_nbzinb_pval/p",perm_label,perm_method,"_jsd_nbzinb_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          jsd_empirical_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_empirical_pval/p",perm_label,perm_method,"_jsd_empirical_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          jsd_direct_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_direct_pval/p",perm_label,perm_method,"_jsd_direct_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          klmean_direct_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_direct_pval/p",perm_label,perm_method,"_klmean_direct_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          klmean_zinb_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_nbzinb_pval/p",perm_label,perm_method,"_klmean_nbzinb_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          klmean_empirical_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_empirical_pval/p",perm_label,perm_method,"_klmean_empirical_",F_method,"_pval_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          deseq2_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_ob_pval_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )

          MAST_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/MAST_pval/p",perm_label,perm_method,"_MAST_pval1_rawcount_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,1,]=deseq2_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,2,]=MAST_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,3,]=jsd_empirical_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,4,]=klmean_empirical_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,5,]=jsd_zinb_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,6,]=klmean_zinb_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,7,]=jsd_direct_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,8,]=klmean_direct_pval
          
          zeros[count,]=c(sum(is.na(deseq2_pval)),sum(is.na(MAST_pval)), sum(is.na(jsd_empirical_pval)),sum(is.na(klmean_empirical_pval)),sum(is.na(jsd_zinb_pval)),sum(is.na(klmean_zinb_pval)),sum(is.na(jsd_direct_pval)),sum(is.na(klmean_direct_pval)))
          rownames_zeros[count]=paste0(perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag)
          count=count+1  
          
          #histogram
          png(paste0("../Data_PRJNA434002/8.Result/fig_final_power/final_power_p",perm_label,perm_method,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 1600,width = 800)
          op=par(mfrow = c(4, 2), pty = "s")
          tryCatch({hist(deseq2_pval,main="pval of deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(MAST_pval,main="pval of MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(jsd_empirical_pval,main="pval of jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(klmean_empirical_pval,main="pval of klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(jsd_zinb_pval,main="pval of jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(klmean_zinb_pval,main="pval of klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(jsd_direct_pval,main="pval of jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(klmean_direct_pval,main="pval of klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          
          par(op)
          
          dev.off()
          
          #record power
          power_matrix=matrix(nrow=8,ncol=1)
          names(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(power_matrix)="pval"
          power_matrix[1]=tryCatch({cal_power(deseq2_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[2]=tryCatch({cal_power(MAST_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[3]=tryCatch({cal_power(jsd_empirical_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[4]=tryCatch({cal_power(klmean_empirical_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[5]=tryCatch({cal_power(jsd_zinb_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[6]=tryCatch({cal_power(klmean_zinb_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[7]=tryCatch({cal_power(jsd_direct_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[8]=tryCatch({cal_power(klmean_direct_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=power_matrix
          
          #record ks test result
          ks_matrix=matrix(nrow=8,ncol=1)
          names(ks_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(ks_matrix)="pval"
          ks_matrix[1]=tryCatch({cal_ks(deseq2_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[2]=tryCatch({cal_ks(MAST_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[3]=tryCatch({cal_ks(jsd_empirical_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[4]=tryCatch({cal_ks(klmean_empirical_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[5]=tryCatch({cal_ks(jsd_zinb_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[6]=tryCatch({cal_ks(klmean_zinb_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[7]=tryCatch({cal_ks(jsd_direct_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[8]=tryCatch({cal_ks(klmean_direct_pval,method = "less")}, error = function(e) {NA} )
          ks_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=ks_matrix
          
          
          #record gene-based cor test result: zero rate ind and expression

          zero_rate_ind=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          zerorate_ind_mean=apply(zero_rate_ind,1,mean)
          nonexpres_ind=apply(zero_rate_ind==10,1,sum)

          zero_rate=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          expression_level=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_gene_read_count_total_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

          zinb_fit=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          overdisp_max=apply(zinb_fit[,,2],1,function(x){return(max(x,na.rm = TRUE))})
          overdisp_median=apply(zinb_fit[,,2],1,function(x){return(median(x,na.rm = TRUE))})
          dropout_max=apply(zinb_fit[,,3],1,function(x){return(max(x,na.rm = TRUE))})
          dropout_median=apply(zinb_fit[,,3],1,function(x){return(median(x,na.rm = TRUE))})
          
          log_deseq2_pval=tryCatch({-log10(deseq2_pval+min(deseq2_pval[deseq2_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_MAST_pval=tryCatch({-log10(MAST_pval+min(MAST_pval[MAST_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_jsd_empirical_pval=tryCatch({-log10(jsd_empirical_pval+min(jsd_empirical_pval[jsd_empirical_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_klmean_empirical_pval=tryCatch({-log10(klmean_empirical_pval+min(klmean_empirical_pval[klmean_empirical_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_jsd_zinb_pval=tryCatch({-log10(jsd_zinb_pval+min(jsd_zinb_pval[jsd_zinb_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_klmean_zinb_pval=tryCatch({-log10(klmean_zinb_pval+min(klmean_zinb_pval[klmean_zinb_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_jsd_direct_pval=tryCatch({-log10(jsd_direct_pval+min(jsd_direct_pval[jsd_direct_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_klmean_direct_pval=tryCatch({-log10(klmean_direct_pval+min(klmean_direct_pval[klmean_direct_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(zerorate_ind_mean,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(zerorate_ind_mean,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(zerorate_ind_mean,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(zerorate_ind_mean,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(zerorate_ind_mean,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(zerorate_ind_mean,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(zerorate_ind_mean,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(zerorate_ind_mean,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_zerorate_ind_mean_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(nonexpres_ind,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(nonexpres_ind,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(nonexpres_ind,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(nonexpres_ind,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(nonexpres_ind,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(nonexpres_ind,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(nonexpres_ind,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(nonexpres_ind,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_nonexpres_ind_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(zero_rate,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(zero_rate,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(zero_rate,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(zero_rate,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(zero_rate,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(zero_rate,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(zero_rate,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(zero_rate,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_zero_rate_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(expression_level,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(expression_level,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(expression_level,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(expression_level,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(expression_level,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(expression_level,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(expression_level,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(expression_level,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_expression_level_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(overdisp_max,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(overdisp_max,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(overdisp_max,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(overdisp_max,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(overdisp_max,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(overdisp_max,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(overdisp_max,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(overdisp_max,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_overdisp_max_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(overdisp_median,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(overdisp_median,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(overdisp_median,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(overdisp_median,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(overdisp_median,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(overdisp_median,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(overdisp_median,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(overdisp_median,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_overdisp_median_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(dropout_max,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(dropout_max,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(dropout_max,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(dropout_max,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(dropout_max,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(dropout_max,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(dropout_max,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(dropout_max,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_dropout_max_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=8,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(dropout_median,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(dropout_median,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(dropout_median,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(dropout_median,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(dropout_median,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(dropout_median,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[7]=tryCatch({cor(dropout_median,log_jsd_direct_pval)}, error = function(e) {NA} )
          cor_matrix[8]=tryCatch({cor(dropout_median,log_klmean_direct_pval)}, error = function(e) {NA} )
          cor_dropout_median_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
        }
        
        #power scatter
        png(paste0("../Data_PRJNA434002/8.Result/fig_power_point/power_point_",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 600,width = 600)
        plot(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,1],power_array[i_file,i_F,i_pre,i_cluster,1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
        abline(v = 0.05, col = "red") 
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,i_perm_label],power_array[i_file,i_F,i_pre,i_cluster,1,i_perm_label],col="blue",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,3],power_array[i_file,i_F,i_pre,i_cluster,1,3],col="pink",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,4],power_array[i_file,i_F,i_pre,i_cluster,1,4],col="brown",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,5],power_array[i_file,i_F,i_pre,i_cluster,1,5],col="orange",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,6],power_array[i_file,i_F,i_pre,i_cluster,1,6],col="green",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,7],power_array[i_file,i_F,i_pre,i_cluster,1,7],col="navy",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,8],power_array[i_file,i_F,i_pre,i_cluster,1,8],col="purple",pch=3,cex=3)
        
        legend("topright",c(names(power_matrix)),pch=rep(3,8),cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
          
        dev.off()
        
        print(paste0("p",perm_label,perm_method,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag))
      }
    }
  }
}

saveRDS(pval_list,paste0("../Data_PRJNA434002/8.Result/final_pval_list.rds"))
saveRDS(power_array,paste0("../Data_PRJNA434002/8.Result/final_power_array.rds"))
saveRDS(ks_array,paste0("../Data_PRJNA434002/8.Result/final_ks_array.rds"))
saveRDS(cor_zerorate_ind_mean_array,paste0("../Data_PRJNA434002/8.Result/final_cor_zerorate_ind_mean_array.rds"))
saveRDS(cor_nonexpres_ind_array,paste0("../Data_PRJNA434002/8.Result/final_cor_nonexpres_ind_array.rds"))
saveRDS(cor_zero_rate_array,paste0("../Data_PRJNA434002/8.Result/final_cor_zero_rate_array.rds"))
saveRDS(cor_expression_level_array,paste0("../Data_PRJNA434002/8.Result/final_cor_expression_level_array.rds"))
saveRDS(cor_overdisp_max_array,paste0("../Data_PRJNA434002/8.Result/final_cor_overdisp_max_array.rds"))
saveRDS(cor_overdisp_median_array,paste0("../Data_PRJNA434002/8.Result/final_cor_overdisp_median_array.rds"))
saveRDS(cor_dropout_max_array,paste0("../Data_PRJNA434002/8.Result/final_cor_dropout_max_array.rds"))
saveRDS(cor_dropout_median_array,paste0("../Data_PRJNA434002/8.Result/final_cor_dropout_median_array.rds"))
rownames(zeros)=rownames_zeros
#View(zeros)
saveRDS(zeros,paste0("../Data_PRJNA434002/10.Result/power_array_NAs_p",perm_label,perm_method,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".rds"))


###############Power array analysis###################
#######one-time plotting functions############
scatter_cell_num_1v8=function(a,...){
  plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,...)
  abline(h = -log10(0.05), col = "red") 
  points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
  points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
  points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
  points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
  points(cell_num,a[,6],pch=8,cex=1.1, col="green")
  points(cell_num,a[,7],pch=9,cex=1.1, col="navy")
  points(cell_num,a[,8],pch=10,cex=1.1, col="purple")
  legend("topright",c(colnames(a)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
}
#usage
#a=-log10(ks_array[i_file,i_F,i_pre,,2,]+min(ks_array[ks_array>0],na.rm = TRUE))
#scatter_cell_num_1v6(a,main="-log10 pvalues from KS test, observed data",ylab="-log10 ks pval",ylim=c(0,max_ab))


scatter_8v8=function(b,a,...){
  plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,...)
  abline(v = 0.05, col = "red") 
  points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
  points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
  points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
  points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
  points(b[,6],a[,6],pch=8,cex=1.1, col="green")
  points(b[,7],a[,7],pch=9,cex=1.1, col="navy")
  points(b[,8],a[,8],pch=10,cex=1.1, col="purple")
  legend("topright",c(colnames(b)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
}
#usage
#a=power_array[i_file,i_F,i_pre,,1,]
#b=power_array[i_file,i_F,i_pre,,2,]
#scatter_cell_num_1v6(b,a,ylab="observed (Power)",xlab="permutated (Type I error)",main="proportion of pval<0.05, observed vs permutated",ylim=c(0,1),xlim=c(0,1))

####################################
power_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_power_array.rds"))
#do boxplot##############
for(i_file in 1:length(file_tag_seq)){
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      file_tag=file_tag_seq[i_file]
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]
      png(paste0("../Data_PRJNA434002/8.Result/fig_boxplot_power/boxplot_power_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 400*(length(perm_label_seq)+1),width = 500)
      op=par(mfrow=c((length(perm_label_seq)+1),1))
      a=power_array[i_file,i_F,i_pre,,1,]
      boxplot(a,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of p-values less than 0.05 of all clusters, observed data",ylab="power")
      abline(h = 0.05, col = "red") 
      
      for(perm_label in perm_label_seq){
        b=power_array[i_file,i_F,i_pre,,perm_label,]
        boxplot(b,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of p-values less than 0.05 of all clusters, permutated data",ylab="type I error")
        abline(h = 0.05, col = "red") 
      }
      
      par(op)
      dev.off()
    }
  }
}


#do scatter plot about cell number vs power##############
for(i_file in 1:length(file_tag_seq)){
  file_tag=file_tag_seq[i_file]
  
  ###calculate cell num
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
  }
  cur_cluster=as.character(unique(tmeta$cluster))
  cell_num=table(tmeta$cluster)
  cell_num=as.numeric(cell_num[match(cur_cluster,names(cell_num))])
  names(cell_num)=cur_cluster
  
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]
      
      png(paste0("../Data_PRJNA434002/8.Result/fig_scatter_power/scatter_power_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500*length(perm_label_seq))
      op=par(mfrow=c(3,length(perm_label_seq)))

      a=power_array[i_file,i_F,i_pre,,1,]
      
      for(perm_label in perm_label_seq){
        b=power_array[i_file,i_F,i_pre,,perm_label,]
        scatter_cell_num_1v8(a,main="proportion of pval<0.05,observed data",ylab="Power",ylim=c(0,1))
        scatter_cell_num_1v8(b,main="proportion of pval<0.05,permutated data",ylab="Power",ylim=c(0,1))
        scatter_8v8(b,a,ylab="observed (Power)",xlab="permutated (Type I error)",main="proportion of pval<0.05, observed vs permutated",ylim=c(0,1),xlim=c(0,1))
      }
      
      par(op)
      dev.off()
    }
  }
}

###############KS test array analysis###################
#######one-time plotting functions############
scatter_cell_num_1v8=function(a,b=cell_num,...){
  plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,...)
  abline(h = -log10(0.05), col = "red") 
  points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
  points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
  points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
  points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
  points(cell_num,a[,6],pch=8,cex=1.1, col="green")
  points(cell_num,a[,7],pch=9,cex=1.1, col="navy")
  points(cell_num,a[,8],pch=10,cex=1.1, col="purple")
  legend("topright",c(colnames(a)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
}

#usage
#a=-log10(ks_array[i_file,i_F,i_pre,,2,]+min(ks_array[ks_array>0],na.rm = TRUE))
#scatter_cell_num_1v6(a,main="-log10 pvalues from KS test, observed data",ylab="-log10 ks pval",ylim=c(0,max_ab))

scatter_8v8=function(b,a,...){
  plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xlab="permutated",main="-log10 pvalues from KS test, observed vs permutated",ylim=c(0,max_ab),xlim=c(0,max_ab))
  abline(v = 0.05, col = "red") 
  points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
  points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
  points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
  points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
  points(b[,6],a[,6],pch=8,cex=1.1, col="green")
  points(b[,7],a[,7],pch=9,cex=1.1, col="navy")
  points(b[,8],a[,8],pch=10,cex=1.1, col="purple")
  legend("topright",c(colnames(b)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
}
#usage
#a=-log10(ks_array[i_file,i_F,i_pre,,1,]+min(ks_array[ks_array>0],na.rm = TRUE))
#b=-log10(ks_array[i_file,i_F,i_pre,,2,]+min(ks_array[ks_array>0],na.rm = TRUE))
#max_ab=max(a,b,na.rm = TRUE)
#scatter_6v6(b,a,ylab="observed",xlab="permutated",main="-log10 pvalues from KS test, observed vs permutated",ylim=c(0,max_ab),xlim=c(0,max_ab))

####################################
ks_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_ks_array.rds"))
#do boxplot##############
for(i_file in 1:length(file_tag_seq)){
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      file_tag=file_tag_seq[i_file]
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]
      png(paste0("../Data_PRJNA434002/8.Result/fig_boxplot_ks/boxplot_ks_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 400*(1+length(perm_label_seq)),width = 500)
      op=par(mfrow=c((length(perm_label_seq)+1),1))
      a=ks_array[i_file,i_F,i_pre,,1,]
      boxplot(a,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="ks test for distribution of pvalues, observed data",ylab="ks pval")
      abline(h = 0.05, col = "red") 
      for(perm_label in 1:length(perm_label_seq)){
        b=ks_array[i_file,i_F,i_pre,,perm_label,]
        boxplot(b,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="ks test for distribution of pvalues, permutated data",ylab="ks pval")
        abline(h = 0.05, col = "red") 
      }
      
      par(op)
      dev.off()
    }
  }
}


#do scatter plot about cell number vs ks##############


for(i_file in 1:length(file_tag_seq)){
  file_tag=file_tag_seq[i_file]
  
  ###calculate cell num
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
  }
  cur_cluster=as.character(unique(tmeta$cluster))
  cell_num=table(tmeta$cluster)
  cell_num=as.numeric(cell_num[match(cur_cluster,names(cell_num))])
  names(cell_num)=cur_cluster
  
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]
      
      png(paste0("../Data_PRJNA434002/8.Result/fig_scatter_ks/scatter_ks_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 500,width = 500*length(perm_label_seq))
      op=par(mfrow=c(3,length(perm_label_seq)))

      a=-log10(ks_array[i_file,i_F,i_pre,,1,]+min(ks_array[ks_array>0],na.rm = TRUE))
      
      for(perm_label in perm_label_seq){
        b=-log10(ks_array[i_file,i_F,i_pre,,perm_label,]+min(ks_array[ks_array>0],na.rm = TRUE))
        max_ab=max(a,b,na.rm = TRUE)
        
        scatter_cell_num_1v8(a,main="-log10 pvalues from KS test, observed data",ylab="-log10 ks pval",ylim=c(0,max_ab))
        scatter_cell_num_1v8(b,main="-log10 pvalues from KS test, permutated data",ylab="-log10 ks pval",ylim=c(0,max_ab))
        scatter_8v8(b,a,ylab="observed",xlab="permutated",main="-log10 pvalues from KS test, observed vs permutated",ylim=c(0,max_ab),xlim=c(0,max_ab))
      }
      par(op)
      dev.off()
    }
  }
}


#do scatter plot about cell number vs cor ##############


#######one-time plotting functions############
cor_scatter_plot=function(cur_array,cur_label){
  for(i_file in 1:length(file_tag_seq)){
    file_tag=file_tag_seq[i_file]
    
    ###calculate cell num
    if(is.na(unlist(strsplit(file_tag,"k"))[2])){
      tmeta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
    }
    if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
      tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
    }
    cur_cluster=as.character(unique(tmeta$cluster))
    cell_num=table(tmeta$cluster)
    cell_num=as.numeric(cell_num[match(cur_cluster,names(cell_num))])
    names(cell_num)=cur_cluster
    
    for(i_F in 1:length(F_method_seq)){
      for(i_pre in 1:length(pre_tag_seq)){
        
        F_method=F_method_seq[i_F]
        pre_tag=pre_tag_seq[i_pre]
        
        png(paste0("../Data_PRJNA434002/8.Result/fig_scatter_",cur_label,"/scatter_",cur_label,"_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500*dim(cur_array)[4])
        op=par(mfrow=c(3,dim(cur_array)[4]))
        
        a=cur_array[i_file,i_F,i_pre,,1,]
        
        for(ip in 2:dim(cur_array)[4]){
          b=cur_array[i_file,i_F,i_pre,,ip,]
          
          plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main=paste0("correlation between pval and ",cur_label,", observed data"),ylim=c(-1,1))
          abline(h = 0.0, col = "red") 
          points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
          points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
          points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
          points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
          points(cell_num,a[,6],pch=8,cex=1.1, col="green")
          points(cell_num,a[,7],pch=9,cex=1.1, col="navy")
          points(cell_num,a[,8],pch=10,cex=1.1, col="purple")
          legend("topright",c(colnames(a)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
          
          plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main=paste0("correlation between pval and ",cur_label,",permutated data"),ylim=c(-1,1))
          abline(h = 0.0, col = "red") 
          points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
          points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
          points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
          points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
          points(cell_num,b[,6],pch=8,cex=1.1, col="green")
          points(cell_num,a[,7],pch=9,cex=1.1, col="navy")
          points(cell_num,a[,8],pch=10,cex=1.1, col="purple")
          legend("topright",c(colnames(a)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
          
          plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main=paste0("correlation between pval and ",cur_label,", observed vs permutated"),ylim=c(-1,1),xlim=c(-1,1))
          abline(v = 0, col = "red") 
          abline(h = 0, col = "red") 
          points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
          points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
          points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
          points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
          points(b[,6],a[,6],pch=8,cex=1.1, col="green")
          points(cell_num,a[,7],pch=9,cex=1.1, col="navy")
          points(cell_num,a[,8],pch=10,cex=1.1, col="purple")
          legend("topright",c(colnames(a)),pch=3:10,cex=1.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
          
          par(op)
          dev.off()
        }
      }
    }
  }
  
}

####################################

cor_zerorate_ind_mean_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_zerorate_ind_mean_array.rds"))
cor_nonexpres_ind_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_nonexpres_ind_array.rds"))
cor_zerorate_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_zerorate_array.rds"))
cor_expression_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_expression_array.rds"))
cor_overdisp_max_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_overdisp_max_array.rds"))
cor_overdisp_median_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_overdisp_median_array.rds"))
cor_dropout_max_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_dropout_max_array.rds"))
cor_dropout_median_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_cor_dropout_median_array.rds"))

cor_scatter_plot(cur_array=cor_zerorate_ind_mean_array,cur_label="zerorate_ind_mean")
cor_scatter_plot(cur_array=cor_nonexpres_ind_array,cur_label="nonexpres_ind")
cor_scatter_plot(cur_array=cor_zerorate_array,cur_label="zerorate")
cor_scatter_plot(cur_array=cor_expression_array,cur_label="expression")
cor_scatter_plot(cur_array=cor_overdisp_max_array,cur_label="overdisp_max")
cor_scatter_plot(cur_array=cor_overdisp_median_array,cur_label="overdisp_median")
cor_scatter_plot(cur_array=cor_dropout_max_array,cur_label="dropout_max")
cor_scatter_plot(cur_array=cor_dropout_median_array,cur_label="dropout_median")





#pval related analysis ##############
pval_list=readRDS(paste0("../Data_PRJNA434002/8.Result/final_pval_list.rds"))

#############one-time functions##############
cur_hist_plot_param = function(xx,i1,i2,i3,...){
  xrange=range(min(xx[i1,,i3],na.rm = TRUE),max(xx[i1,,i3],na.rm = TRUE))
  hist(xx[i1,!i2,i3],xlim=xrange,col=rgb(1,0,0,0.1),...)
  hist(xx[i1,i2,i3],add=T,col=rgb(0,0,1,0.1))
  legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
}
#usage
#cur_hist_plot_param(xx=cur_fit,i1=iplot,i2=cur_phenotype_ind,i3=1,xlab="logmean",main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))


cur_boxplot_param = function(xx,ind_info,i1,i2,i3,...){
  xx0=as.matrix(xx[i1,!i2,i3])
  xx1=as.matrix(xx[i1,i2,i3])
  ind0=ind_info[!i2]
  ind1=ind_info[i2]
  xxc=rbind(xx0,xx1)
  indc=c(ind0,ind1)
  boxplot(xxc~indc,col=c(rep(rgb(1,0,0,0.8),length(unique(ind0))),rep(rgb(0,0,1,0.8),length(unique(ind1)))),...)
  legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.8),rgb(0,0,1,0.8)))
}

library(vioplot)
cur_vioplot_param = function(xx,ind_info,i1,i2,i3,...){
  xx0=as.matrix(xx[i1,!i2,i3])
  xx1=as.matrix(xx[i1,i2,i3])
  ind0=ind_info[!i2]
  ind1=ind_info[i2]
  xxc=rbind(xx0,xx1)
  indc=c(ind0,ind1)
  vioplot(xxc~indc,col=c(rep(rgb(1,0,0,0.8),length(unique(ind0))),rep(rgb(0,0,1,0.8),length(unique(ind1)))),...)
  legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.8),rgb(0,0,1,0.8)))
}

cur_scatter_param = function(xx,ind_info,i1,i2,i3,...){
  xx0=sort(as.numeric(xx[i1,!i2,i3]),decreasing = FALSE)
  xx1=sort(as.numeric(xx[i1,i2,i3]),decreasing = TRUE)
  xxc=c(xx0,xx1)

  plot(1:length(xxc),as.numeric(xxc),pch=16,col=c(rep(rgb(1,0,0,0.8),length(xx0)),rep(rgb(0,0,1,0.8),length(xx1))),...)
  legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.8),rgb(0,0,1,0.8)))
}
#usage

#cur_scatter_param(xx=cur_fit,ind_info=cur_individual_ind,i1=iplot,i2=cur_phenotype_ind,i3=1,xlab="individual",ylab="logmean",main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))

#cur_boxplot_param(xx=sim_data,ind_info=cur_individual,i1=iplot,i2=cur_phenotype,i3=1:dim(sim_data)[[3]],xlab="individual",ylab="log10 read count",main=paste0("sim gene expres,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))



#for plotting distances within/outside groups 
dist_matrix_density_plot = function(mx,i1,i2,...){
  mx=log(mx[i1,,]+1)
  diag(mx)=NA
  #within case
  d1=as.numeric(mx[i2,i2])
  d1=d1[!is.na(d1)]
  #within control
  d2=as.numeric(mx[!i2,!i2])
  d2=d2[!is.na(d2)]
  #outside
  d3=c(as.numeric(mx[i2,!i2]),as.numeric(mx[!i2,i2]))
  d3=d3[!is.na(d3)]
  xrange=range(min(mx,na.rm = TRUE)*0.9,max(mx,na.rm = TRUE)*1.1)
  q=quantile(xrange)
  plot(density(d1),lwd=2,xlim=xrange,ylim=c(0,max(density(d3)$y,density(d2)$y,density(d1)$y)),xaxt = 'n',col=rgb(1,0,0,0.7),...)
  lines(density(d2),lwd=2,col=rgb(0,0,1,0.7))
  lines(density(d3),lwd=2,col=rgb(0,1,0,0.7))
  axis(1, at = q, labels = round(10^q-1,3))
  
  legend("topright",c("within Case","within Control","between Groups"),bty = "n",pch=15,cex=1.1,col=c(rgb(1,0,0,0.7),rgb(0,0,1,0.7),rgb(0,1,0,0.7)))
}
#usage
#dist_matrix_density_plot(mx=klmean_empirical_array,i1=1,i2=cur_phenotype_ind,xlab="Distance",ylab="Density",
#                         main=paste0("Density of distance, klmean_empirical, sig pval of Perm ",iplot))

#############

for(i_file in 1:length(file_tag_seq)){
  file_tag=file_tag_seq[i_file]
  pval_array=pval_list[[file_tag]]
  ###calculate cell num
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
  }
  cur_cluster=as.character(unique(tmeta$cluster))
  cell_num=table(tmeta$cluster)
  cell_num=as.numeric(cell_num[match(cur_cluster,names(cell_num))])
  names(cell_num)=cur_cluster

  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]

      log10_ob_pval=-log10(pval_array[i_F,i_pre,,1,,]+min(pval_array[pval_array>0],na.rm = TRUE))
      
      for(perm_label in 1:length(perm_label_seq)){
        log10_perm_pval=-log10(pval_array[i_F,i_pre,,perm_label,,]+min(pval_array[pval_array>0],na.rm = TRUE))
        
        ##do scatter plot about ob pval vs perm pval
        cor_obperm=array(dim=dim(log10_ob_pval)[1:2],dimnames=list(dimnames(log10_ob_pval)[[1]],dimnames(log10_ob_pval)[[2]]))
        cor_obob=array(dim=dim(log10_ob_pval)[c(1,2,2)],dimnames=list(dimnames(log10_ob_pval)[[1]],dimnames(log10_ob_pval)[[2]],dimnames(log10_ob_pval)[[2]]))
        cor_permperm=array(dim=dim(log10_ob_pval)[c(1,2,2)],dimnames=list(dimnames(log10_ob_pval)[[1]],dimnames(log10_ob_pval)[[2]],dimnames(log10_ob_pval)[[2]]))
        for(ia in 1:dim(log10_ob_pval)[1]){
          for(ib in 1:dim(log10_ob_pval)[2]){
            
            #tryCatch({ }, error = function(e) {NA} )
            cor_obperm[ia,ib]=tryCatch({cor(log10_ob_pval[ia,ib,],log10_perm_pval[ia,ib,],use = "complete.obs")}, error = function(e) {NA} )
            cor_obob[ia,ib,]=tryCatch({apply(log10_ob_pval[ia,,],1, function(x){return(cor(x,log10_ob_pval[ia,ib,],use = "complete.obs"))})}, error = function(e) {NA} )
            cor_permperm[ia,ib,]=tryCatch({apply(log10_perm_pval[ia,,],1, function(x){return(cor(x,log10_perm_pval[ia,ib,],use = "complete.obs"))})}, error = function(e) {NA} )
          }
        }
        
        saveRDS(cor_obperm,paste0("../Data_PRJNA434002/8.Result/cor_pval_obperm_",F_method,"_",pre_tag,"_",file_tag,".rds"))
        saveRDS(cor_obob,paste0("../Data_PRJNA434002/8.Result/cor_pval_obob_",F_method,"_",pre_tag,"_",file_tag,".rds"))
        saveRDS(cor_permperm,paste0("../Data_PRJNA434002/8.Result/cor_pval_permperm_",F_method,"_",pre_tag,"_",file_tag,".rds"))
        #plot cor pval ob vs perm
        png(paste0("../Data_PRJNA434002/8.Result/fig_scatter_cor_pval/scatter_cor_pval_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 400,width = 500)
        plot(cell_num,cor_obperm[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="-log10 pval_cor",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between observed and permutation -log10 pvalues",,ylim=c(-1,1))
        abline(h = -log10(0.05), col = "red")
        points(cell_num,cor_obperm[,2],pch=4,cex=1.1, col="blue")
        points(cell_num,cor_obperm[,3],pch=5,cex=1.1, col="pink")
        points(cell_num,cor_obperm[,4],pch=6,cex=1.1, col="brown")
        points(cell_num,cor_obperm[,5],pch=7,cex=1.1, col="orange")
        points(cell_num,cor_obperm[,6],pch=8,cex=1.1, col="green")
        points(cell_num,cor_obperm[,7],pch=9,cex=1.1, col="navy")
        points(cell_num,cor_obperm[,8],pch=10,cex=1.1, col="purple")
        legend("topright",dimnames(log10_ob_pval)[[2]],pch=3:10,cex=.5,col=c("red","blue","pink","brown","orange","green","navy","purple"))
        
        dev.off()
        
        
        for(cluster_tag in 1:dim(log10_perm_pval)[1]){
          
          fit_data=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          sim_data=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          
          #permutation phenotype
          cur_phenotype=tmeta$diagnosis[tmeta$cluster==cur_cluster[cluster_tag]]!="Control"
          cur_phenotype_ind=tmeta$diagnosis[match(unlist(dimnames(fit_data)[2]),tmeta$individual)]!="Control"
          
          cur_individual=tmeta$individual[tmeta$cluster==cur_cluster[cluster_tag]]
          cur_individual_ind=tmeta$individual[match(unlist(dimnames(fit_data)[2]),tmeta$individual)]
          #count cases and controls
          diag_info=paste0(tmeta$ind,":",tmeta$diagnosis)
          diag_kind=unique(diag_info)
          diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
          if(perm_label>0){
            #permute
            perm_order=readRDS(paste0("../Data_PRJNA434002/7.Result/ind_perm_order.rds"))
            perm_order=as.numeric(perm_order[,perm_label])
            diag_kind[,2]=diag_kind[perm_order,2]
          }
          #match back to each individuals
          diag_kind=diag_kind[match(unlist(dimnames(fit_data)[[2]]),diag_kind[,1]),]
          
          ind_index=match(tmeta$ind,diag_kind[,1])
          ind_index=ind_index[!is.na(ind_index)]
          cur_phenotype=as.factor(diag_kind[ind_index,2])
          cur_phenotype_ind=diag_kind[,2]
          cur_phenotype=cur_phenotype[tmeta$cluster==cur_cluster[cluster_tag]]!="Control"
          cur_phenotype_ind=cur_phenotype_ind!="Control"
          
          
          
          
          #plot first 4 smallest permutated data's pval's gene's re-constructed expression distribution vs the median pvals
          print("sig_gene_count")
          png(paste0("../Data_PRJNA434002/8.Result/fig_sig_gene_count/sig_gene_count_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 2400,width =3200)
          op=par(mfrow = c(6, 8))
          for(i_pval in 1:dim(log10_perm_pval)[2]){
            #sig pval
            cur_index=order(log10_perm_pval[cluster_tag,i_pval,],decreasing = TRUE)[1:4]
            for(iplot in cur_index){
              #cur_hist_plot_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=1,xlab="logmean",
              #                    main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
              #cur_hist_plot_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=2,xlab="dispersion",
              #                    main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
              #cur_hist_plot_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=3,xlab="dropout",
              #                    main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
              
              cur_hist_plot_param(xx=sim_data,i1=iplot,i2=cur_phenotype,i3=1:dim(sim_data)[[3]],xlab="read count",
                                  main=paste0("sim gene expres,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
            }
            
            #median pval
            cur_index=order(log10_perm_pval[cluster_tag,i_pval,],decreasing = TRUE)[1501:1504]
            for(iplot in cur_index){
              #cur_hist_plot_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=1,xlab="logmean",
              #                    main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
              #cur_hist_plot_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=2,xlab="dispersion",
              #                    main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
              #cur_hist_plot_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=3,xlab="dropout",
              #                    main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
              
              cur_hist_plot_param(xx=sim_data,i1=iplot,i2=cur_phenotype,i3=1:dim(sim_data)[[3]],xlab="read count",
                                  main=paste0("sim gene expres,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
            }
          }
          par(op)
          dev.off()
          
          #plot first 4 smallest permutated data's pval's gene's re-constructed expression distribution vs the median pvals with individual info
          
          print("sig_gene_count_ind")
          png(paste0("../Data_PRJNA434002/8.Result/fig_sig_gene_count_ind/sig_gene_count_ind_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 4500,width =6000)
          op=par(mfrow = c(12, 16))
          for(i_pval in 1:dim(log10_perm_pval)[2]){
            #sig pval
            cur_index=order(log10_perm_pval[cluster_tag,i_pval,],decreasing = TRUE)[1:4]
            for(iplot in cur_index){
              cur_scatter_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=1,
                                ind_info=cur_individual_ind,xlab="individual",ylab="logmean",
                                main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
              cur_scatter_param(xx=log10(fit_data),i1=iplot,i2=cur_phenotype_ind,i3=2,
                                ind_info=cur_individual_ind,xlab="individual",ylab="log10 dispersion",
                                main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
              cur_scatter_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=3,
                                ind_info=cur_individual_ind,xlab="individual",ylab="dropout",
                                main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
              cur_boxplot_param(xx=sim_data,i1=iplot,i2=cur_phenotype,i3=1:dim(sim_data)[[3]],
                                ind_info=cur_individual,xlab="individual",ylab="read count",
                                main=paste0("sim gene expres,",dimnames(log10_perm_pval)[[2]][i_pval]," sig pval of Perm ",iplot))
            }
            
            #median pval
            cur_index=order(log10_perm_pval[cluster_tag,i_pval,],decreasing = TRUE)[1501:1504]
            for(iplot in cur_index){
              cur_scatter_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=1,
                                ind_info=cur_individual_ind,xlab="individual",ylab="logmean",
                                main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
              cur_scatter_param(xx=log10(fit_data),i1=iplot,i2=cur_phenotype_ind,i3=2,
                                ind_info=cur_individual_ind,xlab="individual",ylab="log10 dispersion",
                                main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
              cur_scatter_param(xx=fit_data,i1=iplot,i2=cur_phenotype_ind,i3=3,
                                ind_info=cur_individual_ind,xlab="individual",ylab="dropout",
                                main=paste0("ind-zinb param,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
              cur_boxplot_param(xx=sim_data,i1=iplot,i2=cur_phenotype,i3=1:dim(sim_data)[[3]],
                                ind_info=cur_individual,xlab="individual",ylab="read count",
                                main=paste0("sim gene expres,",dimnames(log10_perm_pval)[[2]][i_pval]," median pval of Perm ",iplot))
            }
          }
          par(op)
          dev.off()
          
          
          
          ##plot first 4 smallest permutated data's pval's gene's re-constructed expression distribution vs the median pvals 
          #Draw three density plot of histograms of the distances within cases, within controls and between cases and controls. 
          
          print("ind_dist_density")
          dist_list=list()
          dist_list[["jsd_empirical"]]=readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_empirical_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          dist_list[["klmean_empirical"]]=readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_empirical_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          dist_list[["jsd_zinb"]]=readRDS(paste0("../Data_PRJNA434002/8.Result/jsd_nbzinb_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          dist_list[["klmean_zinb"]]=readRDS(paste0("../Data_PRJNA434002/8.Result/klmean_nbzinb_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          
          
          png(paste0("../Data_PRJNA434002/8.Result/fig_ind_dist_density/ind_dist_density_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 2400,width =1200)
          #png(paste0("../Data_PRJNA434002/8.Result/fig_ind_dist_density/ind_dist_density_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 2400,width =1200)
          op=par(mfrow = c(8, 4))
          for(i_pval in 3:dim(log10_perm_pval)[2]){
            cur_label=dimnames(log10_perm_pval)[[2]][i_pval]
            cur_array=dist_list[[cur_label]]
            #sig pval
            cur_index=order(log10_perm_pval[cluster_tag,i_pval,],decreasing = TRUE)[1:4]
            for(iplot in cur_index){ #plot density
              dist_matrix_density_plot(mx=cur_array,i1=iplot,i2=cur_phenotype_ind,xlab="Distance",ylab="Density",
                                       main=paste0("Density of distance, ",cur_label,", sig pval of Perm ",iplot))
            }
            #median pval
            cur_index=order(log10_perm_pval[cluster_tag,i_pval,],decreasing = TRUE)[1501:1504]
            for(iplot in cur_index){ #plot density
              dist_matrix_density_plot(mx=cur_array,i1=iplot,i2=cur_phenotype_ind,xlab="Distance",ylab="Density",
                                       main=paste0("Density of distance, ",cur_label,", median pval of Perm ",iplot))
            }
          }
          par(op)
          dev.off()
          
          #plot the relationship between log10_perm_pval and a lot of characters
          
          print("scatter_pval_cor_pheno")
          #statistical info
          zero_rate_ind=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          zerorate_ind_mean=apply(zero_rate_ind,1,mean)
          nonexpres_ind=apply(zero_rate_ind==10,1,sum)
          zero_rate=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          expression_level=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_gene_read_count_total_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          zinb_fit=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          logmean_max=apply(zinb_fit[,,1],1,function(x){return(max(x,na.rm = TRUE))})
          logmean_median=apply(zinb_fit[,,1],1,function(x){return(median(x,na.rm = TRUE))})
          overdisp_max=apply(zinb_fit[,,2],1,function(x){return(max(x,na.rm = TRUE))})
          overdisp_median=apply(zinb_fit[,,2],1,function(x){return(median(x,na.rm = TRUE))})
          dropout_max=apply(zinb_fit[,,3],1,function(x){return(max(x,na.rm = TRUE))})
          dropout_median=apply(zinb_fit[,,3],1,function(x){return(median(x,na.rm = TRUE))})
          
          log10_expression_level=log10(expression_level+1)
          log10_overdisp_max=log10(overdisp_max+1)
          log10_overdisp_median=log10(overdisp_median+1)
          #plot
          png(paste0("../Data_PRJNA434002/8.Result/fig_scatter_pval_cor_pheno/scatter_pval_cor_pheno_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 1200,width =2000)
          op=par(mfrow = c(8, 10))
          for(i_pval in 1:dim(log10_perm_pval)[2]){
            cur_plot_pval=NA
            cur_plot_pval=log10_perm_pval[cluster_tag,i_pval,]
            if(sum(!is.na(cur_plot_pval))>0){
              #zerorate_ind_mean
              plot(cur_plot_pval,zerorate_ind_mean,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="zerorate_ind_mean",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x zerorate_ind_mean"))
              abline(v = -log10(0.05), col = "red")
              
              #nonexpres_ind
              plot(cur_plot_pval,nonexpres_ind,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="nonexpres_ind",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x nonexpres_ind"))
              abline(v = -log10(0.05), col = "red")
              
              #zero_rate
              plot(cur_plot_pval,zero_rate,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="zero_rate",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x zero_rate"))
              abline(v = -log10(0.05), col = "red")
              
              #log10_expression_level
              plot(cur_plot_pval,log10_expression_level,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="log10_expression_level",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x log10_expression_level"))
              abline(v = -log10(0.05), col = "red")
              
              #logmean_max
              plot(cur_plot_pval,logmean_max,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="logmean_max",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x logmean_max"))
              abline(v = -log10(0.05), col = "red")
              
              #logmean_median
              plot(cur_plot_pval,logmean_median,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="logmean_median",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x logmean_median"))
              abline(v = -log10(0.05), col = "red")
              
              #log10_overdisp_max
              plot(cur_plot_pval,log10_overdisp_max,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="log10_overdisp_max",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x log10_overdisp_max"))
              abline(v = -log10(0.05), col = "red")
              
              #log10_overdisp_median
              plot(cur_plot_pval,log10_overdisp_median,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="log10_overdisp_median",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x log10_overdisp_median"))
              abline(v = -log10(0.05), col = "red")
              
              #dropout_max
              plot(cur_plot_pval,dropout_max,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="dropout_max",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x dropout_max"))
              abline(v = -log10(0.05), col = "red")
              
              #dropout_median
              plot(cur_plot_pval,dropout_median,type="p",pch=1,cex=0.1, xlab="-log10 pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,ylab="dropout_median",main=paste0(dimnames(log10_perm_pval)[[2]][i_pval]," x dropout_median"))
              abline(v = -log10(0.05), col = "red")
            }
          }
          par(op)
          dev.off()
          
          print(paste0(F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag))
          
        }
      }
    }
  }
}

