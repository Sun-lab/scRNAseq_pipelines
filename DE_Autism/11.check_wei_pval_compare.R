setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")


ref_pval=as.matrix(read.table("~/Desktop/github/ideas/Autism/R/res/step5a_pvals.tsv",header = TRUE,row.names = 1))

cluster_tag=4
file_tag="PFC5k" 
F_method_seq=c("ps","p")
fit_tag=""
resid_flag_seq=c("logresid","","adj")
covariate_flag_seq=c("", "readdepth")
pre_tag="dca" 

ind_covariate_flag="ind"

perm_method=""

method_seq=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")

cur_pval=readRDS("/Volumes/SpecialData/fh_data/Data_PRJNA434002/8.Result/jsd_direct_pval/p0_jsd_direct_p_pval_ind_dca_sim_1_5k.rds")

name_index=match(names(cur_pval),rownames(ref_pval))
#name_index[is.na(name_index)]==0
ref_pval=ref_pval[name_index,]
#rownames(ref_pval)=names(cur_pval)
ref_log_pval=-log10(ref_pval)
pval_array=array(dim=c(
  length(F_method_seq),
  length(resid_flag_seq),
  length(covariate_flag_seq),
  length(method_seq),length(cur_pval)),
  dimnames = list(
    F_method_seq,
    resid_flag_seq,
    covariate_flag_seq,
    method_seq,names(cur_pval)))




for(i_F in 1:length(F_method_seq)){
  for(i_resid in 1:length(resid_flag_seq)){
    for(i_cov in 1:length(covariate_flag_seq)){

      F_method=F_method_seq[i_F]
      resid_flag=resid_flag_seq[i_resid]
      covariate_flag=covariate_flag_seq[i_cov]
      
      pdf(paste0("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/step5a_pvals_compare/",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_",cluster_tag,"_",file_tag,".pdf"),height = 32,width = 16)
      op=par(mfrow=c(8,4))
      
      jsd_zinb_pval=NA
      jsd_empirical_pval=NA
      jsd_direct_pval=NA
      klmean_zinb_pval=NA
      klmean_empirical_pval=NA
      klmean_direct_pval=NA
      deseq2_pval=NA
      MAST_pval=NA
      
      jsd_zinb_pval=tryCatch({readRDS(paste0("./8.Result/jsd_nbzinb_pval/p10",perm_method,"_jsd_nbzinb_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      jsd_empirical_pval=tryCatch({readRDS(paste0("./8.Result/jsd_empirical_pval/p10",perm_method,"_jsd_empirical_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      jsd_direct_pval=tryCatch({readRDS(paste0("./8.Result/jsd_direct_pval/p10",perm_method,"_jsd_direct_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )  
      klmean_zinb_pval=tryCatch({readRDS(paste0("./8.Result/klmean_nbzinb_pval/p10",perm_method,"_klmean_nbzinb_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      klmean_empirical_pval=tryCatch({readRDS(paste0("./8.Result/klmean_empirical_pval/p10",perm_method,"_klmean_empirical_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      klmean_direct_pval=tryCatch({readRDS(paste0("./8.Result/klmean_direct_pval/p10",perm_method,"_klmean_direct_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      
      
      deseq2_pval=tryCatch({readRDS(paste0("./8.Result/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_ob_pval_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      
      MAST_pval=tryCatch({readRDS(paste0("./8.Result/MAST_pval/p",perm_label,perm_method,"_MAST_pval1_rawcount_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
      
      pval_array[i_F,i_resid,i_cov,1,]=deseq2_pval
      pval_array[i_F,i_resid,i_cov,2,]=MAST_pval
      pval_array[i_F,i_resid,i_cov,3,]=jsd_empirical_pval
      pval_array[i_F,i_resid,i_cov,4,]=klmean_empirical_pval
      pval_array[i_F,i_resid,i_cov,5,]=jsd_zinb_pval
      pval_array[i_F,i_resid,i_cov,6,]=klmean_zinb_pval
      pval_array[i_F,i_resid,i_cov,7,]=jsd_direct_pval
      pval_array[i_F,i_resid,i_cov,8,]=klmean_direct_pval
      
      cur_log_pval=-log(t(pval_array[i_F,i_resid,i_cov,,]))
      x_lim=max(ref_log_pval[is.finite(ref_log_pval)],na.rm = TRUE)
      y_lim=max(cur_log_pval[is.finite(ref_log_pval)],na.rm = TRUE)
      xy_lim=max(x_lim,y_lim)
      
      for(i_t in 1:8){
        for(i_ref in 1:ncol(ref_pval)){
          tryCatch({plot(ref_log_pval[,i_ref],cur_log_pval[,i_t],xlim=c(0,xy_lim),ylim=c(0,xy_lim),xlab=colnames(ref_log_pval)[i_ref],ylab=colnames(cur_log_pval)[i_t],cex=.1,sub=cor(ref_log_pval[,i_ref],cur_log_pval[,i_t], use="complete.obs"))}, error = function(e) {NA} )
          tryCatch({lines(c(0,xy_lim),c(0,xy_lim),col="red")}, error = function(e) {NA} )
        }
      }
      
      par(op)
      dev.off()
      print(paste0(F_method,"_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_",cluster_tag,"_",file_tag))
    }
  }
  
}
