#this code plots the comparison between observation and permutation.

setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")

pval_list=readRDS(paste0("./8.Result/final_pval_list.rds"))

#note the 6th dimension are represents the ob-perm situation, 
#they have 11 elements, the 1st is ob pval, the rest 10 are perm pvals.

cluster_tag_seq=1:17
file_tag="PFC3k"
dist_method_seq=c("klmean","jsd")
fit_method_seq=c("empirical","nbzinb")
F_method_seq=c("p","ps")
fit_tag_seq=c("","nb") #"","nb","zinb"
resid_flag_seq=c("", "logresid","adj")
covariate_flag_seq=c("", "readdepth")
method_seq=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")

pdf("~/Desktop/github/scRNAseq_pipelines/DE_Autism/11.check/11.check_ob_perm10_pval_distr_PFC3kL23.pdf",width=32,height=8)
for(i_F in 1:length(F_method_seq)){
  for(i_fit in 1:length(fit_tag_seq)){
    for(i_resid in 1:length(resid_flag_seq)){
      for(i_cov in 1:length(covariate_flag_seq)){
        for(i_cluster in 12){
          F_method=F_method_seq[i_F]
          fit_tag=fit_tag_seq[i_fit]
          resid_flag=resid_flag_seq[i_resid]
          covariate_flag=covariate_flag_seq[i_cov]
          cluster_tag=cluster_tag_seq[i_cluster]
          
          op=par(mfrow = c(2, 8))
          #hist_ob
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,1,],main="ob deseq2_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,2,],main="ob MAST_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,3,],main="ob jsd_empirical_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,4,],main="ob klmean_empirical_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,5,],main="ob jsd_zinb_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,6,],main="ob klmean_zinb_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,7,],main="ob jsd_direct_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,1,8,],main="ob klmean_direct_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )          

          #tryCatch({hist_perm
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,1,],main="perm deseq2_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,2,],main="perm MAST_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,3,],main="perm jsd_empirical_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,4,],main="perm klmean_empirical_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,5,],main="perm jsd_zinb_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,6,],main="perm klmean_zinb_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,7,],main="perm jsd_direct_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          tryCatch({hist(pval_list[[file_tag]][i_F,i_fit,i_cluster,i_resid,i_cov,2:11,8,],main="perm klmean_direct_pval", sub=paste0(F_method,"_",cluster_tag,"_",fit_tag,resid_flag,covariate_flag))}, error = function(e) {NA} )
          par(op)
          
        }
      }
    }
  }
}
dev.off()
