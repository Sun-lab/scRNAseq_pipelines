#this code try to find if here are some genes who has higher expression with the smaller p-values.

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/10.Result/sim_v6/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#complete param
file_tag=1
perm_label=0
perm_method=""
param_tag=1 #c(1,2,3,4)
perm_label_seq=0
param_tag_seq=1:4
perm_method_seq="" #c("","b")
pre_tag_seq=c("", "dca")
fit_tag_seq=c("", "nb") 
covariate_flag_seq=c("","readdepth")
resid_flag_seq=c("","resid")
dist_method_seq=c("mean","JSD")
fit_method_seq=c("empirical","zinb","direct")

#shrink param

param_tag=1
pre_tag="dca"
fit_tag=""
covariate_flag="readdepth"
resid_flag="resid"
perm_method=""
perm_label=0
dist_method="JSD"
fit_method="empirical"

r_mean=1.2
r_var=1.2
r_dp=0.2
r_mult=0.6

n_ind=10
n_cell=400

mean_index=readRDS(paste0("./de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
var_index=readRDS(paste0("./de_label/sim_de.var_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
dp_index=readRDS(paste0("./de_label/sim_de.dp_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
mult_index=readRDS(paste0("./de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))


#readin bulk situation
sim_matrix_bulk=readRDS(paste0("./sim_data/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
bulk_mean=apply(sim_matrix_bulk,1,mean)
bulk_sd=apply(sim_matrix_bulk,1,sd)

#tryCatch({     }, error = function(e) {NA} )

jsd_zinb_pval=NA
jsd_empirical_pval=NA
jsd_direct_pval=NA
klmean_zinb_pval=NA
klmean_empirical_pval=NA
klmean_direct_pval=NA
deseq2_pval=NA
MAST_pval=NA

if(perm_label>0){
  tryCatch({jsd_zinb_pval=readRDS(paste0("./JSD_zinb_pval/p10",perm_method,"_JSD_zinb_perm_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
  tryCatch({jsd_empirical_pval=readRDS(paste0("./JSD_empirical_pval/p10",perm_method,"_JSD_empirical_perm_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
  tryCatch({jsd_direct_pval=readRDS(paste0("./JSD_direct_pval/p10",perm_method,"_JSD_direct_perm_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
  tryCatch({klmean_zinb_pval=readRDS(paste0("./mean_zinb_pval/p10",perm_method,"_mean_zinb_perm_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
  tryCatch({klmean_empirical_pval=readRDS(paste0("./mean_empirical_pval/p10",perm_method,"_mean_empirical_perm_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
  tryCatch({klmean_direct_pval=readRDS(paste0("./mean_direct_pval/p10",perm_method,"_mean_direct_perm_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
  tryCatch({deseq2_pval=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,".rds"))[,perm_label]}, error = function(e) {NA} )
  
  #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
  tryCatch({MAST_pval=readRDS(paste0("./MAST_pval/p",perm_label,perm_method,"_MAST_pval1_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
}
if(perm_label==0){
  tryCatch({jsd_zinb_pval=readRDS(paste0("./JSD_zinb_pval/p",perm_label,perm_method,"_JSD_zinb_raw_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({jsd_empirical_pval=readRDS(paste0("./JSD_empirical_pval/p",perm_label,perm_method,"_JSD_empirical_raw_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({jsd_direct_pval=readRDS(paste0("./JSD_direct_pval/p",perm_label,perm_method,"_JSD_direct_raw_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_zinb_pval=readRDS(paste0("./mean_zinb_pval/p",perm_label,perm_method,"_mean_zinb_raw_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_empirical_pval=readRDS(paste0("./mean_empirical_pval/p",perm_label,perm_method,"_mean_empirical_raw_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_direct_pval=readRDS(paste0("./mean_direct_pval/p",perm_label,perm_method,"_mean_direct_raw_pval_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({deseq2_pval=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,".rds"))}, error = function(e) {NA} )
  
  #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
  tryCatch({MAST_pval=readRDS(paste0("./MAST_pval/p0_MAST_pval1_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
}



#pval vs bulk cor
png(paste0("./fig_pval_bulk_cor/p",perm_label,perm_method,"_pval_bulkmean_cor_",param_tag,"_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 4800,width = 3000)
op=par(mfrow = c(8, 5),cex=2)

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(klmean_zinb_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,klmean_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(klmean_zinb_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,klmean_zinb")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(klmean_zinb_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,klmean_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(klmean_zinb_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,klmean_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,klmean_zinb")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(jsd_zinb_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,jsd_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(jsd_zinb_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,jsd_zinb")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(jsd_zinb_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,jsd_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(jsd_zinb_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,jsd_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,jsd_zinb")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(klmean_empirical_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,klmean_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(klmean_empirical_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,klmean_empirical")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(klmean_empirical_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,klmean_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(klmean_empirical_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,klmean_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,klmean_empirical")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(jsd_empirical_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,jsd_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(jsd_empirical_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,jsd_empirical")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(jsd_empirical_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,jsd_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(jsd_empirical_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,jsd_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,jsd_empirical")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(jsd_direct_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,jsd_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(jsd_direct_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,jsd_direct")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(jsd_direct_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,jsd_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(jsd_direct_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,jsd_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,jsd_direct")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(klmean_direct_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,klmean_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(klmean_direct_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,klmean_direct")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(klmean_direct_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,klmean_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(klmean_direct_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,klmean_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,klmean_direct")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(deseq2_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,deseq2")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(deseq2_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,deseq2")}, error = function(e) {NA} )  
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(deseq2_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,deseq2")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(deseq2_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,deseq2")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,deseq2")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_mean[mean_index==1]),-log10(MAST_pval[mean_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mean-DE genes,MAST")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[var_index==1]),-log10(MAST_pval[var_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk var-DE genes,MAST")}, error = function(e) {NA} )  
tryCatch({plot(log10(bulk_mean[dp_index==1]),-log10(MAST_pval[dp_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk dp-DE genes,MAST")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mult_index==1]),-log10(MAST_pval[mult_index==1]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk mult-DE genes,MAST")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_mean[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_mean",ylab="-log10pval",main="pval vs bulk non-DE genes,MAST")}, error = function(e) {NA} )

par(op)
dev.off()


png(paste0("./fig_pval_bulk_cor/p",perm_label,perm_method,"_pval_bulksd_cor_",param_tag,"_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 4800,width = 3000)
op=par(mfrow = c(8, 5),cex=2)

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(klmean_zinb_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,klmean_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(klmean_zinb_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,klmean_zinb")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(klmean_zinb_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,klmean_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(klmean_zinb_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,klmean_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,klmean_zinb")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(jsd_zinb_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,jsd_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(jsd_zinb_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,jsd_zinb")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(jsd_zinb_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,jsd_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(jsd_zinb_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,jsd_zinb")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,jsd_zinb")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(klmean_empirical_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,klmean_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(klmean_empirical_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,klmean_empirical")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(klmean_empirical_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,klmean_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(klmean_empirical_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,klmean_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,klmean_empirical")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(jsd_empirical_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,jsd_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(jsd_empirical_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,jsd_empirical")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(jsd_empirical_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,jsd_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(jsd_empirical_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,jsd_empirical")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,jsd_empirical")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(jsd_direct_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,jsd_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(jsd_direct_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,jsd_direct")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(jsd_direct_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,jsd_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(jsd_direct_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,jsd_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,jsd_direct")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(klmean_direct_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,klmean_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(klmean_direct_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,klmean_direct")}, error = function(e) {NA} )        
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(klmean_direct_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,klmean_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(klmean_direct_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,klmean_direct")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,klmean_direct")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(deseq2_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,deseq2")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(deseq2_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,deseq2")}, error = function(e) {NA} )  
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(deseq2_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,deseq2")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(deseq2_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,deseq2")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,deseq2")}, error = function(e) {NA} )

tryCatch({plot(log10(bulk_sd[mean_index==1]),-log10(MAST_pval[mean_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mean-DE genes,MAST")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[var_index==1]),-log10(MAST_pval[var_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk var-DE genes,MAST")}, error = function(e) {NA} )  
tryCatch({plot(log10(bulk_sd[dp_index==1]),-log10(MAST_pval[dp_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk dp-DE genes,MAST")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mult_index==1]),-log10(MAST_pval[mult_index==1]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk mult-DE genes,MAST")}, error = function(e) {NA} )
tryCatch({plot(log10(bulk_sd[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),-log10(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0]),cex=.1,xlab="log10_bulk_sd",ylab="-log10pval",main="pval vs bulk non-DE genes,MAST")}, error = function(e) {NA} )

par(op)
dev.off()

