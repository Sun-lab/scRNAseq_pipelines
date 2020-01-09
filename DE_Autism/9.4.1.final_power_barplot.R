
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialPass/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

cluster_tag_seq=1:17
file_tag_seq=c("3k10","5k")
pre_tag_seq=c("dca","scvi")
dist_method_seq=c("klmean","jsd")
fit_method_seq=c("empirical","nbzinb")
F_method_seq=c("p","ps")


perm_label_seq=c(0,1)
ind_covariate_flag=NA


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
  6),
  dimnames = list(
    file_tag_seq,
    F_method_seq,
    pre_tag_seq,
    cluster_tag_seq,
    c("Power","Type I error"),
    c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")))

ks_array=power_array
cor_nonexpres_ind_array=power_array
cor_zerorate_ind_mean_array=power_array
cor_zerorate_array=power_array
cor_expression_array=power_array
cor_overdisp_median_array=power_array
cor_overdisp_max_array=power_array
cor_dropout_median_array=power_array
cor_dropout_max_array=power_array
count=1
zeros=matrix(ncol=6,nrow=length(file_tag_seq)*length(F_method_seq)*length(pre_tag_seq)*length(cluster_tag_seq)*length(perm_label_seq))
rownames_zeros=matrix(ncol=1,nrow=length(file_tag_seq)*length(F_method_seq)*length(pre_tag_seq)*length(cluster_tag_seq)*length(perm_label_seq))
colnames(zeros)=c("deseq2_pval","MAST_pval","jsd_empirical_pval","klmean_empirical_pval","jsd_zinb_pval","klmean_zinb_pval")

pval_list=list()

for(i_file in 1:length(file_tag_seq)){
  file_tag=file_tag_seq[i_file]
  pval_length=as.numeric(unlist(strsplit(file_tag,"k"))[1])*1000
  pval_list[[file_tag]]=array(dim=c(
    length(F_method_seq),
    length(pre_tag_seq),
    length(cluster_tag_seq),
    length(perm_label_seq),
    6,pval_length),
    dimnames = list(
      F_method_seq,
      pre_tag_seq,
      cluster_tag_seq,
      c("Power","Type I error"),
      c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb"),1:pval_length))
  
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
          klmean_zinb_pval=NA
          klmean_empirical_pval=NA
          deseq2_pval=NA
          MAST_pval=NA
          
          jsd_zinb_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_jsd_nbzinb_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          jsd_empirical_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_jsd_empirical_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          klmean_zinb_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_klmean_nbzinb_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          klmean_empirical_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_klmean_empirical_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          deseq2_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_DESeq2_ob_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )

          MAST_pval=tryCatch({readRDS(paste0("../Data_PRJNA434002/8.Result/p",perm_label,"_MAST_pval1_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))}, error = function(e) {NA} )
          
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,1,]=deseq2_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,2,]=MAST_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,3,]=jsd_empirical_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,4,]=klmean_empirical_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,5,]=jsd_zinb_pval
          pval_list[[file_tag]][i_F,i_pre,i_cluster,i_perm_label,6,]=klmean_zinb_pval
          
          zeros[count,]=c(sum(is.na(deseq2_pval)),sum(is.na(MAST_pval)), sum(is.na(jsd_empirical_pval)),sum(is.na(klmean_empirical_pval)),sum(is.na(jsd_zinb_pval)),sum(is.na(klmean_zinb_pval)))
          rownames_zeros[count]=paste0(perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag)
          count=count+1  
          
          #histogram
          png(paste0("../Data_PRJNA434002/8.Result/final_power_p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 1200,width = 800)
          op=par(mfrow = c(3, 2), pty = "s")
          tryCatch({hist(deseq2_pval,main="pval of deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(MAST_pval,main="pval of MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(jsd_empirical_pval,main="pval of jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(klmean_empirical_pval,main="pval of klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(jsd_zinb_pval,main="pval of jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          tryCatch({hist(klmean_zinb_pval,main="pval of klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
          
          par(op)
          
          dev.off()
          
          #record power
          power_matrix=matrix(nrow=6,ncol=1)
          names(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(power_matrix)="pval"
          power_matrix[1]=tryCatch({cal_power(deseq2_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[2]=tryCatch({cal_power(MAST_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[3]=tryCatch({cal_power(jsd_empirical_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[4]=tryCatch({cal_power(klmean_empirical_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[5]=tryCatch({cal_power(jsd_zinb_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_matrix[6]=tryCatch({cal_power(klmean_zinb_pval,threshold = 0.05)}, error = function(e) {NA} )
          power_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=power_matrix
          
          #record ks test result
          ks_matrix=matrix(nrow=6,ncol=1)
          names(ks_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(ks_matrix)="pval"
          ks_matrix[1]=tryCatch({cal_ks(deseq2_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[2]=tryCatch({cal_ks(MAST_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[3]=tryCatch({cal_ks(jsd_empirical_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[4]=tryCatch({cal_ks(klmean_empirical_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[5]=tryCatch({cal_ks(jsd_zinb_pval,method = "less")}, error = function(e) {NA} )
          ks_matrix[6]=tryCatch({cal_ks(klmean_zinb_pval,method = "less")}, error = function(e) {NA} )
          ks_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=ks_matrix
          
          
          #record gene-based cor test result: zero rate ind and expression

          zero_rate_ind=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          zerorate_ind_mean=apply(zero_rate_ind,1,mean)
          nonexpres_ind=apply(zero_rate_ind==10,1,sum)

          zero_rate=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_zero_rate_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          expression_level=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_gene_read_count_total_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

          zinb_fit=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          overdisp_max=apply(log(zinb_fit[,,2]),1,function(x){return(max(x,na.rm = TRUE))})
          overdisp_median=apply(log(zinb_fit[,,2]),1,function(x){return(max(x,na.rm = TRUE))})
          dropout_max=apply(log(zinb_fit[,,3]),1,function(x){return(max(x,na.rm = TRUE))})
          dropout_median=apply(log(zinb_fit[,,3]),1,function(x){return(max(x,na.rm = TRUE))})
          
          log_deseq2_pval=tryCatch({-log10(deseq2_pval+min(deseq2_pval[deseq2_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_MAST_pval=tryCatch({-log10(MAST_pval+min(MAST_pval[MAST_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_jsd_empirical_pval=tryCatch({-log10(jsd_empirical_pval+min(jsd_empirical_pval[jsd_empirical_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_klmean_empirical_pval=tryCatch({-log10(klmean_empirical_pval+min(klmean_empirical_pval[klmean_empirical_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_jsd_zinb_pval=tryCatch({-log10(jsd_zinb_pval+min(jsd_zinb_pval[jsd_zinb_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          log_klmean_zinb_pval=tryCatch({-log10(klmean_zinb_pval+min(klmean_zinb_pval[klmean_zinb_pval>0],na.rm = TRUE))}, error = function(e) {NA} )
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(zerorate_ind_mean,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(zerorate_ind_mean,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(zerorate_ind_mean,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(zerorate_ind_mean,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(zerorate_ind_mean,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(zerorate_ind_mean,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_zerorate_ind_mean_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(nonexpres_ind,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(nonexpres_ind,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(nonexpres_ind,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(nonexpres_ind,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(nonexpres_ind,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(nonexpres_ind,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_nonexpres_ind_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(zero_rate,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(zero_rate,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(zero_rate,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(zero_rate,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(zero_rate,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(zero_rate,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_zerorate_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(expression_level,log_deseq2_pval)}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(expression_level,log_MAST_pval)}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(expression_level,log_jsd_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(expression_level,log_klmean_empirical_pval)}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(expression_level,log_jsd_zinb_pval)}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(expression_level,log_klmean_zinb_pval)}, error = function(e) {NA} )
          cor_expression_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(overdisp_max,log_deseq2_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(overdisp_max,log_MAST_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(overdisp_max,log_jsd_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(overdisp_max,log_klmean_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(overdisp_max,log_jsd_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(overdisp_max,log_klmean_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_overdisp_max_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(overdisp_median,log_deseq2_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(overdisp_median,log_MAST_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(overdisp_median,log_jsd_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(overdisp_median,log_klmean_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(overdisp_median,log_jsd_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(overdisp_median,log_klmean_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_overdisp_median_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(dropout_max,log_deseq2_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(dropout_max,log_MAST_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(dropout_max,log_jsd_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(dropout_max,log_klmean_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(dropout_max,log_jsd_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(dropout_max,log_klmean_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_dropout_max_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
          
          cor_matrix=matrix(nrow=6,ncol=1)
          names(cor_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
          colnames(cor_matrix)="cor"
          cor_matrix[1]=tryCatch({cor(dropout_median,log_deseq2_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[2]=tryCatch({cor(dropout_median,log_MAST_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[3]=tryCatch({cor(dropout_median,log_jsd_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[4]=tryCatch({cor(dropout_median,log_klmean_empirical_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[5]=tryCatch({cor(dropout_median,log_jsd_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_matrix[6]=tryCatch({cor(dropout_median,log_klmean_zinb_pval,use="complete.obs")}, error = function(e) {NA} )
          cor_dropout_median_array[i_file,i_F,i_pre,i_cluster,i_perm_label,]=cor_matrix
        }
        #barplot
        png(paste0("../Data_PRJNA434002/8.Result/barplot_p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 1200,width = 800)
        op=par(mfrow = c(3, 2), pty = "s")
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,1],ylab="power",main=names(power_matrix)[1],ylim=c(0,1))
        abline(h = 0.05, col = "red") 
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,2],ylab="power",main=names(power_matrix)[2],ylim=c(0,1))
        abline(h = 0.05, col = "red") 
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,3],ylab="power",main=names(power_matrix)[3],ylim=c(0,1))
        abline(h = 0.05, col = "red") 
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,4],ylab="power",main=names(power_matrix)[4],ylim=c(0,1))
        abline(h = 0.05, col = "red") 
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,5],ylab="power",main=names(power_matrix)[5],ylim=c(0,1))
        abline(h = 0.05, col = "red") 
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,6],ylab="power",main=names(power_matrix)[6],ylim=c(0,1))
        abline(h = 0.05, col = "red") 
        par(op)
        dev.off()
        
        #power scatter
        png(paste0("../Data_PRJNA434002/8.Result/power_point_",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 600,width = 600)
        plot(power_array[i_file,i_F,i_pre,i_cluster,2,1],power_array[i_file,i_F,i_pre,i_cluster,1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
        abline(v = 0.05, col = "red") 
        points(power_array[i_file,i_F,i_pre,i_cluster,2,2],power_array[i_file,i_F,i_pre,i_cluster,1,2],col="blue",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,3],power_array[i_file,i_F,i_pre,i_cluster,1,3],col="pink",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,4],power_array[i_file,i_F,i_pre,i_cluster,1,4],col="brown",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,5],power_array[i_file,i_F,i_pre,i_cluster,1,5],col="orange",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,6],power_array[i_file,i_F,i_pre,i_cluster,1,6],col="green",pch=3,cex=3)
        
        legend("topright",c(names(power_matrix)),pch=rep(3,6),cex=1.5,col=c("red","blue","pink","brown","orange","green"))
          
        dev.off()
        
        print(paste0("p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag))
      }
    }
  }
}

saveRDS(pval_list,paste0("../Data_PRJNA434002/8.Result/final_pval_list.rds"))
saveRDS(power_array,paste0("../Data_PRJNA434002/8.Result/final_power_array.rds"))
saveRDS(ks_array,paste0("../Data_PRJNA434002/8.Result/final_ks_array.rds"))
saveRDS(cor_zerorate_ind_mean_array,paste0("../Data_PRJNA434002/8.Result/final_cor_zerorate_ind_mean_array.rds"))
saveRDS(cor_nonexpres_ind_array,paste0("../Data_PRJNA434002/8.Result/final_cor_nonexpres_ind_array.rds"))
saveRDS(cor_zerorate_array,paste0("../Data_PRJNA434002/8.Result/final_cor_zerorate_array.rds"))
saveRDS(cor_expression_array,paste0("../Data_PRJNA434002/8.Result/final_cor_expression_array.rds"))
saveRDS(cor_overdisp_max_array,paste0("../Data_PRJNA434002/8.Result/final_cor_overdisp_max_array.rds"))
saveRDS(cor_overdisp_median_array,paste0("../Data_PRJNA434002/8.Result/final_cor_overdisp_median_array.rds"))
saveRDS(cor_dropout_max_array,paste0("../Data_PRJNA434002/8.Result/final_cor_dropout_max_array.rds"))
saveRDS(cor_dropout_median_array,paste0("../Data_PRJNA434002/8.Result/final_cor_dropout_median_array.rds"))
rownames(zeros)=rownames_zeros
View(zeros)
saveRDS(zeros,paste0("../Data_PRJNA434002/10.Result/power_array_NAs_p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".rds"))


###############Power array analysis###################

power_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_power_array.rds"))
#do boxplot##############
for(i_file in 1:length(file_tag_seq)){
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      file_tag=file_tag_seq[i_file]
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]
      png(paste0("../Data_PRJNA434002/8.Result/boxplot_power_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 800,width = 500)
      op=par(mfrow=c(2,1))
      a=power_array[i_file,i_F,i_pre,,1,]
      b=power_array[i_file,i_F,i_pre,,2,]
      boxplot(a,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of p-values less than 0.05 of all clusters, observed data",ylab="power")
      abline(h = 0.05, col = "red") 
      boxplot(b,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of p-values less than 0.05 of all clusters, permutated data",ylab="type I error")
      abline(h = 0.05, col = "red") 
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_power_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))

      a=power_array[i_file,i_F,i_pre,,1,]
      b=power_array[i_file,i_F,i_pre,,2,]

      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="Power",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of pval<0.05,observed data",,ylim=c(0,1))
      abline(h = 0.05, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="Type I error",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of pval<0.05,permutated data",ylim=c(0,1))
      abline(h = 0.05, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed (Power)",xlab="permutated (Type I error)",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="proportion of pval<0.05, observed vs permutated",ylim=c(0,1),xlim=c(0,1))
      abline(v = 0.05, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))

      par(op)
      dev.off()
    }
  }
}

###############KS test array analysis###################

ks_array=readRDS(paste0("../Data_PRJNA434002/8.Result/final_ks_array.rds"))
#do boxplot##############
for(i_file in 1:length(file_tag_seq)){
  for(i_F in 1:length(F_method_seq)){
    for(i_pre in 1:length(pre_tag_seq)){
      file_tag=file_tag_seq[i_file]
      F_method=F_method_seq[i_F]
      pre_tag=pre_tag_seq[i_pre]
      png(paste0("../Data_PRJNA434002/8.Result/boxplot_ks_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 800,width = 500)
      op=par(mfrow=c(2,1))
      a=ks_array[i_file,i_F,i_pre,,1,]
      b=ks_array[i_file,i_F,i_pre,,2,]
      boxplot(a,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="ks test for distribution of pvalues, observed data",ylab="ks pval")
      abline(h = 0.05, col = "red") 
      boxplot(b,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="ks test for distribution of pvalues, permutated data",ylab="ks pval")
      abline(h = 0.05, col = "red") 
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_ks_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))

      a=-log10(ks_array[i_file,i_F,i_pre,,1,]+min(ks_array[ks_array>0],na.rm = TRUE))
      b=-log10(ks_array[i_file,i_F,i_pre,,2,]+min(ks_array[ks_array>0],na.rm = TRUE))
      max_ab=max(a,b,na.rm = TRUE)
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="-log10 ks pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="-log10 pvalues from KS test, observed data",,ylim=c(0,max_ab))
      abline(h = -log10(0.05), col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="-log10 ks pval",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="-log10 pvalues from KS test,permutated data",ylim=c(0,max_ab))
      abline(h = -log10(0.05), col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="-log10 pvalues from KS test, observed vs permutated",ylim=c(0,max_ab),xlim=c(0,max_ab))
      abline(v = -log10(0.05), col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}


#do scatter plot about cell number vs cor ##############

#zeros_mean
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_zerorate_ind_mean_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_zerorate_ind_mean_array[i_file,i_F,i_pre,,1,]
      b=cor_zerorate_ind_mean_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and zerorate_ind_mean, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and zerorate_ind_mean,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and zerorate_ind_mean, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}

#nonexpres_ind
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_nonexpres_ind_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_nonexpres_ind_array[i_file,i_F,i_pre,,1,]
      b=cor_nonexpres_ind_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and nonexpres_ind, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and nonexpres_ind,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and nonexpres_ind, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}



#zero_rate
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_zerorate_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_zerorate_array[i_file,i_F,i_pre,,1,]
      b=cor_zerorate_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and zerorate, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and zerorate,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and zerorate, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}

#expression level
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_expression_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_expression_array[i_file,i_F,i_pre,,1,]
      b=cor_expression_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and expression, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and expression,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and expression, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}

#overdisp_median level
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_overdisp_median_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_overdisp_median_array[i_file,i_F,i_pre,,1,]
      b=cor_overdisp_median_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and overdisp_median, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and overdisp_median,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and overdisp_median, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}


#overdisp_max level
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_overdisp_max_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_overdisp_max_array[i_file,i_F,i_pre,,1,]
      b=cor_overdisp_max_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and overdisp_max, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and overdisp_max,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and overdisp_max, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}

#dropout_median level
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_dropout_median_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_dropout_median_array[i_file,i_F,i_pre,,1,]
      b=cor_dropout_median_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and dropout_median, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and dropout_median,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and dropout_median, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}


#dropout_max level
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
      
      png(paste0("../Data_PRJNA434002/8.Result/scatter_dropout_max_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 1200,width = 500)
      op=par(mfrow=c(3,1))
      
      a=cor_dropout_max_array[i_file,i_F,i_pre,,1,]
      b=cor_dropout_max_array[i_file,i_F,i_pre,,2,]
      
      plot(cell_num,a[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and dropout_max, observed data",,ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,a[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,a[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,a[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,a[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(a)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(cell_num,b[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="correlation",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and dropout_max,permutated data",ylim=c(-1,1))
      abline(h = 0.0, col = "red") 
      points(cell_num,b[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,b[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,b[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,b[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,b[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      plot(b[,1],a[,1],type="p",pch=3,cex=1.1, col="red",ylab="observed",xlab="permutated",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between pval and dropout_max, observed vs permutated",ylim=c(-1,1),xlim=c(-1,1))
      abline(v = 0, col = "red") 
      abline(h = 0, col = "red") 
      points(b[,2],a[,2],pch=4,cex=1.1, col="blue")
      points(b[,3],a[,3],pch=5,cex=1.1, col="pink")
      points(b[,4],a[,4],pch=6,cex=1.1, col="brown")
      points(b[,5],a[,5],pch=7,cex=1.1, col="orange")
      points(b[,6],a[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(b)),pch=3:8,cex=1.5,col=c("red","blue","pink","brown","orange","green"))
      
      par(op)
      dev.off()
    }
  }
}



#do scatter plot about ob pval vs perm pval##############
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

      log_ob_pval=-log10(pval_array[i_F,i_pre,,1,,]+min(pval_array[pval_array>0],na.rm = TRUE))
      log_perm_pval=-log10(pval_array[i_F,i_pre,,2,,]+min(pval_array[pval_array>0],na.rm = TRUE))
      
      cor_obperm=array(dim=dim(log_ob_pval)[1:2],dimnames=list(dimnames(log_ob_pval)[[1]],dimnames(log_ob_pval)[[2]]))
      cor_obob=array(dim=dim(log_ob_pval)[c(1,2,2)],dimnames=list(dimnames(log_ob_pval)[[1]],dimnames(log_ob_pval)[[2]],dimnames(log_ob_pval)[[2]]))
      cor_permperm=array(dim=dim(log_ob_pval)[c(1,2,2)],dimnames=list(dimnames(log_ob_pval)[[1]],dimnames(log_ob_pval)[[2]],dimnames(log_ob_pval)[[2]]))
      for(ia in 1:dim(log_ob_pval)[1]){
        for(ib in 1:dim(log_ob_pval)[2]){
          cor_obperm[ia,ib]=cor(log_ob_pval[ia,ib,],log_perm_pval[ia,ib,],use = "complete.obs")
          cor_obob[ia,ib,]=apply(log_ob_pval[ia,,],1, function(x){return(cor(x,log_ob_pval[ia,ib,],use = "complete.obs"))})
          cor_permperm[ia,ib,]=apply(log_perm_pval[ia,,],1, function(x){return(cor(x,log_perm_pval[ia,ib,],use = "complete.obs"))})
        }
      }
      
      saveRDS(cor_obperm,paste0("../Data_PRJNA434002/8.Result/cor_pval_obperm_",F_method,"_",pre_tag,"_",file_tag,".rds"))
      saveRDS(cor_obob,paste0("../Data_PRJNA434002/8.Result/cor_pval_obob_",F_method,"_",pre_tag,"_",file_tag,".rds"))
      saveRDS(cor_permperm,paste0("../Data_PRJNA434002/8.Result/cor_pval_permperm_",F_method,"_",pre_tag,"_",file_tag,".rds"))
      #plot cor pval ob vs perm
      png(paste0("../Data_PRJNA434002/8.Result/scatter_cor_pval_",F_method,"_",pre_tag,"_",file_tag,".png"),height = 400,width = 500)
      plot(cell_num,cor_obperm[,1],type="p",pch=3,cex=1.1, col="red",xlab="cell number of each cell type",ylab="-log10 pval_cor",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main="correlation between observed and permutation -log10 pvalues",,ylim=c(-1,1))
      abline(h = -log10(0.05), col = "red") 
      points(cell_num,cor_obperm[,2],pch=4,cex=1.1, col="blue")
      points(cell_num,cor_obperm[,3],pch=5,cex=1.1, col="pink")
      points(cell_num,cor_obperm[,4],pch=6,cex=1.1, col="brown")
      points(cell_num,cor_obperm[,5],pch=7,cex=1.1, col="orange")
      points(cell_num,cor_obperm[,6],pch=8,cex=1.1, col="green")
      legend("topright",c(colnames(c)),pch=3:8,cex=.5,col=c("red","blue","pink","brown","orange","green"))
      
      dev.off()

      #plot first 4 smallest pval's gene's re-constructed expression distribution vs the median pvals
      
      
      for(i_cluster in 1:dim(log_ob_pval)[1]){
        cluster_tag=cur_cluster[i_cluster]
        fit_data=readRDS(paste0("../Data_PRJNA434002/7.Result/fit_ind_",pre_tag,"_sim_",i_cluster,"_",file_tag,".rds"))
        sim_data=readRDS(paste0("../Data_PRJNA434002/7.Result/sim_ind_",pre_tag,"_sim_",i_cluster,"_",file_tag,".rds"))
        cur_phenotype=tmeta$diagnosis[tmeta$cluster==cluster_tag]!="Control"
        cur_phenotype_ind=tmeta$diagnosis[match(unlist(dimnames(fit_data)[2]),tmeta$individual)]!="Control"
        png(paste0("../Data_PRJNA434002/8.Result/sig_gene_count_",F_method,"_",pre_tag,"_",i_cluster,"_",file_tag,".png"),height = 4500,width =6000)
        op=par(mfrow = c(12, 16))
        for(i_pval in 1:dim(log_ob_pval)[2]){
          #sig pval
          cur_index=order(log_ob_pval[i_cluster,i_pval,],decreasing = TRUE)[1:4]
          cur_fit=fit_data[cur_index,,]
          cur_sim=sim_data[cur_index,,]
          for(iplot in 1:4){
            
            hist(cur_fit[iplot,!cur_phenotype_ind,1],xlab="logmean",xlim=range(min(cur_fit[iplot,,1],na.rm = TRUE),max(cur_fit[iplot,,1],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("ind-zinb param,",colnames(cor_obperm)[i_pval]," sig pval ",iplot))
            hist(cur_fit[iplot,cur_phenotype_ind,1],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
            hist(cur_fit[iplot,!cur_phenotype_ind,2],xlab="dispersion",xlim=range(min(cur_fit[iplot,,2],na.rm = TRUE),max(cur_fit[iplot,,2],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("ind-zinb param,",colnames(cor_obperm)[i_pval]," sig pval ",iplot))
            hist(cur_fit[iplot,cur_phenotype_ind,2],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
            hist(cur_fit[iplot,!cur_phenotype_ind,3],xlab="dropout",xlim=range(min(cur_fit[iplot,,3],na.rm = TRUE),max(cur_fit[iplot,,3],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("ind-zinb param,",colnames(cor_obperm)[i_pval]," sig pval ",iplot))
            hist(cur_fit[iplot,cur_phenotype_ind,3],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
            hist(cur_sim[iplot,cur_phenotype,],xlab="read count",xlim=range(min(cur_sim[iplot,,],na.rm = TRUE),max(cur_sim[iplot,,],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("sim gene expres,",colnames(cor_obperm)[i_pval]," sig pval ",iplot))
            hist(cur_sim[iplot,!cur_phenotype,],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("case","control"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
          }
          
          #median pval
          cur_index=order(log_ob_pval[i_cluster,i_pval,],decreasing = TRUE)[1501:1504]
          cur_fit=fit_data[cur_index,,]
          cur_sim=sim_data[cur_index,,]
          for(iplot in 1:4){
            hist(cur_fit[iplot,!cur_phenotype_ind,1],xlab="logmean",xlim=range(min(cur_fit[iplot,,1],na.rm = TRUE),max(cur_fit[iplot,,1],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("ind-zinb param,",colnames(cor_obperm)[i_pval]," median pval ",(iplot+1500)))
            hist(cur_fit[iplot,cur_phenotype_ind,1],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
            hist(cur_fit[iplot,!cur_phenotype_ind,2],xlab="dispersion",xlim=range(min(cur_fit[iplot,,2],na.rm = TRUE),max(cur_fit[iplot,,2],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("ind-zinb param,",colnames(cor_obperm)[i_pval]," median pval ",(iplot+1500)))
            hist(cur_fit[iplot,cur_phenotype_ind,2],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
            hist(cur_fit[iplot,!cur_phenotype_ind,3],xlab="dropout",xlim=range(min(cur_fit[iplot,,3],na.rm = TRUE),max(cur_fit[iplot,,3],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("ind-zinb param,",colnames(cor_obperm)[i_pval]," median pval ",(iplot+1500)))
            hist(cur_fit[iplot,cur_phenotype_ind,3],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("control","case"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
            hist(cur_sim[iplot,cur_phenotype,],xlab="read count",xlim=range(min(cur_sim[iplot,,],na.rm = TRUE),max(cur_sim[iplot,,],na.rm = TRUE)),col=rgb(1,0,0,0.1),main=paste0("sim gene expres,",colnames(cor_obperm)[i_pval]," median pval ",(iplot+1500)))
            hist(cur_sim[iplot,!cur_phenotype,],add=T,col=rgb(0,0,1,0.1))
            legend("topright",c("case","control"),pch=15,cex=1.1,col=c(rgb(1,0,0,0.1),rgb(0,0,1,0.1)))
            
          }
        }
        par(op)
        dev.off()
      }
    }
  }
}

