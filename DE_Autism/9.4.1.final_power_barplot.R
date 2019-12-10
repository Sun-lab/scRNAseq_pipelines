
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

count=1
zeros=matrix(ncol=6,nrow=length(file_tag_seq)*length(F_method_seq)*length(pre_tag_seq)*length(cluster_tag_seq)*length(perm_label_seq))
rownames_zeros=matrix(ncol=1,nrow=length(file_tag_seq)*length(F_method_seq)*length(pre_tag_seq)*length(cluster_tag_seq)*length(perm_label_seq))
colnames(zeros)=c("deseq2_pval","MAST_pval","jsd_empirical_pval","klmean_empirical_pval","jsd_zinb_pval","klmean_zinb_pval")

for(i_file in 1:length(file_tag_seq)){
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
          
        }
        #barplot
        png(paste0("../Data_PRJNA434002/8.Result/barplot_p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 1200,width = 800)
        op=par(mfrow = c(3, 2), pty = "s")
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,1],ylab="power",main=names(power_matrix)[1],ylim=c(0,1))
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,2],ylab="power",main=names(power_matrix)[2],ylim=c(0,1))
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,3],ylab="power",main=names(power_matrix)[3],ylim=c(0,1))
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,4],ylab="power",main=names(power_matrix)[4],ylim=c(0,1))
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,5],ylab="power",main=names(power_matrix)[5],ylim=c(0,1))
        barplot(power_array[i_file,i_F,i_pre,i_cluster,,6],ylab="power",main=names(power_matrix)[6],ylim=c(0,1))
        par(op)
        dev.off()
        
        png(paste0("../Data_PRJNA434002/8.Result/power_point_",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".png"),height = 600,width = 600)
        plot(power_array[i_file,i_F,i_pre,i_cluster,2,1],power_array[i_file,i_F,i_pre,i_cluster,1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
        
        points(power_array[i_file,i_F,i_pre,i_cluster,2,2],power_array[i_file,i_F,i_pre,i_cluster,1,2],col="blue",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,3],power_array[i_file,i_F,i_pre,i_cluster,1,3],col="pink",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,4],power_array[i_file,i_F,i_pre,i_cluster,1,4],col="brown",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,5],power_array[i_file,i_F,i_pre,i_cluster,1,5],col="orange",pch=3,cex=3)
        points(power_array[i_file,i_F,i_pre,i_cluster,2,6],power_array[i_file,i_F,i_pre,i_cluster,1,6],col="green",pch=3,cex=3)
        
        legend("topright",c(names(power_matrix)),pch=rep(3,6),cex=1,col=c("red","blue","pink","brown","orange","green"))
          
        dev.off()
          
        print(paste0("p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag))
      }
    }
  }
}


saveRDS(power_array,paste0("../Data_PRJNA434002/8.Result/final_power_array.rds"))

rownames(zeros)=rownames_zeros
View(zeros)
saveRDS(zeros,paste0("../Data_PRJNA434002/10.Result/power_array_NAs_p",perm_label,"_",F_method,"_",pre_tag,"_",cluster_tag,"_",file_tag,".rds"))
