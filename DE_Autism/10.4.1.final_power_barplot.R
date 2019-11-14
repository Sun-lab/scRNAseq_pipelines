
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


param_tag=1

if(param_tag==1){
  r_mean_seq=c(1.2,1.5,2,4,6)
  r_var_seq=2
  
}
if(param_tag==2){
  r_mean_seq=1.5
  r_var_seq=c(1.2,1.5,2,4,6)
}
#r_mean_seq=1.5
#r_var_seq=1.5

sim_method_seq=c("splat.org","zinb.naive") #"splat.mean","splat.var"
file_tag_seq=1:5


# #test
# perm_num=500
# r_mean=1.5  #r_mean/r_var should < 1+mean.shape
# r_var=1.5
# file_tag=1
# sim_method="zinb.naive"


###################functions###################

#power calculation
cal_power=function(x,threshold){
  return(sum(x<threshold,na.rm = TRUE)/length(x))
}
###############################################
power_array=array(dim=c(
                    length(file_tag_seq),
                    length(sim_method_seq),
                    length(r_mean_seq),
                    length(r_var_seq),
                    6,3),
                  dimnames = list(
                    file_tag_seq,
                    sim_method_seq,
                    r_mean_seq,
                    r_var_seq,
                    c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb"),
                    c("mean_diff","var_diff","control(FDR)")
                  ))


for(i_file in 1:length(file_tag_seq)){
  for(i_sim in 1:length(sim_method_seq)){
    for(i_mean in 1:length(r_mean_seq)){
      for(i_var in 1:length(r_var_seq)){
    
        file_tag=file_tag_seq[i_file]
        sim_method=sim_method_seq[i_sim]
        r_mean=r_mean_seq[i_mean]
        r_var=r_var_seq[i_var]
        
        mean_index=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_de.mean_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
        var_index=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_de.var_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
        
        #tryCatch({     }, error = function(e) {NA} )
        
        jsd_zinb_pval=NA
        jsd_empirical_pval=NA
        klmean_zinb_pval=NA
        klmean_empirical_pval=NA
        deseq2_pval=NA
        MAST_pval=NA
        
        tryCatch({jsd_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_zinb_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))}, error = function(e) {NA} )
        tryCatch({jsd_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_empirical_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))}, error = function(e) {NA} )
        tryCatch({klmean_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_zinb_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))}, error = function(e) {NA} )
        tryCatch({klmean_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_empirical_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))}, error = function(e) {NA} )
        tryCatch({deseq2_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/DESeq2_ob_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))}, error = function(e) {NA} )
        
        #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
        tryCatch({MAST_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_0.rds"))}, error = function(e) {NA} )
        
        
        
        
        #histogram
        pdf(paste0("../Data_PRJNA434002/10.Result/final_power_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".pdf"),height = 16,width = 12)
        op=par(mfrow = c(4, 3), pty = "s")
        
        
        tryCatch({hist(jsd_zinb_pval[mean_index==1],main="pval of mean-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(jsd_zinb_pval[var_index==1],main="pval of var-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(jsd_zinb_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        tryCatch({hist(jsd_empirical_pval[mean_index==1],main="pval of mean-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(jsd_empirical_pval[var_index==1],main="pval of var-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(jsd_empirical_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        tryCatch({hist(klmean_zinb_pval[mean_index==1],main="pval of mean-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(klmean_zinb_pval[var_index==1],main="pval of var-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(klmean_zinb_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        tryCatch({hist(klmean_empirical_pval[mean_index==1],main="pval of mean-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(klmean_empirical_pval[var_index==1],main="pval of var-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(klmean_empirical_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        par(op)
        
        #
        op=par(mfrow = c(4, 3), pty = "s")
        tryCatch({hist(jsd_empirical_pval[mean_index==1],main="pval of mean-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(jsd_empirical_pval[var_index==1],main="pval of var-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(jsd_empirical_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        tryCatch({hist(deseq2_pval[mean_index==1],main="pval of mean-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(deseq2_pval[var_index==1],main="pval of var-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(deseq2_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        tryCatch({hist(MAST_pval[mean_index==1],main="pval of mean-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(MAST_pval[var_index==1],main="pval of var-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        tryCatch({hist(MAST_pval[mean_index==0 & var_index==0],main="pval of non-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
        
        par(op)
        
        
        
        power_matrix=matrix(nrow=6,ncol=3)
        
        rownames(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
        colnames(power_matrix)=c("mean_diff","var_diff","control(FDR)")
        
        tryCatch({power_matrix[1,]=c(cal_power(deseq2_pval[mean_index==1],threshold = 0.05),
                                    cal_power(deseq2_pval[var_index==1],threshold = 0.05),
                                    cal_power(deseq2_pval[mean_index==0 & var_index==0],threshold = 0.05))}, error = function(e) {NA} )
        
        tryCatch({power_matrix[2,]=c(cal_power(MAST_pval[mean_index==1],threshold = 0.05),
                                    cal_power(MAST_pval[var_index==1],threshold = 0.05),
                                    cal_power(MAST_pval[mean_index==0 & var_index==0],threshold = 0.05))}, error = function(e) {NA} )
        
        tryCatch({power_matrix[3,]=c(cal_power(jsd_empirical_pval[mean_index==1],threshold = 0.05),
                                    cal_power(jsd_empirical_pval[var_index==1],threshold = 0.05),
                                    cal_power(jsd_empirical_pval[mean_index==0 & var_index==0],threshold = 0.05))}, error = function(e) {NA} )
        
        tryCatch({power_matrix[4,]=c(cal_power(klmean_empirical_pval[mean_index==1],threshold = 0.05),
                                    cal_power(klmean_empirical_pval[var_index==1],threshold = 0.05),
                                    cal_power(klmean_empirical_pval[mean_index==0 & var_index==0],threshold = 0.05))}, error = function(e) {NA} )
        
        tryCatch({power_matrix[5,]=c(cal_power(jsd_zinb_pval[mean_index==1],threshold = 0.05),
                                    cal_power(jsd_zinb_pval[var_index==1],threshold = 0.05),
                                    cal_power(jsd_zinb_pval[mean_index==0 & var_index==0],threshold = 0.05))}, error = function(e) {NA} )
        
        tryCatch({power_matrix[6,]=c(cal_power(klmean_zinb_pval[mean_index==1],threshold = 0.05),
                                    cal_power(klmean_zinb_pval[var_index==1],threshold = 0.05),
                                    cal_power(klmean_zinb_pval[mean_index==0 & var_index==0],threshold = 0.05))}, error = function(e) {NA} )
        
        #barplot
        op=par(mfrow = c(4, 3), pty = "s")
        barplot(power_matrix[1,],ylab="power",main=rownames(power_matrix)[1],ylim=c(0,1))
        barplot(power_matrix[2,],ylab="power",main=rownames(power_matrix)[2],ylim=c(0,1))
        barplot(power_matrix[3,],ylab="power",main=rownames(power_matrix)[3],ylim=c(0,1))
        barplot(power_matrix[4,],ylab="power",main=rownames(power_matrix)[4],ylim=c(0,1))
        barplot(power_matrix[5,],ylab="power",main=rownames(power_matrix)[5],ylim=c(0,1))
        barplot(power_matrix[6,],ylab="power",main=rownames(power_matrix)[6],ylim=c(0,1))
        par(op)
        
        plot(power_matrix[1,3],power_matrix[1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
        points(power_matrix[1,3],power_matrix[1,2],col="red",pch=4,cex=3)
        points(power_matrix[2,3],power_matrix[2,1],col="blue",pch=3,cex=3)
        points(power_matrix[2,3],power_matrix[2,2],col="blue",pch=4,cex=3)
        points(power_matrix[3,3],power_matrix[3,1],col="pink",pch=3,cex=3)
        points(power_matrix[3,3],power_matrix[3,2],col="pink",pch=4,cex=3)
        points(power_matrix[4,3],power_matrix[4,1],col="brown",pch=3,cex=3)
        points(power_matrix[4,3],power_matrix[4,2],col="brown",pch=4,cex=3)
        points(power_matrix[5,3],power_matrix[5,1],col="orange",pch=3,cex=3)
        points(power_matrix[5,3],power_matrix[5,2],col="orange",pch=4,cex=3)
        points(power_matrix[6,3],power_matrix[6,1],col="green",pch=3,cex=3)
        points(power_matrix[6,3],power_matrix[6,2],col="green",pch=4,cex=3)
        
        legend("topright",c(rownames(power_matrix),"mean diff","var diff"),pch=c(rep(15,6),3,4),cex=1,col=c("red","blue","pink","brown","orange","green","black","black"))
        
        dev.off()
        
        
        power_array[i_file,i_sim,i_mean,i_var,,]=power_matrix
        
        print(c(file_tag,sim_method,r_mean,r_var))
        
      }
    }
  }
}

saveRDS(power_array,paste0("../Data_PRJNA434002/10.Result/final_power_array_param",param_tag,".rds"))


#more plot
pdf(paste0("../Data_PRJNA434002/10.Result/final_power_",sim_method,"_param",param_tag,"_",file_tag,".pdf"),height = 16,width = 12)
for(i_sim in length(sim_method_seq):1){
  for(i_file in 1:length(file_tag_seq)){
    file_tag=file_tag_seq[i_file]
    sim_method=sim_method_seq[i_sim]
    cur_power_array=power_array[i_file,i_sim,,,,]
    plot(cur_power_array[,1,3],cur_power_array[,1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3,main=paste0("power scatter: ",file_tag,"_",sim_method))
    points(cur_power_array[,1,3],cur_power_array[,1,2],col="red",pch=4,cex=3)
    points(cur_power_array[,2,3],cur_power_array[,2,1],col="blue",pch=3,cex=3)
    points(cur_power_array[,2,3],cur_power_array[,2,2],col="blue",pch=4,cex=3)
    points(cur_power_array[,3,3],cur_power_array[,3,1],col="pink",pch=3,cex=3)
    points(cur_power_array[,3,3],cur_power_array[,3,2],col="pink",pch=4,cex=3)
    points(cur_power_array[,4,3],cur_power_array[,4,1],col="brown",pch=3,cex=3)
    points(cur_power_array[,4,3],cur_power_array[,4,2],col="brown",pch=4,cex=3)
    points(cur_power_array[,5,3],cur_power_array[,5,1],col="orange",pch=3,cex=3)
    points(cur_power_array[,5,3],cur_power_array[,5,2],col="orange",pch=4,cex=3)
    points(cur_power_array[,6,3],cur_power_array[,6,1],col="green",pch=3,cex=3)
    points(cur_power_array[,6,3],cur_power_array[,6,2],col="green",pch=4,cex=3)
    
    legend("topright",c(dimnames(cur_power_array)[[2]],"mean diff","var diff"),pch=c(rep(15,6),3,4),cex=1,col=c("red","blue","pink","brown","orange","green","black","black"))
    
  }
}
dev.off()


