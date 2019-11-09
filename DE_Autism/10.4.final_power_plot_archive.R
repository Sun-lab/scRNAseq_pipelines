setwd("/Users/mzhang24/Desktop/fh/Data_PRJNA434002/10.Result/")



file_tag=1
sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
# #                         #splat.org--method 4, change the mean.shape and mean.rate originally
# # #zinb.naive--method 5, using naive zinb models to do so.
# # 



#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

perm_num=500
r_mean=1.5  #r_mean/r_var should < 1+mean.shape
r_var=4

mean_index=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_de.mean_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
var_index=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_de.var_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))


jsd_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_zinb_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
jsd_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_empirical_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
klmean_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_zinb_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
klmean_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_empirical_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
deseq2_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/DESeq2_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))



MAST_pval_ob=readRDS(paste0("../Data_PRJNA434002/10.Result/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_0.rds"))


MAST_pval_perm=matrix(ncol=perm_num,nrow=length(MAST_pval_ob))
for(i in 1:perm_num){
  if(file.exists(paste0("../Data_PRJNA434002/10.Result/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",i,".rds"))){
    MAST_pval_perm[,i]=readRDS(paste0("../Data_PRJNA434002/10.Result/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",i,".rds"))
  }
}
  
MAST_pval_perm_flag=apply(MAST_pval_perm,2,function(x){return(x-MAST_pval_ob)})
  
MAST_pval=apply(MAST_pval_perm_flag>=0,1,function(x){return(sum(x,na.rm = TRUE))})
MAST_pval_length=apply(MAST_pval_perm_flag,1,function(x){return(sum(!is.na(x)))})
MAST_pval=MAST_pval/MAST_pval_length



#histogram
pdf(paste0("../Data_PRJNA434002/10.Result/final_power_",sim_method,"_",file_tag,".pdf"),height = 16,width = 12)
op=par(mfrow = c(4, 3), pty = "s")

hist(jsd_zinb_pval[mean_index==1],breaks = 20)
hist(jsd_zinb_pval[var_index==1],breaks = 20)
hist(jsd_zinb_pval[mean_index==0 & var_index==0],breaks = 20)

hist(jsd_empirical_pval[mean_index==1],breaks = 20)
hist(jsd_empirical_pval[var_index==1],breaks = 20)
hist(jsd_empirical_pval[mean_index==0 & var_index==0],breaks = 20)

hist(klmean_zinb_pval[mean_index==1],breaks = 20)
hist(klmean_zinb_pval[var_index==1],breaks = 20)
hist(klmean_zinb_pval[mean_index==0 & var_index==0],breaks = 20)

hist(klmean_empirical_pval[mean_index==1],breaks = 20)
hist(klmean_empirical_pval[var_index==1],breaks = 20)
hist(klmean_empirical_pval[mean_index==0 & var_index==0],breaks = 20)

par(op)

#
op=par(mfrow = c(4, 3), pty = "s")
hist(jsd_empirical_pval[mean_index==1],breaks = 20)
hist(jsd_empirical_pval[var_index==1],breaks = 20)
hist(jsd_empirical_pval[mean_index==0 & var_index==0],breaks = 20)

hist(deseq2_pval[mean_index==1],breaks = 20)
hist(deseq2_pval[var_index==1],breaks = 20)
hist(deseq2_pval[mean_index==0 & var_index==0],breaks = 20)

hist(MAST_pval[mean_index==1],breaks = 20)
hist(MAST_pval[var_index==1],breaks = 20)
hist(MAST_pval[mean_index==0 & var_index==0],breaks = 20)

par(op)

#power calculation
cal_power=function(x,threshold){
  return(sum(x<threshold,na.rm = TRUE)/length(x))
}

power_matrix=matrix(ncol=3,nrow=6)
colnames(power_matrix)=c("mean_diff","var_diff","control(FDR)")
rownames(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb")
cur_pval=deseq2_pval
power_matrix[1,]=c(cal_power(cur_pval[mean_index==1],threshold = 0.05),
              cal_power(cur_pval[var_index==1],threshold = 0.05),
              cal_power(cur_pval[mean_index==0 & var_index==0],threshold = 0.05))

cur_pval=MAST_pval
power_matrix[2,]=c(cal_power(cur_pval[mean_index==1],threshold = 0.05),
              cal_power(cur_pval[var_index==1],threshold = 0.05),
              cal_power(cur_pval[mean_index==0 & var_index==0],threshold = 0.05))

cur_pval=jsd_empirical_pval
power_matrix[3,]=c(cal_power(cur_pval[mean_index==1],threshold = 0.05),
              cal_power(cur_pval[var_index==1],threshold = 0.05),
              cal_power(cur_pval[mean_index==0 & var_index==0],threshold = 0.05))

cur_pval=klmean_empirical_pval
power_matrix[4,]=c(cal_power(cur_pval[mean_index==1],threshold = 0.05),
              cal_power(cur_pval[var_index==1],threshold = 0.05),
              cal_power(cur_pval[mean_index==0 & var_index==0],threshold = 0.05))

cur_pval=jsd_zinb_pval
power_matrix[5,]=c(cal_power(cur_pval[mean_index==1],threshold = 0.05),
              cal_power(cur_pval[var_index==1],threshold = 0.05),
              cal_power(cur_pval[mean_index==0 & var_index==0],threshold = 0.05))

cur_pval=klmean_zinb_pval
power_matrix[6,]=c(cal_power(cur_pval[mean_index==1],threshold = 0.05),
              cal_power(cur_pval[var_index==1],threshold = 0.05),
              cal_power(cur_pval[mean_index==0 & var_index==0],threshold = 0.05))

#barplot
op=par(mfrow = c(4, 3), pty = "s")
barplot(power_matrix[1,],ylab="power",main=rownames(power_matrix)[1])
barplot(power_matrix[2,],ylab="power",main=rownames(power_matrix)[2])
barplot(power_matrix[3,],ylab="power",main=rownames(power_matrix)[3])
barplot(power_matrix[4,],ylab="power",main=rownames(power_matrix)[4])
barplot(power_matrix[5,],ylab="power",main=rownames(power_matrix)[5])
barplot(power_matrix[6,],ylab="power",main=rownames(power_matrix)[6])
par(op)

plot(power_matrix[1,3],power_matrix[1,1],xlim=c(0,0.1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
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

legend("topright",c(rownames(power_matrix),"mean diff","var diff"),pch=c(rep(15,6),4,3),cex=1,col=c("red","blue","pink","brown","orange","green","black","black"))
dev.off()