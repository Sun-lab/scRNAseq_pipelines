#this code try to do some basic character detections on the data.

#Data_PRJNA434002

library("moments")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

exprM=as.matrix(readRDS("../Data_PRJNA434002/exprMatrix.rds"))

exprM_nozero_median=apply(exprM,1,function(x){return(median(x[x>0]))})
exprM_nozero_mean=apply(exprM,1,function(x){return(mean(x[x>0]))})
exprM_nozero_sd=apply(exprM,1,function(x){return(sd(x[x>0]))})
exprM_nozero_skewness=apply(exprM,1,function(x){return(skewness(x[x>0]))})


exprM_zero_rate=apply(exprM==0,1,function(x){return(sum(x)/length(x))})

pdf("histogram_exprMatrix_stat.pdf")
hist(exprM_nozero_median,breaks=50)
hist(exprM_nozero_mean,breaks=50)
hist(exprM_nozero_sd,breaks=50)
hist(exprM_nozero_skewness,breaks=50)
hist(exprM_zero_rate,breaks=50,
     sub=paste0("0rate==1: ",
                round(sum(exprM_zero_rate=1)/length(exprM_zero_rate),3),
                ", 0rate>0.99: ",   
                round(sum(exprM_zero_rate>0.99)/length(exprM_zero_rate),3),
                ", 0rate>0.98: ", 
                round(sum(exprM_zero_rate>0.98)/length(exprM_zero_rate),3)))

dev.off()

min_cell_4_each_gene=min(apply(exprM>0,1,sum))
min_gene_4_each_cell=min(apply(exprM>0,2,sum))
min_cell_4_each_gene
min_gene_4_each_cell

