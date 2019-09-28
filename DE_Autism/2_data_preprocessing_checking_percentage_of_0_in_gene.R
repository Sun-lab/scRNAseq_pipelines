#this code try to do some basic character detections on the data.

#Data_PRJNA434002

library("moments")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

# ##################Section 1:log2 transformed Expression #########################
# exprM=as.matrix(readRDS("../Data_PRJNA434002/exprMatrix.rds"))
# 
# exprM_freq=matrix(ncol=1,nrow=11)
# for(i in 0:10){
#   exprM_freq[i]=sum(exprM>(i-1) & exprM<=(i))
# }
# 
# exprM_freq
# 
# exprM_nozero_median=apply(exprM,1,function(x){return(median(x[x>0]))})
# exprM_nozero_mean=apply(exprM,1,function(x){return(mean(x[x>0]))})
# exprM_nozero_sd=apply(exprM,1,function(x){return(sd(x[x>0]))})
# exprM_nozero_skewness=apply(exprM,1,function(x){return(skewness(x[x>0]))})
# 
# 
# exprM_zero_rate=apply(exprM==0,1,function(x){return(sum(x)/length(x))})
# 
# cell_4_each_gene=apply(exprM>0,1,sum)
# gene_4_each_cell=apply(exprM>0,2,sum)
# min_cell_4_each_gene=min(cell_4_each_gene)
# min_gene_4_each_cell=min(gene_4_each_cell)
# 
# min_cell_4_each_gene
# min_gene_4_each_cell
# 
# pdf("histogram_exprMatrix_stat.pdf")
# hist(log10(cell_4_each_gene),breaks=50)
# hist(log10(gene_4_each_cell),breaks=50)
# hist(exprM_nozero_median,breaks=50)
# hist(exprM_nozero_mean,breaks=50)
# hist(exprM_nozero_sd,breaks=50)
# hist(exprM_nozero_skewness,breaks=50)
# hist(exprM_zero_rate,breaks=50,
#      sub=paste0("0rate==1: ",
#                 round(sum(exprM_zero_rate==1)/length(exprM_zero_rate),3),
#                 ", 0rate>0.99: ",   
#                 round(sum(exprM_zero_rate>0.99)/length(exprM_zero_rate),3),
#                 ", 0rate>0.98: ", 
#                 round(sum(exprM_zero_rate>0.98)/length(exprM_zero_rate),3)))
# 
# 
# dev.off()

##################Section 2:raw Expression #########################
library("GenomicFeatures")
library("GenomicAlignments")
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)


sce   = read10xCounts("../Data_PRJNA434002/rawMatrix/", col.names=TRUE)	

rawM=counts(sce) #matrix to large, so we are not able to convert it as regular matrix

rawM_freq=matrix(ncol=1,nrow=101)
for(i in 0:100){
  rawM_freq[i]=sum(rawM==i)
}

rawM_freq

rawM_nozero_median=matrix(ncol=1,nrow=nrow(rawM))
rawM_nozero_mean=matrix(ncol=1,nrow=nrow(rawM))
rawM_nozero_sd=matrix(ncol=1,nrow=nrow(rawM))
rawM_nozero_skewness=matrix(ncol=1,nrow=nrow(rawM))
rawM_zero_rate=matrix(ncol=1,nrow=nrow(rawM))
cell_4_each_gene=matrix(ncol=1,nrow=nrow(rawM))

for(i in 1:nrow(rawM)){
  cur_count=as.numeric(rawM[i,])
  cur_count=cur_count[cur_count>0]
  rawM_nozero_median[i]=median(cur_count)
  rawM_nozero_mean[i]=mean(cur_count)
  rawM_nozero_sd[i]=sd(cur_count)
  rawM_nozero_skewness[i]=skewness(cur_count)
  cell_4_each_gene[i]=length(cur_count)
  rawM_zero_rate[i]=(ncol(rawM)-cell_4_each_gene[i])/ncol(rawM)
}

for(i in 1:ncol(rawM)){
  cur_count=as.numeric(rawM[,i])
  gene_4_each_cell[i]=sum(cur_count>0)
}


min_cell_4_each_gene=min(cell_4_each_gene)
min_gene_4_each_cell=min(gene_4_each_cell)

min_cell_4_each_gene
min_gene_4_each_cell



pdf("histogram_rawMatrix_stat.pdf")
hist(rawM[rawM>0],breaks=50)
hist(log10(cell_4_each_gene),breaks=50)
hist(log10(gene_4_each_cell),breaks=50)
hist(rawM_nozero_median,breaks=50)
hist(rawM_nozero_mean,breaks=50)
hist(rawM_nozero_sd,breaks=50)
hist(rawM_nozero_skewness,breaks=50)
hist(rawM_zero_rate,breaks=50,
     sub=paste0("0rate==1: ",
                round(sum(rawM_zero_rate==1)/length(rawM_zero_rate),3),
                ", 0rate>0.99: ",   
                round(sum(rawM_zero_rate>0.99)/length(rawM_zero_rate),3),
                ", 0rate>0.98: ", 
                round(sum(rawM_zero_rate>0.98)/length(rawM_zero_rate),3)))


dev.off()



