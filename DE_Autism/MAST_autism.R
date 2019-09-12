#this code try to reproduce the Velmeshev_2019_autism and reproduce their TSNE plots,
#perform the MAST DE analysis and calculate the permutated version of it.

#Data_PRJNA434002
#install.packages("RSpectra")
library("RSpectra")
library("Rtsne")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

exprM=as.matrix(readRDS("../Data_PRJNA434002/exprMatrix.rds"))
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")

# exprM100=as.matrix(exprM[,seq(1,1045)*100])
# saveRDS(exprM100,"../Data_PRJNA434002/exprMatrix100.rds")
# meta100=as.matrix(meta[seq(1,1045)*100,])
# saveRDS(meta100,"../Data_PRJNA434002/meta100.rds")

#exprM100=readRDS("../Data_PRJNA434002/exprMatrix100.rds")
#meta100=readRDS("../Data_PRJNA434002/meta100.rds")

k = 50
svd50=svds(exprM, k)

svd50=readRDS("../Data_PRJNA434002/svd50.rds")
yc=as.matrix(meta[,"Capbatch"])

ys=as.matrix(meta[,"Seqbatch"])

yck=unique(yc)
ysk=unique(ys)
cor_c=matrix(ncol=1,nrow=k)
cor_s=matrix(ncol=1,nrow=k)
for(ik in 1:k){
  x=svd50$v[,ik]
  ##calculate correlation of capbatch
  yck_mean=matrix(ncol=1,nrow=length(yck))
  cur_yc=yc
  for(iyck in 1:length(yck)){
    yck_mean[iyck]=mean(x[yc==yck[iyck]])
  }
  yck_label=order(yck_mean)
  for(iyck in 1:length(yck)){
    cur_yc[yc==yck[iyck]]=as.numeric(yck_label[iyck])
  }
  cor_c[ik]=cor(x,as.numeric(cur_yc),method="pearson")
  
  ##calculate correlation of seqbatch
  ysk_mean=matrix(ncol=1,nrow=length(ysk))
  cur_ys=ys
  for(iysk in 1:length(ysk)){
    ysk_mean[iysk]=mean(x[ys==ysk[iysk]])
  }
  ysk_label=order(ysk_mean)
  for(iysk in 1:length(ysk)){
    cur_ys[ys==ysk[iysk]]=ysk_label[iysk]
  }
  cor_s[ik]=cor(x,as.numeric(cur_ys),method="pearson")
  
}
View(cor_c)
View(cor_s)
flag=(abs(cor_c)<=0.2 & abs(cor_s)<=0.2) 

svd50v=svd50$v[,flag==1]

pdf("tSNE_plots.pdf",height = 8,width = 8)
tsne=Rtsne(svd50v,dims=2, perplexity=50)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=cur_yc)
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=cur_ys)
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()

