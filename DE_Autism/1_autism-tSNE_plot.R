#this code try to reproduce the Velmeshev_2019_autism and reproduce their TSNE plots,


#Data_PRJNA434002
#install.packages("RSpectra")
library("RSpectra")
library("Rtsne")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

exprM=as.matrix(readRDS("../Data_PRJNA434002/exprMatrix.rds"))
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")

k = 50
svd50=svds(exprM, k)
#svd50=readRDS("../Data_PRJNA434002/svd50.rds")

xc=as.matrix(meta[,"Capbatch"])
xs=as.matrix(meta[,"Seqbatch"])
cor_c=matrix(ncol=1,nrow=k)
cor_s=matrix(ncol=1,nrow=k)
for(ik in 1:k){
  y=svd50$v[,ik]
  ##calculate correlation of capbatch
  lmc = lm(y ~ as.factor(xc))
  cor_c[ik]=summary(lmc)$r.square
  ##calculate correlation of seqbatch
  lms = lm(y ~ as.factor(xs))
  cor_s[ik]=summary(lms)$r.square
  
}

#cor 0.05 threshold perplexity 30
flag=(abs(cor_c)<0.05 & abs(cor_s)<0.05) 
svd50v=svd50$v[,flag==1]
tsne=Rtsne(svd50v,dims=2, perplexity=30)

pdf("tSNE_plots_cor0.05xperp30.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()

#cor 0.1 threshold perplexity 40
flag=(abs(cor_c)<0.1 & abs(cor_s)<0.1) 
svd50v=svd50$v[,flag==1]
tsne=Rtsne(svd50v,dims=2, perplexity=40)

pdf("tSNE_plots_cor0.2xperp50.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()

#cor 0.2 threshold perplexity 50
flag=(abs(cor_c)<0.2 & abs(cor_s)<0.2) 
svd50v=svd50$v[,flag==1]
tsne=Rtsne(svd50v,dims=2, perplexity=50)

pdf("tSNE_plots_cor0.2xperp50.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()
