#this code try to reproduce the Velmeshev_2019_autism and reproduce their TSNE plots,


#Data_PRJNA434002
#install.packages("RSpectra")
library("RSpectra")
library("Rtsne")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

exprM=as.matrix(readRDS("../Data_PRJNA434002/exprMatrix5k.rds"))
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")

k = 50
svd50=svds(exprM, k)
#svd50=readRDS("../Data_PRJNA434002/svd50.rds")

pdf("scree_plot.pdf",width = 6,height = 6)
plot(svd50$d^2/sum(svd50$d^2), xlim = c(0, k),ylim=c(0,0.01), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained", main="scree plot")
dev.off()

svd50$d^2/sum(svd50$d^2)

k=17
svd50v=svd50$v[,1:k]

xc=as.matrix(meta[,"Capbatch"])
xs=as.matrix(meta[,"Seqbatch"])
cor_c=matrix(ncol=1,nrow=k)
cor_s=matrix(ncol=1,nrow=k)

for(ik in 1:k){
  y=svd50v[,ik]
  ##calculate correlation of capbatch
  lmc = lm(y ~ as.factor(xc))
  cor_c[ik]=summary(lmc)$r.square
  ##calculate correlation of seqbatch
  lms = lm(y ~ as.factor(xs))
  cor_s[ik]=sqrt(summary(lms)$r.square)
  
}

cor_s
cor_c

# #cor 0.1 threshold 
# flag=(abs(cor_c)<0.1 & abs(cor_s)<0.1) 
# svd50v0.1=svd50v[,flag==1]
# tsne=Rtsne(svd50v0.1,dims=2, perplexity=10)
# 
# pdf("tSNE_plots_cor0.1_10.pdf",height = 8,width = 8)
# plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
# plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
# plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
# plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
# dev.off()

#cor 0.2 threshold 
flag=(abs(cor_c)<0.2 & abs(cor_s)<0.2) 
sum(flag)

svd50v0.2=svd50v[,flag==1]
tsne=Rtsne(svd50v0.2,dims=2, perplexity=10)

pdf("tSNE_plots_cor0.2_5k17_10.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()

tsne=Rtsne(svd50v0.2,dims=2, perplexity=15)

pdf("tSNE_plots_cor0.2_5k17_15.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()

tsne=Rtsne(svd50v0.2,dims=2, perplexity=20)

pdf("tSNE_plots_cor0.2_5k17_20.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()

tsne=Rtsne(svd50v0.2,dims=2, perplexity=30)

pdf("tSNE_plots_cor0.2_5k17_30.pdf",height = 8,width = 8)
plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=as.numeric(as.factor(meta[,"Capbatch"])))
plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
plot(tsne$Y, cex=.2,main="tSNE-cluster",col=as.numeric(as.factor(meta[,"cluster"])))
plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
dev.off()









