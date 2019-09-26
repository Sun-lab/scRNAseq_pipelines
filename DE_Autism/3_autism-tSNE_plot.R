#this code try to reproduce the Velmeshev_2019_autism and reproduce their TSNE plots,

cur_k=17
cur_file="1k"
cor_thres=0.2
#Data_PRJNA434002
#install.packages("RSpectra")
library("RSpectra")
library("Rtsne")

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

exprM=as.matrix(readRDS(paste0("../Data_PRJNA434002/exprMatrix",cur_file,".rds")))
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")


k = 50
svd50=svds(exprM, k)
#svd50=readRDS("../Data_PRJNA434002/svd50.rds")

pdf(paste0("scree_plot.pdf"),width = 6,height = 6)
plot(svd50$d^2/sum(svd50$d^2), xlim = c(0, k),ylim=c(0,0.01), type = "b", pch = 16, xlab = "principal components", 
     ylab = "variance explained", main="scree plot")
dev.off()

svd50$d^2/sum(svd50$d^2)

k=cur_k ###here should depend on the plots
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



#cor 0.2 threshold 
flag=(abs(cor_c)<cor_thres & abs(cor_s)<cor_thres) 
sum(flag)

cur_svd50v=svd50v[,flag==1]



#generate plotting color

library(RColorBrewer)
color_type_num=apply(meta,2,function(x){return(length(table(as.factor(x))))})
n=max(color_type_num[color_type_num<50])
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'seq',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



# ####This part is used for detecting which color provide more distinguish tSNE plots.####
# tsne=Rtsne(cur_svd50v,dims=2, perplexity=15)
# 
# pdf(paste0("tSNE_color_test_",cur_file,"_k=",cur_k,".pdf"),height = 8,width = 8)
# for(i in 1:30){
#   set.seed(i) #2,8,28
#   cust_col=sample(col_vector, n)
#   #pie(rep(1,n), col=cust_col)
#   plot(tsne$Y, cex=.2,main=paste0("tSNE-cluster",i),col=cust_col[as.numeric(as.factor(meta[,"cluster"]))])
# }
# dev.off()

#Then,we choose seed 28 as the best choice for the color generation.
set.seed(28) #2,8,28
cust_col=sample(col_vector, n)

#TSNE_plot is a function to do 2 pdf based tsne plots, designed especially for current situation.
#example
#tsne_obj=tsne
#meta_info=meta
#file_label="k17_cor0.2_10"

TSNE_plot=function(tsne_obj,meta_info,file_label){
  
  png(paste0("tSNE_plots_cluster_",file_label,".png"),height = 800,width = 800)
  plot(tsne$Y, cex=.2,main="tSNE-cluster",col=cust_col[as.numeric(as.factor(meta[,"cluster"]))])
  dev.off()
  
  png(paste0("tSNE_plots_Capbatch_",file_label,".png"),height = 800,width = 800)
  plot(tsne$Y, cex=.1,main="tSNE-Capbatch",col=cust_col[as.numeric(as.factor(meta[,"Capbatch"]))])
  dev.off()
  png(paste0("tSNE_plots_Seqbatch_",file_label,".png"),height = 800,width = 800)
  plot(tsne$Y, cex=.1,main="tSNE-Seqbatch",col=as.numeric(as.factor(meta[,"Seqbatch"])))
  dev.off()
  png(paste0("tSNE_plots_individual_",file_label,".png"),height = 800,width = 800)
  plot(tsne$Y, cex=.2,main="tSNE-individual",col=as.numeric(as.factor(meta[,"individual"])))
  dev.off()
  
  png(paste0("tSNE_sub_cluster_",file_label,".png"),height = 3000,width = 6000)
  op=par(mfrow=c(3,6),mar=c(3, 3, 1, 1), bty="n",cex=0.9)
  for(i in 1:length(as.character(unique(meta$cluster)))){
    cur_cluster=as.character(unique(meta$cluster)[i])
    metaflag=(meta$cluster==cur_cluster)
    #pie(rep(1,n), col=cust_col)
    plot(tsne$Y, cex=.2,main=paste0("tSNE-cluster",i),col=cust_col[as.numeric(as.factor(metaflag))*3-1])
  }
  par(op)
  dev.off()
}


#tsne=Rtsne(cur_svd50v,dims=2, perplexity=15)
#TSNE_plot(tsne,meta,paste0(cur_file,"_k",cur_k,"_cor",cor_thres,"_15"))

tsne=Rtsne(cur_svd50v,dims=2, perplexity=30)
TSNE_plot(tsne,meta,paste0(cur_file,"_k",cur_k,"_cor",cor_thres,"_30"))

tsne=Rtsne(cur_svd50v,dims=2, perplexity=45)
TSNE_plot(tsne,meta,paste0(cur_file,"_k",cur_k,"_cor",cor_thres,"_45"))

sessionInfo()
q(save="no")






