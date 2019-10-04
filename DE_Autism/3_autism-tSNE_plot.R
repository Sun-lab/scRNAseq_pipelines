#this code try to reproduce the Velmeshev_2019_autism and reproduce their TSNE plots,

cur_k=17
cur_file="exprM5k"
cor_thres=0.2
file_label=paste0("k",cur_k,"_cor",cor_thres,"_",cur_file)
#Data_PRJNA434002
#install.packages("RSpectra")
library("RSpectra")
library("Rtsne")




#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

inputM=as.matrix(readRDS(paste0("../Data_PRJNA434002/",cur_file,".rds")))
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")


k = 50
svd50=svds(inputM, k)
#svd50=readRDS("../Data_PRJNA434002/svd50.rds")
pdf(paste0("scree_plot_",file_label,".pdf"),width = 6,height = 6)
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
# set.seed(28) #2,8,28
# cust_col=sample(col_vector, n)
cust_col=c("#023858", "#E31A1C", "#F768A1", "#BD0026", "#D4B9DA", "#7FCDBB", "#CE1256",
           "#88419D", "#FDD0A2", "#4D004B", "#E7298A", "#78C679", "#D9F0A3", "#081D58", 
           "#993404", "#CC4C02", "#FC9272", "#F7FCFD", "#BCBDDC", "#FFEDA0", "#FEE0D2",
           "#D0D1E6", "#7F0000", "#FFF7F3", "#9E9AC8", "#FFFFD9", "#CCEBC5", "#FFFFE5",
           "#014636", "#DADAEB", "#BFD3E6", "#FE9929", "#C994C7", "#FEE8C8", "#FCC5C0",
           "#1D91C0", "#FCFBFD", "#225EA8", "#000000", "#FEC44F", "#41AE76")
  
#TSNE_plot is a function to do 2 pdf based tsne plots, designed especially for current situation.
#example
#tsne_obj=tsne
#meta_info=meta
#file_label="k17_cor0.2_3k10"

TSNE_plot=function(tsne_obj,meta_info,file_label){
  
  ncluster=length(as.character(unique(meta_info$cluster)))
  
  #find label location
  medianY=matrix(ncol=2,nrow=ncluster)
  for(i in 1:ncluster){
    cur_cluster=as.character(unique(meta_info$cluster)[i])
    cur_Y=tsne_obj$Y[(meta_info$cluster==cur_cluster),]
    medianY[i,]=apply(cur_Y,2,median)
  }
  #plot compact tsne
  pdf(paste0("tsne_plots_",file_label,".pdf"),height = 8,width = 8)
  plot(tsne_obj$Y, cex=.2,main="tSNE-cluster",col=cust_col[as.numeric(as.factor(meta_info[,"cluster"]))])
  for(i in 1:ncluster){
    text(medianY[i,1],medianY[i,2],as.character(unique(meta_info$cluster)[i]))
  }
  plot(tsne_obj$Y, cex=.1,main="tSNE-Capbatch",col=cust_col[as.numeric(as.factor(meta_info[,"Capbatch"]))])
  for(i in 1:ncluster){
    text(medianY[i,1],medianY[i,2],as.character(unique(meta_info$cluster)[i]))
  }
  plot(tsne_obj$Y, cex=.1,main="tSNE-Seqbatch",col=cust_col[as.numeric(as.factor(meta_info[,"Seqbatch"]))])
  for(i in 1:ncluster){
    text(medianY[i,1],medianY[i,2],as.character(unique(meta_info$cluster)[i]))
  }
  plot(tsne_obj$Y, cex=.2,main="tSNE-individual",col=cust_col[as.numeric(as.factor(meta_info[,"individual"]))])
  for(i in 1:ncluster){
    text(medianY[i,1],medianY[i,2],as.character(unique(meta_info$cluster)[i]))
  }
  dev.off()
  
  #plot each sub tSNE
  pdf(paste0("tSNE_sub_cluster_",file_label,".pdf"),height = 15,width = 30)
  op=par(mfrow=c(3,6),mar=c(3, 3, 1, 1), bty="n",cex=0.9)
  for(i in 1:ncluster){
    cur_cluster=as.character(unique(meta_info$cluster)[i])
    metaflag=(meta_info$cluster==cur_cluster)
    #pie(rep(1,n), col=cust_col)
    plot(tsne_obj$Y, cex=.2,main=paste0("tSNE-cluster",cur_cluster),col=cust_col[as.numeric(as.factor(metaflag))*13+6])
  }
  par(op)
  dev.off()
}

tsne=Rtsne(cur_svd50v,dims=2, perplexity=30)
saveRDS(tsne,paste0("tsne_",cur_file,"_k",cur_k,"_cor",cor_thres,"_30.rds"))
TSNE_plot(tsne,meta,paste0(cur_file,"_k",cur_k,"_cor",cor_thres,"_30"))


sessionInfo()
q(save="no")






