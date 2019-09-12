#this code convert the data from to .rds format and generate small subsets of them.

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


exprM=read.table("../Data_PRJNA434002/exprMatrix.tsv.gz",header=TRUE,row.names=1)
exprM=as.matrix(exprM)
saveRDS(exprM,"../Data_PRJNA434002/exprMatrix.rds")

exprM100=as.matrix(exprM[,seq(1,1045)*100])
saveRDS(exprM100,"../Data_PRJNA434002/exprMatrix100.rds")

meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
meta100=meta[seq(1,1045)*100,]
saveRDS(meta100,"../Data_PRJNA434002/meta100.tsv")
