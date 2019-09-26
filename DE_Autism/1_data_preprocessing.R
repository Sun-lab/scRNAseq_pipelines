#this code convert the data from to .rds format and generate small subsets of them.

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")

#change the original file into rds for better read-writing-storage
exprM=read.table("../Data_PRJNA434002/exprMatrix.tsv.gz",header=TRUE,row.names=1)

saveRDS(exprM,"../Data_PRJNA434002/exprMatrix.rds")


#generate a 1% percentage cell subfiles, for testing
exprM100=as.matrix(exprM[,seq(1,1045)*100])
saveRDS(exprM100,"../Data_PRJNA434002/exprMatrix100.rds")

meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
meta100=meta[seq(1,1045)*100,]
saveRDS(meta100,"../Data_PRJNA434002/meta100.rds")

meta[1:10,1:10]
exprM[1:10,1:10]

#generate another 10% percentage cells with first 3000 genes, for testing
meta10=meta[seq(1,10455)*10,]
saveRDS(meta100,"../Data_PRJNA434002/meta10.rds")
exprM3k10=exprM[1:3000,seq(1,10455)*10]
saveRDS(exprM3k10,"../Data_PRJNA434002/exprMatrix3k10.rds")
exprM1k10=exprM[1:1000,seq(1,10455)*10]
saveRDS(exprM1k10,"../Data_PRJNA434002/exprMatrix1k10.rds")

#generate another 10% percentage cells with the 5000 lowest 0 rate genes, for tSNE plots
exprM_zero_rate=apply(exprM==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(exprM_zero_rate)

exprM5k=exprM[o_zero_rate[1:5000],]
exprM5k[1:10,1:10]
dim(exprM5k)
saveRDS(exprM5k,"../Data_PRJNA434002/exprMatrix5k.rds")

exprM1k=exprM[o_zero_rate[1:1000],]
exprM1k[1:10,1:10]
dim(exprM1k)
saveRDS(exprM1k,"../Data_PRJNA434002/exprMatrix1k.rds")

#read in raw data directly


library("GenomicFeatures")
library("GenomicAlignments")
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)


sce   = read10xCounts("../Data_PRJNA434002/rawMatrix/", col.names=TRUE)	

rawMatrix=counts(sce) #matrix to large, so we are not able to convert it as regular matrix
saveRDS(rawMatrix,"../Data_PRJNA434002/rawM.rds")

rawM1k10=as.matrix(rawMatrix[1:1000,seq(1,10455)*10])
saveRDS(rawM1k10,"../Data_PRJNA434002/rawM1k10.rds")
write.table(rawM1k10,"../Data_PRJNA434002/rawM1k10.csv",sep=",")


rawM3k10=as.matrix(rawMatrix[1:3000,seq(1,10455)*10])
saveRDS(rawM3k10,"../Data_PRJNA434002/rawM3k10.rds")
write.table(rawM3k10,"../Data_PRJNA434002/rawM3k10.csv",sep=",")


rawM_zero_rate=apply(rawMatrix==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(rawM_zero_rate)

rawM5k=as.matrix(rawMatrix[o_zero_rate[1:5000],])
rawM5k[1:10,1:10]
dim(rawM5k)
saveRDS(rawM5k,"../Data_PRJNA434002/rawM5k.rds")
write.table(rawM5k,"../Data_PRJNA434002/rawM5k.csv",sep=",")

rawM1k=as.matrix(rawMatrix[o_zero_rate[1:1000],])
rawM1k[1:10,1:10]
dim(rawM1k)
saveRDS(rawM1k,"../Data_PRJNA434002/rawM1k.rds")
write.table(rawM1k,"../Data_PRJNA434002/rawM1k.csv",sep=",")

sessionInfo()
q(save="no")
