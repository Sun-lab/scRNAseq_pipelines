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
library(scran)


sce   = read10xCounts("../Data_PRJNA434002/rawMatrix/", col.names=TRUE)	


#gene filter
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")	

attr.string = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name')	
attr.string = c(attr.string, 'start_position', 'end_position', 'strand')	
attr.string = c(attr.string, 'description', 'percentage_gene_gc_content')	
attr.string = c(attr.string, 'gene_biotype')	

rowData(sce)[1:2,]	
gene.annotation = getBM(attributes=attr.string, 	
                        filters =  'ensembl_gene_id', 	
                        values = rowData(sce)$ID, 	
                        mart = ensembl)	

dim(gene.annotation)	
gene.annotation[1:2,]	

t1 = table(gene.annotation$ensembl_gene_id)	
t2 = t1[t1 > 1]	
t2 	

gene.annotation[which(gene.annotation$ensembl_gene_id %in% names(t2)),]	
gene.annotation = distinct(gene.annotation, ensembl_gene_id, 	
                           .keep_all = TRUE)	
dim(gene.annotation)	
gene.annotation[1:2,]	
table(gene.annotation$chromosome_name)	
table(gene.annotation$gene_biotype)	

## some genes do not have annotation because their ids are retired	
gene.missing = dplyr::setdiff(rowData(sce)$ID, gene.annotation$ensembl_gene_id)	
length(gene.missing)	
gene.missing[1:6]	

w2kp = match(gene.annotation$ensembl_gene_id, rowData(sce)$ID)	
sce  = sce[w2kp,]	
dim(sce)	

table(gene.annotation$ensembl_gene_id == rowData(sce)$ID)	
rowData(sce)  = gene.annotation	
rownames(sce) = scater::uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, 	
                                             rowData(sce)$hgnc_symbol)	


###########generate the raw count tables


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


####################generate the normalization tables#############
zero_rate=apply(counts(sce)==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(zero_rate)

#gene normalization
date()
clusters = quickCluster(sce, min.mean=0.1, method="igraph")
date()
sce      = computeSumFactors(sce, cluster=clusters, min.mean=0.1)
date()
summary(sizeFactors(sce))
sce = normalize(sce)

#saveRDS(sce,"../Data_PRJNA434002/rawMnorm.rds")
sce=readRDS("../Data_PRJNA434002/rawMnorm.rds")

#writeout
top5k=as.numeric(o_zero_rate[1:5000])
sce2=counts(sce)[top5k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm5k=as.matrix(sce2)
rawMnorm5k[1:10,1:10]
dim(rawMnorm5k)
saveRDS(rawMnorm5k,"../Data_PRJNA434002/rawMnorm5k.rds")
write.table(rawMnorm5k,"../Data_PRJNA434002/rawMnorm5k.csv",sep=",")


top3k=as.numeric(o_zero_rate[1:3000])
sce2=counts(sce)[top3k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm3k10=as.matrix(sce2)
rawMnorm3k10=as.matrix(rawMnorm3k10[,seq(1,10455)*10])
rawMnorm3k10[1:10,1:10]
dim(rawMnorm3k10)
saveRDS(rawMnorm3k10,"../Data_PRJNA434002/rawMnorm3k10.rds")
write.table(rawMnorm3k10,"../Data_PRJNA434002/rawMnorm3k10.csv",sep=",")




top1k=as.numeric(o_zero_rate[1:1000])
sce2=counts(sce)[top1k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm1k=as.matrix(sce2)
rawMnorm1k[1:10,1:10]
dim(rawMnorm1k)
saveRDS(rawMnorm1k,"../Data_PRJNA434002/rawMnorm1k.rds")
write.table(rawMnorm1k,"../Data_PRJNA434002/rawMnorm1k.csv",sep=",")

rawMnorm1k10=as.matrix(rawMnorm1k[,seq(1,10455)*10])
rawMnorm1k10[1:10,1:10]
dim(rawMnorm1k10)
saveRDS(rawMnorm1k10,"../Data_PRJNA434002/rawMnorm1k10.rds")
write.table(rawMnorm1k10,"../Data_PRJNA434002/rawMnorm1k10.csv",sep=",")







sessionInfo()
q(save="no")
