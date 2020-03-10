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
saveRDS(meta10,"../Data_PRJNA434002/meta10.rds")
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
zero_rate=apply(counts(sce)==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(zero_rate)

rawMatrix=counts(sce) #matrix to large, so we are not able to convert it as regular matrix
saveRDS(rawMatrix,"../Data_PRJNA434002/rawM.rds")

top1k=as.numeric(o_zero_rate[1:1000])
rawM1k10=as.matrix(rawMatrix[top1k,seq(1,10455)*10])
saveRDS(rawM1k10,"../Data_PRJNA434002/rawM1k10.rds")
write.table(rawM1k10,"../Data_PRJNA434002/rawM1k10.csv",sep=",")

top3k=as.numeric(o_zero_rate[1:3000])
rawM3k10=as.matrix(rawMatrix[top3k,seq(1,10455)*10])
saveRDS(rawM3k10,"../Data_PRJNA434002/rawM3k10.rds")
write.table(rawM3k10,"../Data_PRJNA434002/rawM3k10.csv",sep=",")

rawM5k=as.matrix(rawMatrix[o_zero_rate[1:5000],])
rawM5k[1:10,1:10]
dim(rawM5k)
saveRDS(rawM5k,"../Data_PRJNA434002/rawM5k.rds")
write.table(rawM5k,"../Data_PRJNA434002/rawM5k.csv",sep=",")

rawM3k=as.matrix(rawMatrix[o_zero_rate[1:3000],])
rawM3k[1:10,1:10]
dim(rawM3k)
saveRDS(rawM3k,"../Data_PRJNA434002/rawM3k.rds")
write.table(rawM3k,"../Data_PRJNA434002/rawM3k.csv",sep=",")

rawM2k=as.matrix(rawMatrix[o_zero_rate[1:2000],])
rawM2k[1:10,1:10]
dim(rawM2k)
saveRDS(rawM2k,"../Data_PRJNA434002/rawM2k.rds")
write.table(rawM2k,"../Data_PRJNA434002/rawM2k.csv",sep=",")

rawM1k=as.matrix(rawMatrix[o_zero_rate[1:1000],])
rawM1k[1:10,1:10]
dim(rawM1k)
saveRDS(rawM1k,"../Data_PRJNA434002/rawM1k.rds")
write.table(rawM1k,"../Data_PRJNA434002/rawM1k.csv",sep=",")

# %50 version
set.seed(502)
random50=sample.int(104559,52280)
meta2=meta[random50,]
saveRDS(meta2,"../Data_PRJNA434002/meta2.rds")

rawM3k50=as.matrix(rawMatrix[o_zero_rate[1:3000],random50])
dim(rawM3k50)
saveRDS(rawM3k50,"../Data_PRJNA434002/rawM3k50.rds")
write.table(rawM3k50,"../Data_PRJNA434002/rawM3k50.csv",sep=",")

rawM2k50=as.matrix(rawMatrix[o_zero_rate[1:3000],random50])
dim(rawM2k50)
saveRDS(rawM2k50,"../Data_PRJNA434002/rawM2k50.rds")
write.table(rawM2k50,"../Data_PRJNA434002/rawM2k50.csv",sep=",")

rawM1k50=as.matrix(rawMatrix[o_zero_rate[1:3000],random50])
dim(rawM1k50)
saveRDS(rawM1k50,"../Data_PRJNA434002/rawM1k50.rds")
write.table(rawM1k50,"../Data_PRJNA434002/rawM1k50.csv",sep=",")

####################generate the library adjusted(individual level read-depth adjusted) tables#############

#this part using 

rawM3k10=readRDS("../Data_PRJNA434002/rawM3k10.rds")
read_depth_3k10=readRDS("../Data_PRJNA434002/rawM3k10_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_3k10=median(read_depth_3k10)/read_depth_3k10
meta10=readRDS("../Data_PRJNA434002/meta10.rds")
read_depth_ratio_3k10=read_depth_ratio_3k10[match(meta10$individual,rownames(read_depth_3k10))]

rawMrdpadj3k10=t(apply(rawM3k10,1,function(x){round(x*read_depth_ratio_3k10)}))
rawMrdpadj3k10[1:10,1:10]
dim(rawMrdpadj3k10)
saveRDS(rawMrdpadj3k10,"../Data_PRJNA434002/rawMrdpadj3k10.rds")
write.table(rawMrdpadj3k10,"../Data_PRJNA434002/rawMrdpadj3k10.csv",sep=",")



rawM5k=readRDS("../Data_PRJNA434002/rawM5k.rds")
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
read_depth_5k=readRDS("../Data_PRJNA434002/rawM5k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_5k=median(read_depth_5k)/read_depth_5k

read_depth_ratio_5k=read_depth_ratio_5k[match(meta$individual,rownames(read_depth_5k))]

rawMrdpadj5k=t(apply(rawM5k,1,function(x){round(x*read_depth_ratio_5k)}))
rawMrdpadj5k[1:10,1:10]
dim(rawMrdpadj5k)
saveRDS(rawMrdpadj5k,"../Data_PRJNA434002/rawMrdpadj5k.rds")
write.table(rawMrdpadj5k,"../Data_PRJNA434002/rawMrdpadj5k.csv",sep=",")


rawM3k=readRDS("../Data_PRJNA434002/rawM3k.rds")
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
read_depth_3k=readRDS("../Data_PRJNA434002/rawM3k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_3k=median(read_depth_3k)/read_depth_3k

read_depth_ratio_3k=read_depth_ratio_3k[match(meta$individual,rownames(read_depth_3k))]

rawMrdpadj3k=t(apply(rawM3k,1,function(x){round(x*read_depth_ratio_3k)}))
rawMrdpadj3k[1:10,1:10]
dim(rawMrdpadj3k)
saveRDS(rawMrdpadj3k,"../Data_PRJNA434002/rawMrdpadj3k.rds")
write.table(rawMrdpadj3k,"../Data_PRJNA434002/rawMrdpadj3k.csv",sep=",")


rawM2k=readRDS("../Data_PRJNA434002/rawM2k.rds")
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
read_depth_2k=readRDS("../Data_PRJNA434002/rawM2k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_2k=median(read_depth_2k)/read_depth_2k

read_depth_ratio_2k=read_depth_ratio_2k[match(meta$individual,rownames(read_depth_2k))]

rawMrdpadj2k=t(apply(rawM2k,1,function(x){round(x*read_depth_ratio_2k)}))
rawMrdpadj2k[1:10,1:10]
dim(rawMrdpadj2k)
saveRDS(rawMrdpadj2k,"../Data_PRJNA434002/rawMrdpadj2k.rds")
write.table(rawMrdpadj2k,"../Data_PRJNA434002/rawMrdpadj2k.csv",sep=",")
rawM1k=readRDS("../Data_PRJNA434002/rawM1k.rds")
meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
read_depth_1k=readRDS("../Data_PRJNA434002/rawM1k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_1k=median(read_depth_1k)/read_depth_1k

read_depth_ratio_1k=read_depth_ratio_1k[match(meta$individual,rownames(read_depth_1k))]
rawMrdpadj1k=t(apply(rawM1k,1,function(x){round(x*read_depth_ratio_1k)}))
rawMrdpadj1k[1:10,1:10]
dim(rawMrdpadj1k)
saveRDS(rawMrdpadj1k,"../Data_PRJNA434002/rawMrdpadj1k.rds")
write.table(rawMrdpadj1k,"../Data_PRJNA434002/rawMrdpadj1k.csv",sep=",")

#%50 version
rawM3k50=readRDS("../Data_PRJNA434002/rawM3k50.rds")
meta2=readRDS("../Data_PRJNA434002/meta2.rds")
read_depth_3k50=readRDS("../Data_PRJNA434002/rawM3k50_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_3k50=median(read_depth_3k50)/read_depth_3k50

read_depth_ratio_3k50=read_depth_ratio_3k50[match(meta2$individual,rownames(read_depth_3k50))]

rawMrdpadj3k50=t(apply(rawM3k50,1,function(x){round(x*read_depth_ratio_3k50)}))
rawMrdpadj3k50[1:10,1:10]
dim(rawMrdpadj3k50)
saveRDS(rawMrdpadj3k50,"../Data_PRJNA434002/rawMrdpadj3k50.rds")
write.table(rawMrdpadj3k50,"../Data_PRJNA434002/rawMrdpadj3k50.csv",sep=",")


rawM2k50=readRDS("../Data_PRJNA434002/rawM2k50.rds")
meta2=readRDS("../Data_PRJNA434002/meta2.rds")
read_depth_2k50=readRDS("../Data_PRJNA434002/rawM2k50_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_2k50=median(read_depth_2k50)/read_depth_2k50

read_depth_ratio_2k50=read_depth_ratio_2k50[match(meta2$individual,rownames(read_depth_2k50))]

rawMrdpadj2k50=t(apply(rawM2k50,1,function(x){round(x*read_depth_ratio_2k50)}))
rawMrdpadj2k50[1:10,1:10]
dim(rawMrdpadj2k50)
saveRDS(rawMrdpadj2k50,"../Data_PRJNA434002/rawMrdpadj2k50.rds")
write.table(rawMrdpadj2k50,"../Data_PRJNA434002/rawMrdpadj2k50.csv",sep=",")

rawM1k50=readRDS("../Data_PRJNA434002/rawM1k50.rds")
meta2=readRDS("../Data_PRJNA434002/meta2.rds")
read_depth_1k50=readRDS("../Data_PRJNA434002/rawM1k50_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_1k50=median(read_depth_1k50)/read_depth_1k50

read_depth_ratio_1k50=read_depth_ratio_1k50[match(meta2$individual,rownames(read_depth_1k50))]
rawMrdpadj1k50=t(apply(rawM1k50,1,function(x){round(x*read_depth_ratio_1k50)}))
rawMrdpadj1k50[1:10,1:10]
dim(rawMrdpadj1k50)
saveRDS(rawMrdpadj1k50,"../Data_PRJNA434002/rawMrdpadj1k50.rds")
write.table(rawMrdpadj1k50,"../Data_PRJNA434002/rawMrdpadj1k50.csv",sep=",")
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

saveRDS(sce,"../Data_PRJNA434002/rawMnorm.rds")
sce=readRDS("../Data_PRJNA434002/rawMnorm.rds")

#writeout
top5k=as.numeric(o_zero_rate[1:5000])
sce2=logcounts(sce)[top5k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm5k=as.matrix(sce2)
rawMnorm5k[1:10,1:10]
dim(rawMnorm5k)
saveRDS(rawMnorm5k,"../Data_PRJNA434002/rawMnorm5k.rds")
write.table(rawMnorm5k,"../Data_PRJNA434002/rawMnorm5k.csv",sep=",")


top3k=as.numeric(o_zero_rate[1:3000])
sce2=logcounts(sce)[top3k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm3k10=as.matrix(sce2)
rawMnorm3k10=as.matrix(rawMnorm3k10[,seq(1,10455)*10])
rawMnorm3k10[1:10,1:10]
dim(rawMnorm3k10)
saveRDS(rawMnorm3k10,"../Data_PRJNA434002/rawMnorm3k10.rds")
write.table(rawMnorm3k10,"../Data_PRJNA434002/rawMnorm3k10.csv",sep=",")




top1k=as.numeric(o_zero_rate[1:1000])
sce2=logcounts(sce)[top1k,]
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













###############Get PFC info#############################

meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
PFC_index=which(meta$region=="PFC")
meta10=readRDS("../Data_PRJNA434002/meta10.rds")
PFC10_index=which(meta10$region=="PFC")
meta50=readRDS("../Data_PRJNA434002/meta50.rds")
PFC50_index=which(meta50$region=="PFC")


meta_PFC=meta[PFC_index,]
meta10_PFC=meta10[PFC10_index,]
meta50_PFC=meta50[PFC50_index,]

saveRDS(meta_PFC,"../Data_PRJNA434002/meta_PFC.rds")
saveRDS(meta10_PFC,"../Data_PRJNA434002/meta10_PFC.rds")
saveRDS(meta50_PFC,"../Data_PRJNA434002/meta50_PFC.rds")

rawM1k10=readRDS("../Data_PRJNA434002/rawM1k10.rds")
rawMPFC1k10=rawM1k10[,PFC10_index]
saveRDS(rawMPFC1k10,"../Data_PRJNA434002/rawMPFC1k10.rds")
write.table(rawMPFC1k10,"../Data_PRJNA434002/rawMPFC1k10.csv",sep=",")

rawM3k10=readRDS("../Data_PRJNA434002/rawM3k10.rds")
rawMPFC3k10=rawM3k10[,PFC10_index]
saveRDS(rawMPFC3k10,"../Data_PRJNA434002/rawMPFC3k10.rds")
write.table(rawMPFC3k10,"../Data_PRJNA434002/rawMPFC3k10.csv",sep=",")

rawM1k=readRDS("../Data_PRJNA434002/rawM1k.rds")
rawMPFC1k=rawM1k[,PFC_index]
saveRDS(rawMPFC1k,"../Data_PRJNA434002/rawMPFC1k.rds")
write.table(rawMPFC1k,"../Data_PRJNA434002/rawMPFC1k.csv",sep=",")

rawM2k=readRDS("../Data_PRJNA434002/rawM2k.rds")
rawMPFC2k=rawM2k[,PFC_index]
saveRDS(rawMPFC2k,"../Data_PRJNA434002/rawMPFC2k.rds")
write.table(rawMPFC2k,"../Data_PRJNA434002/rawMPFC2k.csv",sep=",")

rawM3k=readRDS("../Data_PRJNA434002/rawM3k.rds")
rawMPFC3k=rawM3k[,PFC_index]
saveRDS(rawMPFC3k,"../Data_PRJNA434002/rawMPFC3k.rds")
write.table(rawMPFC3k,"../Data_PRJNA434002/rawMPFC3k.csv",sep=",")

rawM5k=readRDS("../Data_PRJNA434002/rawM5k.rds")
rawMPFC5k=rawM5k[,PFC_index]
saveRDS(rawMPFC5k,"../Data_PRJNA434002/rawMPFC5k.rds")
write.table(rawMPFC5k,"../Data_PRJNA434002/rawMPFC5k.csv",sep=",")



rawMnorm1k10=readRDS("../Data_PRJNA434002/rawMnorm1k10.rds")
rawMPFCnorm1k10=rawMnorm1k10[,PFC10_index]
saveRDS(rawMPFCnorm1k10,"../Data_PRJNA434002/rawMPFCnorm1k10.rds")
write.table(rawMPFCnorm1k10,"../Data_PRJNA434002/rawMPFCnorm1k10.csv",sep=",")

rawMnorm3k10=readRDS("../Data_PRJNA434002/rawMnorm3k10.rds")
rawMPFCnorm3k10=rawMnorm3k10[,PFC10_index]
saveRDS(rawMPFCnorm3k10,"../Data_PRJNA434002/rawMPFCnorm3k10.rds")
write.table(rawMPFCnorm3k10,"../Data_PRJNA434002/rawMPFCnorm3k10.csv",sep=",")

rawMnorm1k=readRDS("../Data_PRJNA434002/rawMnorm1k.rds")
rawMPFCnorm1k=rawMnorm1k[,PFC_index]
saveRDS(rawMPFCnorm1k,"../Data_PRJNA434002/rawMPFCnorm1k.rds")
write.table(rawMPFCnorm1k,"../Data_PRJNA434002/rawMPFCnorm1k.csv",sep=",")

rawMnorm2k=readRDS("../Data_PRJNA434002/rawMnorm2k.rds")
rawMPFCnorm2k=rawMnorm2k[,PFC_index]
saveRDS(rawMPFCnorm2k,"../Data_PRJNA434002/rawMPFCnorm2k.rds")
write.table(rawMPFCnorm2k,"../Data_PRJNA434002/rawMPFCnorm2k.csv",sep=",")

rawMnorm3k=readRDS("../Data_PRJNA434002/rawMnorm3k.rds")
rawMPFCnorm3k=rawMnorm3k[,PFC_index]
saveRDS(rawMPFCnorm3k,"../Data_PRJNA434002/rawMPFCnorm3k.rds")
write.table(rawMPFCnorm3k,"../Data_PRJNA434002/rawMPFCnorm3k.csv",sep=",")

rawMnorm5k=readRDS("../Data_PRJNA434002/rawMnorm5k.rds")
rawMPFCnorm5k=rawMnorm5k[,PFC_index]
saveRDS(rawMPFCnorm5k,"../Data_PRJNA434002/rawMPFCnorm5k.rds")
write.table(rawMPFCnorm5k,"../Data_PRJNA434002/rawMPFCnorm5k.csv",sep=",")





rawMrdpadj1k10=readRDS("../Data_PRJNA434002/rawMrdpadj1k10.rds")
rawMPFCrdpadj1k10=rawMrdpadj1k10[,PFC10_index]
saveRDS(rawMPFCrdpadj1k10,"../Data_PRJNA434002/rawMPFCrdpadj1k10.rds")
write.table(rawMPFCrdpadj1k10,"../Data_PRJNA434002/rawMPFCrdpadj1k10.csv",sep=",")

rawMrdpadj3k10=readRDS("../Data_PRJNA434002/rawMrdpadj3k10.rds")
rawMPFCrdpadj3k10=rawMrdpadj3k10[,PFC10_index]
saveRDS(rawMPFCrdpadj3k10,"../Data_PRJNA434002/rawMPFCrdpadj3k10.rds")
write.table(rawMPFCrdpadj3k10,"../Data_PRJNA434002/rawMPFCrdpadj3k10.csv",sep=",")

rawMrdpadj1k=readRDS("../Data_PRJNA434002/rawMrdpadj1k.rds")
rawMPFCrdpadj1k=rawMrdpadj1k[,PFC_index]
saveRDS(rawMPFCrdpadj1k,"../Data_PRJNA434002/rawMPFCrdpadj1k.rds")
write.table(rawMPFCrdpadj1k,"../Data_PRJNA434002/rawMPFCrdpadj1k.csv",sep=",")

rawMrdpadj2k=readRDS("../Data_PRJNA434002/rawMrdpadj2k.rds")
rawMPFCrdpadj2k=rawMrdpadj2k[,PFC_index]
saveRDS(rawMPFCrdpadj2k,"../Data_PRJNA434002/rawMPFCrdpadj2k.rds")
write.table(rawMPFCrdpadj2k,"../Data_PRJNA434002/rawMPFCrdpadj2k.csv",sep=",")

rawMrdpadj3k=readRDS("../Data_PRJNA434002/rawMrdpadj3k.rds")
rawMPFCrdpadj3k=rawMrdpadj3k[,PFC_index]
saveRDS(rawMPFCrdpadj3k,"../Data_PRJNA434002/rawMPFCrdpadj3k.rds")
write.table(rawMPFCrdpadj3k,"../Data_PRJNA434002/rawMPFCrdpadj3k.csv",sep=",")

rawMrdpadj5k=readRDS("../Data_PRJNA434002/rawMrdpadj5k.rds")
rawMPFCrdpadj5k=rawMrdpadj5k[,PFC_index]
saveRDS(rawMPFCrdpadj5k,"../Data_PRJNA434002/rawMPFCrdpadj5k.rds")
write.table(rawMPFCrdpadj5k,"../Data_PRJNA434002/rawMPFCrdpadj5k.csv",sep=",")












sessionInfo()
q(save="no")
