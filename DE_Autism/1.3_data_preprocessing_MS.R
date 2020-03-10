#this code convert the data from to .rds format and generate small subsets of them.

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#change the original file into rds for better read-writing-storage
exprM=read.table("../MS/exprMatrix.tsv.gz",header=TRUE,row.names=1)

saveRDS(exprM,"../MS/exprMatrix.rds")


meta[1:10,1:10]
exprM[1:10,1:10]



#generate another 10% percentage cells with the 5000 lowest 0 rate genes, for tSNE plots
exprM_zero_rate=apply(exprM==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(exprM_zero_rate)

exprM5k=exprM[o_zero_rate[1:5000],]
exprM5k[1:10,1:10]
dim(exprM5k)
saveRDS(exprM5k,"../MS/exprMatrix5k.rds")

exprM1k=exprM[o_zero_rate[1:1000],]
exprM1k[1:10,1:10]
dim(exprM1k)
saveRDS(exprM1k,"../MS/exprMatrix1k.rds")

#read in raw data directly


library("GenomicFeatures")
library("GenomicAlignments")
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)
library(scran)


sce   = read10xCounts("../MS/rawMatrix/", col.names=TRUE)	


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
saveRDS(rawMatrix,"../MS/rawM.rds")

rawM5k=as.matrix(rawMatrix[o_zero_rate[1:5000],])
rawM5k[1:10,1:10]
dim(rawM5k)
saveRDS(rawM5k,"../MS/rawM5k.rds")
write.table(rawM5k,"../MS/rawM5k.csv",sep=",")

rawM3k=as.matrix(rawMatrix[o_zero_rate[1:3000],])
rawM3k[1:10,1:10]
dim(rawM3k)
saveRDS(rawM3k,"../MS/rawM3k.rds")
write.table(rawM3k,"../MS/rawM3k.csv",sep=",")

rawM1k=as.matrix(rawMatrix[o_zero_rate[1:1000],])
rawM1k[1:10,1:10]
dim(rawM1k)
saveRDS(rawM1k,"../MS/rawM1k.rds")
write.table(rawM1k,"../MS/rawM1k.csv",sep=",")

####################generate the library adjusted(individual level read-depth adjusted) tables#############
rawM=rawMatrix

cur_individual=unique(meta$sample)
cell_num=matrix(ncol=1,nrow=length(cur_individual))
rownames(cell_num)=cur_individual
colnames(cell_num)="cell_num"
read_depth=matrix(ncol=1,nrow=length(cur_individual))
rownames(read_depth)=cur_individual
colnames(read_depth)="read_depth"
zero_rate_ind=matrix(nrow=nrow(rawM),ncol=length(cur_individual))
rownames(zero_rate_ind)=rownames(rawM)
colnames(zero_rate_ind)=cur_individual

for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=rawM[,meta$sample==cur_ind]
  cell_num[i_ind]=ncol(cur_ind_m)
  read_depth[i_ind]=sum(cur_ind_m,na.rm = TRUE)/cell_num[i_ind]*1000
  zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(x,na.rm = TRUE))})/cell_num[i_ind]
}

saveRDS(cell_num,"../MS/rawM_cell_num_per_ind.rds")
saveRDS(read_depth,"../MS/rawM_read_depth_per_1Kcell_ind.rds")
saveRDS(zero_rate_ind,"../MS/rawM_zero_rate_per_ind.rds")


####################
rawM5k=readRDS("../MS/rawM5k.rds")
meta=read.table("../MS/meta.tsv",header = TRUE, sep = "\t")
read_depth_5k=readRDS("../MS/rawM5k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_5k=median(read_depth_5k)/read_depth_5k

read_depth_ratio_5k=read_depth_ratio_5k[match(meta$sample,rownames(read_depth_5k))]

rawMrdpadj5k=t(apply(rawM5k,1,function(x){round(x*read_depth_ratio_5k)}))
rawMrdpadj5k[1:10,1:10]
dim(rawMrdpadj5k)
saveRDS(rawMrdpadj5k,"../MS/rawMrdpadj5k.rds")
write.table(rawMrdpadj5k,"../MS/rawMrdpadj5k.csv",sep=",")


rawM3k=readRDS("../MS/rawM3k.rds")
meta=read.table("../MS/meta.tsv",header = TRUE, sep = "\t")
read_depth_3k=readRDS("../MS/rawM3k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_3k=median(read_depth_3k)/read_depth_3k

read_depth_ratio_3k=read_depth_ratio_3k[match(meta$sample,rownames(read_depth_3k))]

rawMrdpadj3k=t(apply(rawM3k,1,function(x){round(x*read_depth_ratio_3k)}))
rawMrdpadj3k[1:10,1:10]
dim(rawMrdpadj3k)
saveRDS(rawMrdpadj3k,"../MS/rawMrdpadj3k.rds")
write.table(rawMrdpadj3k,"../MS/rawMrdpadj3k.csv",sep=",")

rawM1k=readRDS("../MS/rawM1k.rds")
meta=read.table("../MS/meta.tsv",header = TRUE, sep = "\t")
read_depth_1k=readRDS("../MS/rawM1k_read_depth_per_1Kcell_ind.rds")
read_depth_ratio_1k=median(read_depth_1k)/read_depth_1k

read_depth_ratio_1k=read_depth_ratio_1k[match(meta$sample,rownames(read_depth_1k))]
rawMrdpadj1k=t(apply(rawM1k,1,function(x){round(x*read_depth_ratio_1k)}))
rawMrdpadj1k[1:10,1:10]
dim(rawMrdpadj1k)
saveRDS(rawMrdpadj1k,"../MS/rawMrdpadj1k.rds")
write.table(rawMrdpadj1k,"../MS/rawMrdpadj1k.csv",sep=",")

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

saveRDS(sce,"../MS/rawMnorm.rds")
sce=readRDS("../MS/rawMnorm.rds")

#writeout
top5k=as.numeric(o_zero_rate[1:5000])
sce2=logcounts(sce)[top5k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm5k=as.matrix(sce2)
rawMnorm5k[1:10,1:10]
dim(rawMnorm5k)
saveRDS(rawMnorm5k,"../MS/rawMnorm5k.rds")
write.table(rawMnorm5k,"../MS/rawMnorm5k.csv",sep=",")

top3k=as.numeric(o_zero_rate[1:3000])
sce2=logcounts(sce)[top3k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm3k=as.matrix(sce2)
rawMnorm3k=as.matrix(rawMnorm3k[,seq(1,4891)*10])
rawMnorm3k[1:10,1:10]
dim(rawMnorm3k)
saveRDS(rawMnorm3k,"../MS/rawMnorm3k.rds")
write.table(rawMnorm3k,"../MS/rawMnorm3k.csv",sep=",")


top1k=as.numeric(o_zero_rate[1:1000])
sce2=logcounts(sce)[top1k,]
dim(sce2)
sce2[1:10,1:10]
rawMnorm1k=as.matrix(sce2)
rawMnorm1k[1:10,1:10]
dim(rawMnorm1k)
saveRDS(rawMnorm1k,"../MS/rawMnorm1k.rds")
write.table(rawMnorm1k,"../MS/rawMnorm1k.csv",sep=",")


sessionInfo()
q(save="no")
