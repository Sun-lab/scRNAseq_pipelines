#this code is an add-on for the code 3_autism-tSNE_plot 
#by generating the file with rawCount data with the normalization



library("GenomicFeatures")
library("GenomicAlignments")
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)
library(scran)


setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

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

rawM_zero_rate=apply(counts(sce)==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(rawM_zero_rate)
top5k=as.numeric(o_zero_rate[1:5000])
sce=sce[top5k,]

table(gene.annotation$ensembl_gene_id == rowData(sce)$ID)	
rowData(sce)  = gene.annotation	
rownames(sce) = scater::uniquifyFeatureNames(rowData(sce)$ensembl_gene_id, 	
                                             rowData(sce)$hgnc_symbol)	

#gene normalization
date()
clusters = quickCluster(sce, min.mean=0.1, method="igraph")
date()
sce      = computeSumFactors(sce, cluster=clusters, min.mean=0.1)
date()
summary(sizeFactors(sce))
sce = normalize(sce)


rawMnorm5k=as.matrix(counts(sce))
rawMnorm5k[1:10,1:10]
dim(rawMnorm5k)

dim(rawMnorm5k)
saveRDS(rawMnorm5k,"../Data_PRJNA434002/rawMnorm5k.rds")
write.table(rawMnorm5k,"../Data_PRJNA434002/rawMnorm5k.csv",sep=",")

