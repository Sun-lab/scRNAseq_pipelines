
setwd("~/research/GitHub/scRNAseq_pipelines/MTG")

genes = readRDS("anno_marker_genes.rds")
dim(genes)
genes[1:2,]

table(genes$cell_type)

q("no")

###
