
setwd("~/research/GitHub/scRNAseq_pipelines/MTG")

genes = readRDS("DE_gene_anno.rds")
dim(genes)
genes[1:2,]

cts = gsub("logFC.", "", colnames(genes)[grep("logFC", colnames(genes))])
cts

for(ct1 in cts){
  
  cname     = paste0("selected.", ct1)
  fdr.col   = paste0("FDR_qvalue.", ct1)
  logfc.col = paste0("logFC.", ct1)
  
  wwi = which(genes[[fdr.col]] < 0.001 & genes[[logfc.col]] > log(2))
  
  genes[[cname]] = rep(0, nrow=(genes))
  genes[[cname]][wwi] = 1
}

genes[["n.cellTypes"]] = rowSums(genes[,31:36])
table(genes$n.cellTypes)
colSums(genes[which(genes$n.cellTypes==1),31:36])

library(ggplot2)
p <- ggplot(genes, aes(x=as.factor(n.cellTypes), y=prop_express)) + 
  geom_boxplot()
ggsave("boxplot_prop_express_vs_n_cell_types.pdf", device="pdf", 
       width=4, height=3, units="in")

tapply(genes$prop_express, genes$n.cellTypes, quantile, 
       probs=c(0.05, 0.1, 0.2,  0.5, 0.8, 0.9, 0.95))

table(genes$prop_express > 0.01, genes$n.cellTypes)

col2use = c(1:6, 31:36)

ww1 = which(genes$prop_express > 0.01 & genes$n.cellTypes==1)
ww2 = which(genes$prop_express > 0.01 & genes$n.cellTypes==0)
genes.ct = genes[ww1,col2use]
genes.control = genes[ww2,col2use]

dim(genes.ct)
genes.ct[1:2,]

dim(genes.control)
genes.control[1:2,]

write.table(genes.ct, file="genes_ct_specific_expression.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

write.table(genes.control, file="genes_control.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

###
