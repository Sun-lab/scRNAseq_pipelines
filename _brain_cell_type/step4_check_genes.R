
# --------------------------------------------------------------
# check the differential expression analysis results from MTG
# and psychENCODE, and then use MTG results to extract genes 
# that have cell type-specific expression
# --------------------------------------------------------------

# --------------------------------------------------------------
# read in results
# --------------------------------------------------------------

genes = readRDS("../MTG/DE_gene_anno.rds")
dim(genes)
genes[1:2,]

genes.Psy = readRDS("../psychENCODE/deconvolution/DE_gene_anno.rds")
dim(genes.Psy)
genes.Psy[1:2,]

# cell types
cts = gsub("logFC.", "", colnames(genes)[grep("logFC", colnames(genes))])
cts

cts2 = gsub("logFC.", "", 
            colnames(genes.Psy)[grep("logFC", colnames(genes.Psy))])
cts2

# --------------------------------------------------------------
# take intersections
# --------------------------------------------------------------

g2use = intersect(genes$gene, genes.Psy$gene)

genes.MTG = genes[match(g2use, genes$gene),]
genes.Psy = genes.Psy[match(g2use, genes.Psy$gene),]

dim(genes.MTG)
dim(genes.Psy)

pdf("figure/compare_MTG_PsychENCODE.pdf", width=8, height=4)
par(mfrow=c(1,2), bty="n")

for(ct1 in cts){
  
  cat(ct1, "\n")
  
  fdr.col   = paste0("FDR_qvalue.", ct1)
  logfc.col = paste0("logFC.", ct1)
  
  sig.MTG = genes.MTG[[fdr.col]] < 0.001 & genes.MTG[[logfc.col]] > log(2)
  sig.Psy = genes.Psy[[fdr.col]] < 0.001 & genes.Psy[[logfc.col]] > log(2)
  
  print(table(sig.MTG, sig.Psy))
  
  cat("\n")
  
  qval.MTG = genes.MTG[[fdr.col]]
  qval.MTG[which(qval.MTG < 1e-10)] = 1e-10
  
  qval.Psy = genes.Psy[[fdr.col]]
  qval.Psy[which(qval.Psy < 1e-10)] = 1e-10
  
  logfc.MTG = genes.MTG[[logfc.col]]
  logfc.MTG = sign(logfc.MTG) * log10(abs(logfc.MTG))

  logfc.Psy = genes.Psy[[logfc.col]]
  logfc.Psy = sign(logfc.Psy) * log10(abs(logfc.Psy))
  
  smoothScatter(-log10(qval.MTG), -log10(qval.Psy), 
                xlab="MGT", ylab="PsychENCODE", 
                main=paste0(ct1, ": -log10(FDR)"))
  
  smoothScatter(logfc.MTG, logfc.Psy, xlab="MGT", ylab="PsychENCODE", 
                main="log10(log fold-change)")
  
}

dev.off()

# --------------------------------------------------------------
# add an indicator for differentially expressed genes
# --------------------------------------------------------------

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

# --------------------------------------------------------------
# plot the relation between prop. of cells expressed vs. 
# number of cell tyeps it is classified as cell type-specific
# --------------------------------------------------------------

library(ggplot2)
p <- ggplot(genes, aes(x=as.factor(n.cellTypes), y=prop_express)) + 
  geom_boxplot()
ggsave("figure/MTG_boxplot_prop_express_vs_n_cell_types.pdf", device="pdf", 
       width=4, height=3, units="in")

tapply(genes$prop_express, genes$n.cellTypes, quantile, 
       probs=c(0.05, 0.1, 0.2,  0.5, 0.8, 0.9, 0.95))

table(genes$prop_express > 0.01, genes$n.cellTypes)

# --------------------------------------------------------------
# write out cell type-specifically expressed genes
# --------------------------------------------------------------

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

# --------------------------------------------------------------
# compare with the cell type-specific genes reported by the 
# psychENCODE project, which has 25 cell types
# --------------------------------------------------------------

marker.psy = "DER-21_Single_cell_markergenes_UMI.txt"
markers = read.table(paste0("../psychENCODE/deconvolution/", marker.psy), 
                     sep="\t", as.is=TRUE, header=TRUE)
dim(markers)
markers[1:2,]
length(unique(markers$Cluster))
table(markers$Cluster)

length(unique(markers$Gene))
tb1 = table(table(markers$Gene))
tb2 = tb1[1:9]
tb2[9] = sum(tb1[9:length(tb1)])
tb2

# --------------------------------------------------------------
# for each gene, check its associaated clusters (cell types)
# similar cell types are often associated
# --------------------------------------------------------------

cluster.group = tapply(markers$Cluster, markers$Gene, paste, collapse=":")
length(cluster.group)
cluster.group[1:2]

tbc = table(cluster.group)
tbc = tbc[grep(":", names(tbc))]
sort(tbc, decreasing=T)[1:20]

# --------------------------------------------------------------
# for each cluster/cell type, check the overlap between 
# our marker genes and PsychENCODE marker genes
# --------------------------------------------------------------

table(genes.ct$gene %in% markers$Gene)
table(markers$Gene %in% genes.ct$gene)
tb3 = table(markers$Gene %in% genes.ct$gene, markers$Cluster)
tb3

pdf("figure/DER-21_Single_cell_markergenes.pdf", width=6, height=4)
par(las=2, cex=0.8)
barplot(tb3)
legend("topright", c("overlap", "no overlap"), 
       fill = c(gray(0.9), gray(0.2)))
dev.off()

sessionInfo()
q(save="no")

