
setwd("~/research/GitHub/scRNAseq_pipelines/dronc/")

# ------------------------------------------------------------------
# here are the cell types 
# ------------------------------------------------------------------

# exPFC, glutamatergic neurons from the PFC; 
# GABA, GABAergic interneurons; 
# exCA1/3, pyramidal neurons from the hip CA region; 
# exDG, granule neurons from the hip dentate gyrus region; 
# ASC, astrocytes; 
# MG, microglia; 
# ODC, oligodendrocytes; 
# OPC, oligodendrocyte precursor cells; 
# NSC, neuronal stem cells; 
# SMC, smooth muscle cells; 
# END, endothelial cells

# ------------------------------------------------------------------
# load clustering results 
# ------------------------------------------------------------------

sce = readRDS("sce.rds")
clusters = readRDS("all_clust_res.rds")

dim(sce)
dim(colData(sce))

dim(clusters)
clusters[1:2,]

t1 = table(clusters$part_cell_id, clusters$Cell_Type)
t2 = table(clusters$part_cell_id, clusters$Cluster.ID)
t1
t2

pdf("part_cell_id_vs_cell_type.pdf", width=5, height=5)
heatmap(t1)
dev.off()

# ------------------------------------------------------------------
# I cannot see lable of samples for each cell, based on supp. Fig 7
# and the part_cell_id, I assume those part_cell_id that have 
# large number of cells from exPFC and GABA. 
# ------------------------------------------------------------------

PFC = c("hCc", "hCd", "hCe", "hCf", "humanPFCa", "humanPFCb", "PFC-CD")

# ------------------------------------------------------------------
# for each cell type, find the largest cluster, except for ODC, use  
# the top 2 clusters
# ------------------------------------------------------------------

t1 = table(clusters$cluster_kmean, clusters$Cell_Type)
t1
clusts = apply(t1, 2, function(v){union(which(v > 1000), which.max(v))})
clusts

# ------------------------------------------------------------------
# align cells of sce object and cells of cluster results
# ------------------------------------------------------------------

colData(sce)[1:2,1:3]
rownames(colData(sce))[1:2]
table(clusters$Cell.ID == rownames(colData(sce)))
setequal(clusters$Cell.ID, rownames(colData(sce)))

mat1 = match(rownames(colData(sce)), clusters$Cell.ID)
clusters = clusters[mat1,]
table(clusters$Cell.ID == rownames(colData(sce)))

# ------------------------------------------------------------------
# collect counts for each cell type
# ------------------------------------------------------------------

celltypes = na.omit(unique(clusters$Cell_Type))
celltypes

zeros  = rep(0,length(celltypes))
nCells = data.frame(Cell_Type=celltypes, nCells_All=zeros, nCells_PFC=zeros)

ct.matrx = ct.matrx.PFC = matrix(NA, nrow=nrow(sce), ncol=length(celltypes))
colnames(ct.matrx) = colnames(ct.matrx.PFC) = celltypes
rownames(ct.matrx) = rownames(ct.matrx.PFC) = rowData(sce)$external_gene_name

geneInfo = rowData(sce)

for(ct1 in celltypes){
  ct.cond    = clusters$Cell_Type == ct1
  clust.cond = clusters$cluster_kmean %in% clusts[[ct1]]
  samp.cond  = clusters$part_cell_id %in% PFC
  
  cells      = which(ct.cond & clust.cond)
  cells.PFC  = which(ct.cond & clust.cond & samp.cond)
  
  nCells[which(nCells$Cell_Type==ct1),2] = length(cells)
  nCells[which(nCells$Cell_Type==ct1),3] = length(cells.PFC)
  
  ct.matrx[,ct1]  = rowSums(counts(sce)[,cells])
  ct.matrx.PFC[,ct1]  = rowSums(counts(sce)[,cells.PFC])
}

dim(ct.matrx)
ct.matrx[1:2,1:3]

dim(ct.matrx.PFC)
ct.matrx.PFC[1:2,1:3]

summary(ct.matrx)
summary(ct.matrx.PFC)

dim(nCells)
nCells

ct.matrx = list(all=ct.matrx, PFC=ct.matrx.PFC)
saveRDS(ct.matrx, "ct_matrix.rds")
saveRDS(nCells, "ct_cells.rds")

q(save="no")


