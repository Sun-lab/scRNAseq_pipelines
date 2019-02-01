
repo_dir  = "~/research/GitHub/scRNAseq_pipelines"
work_dir  = file.path(repo_dir,"dronc")

dronc_dir = "~/research/scRNAseq/data/GTEx_droncseq_hip_pcf"
setwd(dronc_dir)

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

sce      = readRDS("sce.rds")
clusters = readRDS("all_clust_res.rds")

dim(sce)
dim(colData(sce))

dim(clusters)
clusters[1:2,]

t1 = table(clusters$part_cell_id, clusters$Cell_Type)
t2 = table(clusters$part_cell_id, clusters$Cluster.ID)
t1
t2

setwd(work_dir)

pdf("figure/part_cell_id_vs_cell_type.pdf", width=5, height=5)
heatmap(t1)
dev.off()

# ------------------------------------------------------------------
# We want to separate the clusters from PFC versus hippocampus. 
# The paper did not provide such label for each cell. 
# Based on supp. Fig 7d, and the part_cell_id, I identify the 
# following part_cell_id as cells from PFC. 
# ------------------------------------------------------------------

PFC = c("hCc", "hCd", "hCe", "hCf", "humanPFCa", "humanPFCb", "PFC-CD")

# ------------------------------------------------------------------
# Compare the results of kmeans 12 vs. 15 clusters, when 
# there are 12 clusters, they already capture cell type-specific
# infomration. 15 clusters further split clusters for a few 
# cell types. so we choose kmeans 12 clusters. 
# 
# For each cell type, find the largest cluster
# ------------------------------------------------------------------

t1 = table(clusters$KM_15, clusters$Cell_Type)
t1

t1 = table(clusters$KM_12, clusters$Cell_Type)
t1

clusters$cluster_kmean = clusters$KM_12
clusts = apply(t1, 2, function(v){union(which.max(v), which(v > 500))})
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

# ------------------------------------------------------------------
# save count data
# ------------------------------------------------------------------

setwd(work_dir)

geneInfo = as.data.frame(rowData(sce))
dim(geneInfo)
geneInfo[1:2,]

ct.matrx = list(all=ct.matrx, PFC=ct.matrx.PFC)

saveRDS(geneInfo, "gene_info_dronc.rds")
saveRDS(ct.matrx, "ct_matrix_dronc.rds")
saveRDS(nCells,   "ct_cells_dronc.rds")

sessionInfo()
q(save="no")


