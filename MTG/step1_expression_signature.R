
repo_dir  = "~/research/GitHub/scRNAseq_pipelines"
work_dir  = file.path(repo_dir,"MTG")

MTG_dir = "~/research/scRNAseq/data/Allen_BI/human_MTG_gene_expression_matrices_2018-06-14"
setwd(MTG_dir)

library('org.Hs.eg.db')

# ------------------------------------------------------------------
# here are the cell types 
# ------------------------------------------------------------------

# Middle Temporal Gyrus (MTG): This RNA-Seq data set is created from intact 
# nuclei derived from frozen human brain specimens, to survey cell type 
# diversity in the human middle temporal gyrus (MTG). In total, 15,928 nuclei 
# from 8 human tissue donors ranging in age from 24-66 years were analyzed. 
# Analysis of these transcriptional profiles reveals approximately 75 
# transcriptionally distinct cell types, subdivided into 45 inhibitory 
# neuron types, 24 excitatory neuron types, and 6 non-neuronal types.

# Exc: excitatory neurons, or glutamatergic neurons
# Inh: inhibitory neurons, or GABAergic inhibitory interneurons
# Astro: astrocytes
# Endo: endothelial cells
# Micro: microglia
# Oligo: oligodendrocytes
# OPC: oligodendrocyte precursor cells

# The following are the neuron types from dronc paper

# exPFC, glutamatergic neurons from the PFC; 
# GABA, GABAergic inhibitory interneurons; 
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

sce      = readRDS("final_sce.rds")
clusters = readRDS("final_hvg_clust.rds")

dim(sce)
dim(colData(sce))
colData(sce)[1:2,1:5]

table(colData(sce)$cell_type, colData(sce)$class)

dim(clusters)
clusters[1:2,1:5]
names(clusters)

table(clusters$sample_name == colData(sce)$sample_name)
table(clusters$cell_type == colData(sce)$cell_type)

t1 = table(clusters$KM_15, clusters$cell_type)
t1

# based on mannual examinaion of human_MTG.html, we choose to use the 
# clustering result of kmeans with 15 clusters. 
clusters$cluster_kmean = clusters$KM_15
clusts = apply(t1, 2, function(v){union(which.max(v), which(v > 200))})
clusts

# note that for some clusters, some cells belong to one cell type, 
# but other cells belong to another cell type. 
table(unlist(clusts))

# ------------------------------------------------------------------
# process geneInfo
# ------------------------------------------------------------------

geneInfo = as.data.frame(rowData(sce))
dim(geneInfo)
geneInfo[1:2,]
length(unique(geneInfo$gene))

columns(org.Hs.eg.db)
map1 = mapIds(org.Hs.eg.db, keys=as.character(geneInfo$entrez_id), 
              'ENSEMBL', 'ENTREZID')
length(map1)
map1[1:5]

geneInfo$ensembl_gene_id = as.character(map1)
table(names(map1) == geneInfo$entrez_id)

# ------------------------------------------------------------------
# collect counts for each cell type
# ------------------------------------------------------------------

celltypes = setdiff(unique(clusters$cell_type), "unknown")
celltypes

zeros  = rep(0,length(celltypes))
nCells = data.frame(Cell_Type=celltypes, nCells_All=zeros)

ct.matrx = matrix(NA, nrow=nrow(sce), ncol=length(celltypes))
colnames(ct.matrx) = celltypes
rownames(ct.matrx) = rowData(sce)$gene

for(ct1 in celltypes){
  ct.cond    = clusters$cell_type == ct1
  clust.cond = clusters$cluster_kmean %in% clusts[[ct1]]
  cells      = which(ct.cond & clust.cond)

  nCells[which(nCells$Cell_Type==ct1),2] = length(cells)

  ct.matrx[,ct1]      = rowSums(counts(sce)[,cells])
}

dim(ct.matrx)
ct.matrx[1:2,1:3]
summary(ct.matrx)

dim(nCells)
nCells

# ------------------------------------------------------------------
# save count data
# ------------------------------------------------------------------

setwd(work_dir)

saveRDS(geneInfo, "gene_info_human_MTG.rds")
saveRDS(ct.matrx, "ct_matrix_human_MTG.rds")
saveRDS(nCells,   "ct_cells_human_MTG.rds")

sessionInfo()
q(save="no")


