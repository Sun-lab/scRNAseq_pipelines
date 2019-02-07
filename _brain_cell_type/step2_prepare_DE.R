rm(list=ls())

repo_dir	= "~/research/GitHub/scRNAseq_pipelines"
work_dir	= file.path(repo_dir,"MTG")
MTG_dir		= "~/research/scRNAseq/data/Allen_BI/human_MTG_gene_expression_matrices_2018-06-14"

MTG_dir 	= "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/scRNAseq_pipelines/MTG"

# ------------------------------------------------------------------
# load data matrix 
# ------------------------------------------------------------------

setwd(MTG_dir)

sce			= readRDS("final_sce.rds")
clusters	= readRDS("final_hvg_clust.rds")

dim(sce)
dim(colData(sce))
colData(sce)[1:2,1:5]

names(colData(sce))
table(colData(sce)$donor)

dim(rowData(sce))
rowData(sce)[1:2,]

# ------------------------------------------------------------------
# select the subset of cells that are clusted with the cells 
# of the same type
# ------------------------------------------------------------------

dim(clusters)
clusters[1:2,1:5]
table(clusters$sample_name == colData(sce)$sample_name)
table(clusters$cell_type   == colData(sce)$cell_type)

# opt_clust = paste0("KM_",15)
opt_clust = paste0("KM_",18) # our clustering results differ slightly
t1 = table(clusters[,opt_clust], clusters$cell_type)
t1

clusts = apply(t1, 2, function(v){union(which.max(v), which(v > 200))})
clusts

celltypes = setdiff(unique(clusters$cell_type), "unknown")
celltypes

w2kp = NULL

for(ct1 in celltypes){
	ct.cond    = clusters$cell_type == ct1
	clust.cond = clusters[,opt_clust] %in% clusts[[ct1]]
	w2kp = c(w2kp, which(ct.cond & clust.cond))
}
length(w2kp)

dim(sce)
sce = sce[,w2kp]
dim(sce)

saveRDS(sce, file.path(MTG_dir, "final_sce_filtered_by_kmeans.rds"))

sessionInfo()
q(save="no")


