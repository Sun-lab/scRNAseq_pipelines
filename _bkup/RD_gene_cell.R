# Construct gene by cell expression matrix from scRNA clustering

# ----------
# Shortcuts
# ----------
rm(list=ls())
scRNA_dat = c("dronc","MTG")[1]
# scRNA_dat = c("dronc","MTG")[2]
if( scRNA_dat == "dronc" ){
	data_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/scRNAseq_pipelines/dronc"
} else if( scRNA_dat == "MTG" ){
	data_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/scRNAseq_pipelines/MTG"
} else if(FALSE){
	data_dir = "C:/Users/Admin/Desktop/blah"
}
setwd(data_dir)


# ----------
# Functions
# ----------
construct_express = function(data_dir,sce_obj,clust_df,sample_id,
	gene_name,cell_type_name,clust_name,max_count,my_cols){

	if(FALSE){
		sce_obj = sce; clust_df = clusters
		sample_id = "Cell.ID"; gene_name = "external_gene_name"
		cell_type_name = "Cell_Type"; clust_name = c("cluster_kmean","sc3_12_clusters")[1]
		
		sce_obj = sce_kmeans; clust_df = clust_kmeans
		sample_id = "Cell.ID"; gene_name = "gene"
		cell_type_name = "class"; clust_name = "cluster_kmean"
		
		max_count = 1000
		my_cols = c("white",gplots::colorpanel(10,"black","green"))
		
		sce_obj = sce_kmeans; clust_df = clust_all; sample_id = "Cell.ID"
		gene_name = "gene"; cell_type_name = "class"; clust_name = "cluster_MTG"; max_count = 50
		
	}

	pdf(file.path(data_dir,paste0(clust_name,".pdf")), width=5, height=5)
	
	print(dim(sce_obj))
	print(dim(colData(sce_obj[1:5,])))
	# print(colData(sce_obj[1:5,]))
	print(clust_df[1:5,])
	
	t1 = table(clust_df[,c(clust_name,cell_type_name)])
	print(t1)
	heatmap(t1,col=my_cols)
	heatmap(t1,col=my_cols,scale="none")
	heatmap(log10(1+t1),col=my_cols)
	heatmap(log10(1+t1),col=my_cols,scale="none")
	
	# Combine cluster results together to match/map to cell types
	clusts = apply(t1, 2, function(v){union(which(v > max_count), which.max(v))})
	print(clusts)
	
	# align cells of sce object and cells of cluster results
	# print(colData(sce_obj)[1:2,1:3])
	print(rownames(colData(sce_obj))[1:2])
	print(table(clust_df[,sample_id] == rownames(colData(sce_obj))))
	# do the two vectors contain the same set of elements
	setequal(clust_df[,sample_id], rownames(colData(sce_obj))) 

	# Re-arrange clusters rows to match sce column order
	mat1 = match(rownames(colData(sce_obj)), clust_df[,sample_id])
	clust_df = clust_df[mat1,]
	table(clust_df[,sample_id] == rownames(colData(sce_obj)))
	
	# collect counts for each cell type
	celltypes = na.omit(unique(clust_df[,cell_type_name]))
	print(celltypes)
	
	zeros = rep(0,length(celltypes))
	nCells = data.frame(Cell_Type=celltypes,nCells_All=zeros)

	ct.matrx = matrix(NA,nrow=nrow(sce_obj),ncol=length(celltypes))
	colnames(ct.matrx) = celltypes
	rownames(ct.matrx) = rowData(sce_obj)[,gene_name]

	cat("Constructing matrix...\n")
	for(ct1 in celltypes){
		# ct1 = celltypes[1]
		
		# for each cell type ...
		ct.cond = clust_df[,cell_type_name] == ct1
		
		# get the cluster number associated with it
		clust.cond = clust_df[,clust_name] %in% clusts[[ct1]]

		cells = which(ct.cond & clust.cond)

		# Track number of cells for each cell type
		nCells[which(nCells$Cell_Type==ct1),2] = length(cells)

		# producing cell type profile matrix
		ct.matrx[,ct1] = rowSums(counts(sce_obj)[,cells])
	}

	print(dim(ct.matrx))
	print(ct.matrx[1:2,1:3])
	
	print(summary(ct.matrx))

	print(dim(nCells))
	print(nCells)

	ct.matrx = list(all=ct.matrx)
	cat("Save images ...\n")
	saveRDS(ct.matrx,file.path(data_dir,paste0("ct_matrix_",clust_name,".rds")))
	saveRDS(nCells,file.path(data_dir,paste0("ct_cells_",clust_name,".rds")))
	
	cat("Plotting heatmaps ...\n")
	heatmap(ct.matrx$all,col=my_cols,Rowv=NA)
	heatmap(ct.matrx$all,col=my_cols,Rowv=NA,scale="none")
	heatmap(log10(1+ct.matrx$all),col=my_cols,Rowv=NA)
	heatmap(log10(1+ct.matrx$all),col=my_cols,Rowv=NA,scale="none")
	
	dev.off()
}
my_cols = c("white",gplots::colorpanel(10,"black","green"))


# ----------
# Construct expression matrices
# ----------
if( scRNA_dat == "dronc" ){

# Read and Prep Data
sce = readRDS(file.path(data_dir,"sce.rds"))
clusters = readRDS(file.path(data_dir,"all_clust_res.rds"))

construct_express(data_dir = data_dir,
	sce_obj = sce,
	clust_df = clusters,
	sample_id = "Cell.ID",
	gene_name = "external_gene_name",
	cell_type_name = "Cell_Type",
	clust_name = "cluster_kmean",
	max_count = 1000,
	my_cols = my_cols)

construct_express(data_dir = data_dir,
	sce_obj = sce,
	clust_df = clusters,
	sample_id = "Cell.ID",
	gene_name = "external_gene_name",
	cell_type_name = "Cell_Type",
	clust_name = "sc3_12_clusters",
	max_count = 1000,
	my_cols = my_cols)

construct_express(data_dir = data_dir,
	sce_obj = sce,
	clust_df = clusters,
	sample_id = "Cell.ID",
	gene_name = "external_gene_name",
	cell_type_name = "Cell_Type",
	clust_name = "sc3_10_clusters",
	max_count = 1000,
	my_cols = my_cols)

# PFC = c("hCc", "hCd", "hCe", "hCf", "humanPFCa", "humanPFCb", "PFC-CD")

}

if( scRNA_dat == "MTG" ){

# Read and Prep Data
rds = readRDS(file.path(data_dir,"post_cluster.rds"))
sce_kmeans = rds$sce_sub
clust_kmeans = data.frame(Cell.ID = colnames(sce_kmeans),
	rds$df_tsne,stringsAsFactors = FALSE)
clust_kmeans$cluster_kmean = as.numeric(clust_kmeans$cluster_kmean)
clust_kmeans$class = colData(sce_kmeans)$class
clust_kmeans$cluster_MTG = colData(sce_kmeans)$cluster
rm(rds)

# Read and Prep Data
sce_sc3 = readRDS(file.path(data_dir,"post_sc3.rds"))$sce
clust_sc3 = data.frame(colData(sce_sc3)[,c("sample_name","class",
	"cluster","sc3_10_clusters","sc3_15_clusters")],
	stringsAsFactors = FALSE)
rownames(clust_sc3) = NULL
names(clust_sc3)[names(clust_sc3) == "sample_name"] = "Cell.ID"
names(clust_sc3)[names(clust_sc3) == "cluster"] = "cluster_MTG"

# Merge kmeans and sc3 clustering results
clust_all = merge(clust_kmeans,clust_sc3,
	by = intersect(names(clust_kmeans),names(clust_sc3)))
tmp_df = data.frame(cluster_MTG = sort(unique(clust_all$cluster_MTG)),
	stringsAsFactors = FALSE)
tmp_df$cluster_MTG_2 = seq(nrow(tmp_df))
clust_all = merge(clust_all,tmp_df,by=c("cluster_MTG"))

construct_express(data_dir = data_dir,
	sce_obj = sce_kmeans,
	clust_df = clust_all,
	sample_id = "Cell.ID",
	gene_name = "gene",
	cell_type_name = "class",
	clust_name = "cluster_MTG_2",
	max_count = 50,
	my_cols = my_cols)

construct_express(data_dir = data_dir,
	sce_obj = sce_kmeans,
	clust_df = clust_all,
	sample_id = "Cell.ID",
	gene_name = "gene",
	cell_type_name = "class",
	clust_name = "cluster_kmean",
	max_count = 500,
	my_cols = my_cols)

construct_express(data_dir = data_dir,
	sce_obj = sce_sc3,
	clust_df = clust_all,
	sample_id = "Cell.ID",
	gene_name = "gene",
	cell_type_name = "class",
	clust_name = "sc3_10_clusters",
	max_count = 550,
	my_cols = my_cols)

construct_express(data_dir = data_dir,
	sce_obj = sce_sc3,
	clust_df = clust_all,
	sample_id = "Cell.ID",
	gene_name = "gene",
	cell_type_name = "class",
	clust_name = "sc3_15_clusters",
	max_count = 500,
	my_cols = my_cols)

}



###


