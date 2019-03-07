# scRNAseq differential expression (seems to require at least 80 gbs)

# ----------
# Shortcuts
# ----------
rm(list=ls())
repo_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/scRNAseq_pipelines"
MTG_dir = file.path(repo_dir,"MTG")
setwd(MTG_dir)


# ----------
# Libraries/Functions
# ----------
source(file.path(repo_dir,"SOURCE.R"))
if( !("BiocManager" %in% installed.packages()[,"Package"]) ){
	install.packages("BiocManager",repos = "https://mirrors.nics.utk.edu/cran/")
}
if( !("MAST" %in% installed.packages()[,"Package"]) ){
	# MAST = Model-based Analysis of Single Cell Transcriptomics
	BiocManager::install("MAST")
}
library(ggplot2)
library(data.table)
MAST_DEgenes = function(work_dir,num_genes=NULL,sce_obj,one_cell_type,fdr_thres=1e-3,logFC_thres=log(2)){
	
	if(FALSE){
		work_dir = MTG_dir
		num_genes = 100
		sce_obj = sce
		one_cell_type = c("Astro","Exc","Inh","Micro","Oligo","OPC")[1]
		fdr_thres = 1e-3
		logFC_thres = log(2)
	}
	
	# Subset Genes and make SingleCellAssay object
	if( is.null(num_genes) ){
		num_genes = nrow(sce_obj)
	}
	sca = MAST::SceToSingleCellAssay(sce = sce_obj[seq(num_genes),])

	# Calculate CDS
	colData(sca)$cngeneson = scale(colSums(assay(sca) > 0))

	# Calculate Group
	# one_cell_type = "Astro"
	ref_group = paste0("not_",one_cell_type)
	colData(sca)$Group = ifelse(colData(sca)$cell_type == one_cell_type,
		one_cell_type,ref_group)
	colData(sca)$Group = factor(colData(sca)$Group,
		levels = c(ref_group,one_cell_type))

	# zlm analysis: Model log-transformed expression as 
	#	function of clustered cell type and num detected genes
	zlm_fn = file.path(work_dir,paste0("zlm",
		"_nG",num_genes,
		"_cell",one_cell_type,
		".rds"))
	if( !file.exists(zlm_fn) ){
		print(date())
		zlm_output = MAST::zlm(formula = ~ Group + cngeneson,sca = sca)
		print(date())
		cat("Saving image ...\n")
		saveRDS(zlm_output,zlm_fn)
	}
	cat("Reading in image ...\n")
	zlm_output = readRDS(zlm_fn)
	# show(zlm_output)

	# Perform LRT of one cell_type against all others
	cat("Running LRT by excluding the cell_type covariate ...\n")
	print(date())
	ssc = MAST::summary(object = zlm_output,
		doLRT = paste0("Group",one_cell_type))
	print(date())
	ssd0 = smart_df(ssc$datatable)
	ssd = ssd0[which(ssd0$component == "H"),]
	ssd = name_change(ssd,"Pr..Chisq.","pvalue")
	ssd = name_change(ssd,"primerid","gene")
	ssd = ssd[,c("gene","pvalue")]
	tmp_index = which(ssd0$component == "logFC"
					& ssd0$contrast == paste0("Group",one_cell_type))
	ssd$logFC = ssd0$coef[tmp_index]
	ssd$z = ssd0$z[tmp_index]
	ssd = ssd[which(!is.na(ssd$logFC)),]
	ssd$FDR_qvalue = p.adjust(p = ssd$pvalue,method = "fdr")
	ssd = ssd[order(ssd$FDR_qvalue),]
	# Plot
	# smart_hist(ssd$FDR_qvalue,breaks=40,xlab="FDR",main="")
	# smart_hist(ssd$logFC,breaks=40)
	# plot(ssd[,c("FDR_qvalue","logFC")],pch=16,col=rgb(0,0,0,0.5))
	# dev.off()
	
	# Subsetting supposed DE genes
	# ssd = ssd[which(ssd$FDR_qvalue < fdr_thres & ssd$logFC > logFC_thres),]
	
	# dim(ssd); smart_table(ssd$logFC > 0); ssd[1:20,]
	ssd = smart_df(cell_type = one_cell_type,ssd)
	
	ssd_fn = file.path(work_dir,paste0("ssd",
		"_nG",num_genes,
		"_cell",one_cell_type,
		".rds"))
	cat(paste0("Saving image in ",ssd_fn," ...\n"))
	saveRDS(ssd,ssd_fn)

	return(NULL)
}


# ----------
# Analyze
# ----------
sce = readRDS(file.path(MTG_dir,"final_sce_filtered_by_kmeans.rds"))
sce
# sort(names(colData(sce)))

# Subset cell_types with enough cells: Exclude Endo
smart_table(colData(sce)$cell_type)
sce = sce[,colData(sce)$cell_type != "Endo"]
smart_table(colData(sce)$cell_type)
dim(sce)

# Get cell type specific marker genes
one_ct = "blah_one_ct"
# one_ct = "Astro"
# one_ct = "Exc"
# one_ct = "Inh"
# one_ct = "Micro"
# one_ct = "Oligo"
# one_ct = "OPC"
MAST_DEgenes(work_dir = MTG_dir,
			num_genes = nrow(sce),
			# num_genes = 100,
			sce_obj = sce,
			one_cell_type = one_ct,
			fdr_thres = 1e-3,
			logFC_thres = log(2))

# ssd = readRDS(file.path(MTG_dir,paste0("ssd_nG100_cellAstro.rds")))
# dim(ssd)
q("no")


# ----------
# VennDiagram and getting/annotating marker genes
# ----------
if(FALSE){

cell_types = c("Astro","Exc","Inh","Micro","Oligo","OPC")

ct_genes = list()
for(ct in cell_types){
	ct_genes[[ct]] = readRDS(file.path(MTG_dir,paste0("ssd_nG37657_cell",ct,".rds")))$gene
}
saveRDS(ct_genes,file.path(MTG_dir,"ct_genes.rds"))

smart_pack("venn")
pdf(file.path(MTG_dir,"venn.pdf"),width=8,height=8)
venn(x = ct_genes,ilabels = TRUE,zcolor = "style")
dev.off()

# Get marker genes
cts_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	cts_genes[[ct]] = sort(setdiff(ct_genes[[ct]],
		unique(unlist(ct_genes[which(names(ct_genes) != ct)]))))
}
saveRDS(cts_genes,file.path(MTG_dir,"cts_genes.rds"))

# Get gene annotations
sce = readRDS(file.path(MTG_dir,"final_sce_filtered_by_kmeans.rds"))
sce
rowData(sce)

# Import/Prep gtf
rsem_fn = file.path(MTG_dir,"rsem_GRCh38.p2.gtf")
if( !file.exists(rsem_fn) ){
	gtf_link = "http://celltypes.brain-map.org/api/v2/well_known_file_download/502175284"
	gtf_fn0 = strsplit(gtf_link,"/")[[1]]
	gtf_fn0 = gtf_fn0[length(gtf_fn0)]
	gtf_fn = file.path(MTG_dir,gtf_fn0)
	sprintf("cd %s; wget %s; unzip %s",MTG_dir,gtf_link,gtf_fn)
}
rsem = smart_RT(rsem_fn,sep="\t",header=FALSE)
names(rsem) = c("seqname","source","feature","start_position","end_position",
	"score","strand","frame","attributes")
smart_table(rsem$feature)
rsem = rsem[rsem$feature == "gene",]
rsem$gene_id = sapply(rsem$attributes,function(xx) 
	gsub("gene_id ","",strsplit(xx,";")[[1]][1]),USE.NAMES=FALSE)
rsem$symbol = sapply(rsem$attributes,function(xx) 
	gsub(" gene_symbol ","",strsplit(xx,";")[[1]][2]),USE.NAMES=FALSE)
rsem = rsem[,c("symbol","gene_id","start_position","end_position")]
all(rownames(sce) %in% rsem$symbol)
all(rowData(sce)$entrez_id %in% rsem$gene_id)
rsem = rsem[which(rsem$symbol %in% rowData(sce)$gene),]
rsem[1:10,]
rsem = rsem[match(rowData(sce)$gene,rsem$symbol),]
rownames(rsem) = NULL
rsem[1:10,]
all(rowData(sce)$gene == rsem$symbol)
all(rowData(sce)$entrez_id == rsem$gene_id)
rowData(sce)$start_position = rsem$start_position
rowData(sce)$end_position = rsem$end_position

gene_anno = smart_df(rowData(sce)[,c("gene","chromosome",
	"entrez_id","start_position","end_position")])
rownames(gene_anno) = NULL

res_de = c()
for(ct in cell_types){
	# ct = cell_types[1]
	one_ct_genes = cts_genes[[ct]]
	one_res = readRDS(file.path(MTG_dir,paste0("ssd_nG37657_cell",ct,".rds")))
	one_res = one_res[one_res$gene %in% one_ct_genes,]
	# one_res[1:10,]
	res_de = rbind(res_de,one_res)
}

gene_anno = gene_anno[gene_anno$gene %in% res_de$gene,]
gene_anno = gene_anno[match(res_de$gene,gene_anno$gene),]
all(gene_anno$gene == res_de$gene)
gene_anno = cbind(gene_anno,res_de[,names(res_de) != "gene"])
saveRDS(gene_anno,file.path(MTG_dir,"anno_marker_genes.rds"))


}





q("no")

###
