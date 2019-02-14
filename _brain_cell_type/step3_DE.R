# scRNAseq differential expression (seems to require at least 80 gbs)
# srun -c 1 --time=8:00:00 --mem=80000 --job-name="Rint" --pty -p interact /nas/longleaf/home/pllittle/downloads/R-3.5.1/bin/R --vanilla

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
	ssc = MAST::summary(object = zlm_output,
		doLRT = paste0("Group",one_cell_type))
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
	ssd = ssd[which(ssd$FDR_qvalue < fdr_thres 
					& ssd$logFC > logFC_thres),]
	# dim(ssd); smart_table(ssd$logFC > 0); ssd[1:20,]
	ssd = smart_df(cell_type=one_cell_type,ssd)
	
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

###
