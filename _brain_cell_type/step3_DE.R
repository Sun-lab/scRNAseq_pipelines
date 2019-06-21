
# scRNAseq differential expression
# NOTE: the code require large amount of memory and may not work 
#       on a laptop or a desktop. 

# ------------------------------------------------------------
# Shortcuts
# ------------------------------------------------------------

rm(list=ls())
repo_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/scRNAseq_pipelines"
# repo_dir = file.path("..") # relative path of repo directory from this R code

source(file.path(repo_dir,"SOURCE.R"))
MTG_dir = file.path(repo_dir,"MTG")
setwd(MTG_dir)

# ------------------------------------------------------------
# a raw data needed for the analysis is too large ()
# final_sce_filtered_by_kmeans.rds
# ------------------------------------------------------------

rawData_dir = "human_MTG_gene_expression_matrices_2018-06-14"
rawData_dir = file.path("~/research/scRNAseq/data/Allen_BI/", rawData_dir)
rawData_dir

# ------------------------------------------------------------
# Libraries/Functions
# ------------------------------------------------------------

if( !("BiocManager" %in% installed.packages()[,"Package"]) ){
  install.packages("BiocManager",repos = "https://mirrors.nics.utk.edu/cran/")
}

if( !("MAST" %in% installed.packages()[,"Package"]) ){
  # MAST = Model-based Analysis of Single Cell Transcriptomics
  BiocManager::install("MAST")
}

library(ggplot2)
library(data.table)
library(scater)

GTF_calc_gene_length = function(work_dir,rsem_gtf_fn){
	if(FALSE){
		work_dir = MTG_dir
		rsem_gtf_fn = file.path(work_dir,"rsem_GRCh38.p2.gtf")
	}

	# Download GTF file
	if( !file.exists(rsem_gtf_fn) ){
		gtf_link = "http://celltypes.brain-map.org/api/v2/well_known_file_download/502175284"
		gtf_fn = strsplit(gtf_link,"/")[[1]]
		gtf_fn = gtf_fn[length(gtf_fn)]
		gtf_fn = file.path(work_dir,gtf_fn)
		system(sprintf("cd %s; wget %s; unzip %s",work_dir,gtf_link,gtf_fn))
	}
  
	# Get mapping from gene_id to gene_symbol
	gtf = data.table::fread(rsem_gtf_fn,header = FALSE,data.table = FALSE)
	dim(gtf); gtf[1:3,]
	names(gtf) = c("seqname","source","feature",
				"start_position","end_position",
				"score","strand","frame","attributes")
	gtf = gtf[which(gtf$feature %in% c("gene")),]
	gtf$gene_id = sapply(gtf$attributes,function(xx) 
		gsub("\"","",gsub("gene_id ","",strsplit(xx,";")[[1]][1])),
		USE.NAMES = !TRUE)
	gtf$gene_symbol = sapply(gtf$attributes,function(xx) 
		gsub("\"","",gsub(" gene_symbol ","",strsplit(xx,";")[[1]][2])),
		USE.NAMES = !TRUE)
  
	# Import gtf to database
	library(GenomicRanges)
	library(GenomicFeatures)
	cat("Make TxDb...\n")
	exdb = suppressWarnings(GenomicFeatures::makeTxDbFromGFF(
		file = rsem_gtf_fn,format = "gtf"))
	exons_list_per_gene = GenomicFeatures::exonsBy(exdb,by = "gene")
  
	# Get intersection
	gene_intersect = intersect(gtf$gene_id,names(exons_list_per_gene))
	exons_list_per_gene = exons_list_per_gene[names(exons_list_per_gene) %in% gene_intersect]
	gtf = gtf[which(gtf$gene_id %in% gene_intersect),]
  
	# Calculate exonic gene lengths
	tmp_df = smart_df(gene_id = names(exons_list_per_gene),
	gene_length = as.numeric(sum(width(GenomicRanges::reduce(exons_list_per_gene)))))
	tmp_df = tmp_df[match(gtf$gene_id,tmp_df$gene_id),]
	all(gtf$gene_id == tmp_df$gene_id)
	gtf$gene_length = tmp_df$gene_length
  
	# Output
	gtf
}
MAST_DEgenes = function(work_dir,num_genes=NULL,sce_obj,one_cell_type){
  
	if(FALSE){
		work_dir = MTG_dir
		num_genes = NULL
		sce_obj = sce
		one_cell_type = c("Astro","Exc","Inh","Micro","Oligo","OPC")[1]
		# fdr_thres = 1e-3; logFC_thres = log(2)
	}

	setwd(work_dir)
	
	# Subset Genes and make SingleCellAssay object
	if( is.null(num_genes) ){
		num_genes = nrow(sce_obj)
	}
	sca = MAST::SceToSingleCellAssay(sce = sce_obj[seq(num_genes),])
  
	# Exclude any subjects with no gene counts
	sca = sca[,colSums(counts(sca)) > 0]
  
	# Make ET assay log2(TPM+1)
	## Calculate log2(TPM+1)
	cat("Calculating log2(TPM+1)...\n")
	calc_count_to_log2_TPM_1 = function(vec_count,vec_lengths){
		xx = vec_count / vec_lengths
		xx = xx / sum(xx)
		TPM = xx * 1e6
		log2(1 + TPM)
	}
  
	assay(sca,"Et") = apply(counts(sca),2,function(xx) 
		calc_count_to_log2_TPM_1(xx,rowData(sca)$gene_length))
	assay(sca,"counts") = NULL
	assay(sca,"logcounts") = NULL
  
	# Calculate CDS, only after Et assay is created
	colData(sca)$cngeneson = scale(colSums(assay(sca,"Et") > 0))
  
	# Calculate Group
	ref_group = paste0("not_",one_cell_type)
	colData(sca)$Group = ifelse(colData(sca)$cell_type == one_cell_type,
		one_cell_type,ref_group)
	colData(sca)$Group = factor(colData(sca)$Group,
		levels = c(ref_group,one_cell_type))
  
	# zlm analysis: Model log-transformed expression as function 
	#	of clustered cell type and num detected genes
	print(date())
	zlm_output = MAST::zlm(formula = ~ Group + cngeneson,sca = sca)
	print(date())

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
	ssd = smart_df(cell_type = one_cell_type,ssd)
  
	# Output
	ssd_fn = paste0("ssd","_nG",num_genes,"_cell",one_cell_type,".rds")
	cat(paste0("Saving image in ",ssd_fn," ...\n"))
	saveRDS(ssd,ssd_fn)
  
	return(NULL)
}

# ICeDT related Functions
q90 = function(vv){
	qs = quantile(vv, probs=c(0.10, 0.90))
	qs[2] - qs[1]
}
plot_log1p = function(x, y, ...) {
	smoothScatter(log(x+1e-5), log(y+1e-5), xlim=c(-5, 10), ylim=c(-5, 10), ...)
	legend("bottom",bty="n",legend=sprintf("Pearson correlation = %.2f",
		cor(log(x+1e-5), log(y+1e-5))))
}
ICeDT_consistency = function(sig,bulk,ICeDT_out){

	p1 = ICeDT_out$cProb
	prop_icedt = t(ICeDT_out$rho)[,-1]
	# dim(p1)
	p1 = data.matrix(p1)
	
	par(mfrow=c(2,3),mar=c(4,4,1,0.5),bty="n")
	# par(mfrow=c(1,3))
	plot(density(c(p1))$y,main="",xlim=c(0,1),xlab="probability consistent",
		ylab="density",type="n")
	lines(density(c(p1)), lty=1, col="black")
	legend("topright", c("no weight"), lty=c(1,2),col=c("black"), bty="n")
	plot(apply(p1,1,median),apply(p1,1,q90),xlab="median prob. consistent",
		ylab="90 percentile - 10 percentile",main="Gene Consistency across subjects",
		xlim=c(0,1),ylim=c(0,1))
	boxplot(prop_icedt,main = "Simulation - ICeDT")
	# par(mfrow=c(1,1))
  
	predicted_bulk_w0 = sig %*% t(prop_icedt)
	p1_cutoffs = quantile(p1,c(1/3,2/3)); p1_cutoffs
  
	# par(mfrow=c(1,3))
	plot_log1p(x = c(predicted_bulk_w0)[p1 < p1_cutoffs[1]],
		y = c(bulk)[p1 < p1_cutoffs[1]], 
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "low prob of being consistent")
	plot_log1p(x = c(predicted_bulk_w0)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]],
		y = c(bulk)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]], 
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "med prob of being consistent")
	plot_log1p(x = c(predicted_bulk_w0)[p1 > p1_cutoffs[2]],
		y = c(bulk)[p1 > p1_cutoffs[2]],
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "high prob of being consistent")
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,bty="o")
}



# ------------------------------------------------------------
# Import and Prep Data
# ------------------------------------------------------------
# Import gtf with gene lengths
rsem_fn = "rsem_GRCh38.p2.gtf"
gtf = GTF_calc_gene_length(work_dir = MTG_dir,rsem_gtf_fn = rsem_fn)
head(gtf)

# sce = readRDS(file.path(rawData_dir,"final_sce_filtered_by_kmeans.rds"))
sce = readRDS(file.path(MTG_dir,"final_sce_filtered_by_kmeans.rds"))
sce

# Append gene_lengths to sce
inter_genes = intersect(rownames(sce),gtf$gene_symbol)
length(inter_genes)

sce = sce[inter_genes,]
gtf = gtf[which(gtf$gene_symbol %in% inter_genes),]
gtf = gtf[match(rownames(sce),gtf$gene_symbol),]
rowData(sce)$gene_length = gtf$gene_length

# Subset cell_types with enough cells: Exclude Endo
smart_table(colData(sce)$cell_type)
sce = sce[,colData(sce)$cell_type != "Endo"]
smart_table(colData(sce)$cell_type)
dim(sce)
cell_types 	= c("Astro","Exc","Inh","Micro","Oligo","OPC")
num_genes = nrow(sce)


# ------------------------------------------------------------
# Outlined Steps
# ------------------------------------------------------------
# Run MAST to get differentially expressed genes
# Get disjoint set of marker genes per cell type (filter on logFC and qvalue), 
#    subset approximately 100 genes per cell type
# Use sce counts and exonic gene lengths to calculate TPM counts
# Bulk RNAseq deconvolution: Run CIBERSORT, ICEDT


# ------------------------------------------------------------
# Run MAST: Run each cell type as it's own job
# ------------------------------------------------------------
for(ct in cell_types){
	print(ct)
	if(FALSE){
	MAST_DEgenes(work_dir = MTG_dir,
				num_genes = num_genes,
				sce_obj = sce,
				one_cell_type = one_ct)
	}
}

# ------------------------------------------------------------
# Combine all genes with annotations and proportion of 
# expressed cells per cell type
# ------------------------------------------------------------
# Get gene annotations
all(rowData(sce)$gene == gtf$gene_symbol)
gtf$chromosome = rowData(sce)$chromosome
gtf = name_change(gtf,"gene_id","entrez_id")
gtf$prop_express = as.numeric(apply(assay(sce),1,function(xx) mean(xx > 0)))

# For each cell type, gather logFC, qvalues, prop cells expressed in genes
for(ct in cell_types){
	# ct = cell_types[1]
	cat(paste0(ct," "))
	tmp_df = readRDS(paste0("ssd_nG",num_genes,"_cell",ct,".rds"))
	tmp_df = tmp_df[,c("gene","pvalue","logFC","FDR_qvalue")]
	tmp_df = name_change(tmp_df,"gene","gene_symbol")
	names(tmp_df)[-1] = paste0(names(tmp_df)[-1],".",ct)
	# Calculating proportion of cells expressed in a gene and cell_type
	bb = apply(assay(sce[,colData(sce)$cell_type == ct]),1,
		function(xx) mean(xx > 0))
	bb = smart_df(gene_symbol = names(bb),prop_express = as.numeric(bb))
	names(bb)[2] = paste0(names(bb)[2],".",ct)
	bb = smart_merge(bb,tmp_df,all.x=TRUE)
	all(bb$gene_symbol == gtf$gene_symbol)
	bb = bb[match(gtf$gene_symbol,bb$gene_symbol),]
	all(bb$gene_symbol == gtf$gene_symbol)
	gtf = cbind(gtf,bb[,names(bb) != "gene_symbol"])
	rm(tmp_df,bb)
}
saveRDS(gtf,"DE_gene_anno.rds")


# ------------------------------------------------------------
# Get marker genes per cell type
# ------------------------------------------------------------
dat = smart_df(rowData(sce)[,c("gene","chromosome","entrez_id")])
fdr_thres 		= 1e-3
logFC_thres	 	= log(2)

# Get each cell type's differentially expressed genes
ct_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	print(ct)
	rds_fn = paste0("ssd_nG",num_genes,"_cell",ct,".rds")
	rds = readRDS(rds_fn)
	rds = smart_merge(rds,dat,all.x=TRUE)
	rds = rds[which(rds$logFC > logFC_thres 
					& rds$FDR_qvalue < fdr_thres
					& rds$chromosome %in% c(1:22,"X","Y")
					& !grepl("^LOC10",rds$gene)),]
	ct_genes[[ct]] = rds$gene; rm(rds)
}
sapply(ct_genes,length)
saveRDS(ct_genes,"ct_genes.rds")

# Get each cell type's disjoint set of differentially expressed genes
cts_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	print(ct)
	cts_genes[[ct]] = sort(setdiff(ct_genes[[ct]],
		unique(unlist(ct_genes[which(names(ct_genes) != ct)]))))
}
sapply(cts_genes,length)
saveRDS(cts_genes,"cts_genes.rds")

# Get approximate top 120 marker genes per cell type
## Rather than arbitrarily selecting 100 genes, I apply 
##	a quantile threshold on both logFC and FDR_qvalue to 
##	select approximately 100 genes per cell type.
set.seed(1)
mark_genes = list()
num_mgene_per_ct = 120 # number of marker genes selected per cell type
for(ct in cell_types){
	# ct = cell_types[1]
	cat(paste0(ct,": "))
	rds_fn = paste0("ssd_nG",num_genes,"_cell",ct,".rds")
	rds = readRDS(rds_fn)
	rds = rds[which(rds$logFC > logFC_thres 
					& rds$FDR_qvalue < fdr_thres
					& rds$gene %in% cts_genes[[ct]]),]
	# dim(rds)
	low_bound = 0
	while(TRUE){
		# qu = 0.25
		qu = runif(1,low_bound,1)
		num_rows = nrow(rds[which(rds$logFC > quantile(rds$logFC,qu) 
			& rds$FDR_qvalue < quantile(rds$FDR_qvalue,1-qu)),])
		if( num_rows > num_mgene_per_ct ){
			# Tighten bounds
			low_bound = qu # mean(c(low_bound[1],qu))
			cat(paste0(round(low_bound,3)," "))
			if( num_rows <= num_mgene_per_ct+10 ) break
		}
	}
	cat("\n")
	opt_qu = qu
	rds = rds[which(rds$logFC > quantile(rds$logFC,opt_qu) 
					& rds$FDR_qvalue < quantile(rds$FDR_qvalue,1-opt_qu)),]
	mark_genes[[ct]] = rds$gene; rm(rds)
}
sapply(mark_genes,length)
saveRDS(mark_genes,"mark_genes.rds")

smart_pack("venn")
pdf("ct_genes_MTG.pdf",width=8,height=8)
venn(x = ct_genes,ilabels = TRUE,zcolor = "style")
dev.off()


# ------------------------------------------------------------
# Calculate TPM and get signature matrix: Based on Chong Jin's code
# ------------------------------------------------------------
table(colData(sce)$cell_type)

# Construct signature matrix
SIG = matrix(nrow = nrow(sce),ncol = length(cell_types))
colnames(SIG) = cell_types
rownames(SIG) = rownames(sce)
for(ct in cell_types){
	print(ct)
	# SIG[,ct] = rowSums(counts(sce)[,colData(sce)$cell_type == ct])
	SIG[,ct] = rowMeans(counts(sce)[,colData(sce)$cell_type == ct])
}
dim(SIG)
head(SIG)

# Subset marker genes
SIG = SIG[sort(as.character(unlist(mark_genes))),]
dim(SIG)
head(SIG)
summary(SIG)

# Get ENSG id
dat	= gtf[,c("gene_symbol","chromosome","entrez_id")]
rownames(dat) = NULL
dat = dat[match(rownames(SIG),dat$gene_symbol),]
biomaRt::listEnsembl()
ensembl = biomaRt::useMart("ensembl")
dd = biomaRt::listDatasets(ensembl); dim(dd); dd[1:5,]
dd[grep("hsap",dd$dataset),]
ee = biomaRt::useDataset("hsapiens_gene_ensembl",mart = ensembl)
# ff = biomaRt::listFilters(ee); dim(ff); ff[1:5,]
# aa = biomaRt::listAttributes(ee); dim(aa); aa[1:5,]

attr_string = c("hgnc_symbol","ensembl_gene_id","external_gene_name",
				"entrezgene","chromosome_name","start_position",
				"end_position","strand","gene_biotype")
filters = c("hgnc_symbol","chromosome_name","entrezgene")
values = list(dat$gene_symbol,dat$chromosome,dat$entrez_id)
gene_anno = biomaRt::getBM(attributes = attr_string,
						filters = filters, 
						values = values,
						mart = ee)
gene_anno = unique(gene_anno)
dim(gene_anno); gene_anno[1:5,]
smart_table(rownames(SIG) %in% gene_anno$hgnc_symbol)

# What are the missing genes?
miss_genes = setdiff(rownames(SIG),gene_anno$hgnc_symbol)
miss_genes

# Subset and order genes
inter_genes = intersect(rownames(SIG),gene_anno$hgnc_symbol)
gene_anno = gene_anno[which(gene_anno$hgnc_symbol %in% inter_genes),]
gtf = gtf[which(gtf$gene_symbol %in% inter_genes),]
SIG = SIG[which(rownames(SIG) %in% inter_genes),]
gtf2 = gtf[,c("gene_symbol","gene_length")]
	gtf2 = name_change(gtf2,"gene_symbol","hgnc_symbol")
gene_anno = smart_merge(gene_anno,gtf2[,c("hgnc_symbol","gene_length")])
gene_anno = gene_anno[match(rownames(SIG),gene_anno$hgnc_symbol),]
dim(gene_anno)
head(gene_anno)

# Output
saveRDS(list(anno = gene_anno,SIG = SIG),"signature_MTG.rds")


# ------------------------------------------------------------
# Import CMC Bulk RNA, get gene lengths
# ------------------------------------------------------------
bulk = readRDS("CMC_MSSM-Penn-Pitt_Paul_geneExpressionRaw.rds")$so1
dim(bulk); bulk[1:5,1:5]

gtf_fn = "Homo_sapiens.GRCh37.70.processed.gtf"
exdb = GenomicFeatures::makeTxDbFromGFF(file = gtf_fn,format = "gtf")
exons_list_per_gene = GenomicFeatures::exonsBy(exdb,by = "gene")
tmp_df = smart_df(ensembl_gene_id = names(exons_list_per_gene),
				gene_length = as.numeric(sum(width(GenomicRanges::reduce(exons_list_per_gene)))))
tmp_df[1:5,]
all(tmp_df$ensembl_gene_id %in% rownames(bulk))
all(tmp_df$ensembl_gene_id == rownames(bulk))

# Subset/Order marker genes
gene_anno 		= readRDS("signature_MTG.rds")$anno
inter_genes2 	= intersect(tmp_df$ensembl_gene_id,gene_anno$ensembl_gene_id)
gene_anno2 		= gene_anno[which(gene_anno$ensembl_gene_id %in% inter_genes2),]
bulk 					= bulk[which(rownames(bulk) %in% inter_genes2),]
bulk 					= bulk[match(gene_anno2$ensembl_gene_id,rownames(bulk)),]
tmp_df 				= tmp_df[which(tmp_df$ensembl_gene_id %in% inter_genes2),]
tmp_df 				= tmp_df[match(gene_anno2$ensembl_gene_id,tmp_df$ensembl_gene_id),]

saveRDS(list(anno = tmp_df,bulk = bulk),"bulk_marker_MTG.rds")


# ------------------------------------------------------------
# Deconvolution: Based on Chong Jin's deconvolution.Rmd code
# ------------------------------------------------------------
# Refer to EPIC paper for variable definitions (bb = bulk,cc = signature,pp = cell type proportions)
bulk_rds 	= readRDS("bulk_marker_MTG.rds")
sig_rds 	= readRDS("signature_MTG.rds")
inter_genes3 	= intersect(bulk_rds$anno$ensembl_gene_id,sig_rds$anno$ensembl_gene_id)
bulk_rds$anno = bulk_rds$anno[which(bulk_rds$anno$ensembl_gene_id %in% inter_genes3),]
bulk_rds$bulk = bulk_rds$bulk[which(rownames(bulk_rds$bulk) %in% inter_genes3),]
sig_rds$anno 	= sig_rds$anno[which(sig_rds$anno$ensembl_gene_id %in% inter_genes3),]
sig_rds$SIG 	= sig_rds$SIG[which(rownames(sig_rds$SIG) %in% sig_rds$anno$hgnc_symbol),]
rownames(sig_rds$SIG) = sig_rds$anno$ensembl_gene_id

bb_tpm = apply(bulk_rds$bulk,2,function(xx) xx / bulk_rds$anno$gene_length)
bb_tpm = 1e6 * apply(bb_tpm,2,function(xx) xx / sum(xx))
cc_tpm = apply(sig_rds$SIG,2,function(xx) xx / sig_rds$anno$gene_length)
cc_tpm = 1e6 * apply(cc_tpm,2,function(xx) xx / sum(xx))

cell_sizes = colSums(apply(sig_rds$SIG,2,function(xx) xx / sig_rds$anno$gene_length))
cell_sizes

# Run ICeDT
pack = "ICeDT"
if( !(pack %in% installed.packages()[,"Package"]) ){
	devtools::install_github("Sun-lab/ICeDT")
}
date()
fit = ICeDT::ICeDT(Y = bb_tpm,Z = cc_tpm,tumorPurity = rep(0,ncol(bb_tpm)),
	refVar = NULL,rhoInit = NULL,maxIter_prop = 4e3,maxIter_PP = 4e3,
	rhoConverge = 1e-2)
date()
pp_bar_icedt = t(fit$rho)[,-1]
# calculate p_hat using cell_sizes
pp_hat_icedt = t(apply(pp_bar_icedt,1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
# Look at distribution of pp_bar and pp_hat
summary(pp_bar_icedt)
summary(pp_hat_icedt)

# Run CIBERSORT
sig_fn = file.path(MTG_dir,"signature_MTG.txt")
mix_fn = file.path(MTG_dir,"mixture_CMC.txt")
write.table(cbind(rowname=rownames(cc_tpm),cc_tpm),
	file = sig_fn,sep = "\t",quote = FALSE,row.names = FALSE)
write.table(cbind(rowname=rownames(bb_tpm),bb_tpm),
	file = mix_fn,sep = "\t",quote = FALSE,row.names = FALSE)
# Login to CIBERSORT website, uploaded above two files, 
#	specified no quantile normalization, ran 1000 permutations
# OR Request CIBERSORT.R file from website, install dependent packages
cibersort_src_fn = "~/github/CSeQTL/R/CIBERSORT.R"
# cibersort_src_fn = "." # whichever directory contains CIBERSORT.R
source(cibersort_src_fn)
print(date())
results = CIBERSORT(sig_matrix = sig_fn,mixture_file = mix_fn,
	perm = 0,QN = FALSE,absolute = FALSE,abs_method = 'sig.score',
	filename = "MTG") # added filename argument for function
print(date())
unlink(sig_fn)
unlink(mix_fn)
ciber_fn = sprintf("CIBERSORT-Results_%s.txt","MTG")
unlink(ciber_fn)
QQ = ncol(cc_tpm)
pp_bar_ciber = results[,seq(QQ)]
# calculate p_hat using cell sizes
pp_hat_ciber = t(apply(pp_bar_ciber,1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
summary(pp_bar_ciber)
summary(pp_hat_ciber)

# Output MTG deconvolution results
pdf("MTG_bulk_deconvolution_summary.pdf",height=8,width=12)
	
	ICeDT_consistency(sig = cc_tpm,bulk = bb_tpm,ICeDT_out = fit)
	
	par(mfrow=c(1,2),mar=c(4,4,1,0.5),oma=c(0,0,2,0),bty="n")
	boxplot(pp_bar_icedt,main="Expression Proportions",ylim=c(0,1))
	boxplot(pp_hat_icedt,main="Cell type Proportions",ylim=c(0,1))
	mtext("ICeDT Results",outer=TRUE,cex=1.3)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
	
	par(mfrow=c(1,2),mar=c(4,4,1,0.5),oma=c(0,0,2,0),bty="n")
	boxplot(pp_bar_ciber,main="Expression Proportions",ylim=c(0,1))
	boxplot(pp_hat_ciber,main="Cell type Proportions",ylim=c(0,1))
	mtext("CIBERSORT Results",outer=TRUE,cex=1.3)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
	
	par(mfrow=c(2,3),mar=c(4,4,1,0.5),oma=c(0,0,2,0),bty="n")
	for(ct in seq(QQ)){
		v1 = pp_hat_ciber[,ct]
		v2 = pp_hat_icedt[,ct]
		tmp_range = range(c(v1,v2))
		pcor = round(cor(v1,v2),4)
		scor = round(cor(v1,v2,method="spear"),4)
		plot(v1,v2,xlim=tmp_range,ylim=tmp_range,
			xlab="CIBERSORT",ylab="ICeDT",
			main=sprintf("%s:Pear=%s;Spear=%s",colnames(pp_hat_ciber)[ct],
				pcor,scor))
		abline(a=0,b=1,lty=2,lwd=2,col="red")
	}
	mtext("CIBERSORT v. ICeDT",outer=TRUE,cex=1.4)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
	
dev.off()

saveRDS(list(ICeDT = pp_hat_icedt,CIBERSORT = pp_hat_ciber),"prop_MTG.rds")


# ------------------------------------------------------------
# Run deconvolution using Chong's psychENCODE scRNA clustering results (copying small lines of code from signature_genes.Rmd)
# ------------------------------------------------------------
options(width=150)
psychENCODE_dir = "/pine/scr/c/h/chongjin/psychENCODE_data"

# Import psychENCODE data
sce = readRDS(file.path(psychENCODE_dir,"sce.rds"))
sce
rowData(sce)[1:3,] 
	# Columns ensembl_gene_id, hgnc_symbol, 
	#		chromosome_name, gene_length available
colData(sce)[1:3,]
	table(colData(sce)$part_cell_id)
	colData(sce)$cell_type = colData(sce)$part_cell_id
	colData(sce)$cell_type[grep("In", colData(sce)$part_cell_id)] = "Inh"
	colData(sce)$cell_type[grep("Ex", colData(sce)$part_cell_id)] = "Exc"
	colData(sce)$cell_type[grep("Oligo", colData(sce)$part_cell_id)] = "Oligo"  # no change
	colData(sce)$cell_type[grep("OPC", colData(sce)$part_cell_id)] = "OPC"  # no change
	colData(sce)$cell_type[grep("Astro", colData(sce)$part_cell_id)] = "Astro"  # no change
	colData(sce)$cell_type[grep("Microglia", colData(sce)$part_cell_id)] = "Micro"  #   no change
	colData(sce)$cell_type[grep("Endo", colData(sce)$part_cell_id)] = "Endo"  #   no change
	table(colData(sce)$cell_type)

# Import cluster results
clusters = readRDS(file.path(psychENCODE_dir,"k25_50_pcs.rds"))
clusters = data.frame(KM25=clusters$cluster)
clusters$cell_type = colData(sce)$cell_type
opt_clust = paste0("KM",25) # our clustering results differ slightly
t1 = table(clusters[,opt_clust], clusters$cell_type)
t1
t1[c(15,23),] = 0
clusts = apply(t1, 2, function(v){union(which.max(v), which(v > 100))})
clusts
cell_types = setdiff(unique(clusters$cell_type), "unknown")
cell_types
w2kp = NULL
for(ct1 in cell_types){
	# ct1 = celltypes[1]
	ct.cond    = clusters$cell_type == ct1
	clust.cond = clusters[,opt_clust] %in% clusts[[ct1]]
	# Intersect previous clustering result with ours
	w2kp = c(w2kp, which(ct.cond & clust.cond))
}
length(w2kp)
sce = sce[,w2kp]
dim(sce)

# Import gene annotation with differential expression results 
gene_anno = as.data.frame(rowData(sce)[,c("ensembl_gene_id",
	"external_gene_name","chromosome_name","start_position","end_position",
	"gene_length")])
gene_anno = dplyr::rename(gene_anno, gene=external_gene_name)  
gene_anno$prop_express = as.numeric(apply(assay(sce),1,function(xx) mean(xx > 0)))

# For each cell type, gather logFC, qvalues, prop cells expressed in genes
cell_types = c("Astro","Exc","Inh","Micro","Oligo","OPC")
num_genes = nrow(sce)
for(ct in cell_types){
	# ct = "Astro"
	cat(paste0(ct," "))
	tmp_df = readRDS(file.path(psychENCODE_dir,paste0("ssd_nG",num_genes,"_cell",ct,".rds")))
	tmp_df = tmp_df[,c("gene","pvalue","logFC","FDR_qvalue")]
	names(tmp_df)[-1] = paste0(names(tmp_df)[-1],".",ct)
	# Calculating proportion of cells expressed in a gene and cell_type
	bb = apply(assay(sce[,colData(sce)$cell_type %in% ct]),1,
			function(xx) mean(xx > 0))
	bb = data.frame(gene = names(bb),prop_express = as.numeric(bb), stringsAsFactors = FALSE)
	names(bb)[2] = paste0(names(bb)[2],".",ct)
	bb = merge(bb,tmp_df,by="gene",all.x=TRUE,sort=FALSE)
	all(bb$gene == gene_anno$gene)
	bb = bb[match(gene_anno$gene,bb$gene),]
	all(bb$gene == gene_anno$gene)
	gene_anno = cbind(gene_anno, bb[,names(bb) != "gene"])
	rm(tmp_df,bb)
}
dim(gene_anno)
head(gene_anno)
saveRDS(gene_anno,file.path(MTG_dir,"DE_gene_anno_psychENCODE.rds"))

# Get each cell type's differentially expressed genes
fdr_thres 		= 1e-3
logFC_thres	 	= log(2)
ct_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	print(ct)
	ct_genes[[ct]] = gene_anno$gene[which(
		gene_anno[,sprintf("logFC.%s",ct)] > logFC_thres
		& gene_anno[,sprintf("FDR_qvalue.%s",ct)] < fdr_thres
		& gene_anno$chromosome_name %in% c(1:22,"X","Y")
		)]
}
sapply(ct_genes,length)

# Get each cell type's disjoint set of differentially expressed genes
cts_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	print(ct)
	cts_genes[[ct]] = sort(setdiff(ct_genes[[ct]],
		unique(unlist(ct_genes[which(names(ct_genes) != ct)]))))
}
sapply(cts_genes,length)

# Get approximate top 120 marker genes per cell type
## Rather than arbitrarily selecting 120 genes, I apply 
##	a quantile threshold on both logFC and FDR_qvalue to 
##	select approximately 120 genes per cell type.
set.seed(1)
mark_genes = list()
num_mgene_per_ct = 120 # number of marker genes selected per cell type
for(ct in cell_types){
	# ct = cell_types[1]
	cat(paste0(ct,": "))
	rds_fn = file.path(psychENCODE_dir,paste0("ssd_nG",num_genes,"_cell",ct,".rds"))
	rds = readRDS(rds_fn)
	rds = rds[which(rds$logFC > logFC_thres 
					& rds$FDR_qvalue < fdr_thres
					& rds$gene %in% cts_genes[[ct]]),]
	# dim(rds)
	low_bound = 0
	while(TRUE){
		# qu = 0.25
		qu = runif(1,low_bound,1)
		num_rows = nrow(rds[which(rds$logFC > quantile(rds$logFC,qu) 
			& rds$FDR_qvalue < quantile(rds$FDR_qvalue,1-qu)),])
		if( num_rows > num_mgene_per_ct ){
			# Tighten bounds
			low_bound = qu # mean(c(low_bound[1],qu))
			cat(paste0(round(low_bound,3)," "))
			if( num_rows <= num_mgene_per_ct+10 ) break
		}
	}
	cat("\n")
	opt_qu = qu
	rds = rds[which(rds$logFC > quantile(rds$logFC,opt_qu) 
					& rds$FDR_qvalue < quantile(rds$FDR_qvalue,1-opt_qu)),]
	mark_genes[[ct]] = rds$gene; rm(rds)
}
sapply(mark_genes,length)

smart_pack("venn")
pdf("ct_genes_psychENCODE.pdf",width=8,height=8)
venn(x = ct_genes,ilabels = TRUE,zcolor = "style")
dev.off()

# Construct signature matrix
table(colData(sce)$cell_type)
SIG = matrix(NA,nrow=nrow(sce),ncol=length(cell_types))
colnames(SIG) = cell_types
rownames(SIG) = rownames(sce)
for(ct in cell_types){
	cat(paste0(ct,"\n"))
  # SIG[,ct] = rowSums(counts(sce)[, colData(sce)$cell_type == ct])
	SIG[,ct] = rowMeans(counts(sce)[, colData(sce)$cell_type == ct])
}
dim(SIG)
head(SIG)
summary(SIG)

# Subset marker genes
SIG = SIG[sort(as.character(unlist(mark_genes))),]
dim(SIG)
head(SIG)
summary(SIG)

# Subset and order genes
inter_genes = intersect(rownames(SIG),gene_anno$gene)
gene_anno = gene_anno[which(gene_anno$gene %in% inter_genes),]
SIG = SIG[which(rownames(SIG) %in% inter_genes),]
gene_anno = gene_anno[match(rownames(SIG),gene_anno$gene),]

# Output
saveRDS(list(anno = gene_anno,SIG = SIG),"signature_psychENCODE.rds")

# Import CMC Bulk RNA, get gene lengths
bulk = readRDS("CMC_MSSM-Penn-Pitt_Paul_geneExpressionRaw.rds")$so1
dim(bulk); bulk[1:5,1:5]
gtf_fn = "Homo_sapiens.GRCh37.70.processed.gtf"
exdb = GenomicFeatures::makeTxDbFromGFF(file = gtf_fn,format = "gtf")
exons_list_per_gene = GenomicFeatures::exonsBy(exdb,by = "gene")
tmp_df = smart_df(ensembl_gene_id = names(exons_list_per_gene),
	gene_length = as.numeric(sum(width(GenomicRanges::reduce(exons_list_per_gene)))))
tmp_df[1:5,]
all(tmp_df$ensembl_gene_id %in% rownames(bulk))
all(tmp_df$ensembl_gene_id == rownames(bulk))

# Subset/Order marker genes
inter_genes2 	= intersect(tmp_df$ensembl_gene_id,gene_anno$ensembl_gene_id)
gene_anno2 		= gene_anno[which(gene_anno$ensembl_gene_id %in% inter_genes2),]
bulk 					= bulk[which(rownames(bulk) %in% inter_genes2),]
bulk 					= bulk[match(gene_anno2$ensembl_gene_id,rownames(bulk)),]
tmp_df 				= tmp_df[which(tmp_df$ensembl_gene_id %in% inter_genes2),]
tmp_df 				= tmp_df[match(gene_anno2$ensembl_gene_id,tmp_df$ensembl_gene_id),]
saveRDS(list(anno = tmp_df,bulk = bulk),"bulk_marker_psychENCODE.rds")

# Refer to EPIC paper for variable definitions (bb = bulk,cc = signature,pp = cell type proportions)
bulk_rds 	= readRDS("bulk_marker_psychENCODE.rds")
sig_rds 	= readRDS("signature_psychENCODE.rds")
inter_genes3 	= intersect(bulk_rds$anno$ensembl_gene_id,sig_rds$anno$ensembl_gene_id)
bulk_rds$anno = bulk_rds$anno[which(bulk_rds$anno$ensembl_gene_id %in% inter_genes3),]
bulk_rds$bulk = bulk_rds$bulk[which(rownames(bulk_rds$bulk) %in% inter_genes3),]
sig_rds$anno 	= sig_rds$anno[which(sig_rds$anno$ensembl_gene_id %in% inter_genes3),]
sig_rds$SIG 	= sig_rds$SIG[which(rownames(sig_rds$SIG) %in% sig_rds$anno$gene),]
rownames(sig_rds$SIG) = sig_rds$anno$ensembl_gene_id

bb_tpm = apply(bulk_rds$bulk,2,function(xx) xx / bulk_rds$anno$gene_length)
bb_tpm = 1e6 * apply(bb_tpm,2,function(xx) xx / sum(xx))
cc_tpm = apply(sig_rds$SIG,2,function(xx) xx / sig_rds$anno$gene_length)
cc_tpm = 1e6 * apply(cc_tpm,2,function(xx) xx / sum(xx))

cell_sizes = colSums(apply(sig_rds$SIG,2,function(xx) xx / sig_rds$anno$gene_length))
cell_sizes

# Run ICeDT
pack = "ICeDT"
if( !(pack %in% installed.packages()[,"Package"]) ){
	devtools::install_github("Sun-lab/ICeDT")
}
date()
fit = ICeDT::ICeDT(Y = bb_tpm,Z = cc_tpm,tumorPurity = rep(0,ncol(bb_tpm)),
	refVar = NULL,rhoInit = NULL,maxIter_prop = 4e3,maxIter_PP = 4e3,
	rhoConverge = 1e-2)
date()
pp_bar_icedt = t(fit$rho)[,-1]
# calculate p_hat using cell_sizes
pp_hat_icedt = t(apply(pp_bar_icedt,1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
# Look at distribution of pp_bar and pp_hat
summary(pp_bar_icedt)
summary(pp_hat_icedt)

# Run CIBERSORT
sig_fn = file.path(MTG_dir,"signature_psychENCODE.txt")
mix_fn = file.path(MTG_dir,"mixture_CMC.txt")
write.table(cbind(rowname=rownames(cc_tpm),cc_tpm),
	file = sig_fn,sep = "\t",quote = FALSE,row.names = FALSE)
write.table(cbind(rowname=rownames(bb_tpm),bb_tpm),
	file = mix_fn,sep = "\t",quote = FALSE,row.names = FALSE)
cibersort_src_fn = "~/github/CSeQTL/R/CIBERSORT.R"
# cibersort_src_fn = "." # whichever directory contains CIBERSORT.R
source(cibersort_src_fn)
print(date())
results = CIBERSORT(sig_matrix = sig_fn,mixture_file = mix_fn,
	perm = 0,QN = FALSE,absolute = FALSE,abs_method = 'sig.score',
	filename = "psychENCODE")
print(date())
unlink(sig_fn)
unlink(mix_fn)
ciber_fn = sprintf("CIBERSORT-Results_%s.txt","psychENCODE")
unlink(ciber_fn)
QQ = ncol(cc_tpm)
pp_bar_ciber = results[,seq(QQ)]
# calculate p_hat using cell sizes
pp_hat_ciber = t(apply(pp_bar_ciber,1,function(xx){yy = xx / cell_sizes; yy / sum(yy)}))
summary(pp_bar_ciber)
summary(pp_hat_ciber)

# Output MTG deconvolution results
pdf("psychENCODE_bulk_deconvolution_summary.pdf",height=8,width=12)
	
	ICeDT_consistency(sig = cc_tpm,bulk = bb_tpm,ICeDT_out = fit)
	
	par(mfrow=c(1,2),mar=c(4,4,1,0.5),oma=c(0,0,2,0),bty="n")
	boxplot(pp_bar_icedt,main="Expression Proportions",ylim=c(0,1))
	boxplot(pp_hat_icedt,main="Cell type Proportions",ylim=c(0,1))
	mtext("ICeDT Results",outer=TRUE,cex=1.3)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
	
	par(mfrow=c(1,2),mar=c(4,4,1,0.5),oma=c(0,0,2,0),bty="n")
	boxplot(pp_bar_ciber,main="Expression Proportions",ylim=c(0,1))
	boxplot(pp_hat_ciber,main="Cell type Proportions",ylim=c(0,1))
	mtext("CIBERSORT Results",outer=TRUE,cex=1.3)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
	
	par(mfrow=c(2,3),mar=c(4,4,1,0.5),oma=c(0,0,2,0),bty="n")
	for(ct in seq(QQ)){
		v1 = pp_hat_ciber[,ct]
		v2 = pp_hat_icedt[,ct]
		tmp_range = range(c(v1,v2))
		pcor = round(cor(v1,v2),4)
		scor = round(cor(v1,v2,method="spear"),4)
		plot(v1,v2,xlim=tmp_range,ylim=tmp_range,
			xlab="CIBERSORT",ylab="ICeDT",
			main=sprintf("%s:Pear=%s;Spear=%s",colnames(pp_hat_ciber)[ct],
				pcor,scor))
		abline(a=0,b=1,lty=2,lwd=2,col="red")
	}
	mtext("CIBERSORT v. ICeDT",outer=TRUE,cex=1.4)
	par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,oma=rep(0,4),bty="o")
	
dev.off()

saveRDS(list(ICeDT = pp_hat_icedt,CIBERSORT = pp_hat_ciber),"prop_psychENCODE.rds")











# ------------------------------------------------------------
# Re-run deconvolution with intersection of MTG and psychENCODE signature genes
# ------------------------------------------------------------
if(FALSE){

fdr_thres = 1e-3
logFC_thres = log(2)

# Import MTG cts_genes (cts = cell type specific genes)
mtg_genes = readRDS("cts_genes.rds")
sapply(mtg_genes,length)

# Get psychENCODE cts_genes
DE_gene_anno = readRDS("../psychENCODE/deconvolution/DE_gene_anno.rds")
ct_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	print(ct)
	tmp_df = DE_gene_anno[which(DE_gene_anno[,paste0("logFC.",ct)] > logFC_thres 
					& DE_gene_anno[,paste0("FDR_qvalue.",ct)] < fdr_thres
					& DE_gene_anno$chromosome %in% c(1:22,"X","Y")
					# & !grepl("^LOC10",DE_gene_anno$gene)
					),]
	ct_genes[[ct]] = tmp_df$gene
}
sapply(ct_genes,length)

# smart_pack("venn")
# venn(x = ct_genes,ilabels = TRUE,zcolor = "style")
# Able to replicate the psychENCODE venndiagram

# Get the disjoint set of psychENCODE genes
psy_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	print(ct)
	psy_genes[[ct]] = sort(setdiff(ct_genes[[ct]],
		unique(unlist(ct_genes[which(names(ct_genes) != ct)]))))
}
sapply(psy_genes,length)

# Intersect MTG/psychENCODE genes
int_genes = list()
for(ct in cell_types){
	int_genes[[ct]] = intersect(mtg_genes[[ct]],psy_genes[[ct]])
}
sapply(int_genes,length)
saveRDS(int_genes,"int_MTG_psychENCODE_genes.rds")

# Select top 100 marker genes (using MTG's logFC and FDR values to filter)
set.seed(1)
num_genes = nrow(sce)
int_mark_genes = list()
for(ct in cell_types){
	# ct = cell_types[1]
	cat(paste0(ct,": "))
	rds_fn = paste0("ssd_nG",num_genes,"_cell",ct,".rds")
	rds = readRDS(rds_fn)
	rds = rds[which(rds$logFC > logFC_thres 
					& rds$FDR_qvalue < fdr_thres
					& rds$gene %in% int_genes[[ct]]),]
	# dim(rds)
	low_bound = 0
	while(TRUE){
		# qu = 0.25
		qu = runif(1,low_bound,1)
		num_rows = nrow(rds[which(rds$logFC > quantile(rds$logFC,qu) 
			& rds$FDR_qvalue < quantile(rds$FDR_qvalue,1-qu)),])
		if( num_rows > 100 ){
			# Tighten bounds
			low_bound = qu
			cat(paste0(round(low_bound,3)," "))
			if( num_rows <= 110 ) break
		}
	}
	cat("\n")
	opt_qu = qu
	rds = rds[which(rds$logFC > quantile(rds$logFC,opt_qu) 
					& rds$FDR_qvalue < quantile(rds$FDR_qvalue,1-opt_qu)),]
	int_mark_genes[[ct]] = rds$gene; rm(rds)
}
sapply(int_mark_genes,length)
saveRDS(int_mark_genes,"int_mark_genes.rds")

# Make intersected genes signature matrix
table(colData(sce)$cell_type)
sig_cts = matrix(nrow = nrow(sce),ncol = length(cell_types))
colnames(sig_cts) = cell_types
rownames(sig_cts) = rownames(sce)
for(ct in cell_types){
	print(ct)
	sig_cts[,ct] = rowSums(counts(sce)[,colData(sce)$cell_type == ct])
}
dim(sig_cts)
sig_cts = sig_cts/rowData(sce)$gene_length
sig_cts = 1e6 * apply(sig_cts,2,function(xx) xx/sum(xx))
sig_cts = sig_cts[sort(as.character(unlist(int_mark_genes))),]

# Get ENSG id
dat	= smart_df(rowData(sce)[,c("gene","chromosome","entrez_id")])
dat = dat[which(dat$gene %in% rownames(sig_cts)),]
rownames(dat) = NULL
dat = dat[match(rownames(sig_cts),dat$gene),]

biomaRt::listEnsembl()
ensembl = biomaRt::useMart("ensembl")
dd = biomaRt::listDatasets(ensembl)
	dim(dd); dd[1:5,]
	dd[grep("hsap",dd$dataset),]
ee = biomaRt::useDataset("hsapiens_gene_ensembl",mart = ensembl)
ff = biomaRt::listFilters(ee)
	dim(ff); ff[1:5,]
aa = biomaRt::listAttributes(ee)
	dim(aa); aa[1:5,]

attr_string = c("hgnc_symbol","ensembl_gene_id","external_gene_name",
				"entrezgene","chromosome_name","start_position",
				"end_position","strand","gene_biotype")
filters = c("hgnc_symbol","chromosome_name","entrezgene")
values = list(dat$gene,dat$chromosome,dat$entrez_id)
gene_anno = biomaRt::getBM(attributes = attr_string,
						filters = filters, 
						values = values,
						mart = ee)
gene_anno = unique(gene_anno)
dim(gene_anno); gene_anno[1:5,]
smart_table(rownames(sig_cts) %in% gene_anno$hgnc_symbol)

# What are the missing genes?
miss_genes = setdiff(rownames(sig_cts),gene_anno$hgnc_symbol)
miss_genes

# Output
saveRDS(list(anno = gene_anno,sig_cts = sig_cts),"int_signature.rds")

# Import CMC Bulk RNA, get gene lengths, calculate TPM
bulk = readRDS("CMC_MSSM-Penn-Pitt_Paul_geneExpressionRaw.rds")$so1
dim(bulk); bulk[1:5,1:5]

rds_fn = "Homo_sapiens.GRCh37.70.processed.rds"
if( !file.exists(rds_fn) ){
	gtf_fn = "Homo_sapiens.GRCh37.70.processed.gtf"
	exdb = GenomicFeatures::makeTxDbFromGFF(file = gtf_fn,format = "gtf")
	exons_list_per_gene = GenomicFeatures::exonsBy(exdb,by = "gene")
	rds = smart_df(ensembl_gene_id = names(exons_list_per_gene),
					gene_length = as.numeric(sum(width(GenomicRanges::reduce(exons_list_per_gene)))))
	saveRDS(rds,rds_fn)
}
rds = readRDS(rds_fn)

all(rds$ensembl_gene_id %in% rownames(bulk))
all(rds$ensembl_gene_id == rownames(bulk))
bulk = bulk / rds$gene_length
bulk = 1e6 * apply(bulk,2,function(xx) xx/sum(xx))
bulk[1:5,1:5] # tpm unit

# Import signature matrix and annotation
rds = readRDS("int_signature.rds")
gene_anno = rds$anno
sig_cts = rds$sig_cts
rm(rds)

# Get gene intersection, subset bulk and sig_cts, matching gene order
smart_table(rownames(bulk) %in% gene_anno$ensembl_gene_id)
inter_genes = intersect(rownames(bulk),gene_anno$ensembl_gene_id)
gene_anno = gene_anno[which(gene_anno$ensembl_gene_id %in% inter_genes),]

bulk = bulk[which(rownames(bulk) %in% inter_genes),]
bulk = bulk[match(gene_anno$ensembl_gene_id,rownames(bulk)),]

sig_cts = sig_cts[which(rownames(sig_cts) %in% gene_anno$hgnc_symbol),]
sig_cts = sig_cts[match(gene_anno$hgnc_symbol,rownames(sig_cts)),]

all(rownames(sig_cts) == gene_anno$hgnc_symbol)
all(rownames(bulk) == gene_anno$ensembl_gene_id)

# Rename rownames of sig_cts to match bulk
rownames(sig_cts) = gene_anno$ensembl_gene_id

# Run ICeDT
pack = "ICeDT"
if( !(pack %in% installed.packages()[,"Package"]) ){
	devtools::install_github("Sun-lab/ICeDT")
}
icedt_fn = "int_icedt.rds"
date()
fitw0 = ICeDT::ICeDT(Y = bulk,Z = sig_cts,tumorPurity = rep(0,ncol(bulk)),
	refVar = NULL,rhoInit = NULL,maxIter_prop = 500,maxIter_PP = 250,
	rhoConverge = 1e-2)
saveRDS(fitw0,icedt_fn)
date()
fitw0 = readRDS(icedt_fn)
p1 = fitw0$cProb
prop_icedt = t(fitw0$rho)[,-1]
dim(p1)
p1[1:2,1:5]
p1 = data.matrix(p1)
q90 <- function(v){
	qs = quantile(v, probs=c(0.10, 0.90))
	qs[2] - qs[1]
}

pdf("probConsistent_int_GeneSet.pdf",width=11,height=4)

par(mar=c(5,4,1,1),bty="n",mfrow=c(1,3),cex=0.8)
plot(density(c(p1))$y,main="",xlim=c(0,1),xlab="probability consistent",
	ylab="density",type="n")
lines(density(c(p1)), lty=1, col="black")
legend("topright", c("no weight"), lty=c(1,2),col=c("black"), bty="n")
plot(apply(p1,1,median),apply(p1,1,q90),xlab="median prob. consistent",
	ylab="90 percentile - 10 percentile",main="Gene Consistency across subjects")
boxplot(prop_icedt,main = "MTG - ICeDT")
par(mfrow=c(1,1))

# Predicted vs Observed gene expression
predicted_bulk_w0 = sig_cts %*% t(prop_icedt)
p1_cutoffs = quantile(p1,c(1/3,2/3)); p1_cutoffs
plot_log1p = function(x, y, ...) {
	smoothScatter(log(x+1e-5), log(y+1e-5), xlim=c(-5, 10), ylim=c(-5, 10), ...)
	legend("bottomright",bty="n",legend=sprintf("Pearson correlation = %.2f",
		cor(log(x+1e-5), log(y+1e-5))))
}
par(mar=c(5,4,1,1),bty="n",mfrow=c(1,3),cex=0.6)
plot_log1p(x = c(predicted_bulk_w0)[p1 < p1_cutoffs[1]],
		y = c(bulk)[p1 < p1_cutoffs[1]], 
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "low prob of being consistent")
plot_log1p(x = c(predicted_bulk_w0)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]],
		y = c(bulk)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]], 
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "med prob of being consistent")
plot_log1p(x = c(predicted_bulk_w0)[p1 > p1_cutoffs[2]],
		y = c(bulk)[p1 > p1_cutoffs[2]],
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "high prob of being consistent")

dev.off()

# Run CIBERSORT
write.table(cbind(rowname=rownames(sig_cts),sig_cts),
			file = "signature_INT.txt",
			sep = "\t",quote = FALSE,row.names = FALSE)
write.table(cbind(rowname=rownames(bulk),bulk),
			file = "mixture_CMC_INT.txt",
			sep = "\t",quote = FALSE,row.names = FALSE)
# Login to CIBERSORT website, uploaded above two files, 
#	specified no quantile normalization, ran 1000 permutations
cib = smart_RT(file.path(MTG_dir,"CIBERSORT.Output_Job3.txt"),
			sep = "\t",header = TRUE)
dim(cib); cib[1:5,]
prop_cib = as.matrix(cib[,2:7])

# Compare proportions
pdf("compare_props_int.pdf",width=9,height=5)

# Proportions across cell types
par(mfrow=c(1,2),mar=c(5,4,1,1))
boxplot(prop_icedt,main = "MTG - ICeDT")
boxplot(prop_cib,main = "MTG - CIBERSORT")

# Proportions by cell type
cell_types = colnames(prop_cib)
par(mar=c(5,4,1,1),bty="n",mfrow=c(2,3),cex=0.8)
for(ct in cell_types){
	tmp_range = range(c(prop_icedt[,ct],prop_cib[,ct]))
	plot(prop_icedt[,ct],prop_cib[,ct],xlab="prop est in ICeDT",
		ylab="prop est in CIBERSORT",pch=16,cex=1.1,cex.lab=1.1,
		col=rgb(0,0,0,0.25),xlim=tmp_range,ylim=tmp_range,main=ct)
	abline(a=0,b=1,lty=2,col="blue",lwd=2)
}

saveRDS(list(ICeDT = prop_icedt,CIBERSORT = prop_cib),"prop_int.rds")

dev.off()

}

sessionInfo()

q("no")

###
