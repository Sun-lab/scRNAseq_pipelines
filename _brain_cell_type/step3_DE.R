# scRNAseq differential expression

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

library(ggplot2); library(data.table)
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


# ----------
# Import and Prep Data
# ----------
# Import gtf with gene lengths
rsem_fn = "rsem_GRCh38.p2.gtf"
gtf = GTF_calc_gene_length(work_dir = MTG_dir,rsem_gtf_fn = rsem_fn)

sce = readRDS(file.path(MTG_dir,"final_sce_filtered_by_kmeans.rds"))
sce

# Append gene_lengths to sce
inter_genes = intersect(rownames(sce),gtf$gene_symbol)
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


# ----------
# Outlined Steps
# ----------
# Run MAST to get differentially expressed genes
# Get disjoint set of marker genes per cell type (filter on logFC and qvalue), 
#    subset approximately 100 genes per cell type
# Use sce counts and exonic gene lengths to calculate TPM counts
# Bulk RNAseq deconvolution: Run CIBERSORT, ICEDT


# ----------
# Run MAST: Run each cell type as it's own job
# ----------
num_genes = nrow(sce)
for(ct in cell_types){
	print(ct)
	MAST_DEgenes(work_dir = MTG_dir,
				num_genes = num_genes,
				sce_obj = sce,
				one_cell_type = one_ct)
}


# ----------
# Combine all genes with annotations and proportion of expressed cells per cell type
# ----------
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


# ----------
# Get marker genes per cell type
# ----------
dat = smart_df(rowData(sce)[,c("gene","chromosome","entrez_id")])
fdr_thres = 1e-3
logFC_thres = log(2)

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

# Get approximate top 100 marker genes per cell type
## Rather than arbitrarily selecting 100 genes, I apply 
##	a quantile threshold on both logFC and FDR_qvalue to 
##	select approximately 100 genes per cell type.
mark_genes = list()
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
		if( num_rows > 100 ){
			# Tighten bounds
			low_bound = qu # mean(c(low_bound[1],qu))
			cat(paste0(round(low_bound,3)," "))
			if( num_rows <= 110 ) break
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
pdf("ct_genes.pdf",width=8,height=8)
venn(x = ct_genes,ilabels = TRUE,zcolor = "style")
dev.off()


# ----------
# Calculate TPM and get signature matrix: From Chong Jin's code
# ----------
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
sig_cts = 1e6 * t(t(sig_cts)/colSums(sig_cts))
sig_cts = sig_cts[sort(as.character(unlist(mark_genes))),]

# Get ENSG id
dat	= smart_df(rowData(sce)[,c("gene","chromosome","entrez_id")])
dat = dat[which(dat$gene %in% rownames(sig_cts)),]
rownames(dat) = NULL
dat = dat[match(rownames(sig_cts),dat$gene),]
  
biomaRt::listEnsembl()
ensembl = biomaRt::useMart("ensembl")
dd = biomaRt::listDatasets(ensembl); dim(dd); dd[1:5,]
dd[grep("hsap",dd$dataset),]
ee = biomaRt::useDataset("hsapiens_gene_ensembl",mart = ensembl)
ff = biomaRt::listFilters(ee); dim(ff); ff[1:5,]
aa = biomaRt::listAttributes(ee); dim(aa); aa[1:5,]

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
saveRDS(list(anno = gene_anno,sig_cts = sig_cts),"signature.rds")


# ----------
# Deconvolution: Based on Chong Jin's deconvolution.Rmd code
# ----------
# Import CMC Bulk RNA, get gene lengths, calculate TPM
bulk = readRDS("CMC_MSSM-Penn-Pitt_Paul_geneExpressionRaw.rds")$so1
dim(bulk); bulk[1:5,1:5]
  
gtf_fn = "Homo_sapiens.GRCh37.70.processed.gtf"
exdb = GenomicFeatures::makeTxDbFromGFF(file = gtf_fn,format = "gtf")
exons_list_per_gene = GenomicFeatures::exonsBy(exdb,by = "gene")
tmp_df = smart_df(ensembl_gene_id = names(exons_list_per_gene),
				gene_length = as.numeric(sum(width(GenomicRanges::reduce(exons_list_per_gene)))))
  
all(tmp_df$ensembl_gene_id %in% rownames(bulk))
all(tmp_df$ensembl_gene_id == rownames(bulk))
bulk = bulk / tmp_df$gene_length
bulk = 1e6 * t(t(bulk/colSums(bulk)))
bulk[1:5,1:5] # tpm unit
  
# Import signature matrix and annotation
rds = readRDS("signature.rds")
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
icedt_fn = "icedt.rds"
date()
if( !file.exists(icedt_fn) ){
	fitw0 = ICeDT::ICeDT(Y = bulk,Z = sig_cts,tumorPurity = rep(0,ncol(bulk)),
		refVar = NULL,rhoInit = NULL,maxIter_prop = 500,maxIter_PP = 250,
		rhoConverge = 1e-2)
	saveRDS(fitw0,icedt_fn)
}
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

pdf("probConsistent_GeneSet.pdf",width=9,height=4)

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
		y = bulk[p1 < p1_cutoffs[1]], 
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "low prob of being consistent")
plot_log1p(x = c(predicted_bulk_w0)[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]],
		y = bulk[p1 >= p1_cutoffs[1] & p1 <= p1_cutoffs[2]], 
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "med prob of being consistent")
plot_log1p(x = c(predicted_bulk_w0)[p1 > p1_cutoffs[2]],
		y = bulk[p1 > p1_cutoffs[2]],
		xlab = "Predicted gene expression",ylab = "Observed gene expression",
		sub = "model w/ weight",main = "high prob of being consistent")

dev.off()
  
# Run CIBERSORT
write.table(cbind(rowname=rownames(sig_cts),sig_cts),
			file = file.path(MTG_dir,"signature_MTG.txt"),
			sep = "\t",quote = FALSE,row.names = FALSE)
write.table(cbind(rowname=rownames(bulk),bulk),
			file = file.path(MTG_dir,"mixture_CMC.txt"),
			sep = "\t",quote = FALSE,row.names = FALSE)
# Login to CIBERSORT website, uploaded above two files, 
#	specified no quantile normalization, ran 1000 permutations
cib = smart_RT(file.path(MTG_dir,"CIBERSORT.Output_Job2.txt"),
			sep = "\t",header = TRUE)
dim(cib); cib[1:5,]
prop_cib = as.matrix(cib[,2:7])

# Compare proportions
pdf("compare_props.pdf",width=9,height=5)

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

saveRDS(list(ICeDT = prop_icedt,CIBERSORT = prop_cib),"prop.rds")

dev.off()

sessionInfo()

q("no")

###
