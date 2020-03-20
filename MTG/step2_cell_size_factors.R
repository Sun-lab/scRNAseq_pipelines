
gtf_dir = "GTF2" # wherever the gtf file is located
git_dir = "." # wherever the scRNAseq_pipelines repo is located

gtf_url = "https://personal.broadinstitute.org/francois/topmed/gencode.v26.GRCh38.ERCC.genes.gtf.gz"
gtf_fn = file.path(gtf_dir,"gene_info_v26.rds") # this is from gunzipping the above url file
gtf = readRDS(gtf_fn)

# Calculate gene length
cat(sprintf("%s: Calculate gene length ...\n",date()))
exon_fn = file.path(gtf_dir,"exon_by_genes_gencode_v26.rds")

if( !file.exists(exon_fn) ){
  cat("Make TxDb...\n")
  exdb = suppressWarnings(GenomicFeatures::makeTxDbFromGFF(file = gtf_fn,format = "gtf"))
  exons_list_per_gene = GenomicFeatures::exonsBy(exdb,by = "gene")
  saveRDS(exons_list_per_gene,exon_fn)
}

exon = readRDS(exon_fn)
all(names(exon) == gtf$ENSG)

gtf = gtf[match(names(exon),gtf$ENSG),]
all(names(exon) == gtf$ENSG)
gtf$gene_length = as.numeric(sum(BiocGenerics::width(GenomicRanges::reduce(exon))))
rm(exon)

# Prepare annotation
gtf = gtf[which(gtf$Chr %in% paste0("chr",1:22)),]
gtf$ChrNum = as.integer(gsub("chr","",gtf$Chr))
gtf = gtf[order(gtf$ChrNum,gtf$Start),]

# Import Signature Matrices and Marker genes: 
cat(sprintf("%s: Import signature matrices and marker genes ...\n",date()))
rds = readRDS("all_genes_MTG.rds")
names(rds)
sig = rds$SIG; anno = rds$anno; rm(rds)
dim(sig)
sig[1:5,]
colSums(sig)

rds = readRDS("signature_MTG.rds")
mark = rownames(rds$SIG); rm(rds)

# Intersect sig genes with gtf genes for cell size calculation after gene length adjustment
cat(sprintf("%s: Int sig & gtf genes ...\n",date()))
gene_var = ifelse(dataset %in% c("EGA","CMC","gtexBlood","gtexBrain"),"GENE","ENSG")

## Removing genes with multiple ids
tt = table(gtf[[gene_var]]); mult_gene = which(tt > 1)
if( length(mult_gene) > 0 ) gtf = gtf[which(!(gtf[[gene_var]] %in% names(mult_gene))),]

int_gene_gtf_sig = intersect(gtf[[gene_var]],rownames(sig)); length(int_gene_gtf_sig)
sig = sig[rownames(sig) %in% int_gene_gtf_sig,]
gtf = gtf[which(gtf[[gene_var]] %in% int_gene_gtf_sig),]
sig = sig[match(gtf[[gene_var]],rownames(sig)),]
all(rownames(sig) == gtf[[gene_var]])

# Calculate cell size from colSums() from available genes
cat(sprintf("%s: Calculate cell sizes ...\n",date()))
cell_sizes = apply(sig,2,function(xx) sum(xx / gtf$gene_length))
cell_sizes

sessionInfo()
q(save="no")

