# scRNAseq differential expression (seems to require at least 80 gbs)

# ----------
# Shortcuts
# ----------
rm(list=ls())
repo_dir = "/pine/scr/p/l/pllittle/CS_eQTL/s3_Real/scRNAseq_pipelines"
MTG_dir = file.path(repo_dir,"MTG")
setwd(MTG_dir)


# ----------
# Libraries
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


# ----------
# Analyze
# ----------
rds = readRDS(file.path(MTG_dir,"old_new_clusters.rds"))
	names(rds)
	opt_clust = rds$opt_clust
	clusts = rds$clusts; clusts
	t1 = rds$t1; t1
sce = readRDS(file.path(MTG_dir,"final_sce_filtered_by_kmeans.rds"))
	sce
	sort(names(colData(sce)))
	smart_table(colData(sce)$donor)
if(FALSE){
clusters = readRDS(file.path(MTG_dir,"final_hvg_clust.rds"))
	clusters[1:4,]
	length(clusters$sample_name)
	length(colData(sce)$sample_name)
	clusters = clusters[match(colData(sce)$sample_name,clusters$sample_name),]
	table(clusters$sample_name == colData(sce)$sample_name)
	table(clusters[,c(opt_clust,"cell_type")])
	colData(sce)[[opt_clust]] = clusters[,opt_clust]
}

# zlm analysis
# sca = as(sce,"SingleCellAssay")
sca = MAST::SceToSingleCellAssay(sce = sce)
zlm_output = MAST::zlm(formula = ~ cell_type,sca = sca)
show(zlm_output)

coefAndCI = summary(zlm_output,logFC = FALSE)$datatable
coefAndCI
dim(coefAndCI)
smart_table(coefAndCI$component) # C = continuous betas, D = discrete betas
smart_table(coefAndCI$contrast)
coefAndCI = coefAndCI[contrast != '(Intercept)',]
dim(coefAndCI)

ggplot(coefAndCI,aes(x = contrast,y = coef,ymin = ci.lo,ymax = ci.hi,col = component)) +
	geom_pointrange(position = position_dodge(width=.5)) +
	facet_wrap(~ primerid) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
	coord_cartesian(ylim = 12*c(-1,1))
dev.off()






# zlm Practice
sce2 = sce[1:20,] # first 20 genes
sca2 = as(sce2,"SingleCellAssay")
zlm_output = MAST::zlm(formula = ~ donor,sca = sca2)
show(zlm_output)

coefAndCI = summary(zlm_output,logFC = FALSE)$datatable
coefAndCI
dim(coefAndCI)
smart_table(coefAndCI$component) # C = continuous betas, D = discrete betas
smart_table(coefAndCI$contrast)
coefAndCI = coefAndCI[contrast != '(Intercept)',]
dim(coefAndCI)
# coefAndCI[,contrast:=abbreviate(contrast)]

ggplot(coefAndCI,aes(x = contrast,y = coef,ymin = ci.lo,ymax = ci.hi,col = component)) +
	geom_pointrange(position = position_dodge(width=.5)) +
	facet_wrap(~ primerid) + theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
	coord_cartesian(ylim = c(-3,3))

zlm_lr = MAST::lrTest(zlm_output, 'donor')
dimnames(zlm_lr)

pvalue = ggplot(melt(zlm.lr[,,'Pr(>Chisq)']),aes(x=primerid, y=-log10(value)))+
    geom_bar(stat='identity')+facet_wrap(~test.type) + coord_flip()
print(pvalue)


q("no")

###
