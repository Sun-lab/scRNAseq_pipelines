# 80GB of memory recommended.
library(data.table)
library(SingleCellExperiment)
# library(DropletUtils)
# library(biomaRt)
# library(scater)
# library(scran)
library(limma)
library(ggplot2)


work_dir = psychENCODE_dir = "~/scr/psychENCODE_data"
sce = readRDS(file.path(psychENCODE_dir,"sce.rds"))
 
# calculate TPM
tpm(sce) = counts(sce)/rowData(sce)$gene_length
tpm(sce) = (1e6)*t(t(tpm(sce))/colSums(tpm(sce)))
logcounts(sce) = log2(1 + tpm(sce))

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

clusters        =  readRDS(file.path(psychENCODE_dir, "k25_50_pcs.rds"))

dim(sce)
dim(colData(sce))
colData(sce)[1:2,1:5]

names(colData(sce))
table(colData(sce)$donor)

dim(rowData(sce))
rowData(sce)[1:2,]

clusters = data.frame(KM25=clusters$cluster)
clusters$cell_type = colData(sce)$cell_type

opt_clust = paste0("KM",25) # our clustering results differ slightly
t1 = table(clusters[,opt_clust], clusters$cell_type)
t1

# remove cluster 15 & 23 since they belong to multiple clusters
t1[c(15,23), ] = 0

clusts = apply(t1, 2, function(v){union(which.max(v), which(v > 100))})
clusts

celltypes = setdiff(unique(clusters$cell_type), "unknown")
celltypes

w2kp = NULL

for(ct1 in celltypes){
  # ct1 = celltypes[1]
  ct.cond    = clusters$cell_type == ct1
  clust.cond = clusters[,opt_clust] %in% clusts[[ct1]]
  # Intersect previous clustering result with ours
  w2kp = c(w2kp, which(ct.cond & clust.cond))
}
length(w2kp)

dim(sce)
sce = sce[,w2kp]
dim(sce)

name_change = function(dat,orig_name,new_name){
  index = which(names(dat) == orig_name)
  if( length(index) > 0 ){
    names(dat)[index] = new_name
  }
  dat
}


MAST_DEgenes = function(work_dir,num_genes=NULL,sce_obj,one_cell_type,fdr_thres=1e-3,logFC_thres=log(2)){
  
  if(FALSE){
    work_dir = psychENCODE_dir
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
  colData(sca)$cngeneson = scale(colSums(logcounts(sca) > 0))
  
  # Calculate Group
  # one_cell_type = "Astro"
  ref_group = paste0("not_",one_cell_type)
  colData(sca)$Group = ifelse(colData(sca)$cell_type == one_cell_type,
                              one_cell_type,ref_group)
  colData(sca)$Group = factor(colData(sca)$Group,
                              levels = c(ref_group,one_cell_type))
  
  # zlm analysis: Model log-transformed expression as 
  #     function of clustered cell type and num detected genes
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
  ssc_fn = file.path(work_dir,paste0("ssc",
                                     "_nG",num_genes,
                                     "_cell",one_cell_type,
                                     ".rds"))
  if ( !file.exists(ssc_fn) ){
    print(date())
    ssc = MAST::summary(object = zlm_output,
                        doLRT = paste0("Group",one_cell_type))
    print(date())
    saveRDS(ssc, ssc_fn)
  }
  ssc = readRDS(ssc_fn)
  ssd0 = data.frame(ssc$datatable, stringsAsFactors = FALSE)
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
  # ssd = ssd[which(ssd$FDR_qvalue < fdr_thres 
  #                 & ssd$logFC > logFC_thres),]
  # dim(ssd); smart_table(ssd$logFC > 0); ssd[1:20,]
  ssd = data.frame(cell_type=one_cell_type,ssd,stringsAsFactors = FALSE)
  
  ssd_fn = file.path(work_dir,paste0("ssd",
                                     "_nG",num_genes,
                                     "_cell",one_cell_type,
                                     ".rds"))
  cat(paste0("Saving image in ",ssd_fn," ...\n"))
  saveRDS(ssd,ssd_fn)
  
  return(NULL)
}

# Subset cell_types with enough cells: Exclude Endo
table(colData(sce)$cell_type, useNA="ifany")
sce = sce[,!(colData(sce)$cell_type %in% c("Per","Endo"))]
table(colData(sce)$cell_type, useNA="ifany")
dim(sce)

# Get cell type specific marker genes
# one_ct = "blah_one_ct"
cts = c("Astro", "Exc", "Inh", "Micro", "Oligo", "OPC")
# one_ct = "Astro"
# one_ct = "Exc"
# one_ct = "Inh"
# one_ct = "Micro"
# one_ct = "Oligo"
# one_ct = "OPC"
for (one_ct in cts) {
  MAST_DEgenes(work_dir = psychENCODE_dir,
               num_genes = nrow(sce),
               # num_genes = 100,
               sce_obj = sce,
               one_cell_type = one_ct,
               fdr_thres = 1e-2,
               logFC_thres = log(2))
  gc()
}



