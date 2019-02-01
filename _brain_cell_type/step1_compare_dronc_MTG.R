
library(org.Hs.eg.db)
library(reshape2)
library(ggplot2)

repo_dir  = "~/research/GitHub/scRNAseq_pipelines"
setwd(repo_dir)

# ------------------------------------------------------------------
# load data matrix 
# ------------------------------------------------------------------

dronc_n_cell = readRDS("dronc/ct_cells_dronc.rds")
dronc_g_info = readRDS("dronc/gene_info_dronc.rds")
dronc_ct_mat = readRDS("dronc/ct_matrix_dronc.rds")

dronc_n_cell
dim(dronc_g_info)
dronc_g_info[1:2,]

names(dronc_ct_mat)
dronc_ct_mat$PFC[1:2,1:3]

MTG_n_cell = readRDS("MTG/ct_cells_human_MTG.rds")
MTG_g_info = readRDS("MTG/gene_info_human_MTG.rds")
MTG_ct_mat = readRDS("MTG/ct_matrix_human_MTG.rds")

MTG_n_cell
dim(MTG_g_info)
MTG_g_info[1:2,]

# ------------------------------------------------------------------
# check gene information 
# ------------------------------------------------------------------

length(unique(dronc_g_info$ensembl_gene_id))
table(is.na(dronc_g_info$ensembl_gene_id))
table(is.na(MTG_g_info$ensembl_gene_id))
table(is.na(MTG_g_info$entrez_id))

columns(org.Hs.eg.db)
map1 = mapIds(org.Hs.eg.db, keys=as.character(dronc_g_info$ensembl_gene_id), 
              'ENTREZID', 'ENSEMBL')
length(map1)
map1[1:5]
table(is.na(map1))
dronc_g_info$entrez_id = as.character(map1)

length(intersect(MTG_g_info$ensembl_gene_id, dronc_g_info$ensembl_gene_id))
length(intersect(MTG_g_info$entrez_id, dronc_g_info$entrez_id))

table(MTG_g_info$pct_dropout_counts<90,!is.na(MTG_g_info$ensembl_gene_id))
table(MTG_g_info$pct_dropout_counts<80,!is.na(MTG_g_info$ensembl_gene_id))

wNA = which(is.na(MTG_g_info$ensembl_gene_id))

summary(MTG_g_info$pct_dropout_counts)
log10_p_dropout = log10(1 - MTG_g_info$pct_dropout_counts/100)

pdf("_brain_cell_type/figure/MTG_with_ensembl_id_or_not.pdf", width=7, 
    height=3.5)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
plot(density(MTG_g_info$log10_mean_counts[wNA]), main="", 
     xlab="log10(mean counts)")
lines(density(MTG_g_info$log10_mean_counts[-wNA]), 
      lty=2, col="red")
text(3.5, 5, labels="Ensembl ID", pos=2)
legend("topright", legend=c("", "No", "Yes"), 
       lty=c(1, 1,2), col=c("white", "black", "red"), bty="n")
plot(density(log10_p_dropout[wNA]), main="", ylim=c(0,0.9), xaxt="n", 
     xlab="percent dropout")
lines(density(log10_p_dropout[-wNA]), lty=2, col="red")
pts = c(0, 70, 90, 99, 99.9)
axis(side=1, at=log10(1-pts/100), labels=pts)
dev.off()

# ------------------------------------------------------------------
# take intersection 
# ------------------------------------------------------------------

ens_id = intersect(MTG_g_info$ensembl_gene_id, dronc_g_info$ensembl_gene_id)
dronc_all = dronc_ct_mat$all[match(ens_id, dronc_g_info$ensembl_gene_id),]
dronc_pfc = dronc_ct_mat$PFC[match(ens_id, dronc_g_info$ensembl_gene_id),]
MTG = MTG_ct_mat[match(ens_id, MTG_g_info$ensembl_gene_id),]

dim(MTG)
dim(dronc_all)
dim(dronc_pfc)

dronc_all[1:2,]
dronc_pfc[1:2,]
MTG[1:2,]

col_names = c("GABA", "exPFC", "ODC", "OPC", "ASC", "MG", "END")
dronc_all = dronc_all[,col_names]
dronc_pfc = dronc_pfc[,col_names]

col_names = c("GABA", "Gluta", "ODC", "OPC", "ASC", "MG", "END")
colnames(MTG) = colnames(dronc_all) = colnames(dronc_pfc) = col_names

# ------------------------------------------------------------------
# comparison for each cell type 
# ------------------------------------------------------------------

for(i in 1:length(col_names)){
  ci = col_names[i]
  png(sprintf("_brain_cell_type/figure/ct_%s.png", ci), 
      width=9, height=6, units="in", res=300)
  
  x = log10(dronc_all[,ci]+1)
  y = log10(dronc_pfc[,ci]+1)
  z = log10(MTG[,ci]+1)
  
  par(mfrow=c(2,3), mar=c(5,4,2,1), bty="n")
  hist(x, xlab="log10(count+1)", main="dronc all", breaks=50)
  hist(y, xlab="log10(count+1)", main="dronc PFC", breaks=50)
  hist(z, xlab="log10(count+1)", main="MTG", breaks=50)
  plot(x, y, cex=0.2, xlab="dronc all", ylab="dronc PFC", 
       main=round(cor(x, y, method="spearman"),2))
  plot(x, z, cex=0.2, xlab="dronc all", ylab="MTG", 
       main=round(cor(x, z, method="spearman"),2))
  plot(y, z, cex=0.2, xlab="dronc all", ylab="MTG", 
       main=round(cor(y, z, method="spearman"),2))
  dev.off()
}

# ------------------------------------------------------------------
# comparison across cell types 
# ------------------------------------------------------------------

png("_brain_cell_type/figure/pairs_MTG.png", 
    width=18, height=18, units="in", res=300)
pairs(log10(MTG+1))
dev.off()

png("_brain_cell_type/figure/pairs_dronc_all.png", 
    width=18, height=18, units="in", res=300)
pairs(log10(dronc_all+1))
dev.off()

png("_brain_cell_type/figure/pairs_dronc_pfc.png", 
    width=18, height=18, units="in", res=300)
pairs(log10(dronc_pfc+1))
dev.off()

# ------------------------------------------------------------------
# correlations 
# ------------------------------------------------------------------

cr.MTG = cor(log10(MTG+1), method="spearman")
cr.dronc_all = cor(log10(dronc_all+1), method="spearman")
cr.dronc_pfc = cor(log10(dronc_pfc+1), method="spearman")

diag(cr.MTG) = diag(cor(log10(MTG+1), log10(dronc_all+1), 
                        method="spearman"))
diag(cr.dronc_all) = diag(cor(log10(MTG+1), log10(dronc_all+1), 
                              method="spearman"))
diag(cr.dronc_pfc) = diag(cor(log10(MTG+1), log10(dronc_pfc+1), 
                              method="spearman"))

gplot1 <- function(melted_cormat){
  g1 = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + labs(x="",y="") + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0.5, limit = c(0,1)) + 
    geom_text(aes(label = round(value, 2)))
  g1
}

g1 = gplot1(melt(cr.MTG))
g2 = gplot1(melt(cr.dronc_all))
g3 = gplot1(melt(cr.dronc_pfc))

ggsave("_brain_cell_type/figure/cr_MTG.png", g1, width=5, 
       height=4, units="in")

ggsave("_brain_cell_type/figure/cr_dronc_all.png", g2, width=5, 
       height=4, units="in")

ggsave("_brain_cell_type/figure/cr_dronc_pfc.png", g3, width=5, 
       height=4, units="in")

if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

sessionInfo()
q(save="no")


