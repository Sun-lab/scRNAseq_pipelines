
library(Matrix)
library(data.table)
library(dplyr)
library(DESeq2)
library(svd)

data.dir = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data"


args=(commandArgs(TRUE))
args

if (length(args) != 1) {
  message("one argument is expected, use 'PFC_L2_3' as default.\n")
  grp = "PFC_L2_3"
}else{
  eval(parse(text=args[[1]]))
}

grp

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = read.table(paste0("/fh/fast/sun_w/mengqi/Data_PRJNA434002/meta.tsv"),header = TRUE, sep = "\t")
dim(cell_info)
cell_info[1:2,]

table(cell_info$region)

cell_info = cell_info[which(cell_info$region=="PFC"),]

sort(table(paste(cell_info$diagnosis, cell_info$sample, sep=":")))

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

dat1 = readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data/ct_mtx/",grp,".rds"))
dim(dat1)
class(dat1)
dat1[1:5,1:4]

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat1) %in% cell_info$cell)
meta = cell_info[match(colnames(dat1), cell_info$cell),]
dim(meta)
meta[1:2,]

summary(meta$UMIs/meta$genes)

# check each individual has a unique sample
table(tapply(meta$sample, meta$individual, function(v){length(unique(v))}))

# check each individual has a unique Capbatch
table(tapply(meta$Capbatch, meta$individual, function(v){length(unique(v))}))

table(meta$cluster)
table(meta$region)

sort(table(paste(meta$diagnosis, meta$sample, sep=":")))

# ------------------------------------------------------------------------
# generate individual level information
# ------------------------------------------------------------------------

meta_ind = distinct(meta[,3:12])
dim(meta_ind)
meta_ind[1:2,]
names(meta_ind)[9:10] = c("PMI", "RIN")

length(unique(meta$individual))

if(nrow(meta_ind) != length(unique(meta$individual))){
  stop("there is non-unique information\n")
}

table(meta_ind$Seqbatch, meta_ind$Capbatch)

# ------------------------------------------------------------------------
# collect count data
# ------------------------------------------------------------------------

trec1 = matrix(NA, nrow=nrow(dat1), ncol=nrow(meta_ind))
colnames(trec1) = meta_ind$sample
rownames(trec1) = rownames(dat1)
dim(trec1)
trec1[1:2,1:3]

for(i in 1:ncol(trec1)){
  wi = which(meta$sample == meta_ind$sample[i])
  trec1[,i] = rowSums(dat1[,wi])
}

dim(trec1)
trec1[1:2,1:3]

summary(apply(trec1, 1, median))
q75 = apply(trec1, 1, quantile, probs=0.75)
summary(q75)
table(q75 >= 20)

# ------------------------------------------------------------------------
# run DESeq2
# ------------------------------------------------------------------------

colData = meta_ind
for(i in 1:ncol(colData)){
  if(is.character(colData[[i]])){
    colData[[i]] = as.factor(colData[[i]])
  }
}
dim(colData)
colData[1:2,]
summary(colData)

colData$diagnosis = factor(colData$diagnosis, levels=c("Control", "ASD"))

dd0 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ diagnosis)
dd0  = DESeq(dd0)
res0 = results(dd0)
dim(res0)
res0[1:2,]

dd1 = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ age + sex + Capbatch + PMI + 
                               RIN + diagnosis)
dd1 = DESeq(dd1)

res1 = results(dd1)
dim(res1)
res1[1:2,]

summary(res0)
summary(res1)

# ------------------------------------------------------------------------
# summarize p-value distribution
# ------------------------------------------------------------------------

pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step1_DESeq2_",grp,"_compare_pval_org.pdf"),width=9, height=3)

par(mfrow=c(1,3), bty="n", mar=c(5,4,1,1))
hist(res0$pvalue, main="without covariates", xlab="p-value")
hist(res1$pvalue, main="with covariates", xlab="p-value")
plot(-log10(res0$pvalue), -log10(res1$pvalue), main="-log10(p-value)", 
     xlab="without covariates", ylab="with covariates", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")
dev.off()

nms = resultsNames(dd1)
nms
nms = nms[-1]

pvals = matrix(NA, nrow=nrow(trec1), ncol=length(nms))

for(k in 1:length(nms)){
  rk = results(dd1, name=nms[k])
  pvals[,k] = rk$pvalue
}

colnames(pvals) = nms
dim(pvals)
summary(pvals)

pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step1_DESeq2_",grp,"_pval_hist_org.pdf"), width=9, height=9)
par(mfrow=c(3,3), bty="n", mar=c(5,4,2,1))
for(k in 1:length(nms)){
  hist(pvals[,k], main=nms[k], xlab="p-value", breaks=50)
}

rbatch = results(dd1, contrast=c("Capbatch","CB6","CB7"))
dim(rbatch)
rbatch[1:2,]
hist(rbatch$pvalue, main="Capbatch_CB7_vs_CB6", xlab="p-value", breaks=50)
dev.off()

# ------------------------------------------------------------------------
# summarize p-value distribution after removing lowly expressed genes
# ------------------------------------------------------------------------
# 
# table(q75 >= 20)
# w2kp = which(q75 >= 20)
# 
# pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step1_DESeq2_",grp,"_pval_hist_q75_ge20_org_org.pdf"), width=9, height=9)
# par(mfrow=c(3,3), bty="n", mar=c(5,4,2,1))
# for(k in 1:length(nms)){
#   hist(pvals[w2kp,k], main=nms[k], xlab="p-value", breaks=50)
# }
# 
# hist(rbatch$pvalue[w2kp], main="Capbatch_CB7_vs_CB6", xlab="p-value", breaks=50)
# dev.off()
# 

# ------------------------------------------------------------------------
# now we conclude that RMI is not significantly associated with 
# gene expression, and Capbatch CB1 and CB2 are similar, so does 
# Capbatch CB6 and CB7. so we can replace them by sequencing batch
# ------------------------------------------------------------------------

table(meta_ind$Seqbatch, meta_ind$Capbatch)
table(meta_ind$Seqbatch, meta_ind$diagnosis)

dds = DESeqDataSetFromMatrix(countData = trec1, 
                             colData = colData,
                             design = ~ age + sex + Seqbatch + RIN + diagnosis)
dds = DESeq(dds)

res = results(dds)
dim(res)
res[1:2,]
summary(res)

nms = resultsNames(dds)
nms
nms = nms[-1]

pvals2 = matrix(NA, nrow=nrow(trec1), ncol=length(nms))

for(k in 1:length(nms)){
  rk = results(dds, name=nms[k])
  pvals2[,k] = rk$pvalue
}

colnames(pvals2) = nms
dim(pvals2)
summary(pvals2)

pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step1_DESeq2_",grp,"_pval_hist_final_org.pdf"), width=9, height=6)
par(mfrow=c(2,3), bty="n", mar=c(5,4,2,1))
for(k in 1:length(nms)){
  hist(pvals2[,k], main=nms[k], xlab="p-value", breaks=50)
}

plot(-log10(res1$pvalue), -log10(res$pvalue), main="-log10(p-value)", 
     xlab="with all covariates", ylab="with selected covariates", 
     pch=20, cex=0.2, col=rgb(0.8, 0.1, 0.1, 0.5))
abline(0, 1, col="darkblue")

dev.off()

# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

dim(res0)
res0[1:2,]

dim(res)
res[1:2,]

res0 = as.data.frame(res0)
res  = as.data.frame(res)

write.table(res0, paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_no_covariates_org.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(res, paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_adj_covariates_org.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(trec1,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_edata_org.txt"),quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

write.table(meta, paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_meta_cell_org.txt"), sep="\t")
write.table(meta_ind, paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_meta_ind_org.txt"), sep="\t")

gc()

# ------------------------------------------------------------------------
# re-run analysis after permuting case/control labels
# ------------------------------------------------------------------------

set.seed(2020)
nrep  = 9
pvals = matrix(NA, nrow=nrow(trec1), ncol=nrep)
permuted_x = list()

for(k in 1:nrep){
  cat(k, date(), "\n")
  
  colData$diagnosisP = colData$diagnosis
  wControl = which(colData$diagnosis == "Control")
  wCase    = which(colData$diagnosis == "ASD")
  
  ww1 = sample(wControl, size=round(length(wControl)/2))
  ww2 = sample(wCase, size=round(length(wControl)/2))
  
  colData$diagnosisP[ww1] = "ASD"
  colData$diagnosisP[ww2] = "Control"
  
  colData$diagnosisP = factor(colData$diagnosisP, levels=c("Control", "ASD"))
  table(colData$diagnosis, colData$diagnosisP)
  
  permuted_x[[k]] = colData$diagnosisP
  
  ddP = DESeqDataSetFromMatrix(countData = trec1, 
                               colData = colData,
                               design = ~ age + sex + Seqbatch + RIN + diagnosisP)
  ddP  = DESeq(ddP)
  resP = results(ddP)
  
  pvals[,k] = resP$pvalue
}

pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step1_DESeq2_",grp,"_pval_hist_final_permuted_org.pdf"),width=9, height=9)
par(mfrow=c(3,3), bty="n", mar=c(5,4,2,1))
for(k in 1:nrep){
  hist(pvals[,k], xlab="p-value", main="", breaks=50)
}

dev.off()

# ------------------------------------------------------------------------
# occasionally, permuted lable can be significant, check PCA of trec1
# ------------------------------------------------------------------------

edat = t(trec1)
dim(edat)
edat[1:2,1:5]
table(edat==0)

edat = log((edat + 0.5)/(rowSums(edat)/1e6))
dim(edat)
edat[1:2,1:4]

dim(meta_ind)
meta_ind[1:2,]

table(meta_ind$sample == rownames(edat))
meta_ind$rd = colSums(trec1)/1e6
summary(meta_ind$rd)

edat.r1 = edat.r2 = edat
r2s = matrix(NA, nrow=ncol(edat), ncol=2)

for(i in 1:ncol(edat)){
  lmi1 = lm(edat[,i] ~ log(rd), data=meta_ind)
  lmi2 = lm(edat[,i] ~ log(rd) + age + sex + Seqbatch + RIN, data=meta_ind)
  
  edat.r1[,i] = lmi1$residuals
  edat.r2[,i] = lmi2$residuals
  
  s1 = summary(lmi1)
  s2 = summary(lmi2)
  
  r2s[i,] = c(s1$r.squared, s2$r.squared)
}

summary(r2s)

ppk  = propack.svd(edat, neig=50)
lapply(ppk, dim)
ppk$d[1:10]

ppk1  = propack.svd(edat.r1, neig=50)
ppk1$d[1:10]
pca1 = t(ppk1$d*t(ppk1$u))

ppk2  = propack.svd(edat.r2, neig=50)
ppk2$d[1:10]
pca2 = t(ppk2$d*t(ppk2$u))

z = meta_ind[,c("age", "sex",  "Seqbatch", "RIN")]
dim(z)
z[1:2,]

z = model.matrix(~-1 + ., data=z)
dim(z)
z[1:2,]

y0 = as.numeric(colData$diagnosis)
y1 = as.numeric(permuted_x[[9]])

table(y0, y1)

cor(pca1[,1:3], y1)
cor(pca1[,1:3], cbind(z, y0))

cor(pca2[,1:3], y1)
cor(pca2[,1:3], cbind(z, y0))

pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step1_DESeq2_",grp,"_resid_PCs_vs_permuted_case_control_org.pdf"), width=9, height=9)
par(mfrow=c(3,3), mar=c(5,4,1,1), bty="n")
for(kk in 1:3){
  boxplot(pca2[,kk] ~ y0, xlab="case/control", 
          ylab=paste("PC", kk), outline=FALSE)
  stripchart(pca2[,kk] ~ y0, method = "jitter", add=TRUE, vertical=TRUE)
}
for(kk in 1:3){
  boxplot(pca2[,kk] ~ y1, xlab="permuted case/control", 
          ylab=paste("PC", kk), outline=FALSE)
  stripchart(pca2[,kk] ~ y1, method = "jitter", add=TRUE, vertical=TRUE)
}
barplot(ppk2$d^2)
dev.off()


sessionInfo()
#q(save="no")
