
library(data.table)
library(Matrix)
library(MASS)
library(dplyr)
library(doParallel)

library(ggplot2)
library(ggpubr)
library(ggpointdensity)
theme_set(theme_bw())

dist_method = "JSD"
fit_method  = "zinb"

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/functions/ZINB_fit_functions.R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/functions/kl_divergence_functions.R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/functions/Fstat_functions.R")

# number of cores for multi-core computation
nCore = 5

registerDoParallel(cores=nCore)
options(mc.cores=nCore)

# directory of input data. These files are too large to save in GitHubs
data.dir = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data"

args = commandArgs(TRUE)
args

if (length(args) != 1) {
  message("one argument is expected, use 'PFC_L2_3' as default.\n")
  grp = "PFC_L2_3"
}else{
  eval(parse(text=args[[1]]))
}

grp
grp1 = gsub("PFC_", "", grp)
grp1

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

meta = read.table("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step4_meta_org.tsv",header = TRUE, sep = "\t")
dim(meta)
meta[1:2,]

meta_ind = read.table("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step4_meta_ind_org.tsv",header = TRUE, sep = "\t")
dim(meta_ind)
meta_ind[1:2,]

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

dat1 = readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data/ct_mtx/",grp,".rds"))
class(dat1)

dim(dat1)
dat1[1:5,1:4]


dat1 = dat1[,match(meta$cell, colnames(dat1))]
dim(dat1)

n.zeros = rowSums(dat1 == 0)
summary(n.zeros)

table(n.zeros < 0.8*ncol(dat1))

w2kp = which(n.zeros < 0.8*ncol(dat1))
dat1 = dat1[w2kp,]

dim(dat1)
dat1[1:5,1:4]

# ------------------------------------------------------------------------
# Calculate distance across individuals
# ------------------------------------------------------------------------
# 
#   suppose we have distribution P and Q
#   Define: M=(P+Q)/2 and D(P||Q)=-sum((P * (log(Q/P)))), then
#
#   the Jensen-Shannon divergence:
#           JSD=1/2 * (D(P||M)+D(Q||M))

zinb_fit_dca = readRDS("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step4a_zinb_fit_dca_org.rds")

dist_array = array(
  dim = c(
    length(zinb_fit_dca),
    nrow(meta_ind),
    nrow(meta_ind)
  ),
  dimnames = list(names(zinb_fit_dca), meta_ind$individual, meta_ind$individual)
)

dim(dist_array)

table(names(zinb_fit_dca[[1]]) == meta_ind$individual)

log10.rd = 4

extrac_zinb_par <- function(zb_fit, log10.rd){
  mean_a = exp(t(zb_fit$logmean) %*% c(1, log10.rd))
  disp_a = zb_fit$dispersion
  drop_a = exp(t(zb_fit$logitdropout) %*% c(1, log10.rd))
  drop_a = drop_a/(1 + drop_a)
  c(mean_a, disp_a, drop_a)
}

date()
dist_array_list=foreach (i = 1:length(zinb_fit_dca)) %dopar% {
  gene_i_fit  = zinb_fit_dca[[i]]
  dist_array1 = array(dim=rep(nrow(meta_ind), 2))
  rownames(dist_array1) = meta_ind$individual
  colnames(dist_array1) = meta_ind$individual
  
  for (j_a in 1:nrow(meta_ind)) {
    for (j_b in 1:nrow(meta_ind)) {
      par_a = extrac_zinb_par(gene_i_fit[[j_a]], log10.rd)
      par_b = extrac_zinb_par(gene_i_fit[[j_b]], log10.rd)
      
      dist_array1[j_a, j_b] = tryCatch(
        divergence(par_a, par_b, method = dist_method,
                   zinb.quantile = 0.975, fit_model = fit_method), 
        error = function(e) { NA }
      )
    }
  }
  
  dist_array1
}
date()

for (i in 1:length(zinb_fit_dca)){
  dist_array[i,,] = dist_array_list[[i]]
}

dim(dist_array)
dist_array[1,1:2,1:2]

# ------------------------------------------------------------------------
# compare with results from DESeq2
# ------------------------------------------------------------------------

deseq2_no_cov = read.table(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_no_covariates_org.txt"), sep="\t", header=TRUE)
dim(deseq2_no_cov)
deseq2_no_cov[1:2,]

deseq2 = read.table(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step1_DESeq2_",grp,"_adj_covariates_org.txt"), sep="\t", header=TRUE)
dim(deseq2)
deseq2[1:2,]

df = read.table("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step5a_pvals_org.tsv",header = TRUE, sep = "\t")
dim(df)
df[1:2,]

table(df$gene %in% rownames(deseq2_no_cov))
deseq2_no_cov = deseq2_no_cov[match(df$gene, rownames(deseq2_no_cov)),]
dim(deseq2_no_cov)
deseq2_no_cov[1:2,]

deseq2 = deseq2[match(df$gene, rownames(deseq2)),]
dim(deseq2)
deseq2[1:2,]

dim(df)
df[1:2,]

df$deseq2_rd = deseq2_no_cov$pvalue
df$deseq2    = deseq2$pvalue

dim(df)
df[1:2,]

minP = min(df$pval_no_cov[which(df$pval_no_cov > 0)])
minP

for(k in 2:5){
  df[[k]][which(df[[k]] == 0)] = minP/2
}

gs1 = ggplot(df,aes(x=-log10(deseq2_rd),y=-log10(pval_no_cov))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, rd only)", 
       x = "-log10(pval DESeq2, rd only)")

gs2 = ggplot(df,aes(x=-log10(deseq2_rd),y=-log10(pval_Z))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute G)", 
       x = "-log10(pval DESeq2, rd only)")

gs3 = ggplot(df,aes(x=-log10(deseq2_rd),y=-log10(pval_SZ))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute x, adj z)", 
       x = "-log10(pval DESeq2, rd only)")

gs4 = ggplot(df,aes(x=-log10(deseq2),y=-log10(pval_Z))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute G)", x = "-log10(pval DESeq2)")

gs5 = ggplot(df,aes(x=-log10(deseq2),y=-log10(pval_S))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute x)", x = "-log10(pval DESeq2)")

gs6 = ggplot(df,aes(x=-log10(deseq2),y=-log10(pval_SZ))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute x, adj z)", 
       x = "-log10(pval DESeq2)")

gs7 = ggplot(df,aes(x=-log10(pval_no_cov),y=-log10(pval_SZ))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute x, adj z)", 
       x = "-log10(pval IDEAS, rd only)")

gs8 = ggplot(df,aes(x=-log10(pval_Z),y=-log10(pval_SZ))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(x="-log10(pval IDEAS, permute G)", 
       y="pval IDEAS, permute x, adj z)")

gs9 = ggplot(df,aes(x=-log10(pval_S),y=-log10(pval_SZ))) +
  geom_pointdensity(size = 0.6) + scale_color_viridis_c() +
  labs(y= "-log10(pval IDEAS, permute x, adj z)", 
       x = "-log10(pval IDEAS, permute x)")

pdf(file="/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_DESeq2_vs_IDEAS_org.pdf", width=13.5, height=9)
ggarrange(gs1, gs2, gs3, gs4, gs5, gs6, gs7, gs8, gs9, 
          ncol = 3, nrow = 3)
dev.off()

round(cor(-log10(df[,-1]), method="spearman"),3)

# ------------------------------------------------------------------------
# function to produce summary and plots
# ------------------------------------------------------------------------

meta_ind$rd = meta_ind$rd/1e6
meta_ind$diagnosis = factor(meta_ind$diagnosis, levels=c("Control", "ASD"))

plot.deseq2.vs.ideas <- function(meta_ind, gnm, gdata, gene_i_fit, dist1, 
                                 fnm, max.x=10){
  meta_ind$gene1 = gdata
  
  nb1 = glm.nb(gene1 ~ log10(rd) + diagnosis, data=meta_ind)
  print(summary(nb1))
  
  nb2 = glm.nb(gene1 ~ log10(rd) + diagnosis + age + sex + Seqbatch + RIN, 
               data=meta_ind)
  print(summary(nb2))
  
  # ------------------------------------------------------------------------
  # boxplot of gene expression vs. diagnosis
  # ------------------------------------------------------------------------
  
  p1 = ggplot(meta_ind, aes(x=diagnosis, y=log10(gene1), col=diagnosis)) + 
    geom_boxplot() + labs(y=sprintf("log10(%s)", gnm))
  p1 = p1 + geom_jitter(shape=16, position=position_jitter(0.2))
  
  p2 = ggplot(meta_ind, aes(x=diagnosis, y=log10(gene1/rd), col=diagnosis)) + 
    geom_boxplot() + labs(y=sprintf("log10(%s/rd)", gnm))
  p2 = p2 + geom_jitter(shape=16, position=position_jitter(0.2))
  
  # ------------------------------------------------------------------------
  # plot density for each sample
  # ------------------------------------------------------------------------
  
  log10.rd = 4
  
  xx = yy = NULL
  for(i in 1:nrow(meta_ind)){
    xx = c(xx, 0:max.x)
    pa = extrac_zinb_par(gene_i_fit[[i]], log10.rd)
    yy = c(yy, emdbook::dzinbinom(0:max.x, pa[1], pa[2], pa[3]))
  }
  
  df.den = data.frame(ind=rep(meta_ind$individual, each=(max.x+1)), 
                      diagnosis=rep(meta_ind$diagnosis, each=(max.x+1)), 
                      x=xx, density=yy)
  dim(df.den)
  df.den[1:2,]
  
  p3 = ggplot(df.den, aes(x=x, y=density, group=ind, col=diagnosis)) +
    geom_line(aes(linetype=diagnosis))
  
  # ------------------------------------------------------------------------
  # boxplot of distance, within group vs. between group
  # ------------------------------------------------------------------------
  
  n = nrow(meta_ind)
  row.status = matrix(rep(meta_ind$diagnosis, times = n), ncol=n)
  col.status = matrix(rep(meta_ind$diagnosis, each = n),  ncol=n)
  row.status[1:2,1:2]
  col.status[1:2,1:2]
  
  same.group = (row.status == col.status)
  same.group[1:2,1:2]
  
  d.within.group  = dist1[same.group & upper.tri(dist1)]
  d.between.group = dist1[(!same.group) & upper.tri(dist1)]
  
  dist2 = data.frame(dist=c(d.within.group, d.between.group))
  dist2$label = rep(c("within", "between"), 
                    times=c(length(d.within.group), length(d.between.group)))
  
  p4 = ggplot(dist2, aes(x=label, y=dist, col=label)) + 
    geom_boxplot() + scale_color_manual(values=c("#E69F00", "#56B4E9"))
  p4 = p4 + geom_jitter(shape=16, position=position_jitter(0.2))
  
  pdf(fnm, width=7.5, height=5)
  print(ggarrange(p1, p2, p3, p4, labels = c("A","B","C","D"), ncol=2, nrow=2))
  dev.off()
}

# ------------------------------------------------------------------------
# compare DESeq2 means missed and captured by ideas
# ------------------------------------------------------------------------

table(df$deseq2 < 1e-2)

c1 = chisq.test(df$deseq2 < 1e-2, df$pval_SZ > 0.1)
c1
c1$observed
c1$expected

grp1 = which(df$deseq2 < 1e-2 & df$pval_SZ > 0.1)
grp2 = which(df$deseq2 < 1e-2 & df$pval_SZ > 0.02 & df$pval_SZ <= 0.1)
grp3 = which(df$deseq2 < 1e-2 & df$pval_SZ <= 0.02)

length(grp1)
length(grp2)
length(grp3)

summary(df$deseq2_rd[grp1])
summary(df$deseq2_rd[grp2])
summary(df$deseq2_rd[grp3])

pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_DESeq2_rd_only_pval_",grp,"_by_DESeq2_and_ideas_pval_org.pdf"), width=5, height=3)
par(mar=c(5,4,1,1), bty="n")
plot(density(-log10(df$deseq2_rd[grp1]), adjust=0.5), 
     xlab="DESeq2 p-value, rd only", ylim=c(0,1), main="")
lines(density(-log10(df$deseq2_rd[grp2]), adjust=0.5), lty=2, 
      col="blue", lwd=1.5)
lines(density(-log10(df$deseq2_rd[grp3]), adjust=0.5), lty=4, 
      col="red", lwd=1.5)
legend("topleft", lty=c(1,2,4), col=c("black", "blue", "red"), bty="n", 
       legend=c("p > 0.1", "0.02 < p <= 0.1", "p <= 0.02"), 
       lwd=c(1,1.5,1.5))
dev.off()

# ------------------------------------------------------------------------
# pick one gene missed by ideas
# ------------------------------------------------------------------------

w2check = which(df$deseq2 < 1e-3 & df$pval_SZ > 0.1)
df[w2check,]

dim(dat1)
dat1[1:2,1:4]

gene1 = dat1[which(rownames(dat1)=="CNOT2-DT"),]
length(gene1)
table(gene1)

gdata = tapply(gene1, as.character(meta$sample), sum)
gdata
table(names(gdata) == meta_ind$sample)

gnm = "CNOT2.DT"
gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == "CNOT2-DT")]]
table(names(gene_i_fit) == meta_ind$individual)

dist1 = dist_array[which(dimnames(dist_array)[[1]]=="CNOT2-DT"),,]
dim(dist1)
dist1[1:5,1:5]
table(rownames(dist1) == meta_ind$individual)

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=10, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2+_ideas-CNOT2-DT_org.pdf")
  
# ------------------------------------------------------------------------
# pick one gene missed by DESeq2
# ------------------------------------------------------------------------

w2check = which(df$deseq2 > 0.1 & df$pval_SZ < 0.001)
df[w2check,]

gnm = "MICU2"

gene1 = dat1[which(rownames(dat1)==gnm),]
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=7, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2-_ideas+_MICU2_org.pdf")


gnm = "MBP"

gene1 = dat1[which(rownames(dat1)==gnm),]
table(gene1)
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=7, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2-_ideas+_MBP_org.pdf")


gnm = "SLC38A10"

gene1 = dat1[which(rownames(dat1)==gnm),]
table(gene1)
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=4, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2-_ideas+_SLC38A10_org.pdf")

gnm = "GNG7"

gene1 = dat1[which(rownames(dat1)==gnm),]
table(gene1)
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=7, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2-_ideas+_GNG7_org.pdf")

# ------------------------------------------------------------------------
# pick one gene missed by DESeq2
# ------------------------------------------------------------------------

w2check = which(df$deseq2 <= 0.001 & df$pval_SZ <= 0.001)
df[w2check,]

gnm = "CDK11B"

gene1 = dat1[which(rownames(dat1)==gnm),]
table(gene1)
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=5, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2+_ideas+_CDK11B_org.pdf")

gnm = "TCF25"

gene1 = dat1[which(rownames(dat1)==gnm),]
table(gene1)
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=5, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2+_ideas+_TCF25_org.pdf")

gnm = "PNCK"

gene1 = dat1[which(rownames(dat1)==gnm),]
table(gene1)
gdata = tapply(gene1, as.character(meta$sample), sum)

gene_i_fit  = zinb_fit_dca[[which(names(zinb_fit_dca) == gnm)]]
dist1 = dist_array[which(dimnames(dist_array)[[1]]==gnm),,]

plot.deseq2.vs.ideas(meta_ind, gnm, gdata, gene_i_fit, dist1, max.x=5, 
                     fnm = "/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step6_ex_DESeq2+_ideas+_PNCK_org.pdf")


gc()

sessionInfo()
#q(save="no")
