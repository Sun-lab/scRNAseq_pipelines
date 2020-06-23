
library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(doParallel)
library(svd)

library(ggplot2)
library(ggpubr)
library(ggpointdensity)

dist_method = "JSD"
fit_method  = "zinb"

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/functions/ZINB_fit_functions.R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/functions/kl_divergence_functions.R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/functions/Fstat_functions.R")

# number of cores for multi-core computation
nCore = 10

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

# ========================================================================
# read input
# ========================================================================

# ------------------------------------------------------------------------
# read in cell information
# ------------------------------------------------------------------------

cell_info = read.table(paste0("/fh/fast/sun_w/mengqi/Data_PRJNA434002/meta.tsv"),header = TRUE, sep = "\t")
dim(cell_info)
cell_info[1:2,]

# ------------------------------------------------------------------------
# read in count data of one region and one cluster
# ------------------------------------------------------------------------

dat1 = readRDS(file.path(sprintf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data/ct_mtx/%s.rds", grp)))
class(dat1)

dim(dat1)
dat1[1:5,1:4]

# ------------------------------------------------------------------------
# subset cell information
# ------------------------------------------------------------------------

table(colnames(dat1) %in% cell_info$cell)
meta = cell_info[match(colnames(dat1), cell_info$cell),]
dim(meta)
meta[1:2,]

names(meta)[11:12] = c("PMI", "RIN")
names(meta)[15:16] = c("mitoPercent", "riboPercent")
dim(meta)
meta[1:2,]

summary(meta)
meta$age = scale(meta$age)
meta$PMI = scale(meta$PMI)
meta$RIN = scale(meta$RIN)
meta$Capbatch = meta$Capbatch
meta$Seqbatch = meta$Seqbatch
#meta$Capbatch = droplevels(meta$Capbatch)
#meta$Seqbatch = droplevels(meta$Seqbatch)
meta$individual = as.factor(meta$individual)
summary(meta)

table(meta$Capbatch, meta$Seqbatch)

summary(meta$UMIs/meta$genes)

# check each individual has a unique sample
table(tapply(meta$sample, meta$individual, function(v){length(unique(v))}))

# check each individual has a unique Capbatch
table(tapply(meta$Capbatch, meta$individual, function(v){length(unique(v))}))

table(meta$cluster)
table(meta$region)
table(meta$diagnosis)
sort(table(paste(meta$individual, meta$diagnosis, sep=":")))
tt1 = table(meta$individual)
sort(tt1)

n.zeros = rowSums(dat1 == 0)
summary(n.zeros)

0.6*ncol(dat1)
0.8*ncol(dat1)

table(n.zeros < 0.6*ncol(dat1))
table(n.zeros < 0.8*ncol(dat1))

w2kp = which(n.zeros < 0.8*ncol(dat1))
dat1 = dat1[w2kp,]

dim(dat1)
dat1[1:5,1:4]

# ------------------------------------------------------------------------
# generate individual level information
# ------------------------------------------------------------------------

length(unique(meta$individual))

meta_ind = distinct(meta[,3:12])
dim(meta_ind)
meta_ind[1:2,]
#meta_ind$diagnosis = relevel(meta_ind$diagnosis, ref="Control")
table(meta_ind$diagnosis)

if(nrow(meta_ind) != length(unique(meta$individual))){
  stop("there is non-unique information\n")
}

table(meta_ind$Seqbatch, meta_ind$Capbatch)

# ------------------------------------------------------------------------
# add read-depth information
# ------------------------------------------------------------------------

dim(dat1)
dat1[1:2,1:5]

dim(meta)
meta[1:2,]

table(meta$cell == colnames(dat1))

rd_cell = colSums(dat1)
summary(rd_cell)

meta$rd = rd_cell

med_rd_cell = tapply(rd_cell, meta$individual, median)
med_rd_cell

rd = tapply(rd_cell, meta$individual, sum)
rd

table(names(med_rd_cell) == meta_ind$individual)
meta_ind$med_rd_cell = med_rd_cell
meta_ind$rd = rd
meta_ind[1:2,]

gh = ggplot(meta, aes(x=log10(rd))) +
  geom_histogram(color="darkblue", fill="lightblue") + 
  geom_vline(aes(xintercept=median(log10(rd))),
             color="blue", linetype="dashed", size=1)

pdf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step4_rd_per_cell_org.pdf", width=4, height=3)
gh
dev.off()

# ------------------------------------------------------------------------
# save cell level and individual level information
# ------------------------------------------------------------------------

gc()

write.table(meta_ind, file="/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step4_meta_ind_org.tsv", sep="\t")
write.table(meta, file="/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step4_meta_org.tsv", sep="\t")

# ------------------------------------------------------------------------
# read in DCA estimates
# ------------------------------------------------------------------------

dca_mean = read.table("/fh/scratch/delete90/sun_w/autism/dca_PFC_all/L2_3_mean.tsv.gz",header = TRUE, sep = "\t")
dca_disp = read.table("/fh/scratch/delete90/sun_w/autism/dca_PFC_all/L2_3_dispersion.tsv.gz",header = TRUE, sep = "\t")
dca_pi   = read.table("/fh/scratch/delete90/sun_w/autism/dca_PFC_all/L2_3_pi.tsv.gz",header = TRUE, sep = "\t")

dim(dca_mean)
dim(dca_disp)
dim(dca_pi)

dca_mean[1:2,1:5]
dca_disp[1:2,1:5]
dca_pi[1:2,1:5]

table(meta$cell == colnames(dca_mean)[-1])

dca_mean = dca_mean[w2kp,]
dca_disp = dca_disp[w2kp,]
dca_pi   = dca_pi[w2kp,]

table(rownames(dat1) == dca_mean$V1)
table(rownames(dat1) == dca_disp$V1)
table(rownames(dat1) == dca_pi$V1)

dca_mean = data.matrix(dca_mean[,-1])
dca_disp = data.matrix(dca_disp[,-1])
dca_pi   = data.matrix(dca_pi[,-1])

colnames(dca_mean)=colnames(dat1)
colnames(dca_disp)=colnames(dat1)
colnames(dca_disp)=colnames(dat1)


gc()

# ========================================================================
# perform testing
# ========================================================================

# ------------------------------------------------------------------------
# for all the genes, estimate ZINB estimate using DCA output
# given read-depth per cell and output the estimate when 
# read-depth per cell being 10,000
# ------------------------------------------------------------------------

set.seed(2020)

# number of counts to sample for each gene and each cell
sim_n = 10

date()

# parallele for each gene
zinb_fit_dca=foreach (i = 1:nrow(dat1)) %dopar% {
  
  zinb_fit1 = list()
  
  # each row of dat_simu correspond to a cell
  dat_simu = matrix(nrow=ncol(dat1), ncol=sim_n)
  
  # simulate data per cell
  for(k in 1:ncol(dat1)){
    mean_ik = dca_mean[i,k]
    disp_ik = dca_disp[i,k]
    pi_ik   = dca_pi[i,k]
    dat_simu[k,] = emdbook::rzinbinom(sim_n, mean_ik, disp_ik, pi_ik)
  }
  
  # collapose counts per individual and estimate ZINB
  for (j in 1:nrow(meta_ind)) {
    cur_ind = meta_ind$individual[j]
    w2use   = which(meta$individual == cur_ind)
    dat_ind = c(dat_simu[w2use,])
    rd_cell = rep(log10(meta$rd[w2use]), times=sim_n)
    zinb_fit1[[j]] = fit_zinb(x=dat_ind, z=rd_cell)
  }
  names(zinb_fit1) = as.character(meta_ind$individual)
  
  # the return is logmean, dispersion and logit_dropout_rate
  zinb_fit1
}
date()

length(zinb_fit_dca)
table(sapply(zinb_fit_dca, length))
zinb_fit_dca[[1]][[1]]

names(zinb_fit_dca) = rownames(dat1)

saveRDS(zinb_fit_dca, file="/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step4a_zinb_fit_dca_org.rds")

gc()

sessionInfo()
#q(save="no")
