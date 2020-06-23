#The following codes are modified from 


#the change includes
#1.folder directory 
#2. input/output data
#3. input/output data rename  
 ###"_org" labels original outputs from wei's data
 ###"_mengqi" labels

#require
#step4a_ideas_zinb_fit.R
#####################
library(MASS)
library(Matrix)
library(data.table)
library(dplyr)
library(doParallel)
library(DESeq2)

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

# ========================================================================
# read input
# ========================================================================
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

dat1 = readRDS(file.path(sprintf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/data/ct_mtx/%s.rds", grp)))
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

# zinb_fit_dca is a list of all cells, each list include a list of all individuals, while each list include 3 parameters from dca
# > class(zinb_fit_dca)
# [1] "list"
# > length(zinb_fit_dca)
# [1] 8260
# > class(zinb_fit_dca[[1]])
# [1] "list"
# > length(zinb_fit_dca[[1]])
# [1] 23
# > class(zinb_fit_dca[[1]][[1]])
# [1] "list"
# > length(zinb_fit_dca[[1]][[1]])
# [1] 3
# > zinb_fit_dca[[1]][[1]]
# $logmean
# (Intercept)           z
# -9.132280    2.117404
# 
# $dispersion
# [1] 7.948641
# 
# $logitdropout
# (Intercept)           z
# 2.742727   -1.233726



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
saveRDS(dist_array,"/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step4a_zinb_fit_dca_dist_array_org.rds")
# ------------------------------------------------------------------------
# check the assocaition between distance matrix and each confoudner
# ------------------------------------------------------------------------

z = meta_ind[,c("age", "sex",  "Seqbatch", "RIN")]
dim(z)
z[1:2,]

z = model.matrix(~., data=z)
dim(z)
z[1:2,]

pval.z = matrix(NA, nrow=dim(dist_array)[1], ncol=ncol(z)-1)

for(k in 2:ncol(z)){
  set.seed(2020)
  
  cat(k, date(), "\n")
  pval.z[,k-1] = permanova(dist_array, z[,k])$pval
}

pval.z = data.frame(pval.z)
names(pval.z) = colnames(z)[-1]
dim(pval.z)
pval.z[1:2,]

glist = list()
for(k in 1:ncol(pval.z)){
  glist[[k]] = ggplot(pval.z, aes_string(x=names(pval.z)[k]))+
    geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02))
}

pdf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step5a_pval_hist_confounders2_org.pdf", width=6, height=6)
ggarrange(glist[[1]], glist[[2]], glist[[3]], glist[[4]], 
          ncol = 2, nrow = 2)
dev.off()

# ------------------------------------------------------------------------
# pval calculation by Permanova
# ------------------------------------------------------------------------

set.seed(2020)
nperm = 1000

x = model.matrix(~diagnosis, data=meta_ind)
dim(x)
x[1:2,]

# first do not use any covariate
date()
ideas_no_cov  = permanova_binary(dist_array, x=x, z=NULL, n_perm=nperm)
date()

# next consider all the covariates, without rd since we have accounted
# for rd in the way to generate distance array

date()
ideas_Z  = permanova_binary(dist_array, x=x, z=z, n_perm=nperm, 
                            permuteDistance=TRUE)
date()
ideas_S = permanova_binary(dist_array, x=x, z=z, n_perm=nperm, 
                           permuteDistance=FALSE, adjust_z=FALSE)
date()
ideas_SZ = permanova_binary(dist_array, x=x, z=z, n_perm=nperm, 
                            permuteDistance=FALSE, adjust_z=TRUE)
date()

df = data.frame(gene=rownames(dat1), pval_no_cov=ideas_no_cov$pval, 
                pval_Z=ideas_Z$pval,  pval_S=ideas_S$pval, 
                pval_SZ=ideas_SZ$pval)
dim(df)
df[1:5,]

gh0 = ggplot(df, aes(x=pval_no_cov)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) +
  labs(title="no covariate")

gh1 = ggplot(df, aes(x=pval_Z)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) +
  labs(title="w/ covariates, permute G")

gh2 = ggplot(df, aes(x=pval_S)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) + 
  labs(title="w/ covariates, permute x")

gh3 = ggplot(df, aes(x=pval_SZ)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) + 
  labs(title="w/ covariates, permute x, adj z")

pdf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step5a_pval_hist2_org.pdf", width=6, height=6)
ggarrange(gh0, gh1, gh2, gh3, ncol = 2, nrow = 2)
dev.off()

# ------------------------------------------------------------------------
# save the results
# ------------------------------------------------------------------------

fwrite(df, file="/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/github/ideas/Autism/R/res/step5a_pvals2_org.tsv", sep="\t")

# ------------------------------------------------------------------------
# it seems the results for a permuted case/control lable can be unstable
# now try to do it 10 times to check, using balanced permuation or not
# ------------------------------------------------------------------------

set.seed(2020)

meta_ind$diagnosis_n = as.numeric(as.factor(meta_ind$diagnosis)) - 1
table(meta_ind$diagnosis, meta_ind$diagnosis_n)

glm1 = glm(meta_ind$diagnosis_n ~ -1 + z, family=binomial(link="logit"))
summary(glm1)

permuted_x = list()
nrep = 9
devs = rep(NA, nrep)

for(k in 1:nrep){
  
  while(1){
    rand_diagnosis  = sample(meta_ind$diagnosis_n)
    permuted_x[[k]] = rand_diagnosis
    
    glmk = tryCatch(glm(rand_diagnosis ~ -1 + z, 
                        family=binomial(link="logit")), 
                    warning = function(w) { NULL }, 
                    error = function(e) { NULL},
                    finally = {})
    
    if(! is.null(glmk)){
      sk   = summary(glmk)
      if(all(sk$coefficients[,4] > 0.25)){ 
        devs[k] = sk$deviance
        break 
      }
    }
  }
}

devs

p0.null = p1.null = p2.null = p3.null = NULL

for(k in 1:nrep){
  cat(k, date(), "\n")
  
  rand_diagnosis = permuted_x[[k]]
  
  dist_test0p = permanova_binary(dist_array, x=rand_diagnosis, 
                                 z=NULL, n_perm=500)
  
  dist_test1p = permanova_binary(dist_array, x=rand_diagnosis, z = z, 
                                 n_perm=500, permuteDistance=TRUE)
  
  dist_test2p = permanova_binary(dist_array, x=rand_diagnosis, z = z, 
                                 n_perm=500, permuteDistance=FALSE,
                                 adjust_z=FALSE)
  
  dist_test3p = permanova_binary(dist_array, x=rand_diagnosis, z = z, 
                                 n_perm=500, permuteDistance=FALSE,
                                 adjust_z=TRUE)
  
  p0.null = c(p0.null, dist_test0p$pval)
  p1.null = c(p1.null, dist_test1p$pval)
  p2.null = c(p2.null, dist_test2p$pval)
  p3.null = c(p3.null, dist_test3p$pval)
}

# ------------------------------------------------------------------------
# illustrate the results for permuted data
# ------------------------------------------------------------------------

df0 = data.frame(p0.null, p1.null, p2.null, p3.null)
dim(df0)
df0[1:2,]

gh0.n = ggplot(df0, aes(x=p0.null)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) +
  labs(title="null, no covariate")

gh1.n = ggplot(df0, aes(x=p1.null)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) +
  labs(title="null, permute G")

gh2.n = ggplot(df0, aes(x=p2.null)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) + 
  labs(title="null, permute x")

gh3.n = ggplot(df0, aes(x=p3.null)) +
  geom_histogram(color="darkblue", fill="lightblue", breaks=seq(0,1,by=0.02)) + 
  labs(title="null, permute x, adj z")

pdf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step5a_pval_hist_permuted_x2_org.pdf", width=6, height=6)
ggarrange(gh0.n, gh1.n, gh2.n, gh3.n, ncol = 2, nrow = 2)
dev.off()

pdf("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.wei_code_check/step5a_pval_hist_permuted_x_details2_org.pdf", width=12, height=9)
par(mfrow=c(3,4))
n_genes = dim(dist_array)[1]
for(k in 1:nrep){
  idx1 = n_genes*(k-1)+1
  idx2 = n_genes*k
  dfk  = df0[idx1:idx2,]
  hist(dfk$p0.null, main="null, no covariate", xlab="p-value", breaks=50)
  hist(dfk$p1.null, main="null, permute G", xlab="p-value", breaks=50)
  hist(dfk$p2.null, main="null, permute x", xlab="p-value", breaks=50)
  hist(dfk$p3.null, main="null, permute x, adj z", xlab="p-value", breaks=50)
}
dev.off()

gc()

sessionInfo()
#q(save="no")

