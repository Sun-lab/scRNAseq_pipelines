#this code used for the

# Simulate scRNA-seq data using splatter
# Using their default simulator for now. We can consider simulate 5000 genes in 40 individuals, with 20 cases vs. 20 controls. We can set DE for 500 genes in terms of mean expression difference, and another 500 genes DE in terms of variance difference. In the simulation, try to set the total read-depth across samples to be the same. 
# As an initial analysis, do not add any covariates. 
# For each gene, calculate density, and then JSD across samples, and then use PERMANOVA to calculate p-value for each gene. 
# Also calculate p-value for each gene using MAST-mixed effect model (less priority for now). 
# Collapse gene expression across cells per individual, then run DESeq2 for differential expression testing. 
# Type I error is the proportion of those 4000 equivalently expressed genes with p-value smaller than 0.05. 
# Power1, the proportion of the first 500 genest with p-value < 0.05
# Power2, the proportion of the next 500 genest with p-value < 0.05






library("GenomicFeatures")
library("GenomicAlignments")
library(DropletUtils)
library(biomaRt)
library(dplyr)
library(scater)


library("splatter")
library("scater")
library("ggplot2")

rawM3k10 = readRDS("../Data_PRJNA434002/rawM3k10.rds")
params=splatEstimate(rawM3k10)



#splat examples
#param default
params=newSplatParams()
params
#genes

params=setParams(params, update = list(nGenes = 5000, 
                                       mean.rate = 0.3,
                                       mean.shape= 0.6
                                       batchCells = rep(100, 40),
                                       ))


# Extract multiple parameters as a list

sim_splat=splatSimulate(params, batchCells=100,verbose = FALSE)
counts(sim_splat)[1:10,1:10]

params=newSplatParams(batchCells=rep(100,40),batch.facLoc = 0.001, batch.facScale = 0.001, nGenes = 5000)

sim1 <- splatSimulate(params, 
                            de.prob = 0.1, verbose = FALSE)









#for case

params=setParams(params, update = list(nGenes = 5000, 
                                       mean.rate = 0.3,
                                       mean.shape= 0.6,
                                       batchCells = rep(100,20),
                                       batch.facLoc = 0.001, batch.facScale = 0.001
                                        ))
#control
control=splatSimulate(params,verbose = FALSE)
#case_mean
case_mean=splatSimulate(params,de.prob= 0.1,de.facLoc= 0.1,verbose = FALSE)
#case_variant
case_var=splatSimulate(params,de.prob= 0.1,de.facScale= 0.1,verbose = FALSE)
