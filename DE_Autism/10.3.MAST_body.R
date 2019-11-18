#this code care with the results of simulation data, and calculate the MAST

#note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.

#file_tag=1
#sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
#splat.org--method 4, change the mean.shape and mean.rate originally
#zinb.naive--method 5, using naive zinb models to do so.

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#r_mean=1.5  #r_mean/r_var should < 1+mean.shape
#r_var=4

sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

######################Other method comparison: MAST #################################
print("start MAST calculation: Part II: ZINB KLmean and JSD")

library("moments")
library("MAST")
library("lme4")
#######################  3. MAST analysis  #############################

# the input of MAST analysis can be many format including matrix and 
# SingleCellAssay.

# input:
# (1)the log-transformed expression count matrix,
#     with each column represents the cells and each rows represents the genes.
# (2)the meta data, including cell information and individual information.

# we get the p-values based on the Hurdle model ("H" model)

# due to the speed of MAST analysis, we run each permutation on seperate jobs
#  and then gather all results to calculate permutation p-values.

sim_matrix_log = log2(1 + sim_matrix) #log transformed data

dim(sim_matrix_log)
sim_matrix_log[1:10, 1:10]
cell_id = colnames(sim_matrix_log)   #get the cell id from the data
gene_id = rownames(sim_matrix_log)   #get the gene id from the data

diagnosis = as.character(meta$phenotype) #
diagnosis[diagnosis == 1] = "Case"
diagnosis[diagnosis == 0] = "Control"

fData = data.frame(primerid = gene_id)
cData = data.frame(wellKey = cell_id)
colnames(meta)
length(fData)
length(cData)

sca = MAST::FromMatrix(sim_matrix_log, cData, fData)
colData(sca)$cngeneson = as.numeric(meta$CDR)
colData(sca)$diagnosis = as.factor(diagnosis)
colData(sca)$ind = as.factor(meta$individual)

colData(sca)

date()
b0 = MAST::zlm(formula = ~ diagnosis, sca = sca, parallel = TRUE)
date()
b1 = MAST::zlm(formula = ~ diagnosis + ( 1 | ind ), sca = sca, method = 'glmer', 
         ebayes = FALSE, parallel = TRUE)
date()

b0
b1

lrt0 = MAST::lrTest(b0, "diagnosis")
lrt1 = MAST::lrTest(b1, "diagnosis")

dim(lrt1)
lrt1[1,,]

MAST_pval0 = apply(lrt0, 1, function(x){x[3,3]})
length(MAST_pval0)
MAST_pval0[1:4]

MAST_pval1 = apply(lrt1, 1, function(x){x[3,3]})
length(MAST_pval1)
MAST_pval1[1:4]


saveRDS(MAST_pval0,paste0("../Data_PRJNA434002/10.Result/MAST_pval0_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
saveRDS(MAST_pval1,paste0("../Data_PRJNA434002/10.Result/MAST_pval1_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

sessionInfo()
q(save="no")
