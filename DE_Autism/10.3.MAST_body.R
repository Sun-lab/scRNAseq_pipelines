#this code care with the results of simulation data, and calculate the MAST

#note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
# 
# file_tag=1
# r_mean=1.5
# r_var=1.5
# r_disp=1.5
# r_change_prop=0.75
# 
# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

n_seq=c(20,15,10,5)
ncell_seq=c(100,80,60,40,20)

sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

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




for(ncell in ncell_seq){
  for(n in n_seq){
    #set labels
    selected_index=sample.int(20,n)
    total_cell_index=matrix(ncol=1,nrow=0)
    for(i_s in c(selected_index,(20+selected_index))){
      cell_index=(100*i_s-ncell+1):(100*i_s)
      total_cell_index=c(total_cell_index,cell_index)
    }
    
    #calculation
    cur_sim_matrix_log=sim_matrix_log[,total_cell_index]
    
    cell_id = colnames(cur_sim_matrix_log)   #get the cell id from the data
    gene_id = rownames(cur_sim_matrix_log)   #get the gene id from the data
    
    fData = data.frame(primerid = gene_id)
    cData = data.frame(wellKey = cell_id)
    
    meta=t_meta[total_cell_index,]
    
    
    diagnosis = as.character(meta$phenotype) #
    diagnosis[diagnosis == 1] = "Case"
    diagnosis[diagnosis == 0] = "Control"
    
    sca = MAST::FromMatrix(cur_sim_matrix_log, cData, fData)
    colData(sca)$cngeneson = as.numeric(meta$CDR)
    colData(sca)$diagnosis = as.factor(diagnosis)
    colData(sca)$ind = as.factor(meta$individual)
    
    colData(sca)
    
    b0=NA
    b1=NA
    lrt0=NA
    lrt1=NA
    MAST_pval0=NA
    MAST_pval1=NA
    
    date()
    b0 = tryCatch(MAST::zlm(formula = ~ diagnosis, sca = sca, parallel = TRUE), error = function(e) {NA} )
    date()
    b1 = tryCatch(MAST::zlm(formula = ~ diagnosis + ( 1 | ind ), sca = sca, method = 'glmer', ebayes = FALSE, parallel = TRUE), error = function(e) {NA} )
    date()
    
    lrt0 = tryCatch(MAST::lrTest(b0, "diagnosis"), error = function(e) {NA} )
    lrt1 = tryCatch(MAST::lrTest(b1, "diagnosis"), error = function(e) {NA} )
    date()
    
    MAST_pval0 = tryCatch(apply(lrt0, 1, function(x){x[3,3]}), error = function(e) {NA} )
    
    MAST_pval1 = tryCatch(apply(lrt1, 1, function(x){x[3,3]}), error = function(e) {NA} )
    
    tryCatch(saveRDS(MAST_pval0,paste0("../Data_PRJNA434002/10.Result/MAST_pval0_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds")), error = function(e) {NA} )
    tryCatch(saveRDS(MAST_pval1,paste0("../Data_PRJNA434002/10.Result/MAST_pval1_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds")), error = function(e) {NA} )
    print(c(n,ncell))
  }
}



sessionInfo()
q(save="no")
