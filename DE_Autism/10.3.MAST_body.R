#this code care with the results of simulation data, and calculate the MAST

#note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
# 
# file_tag=1
# r_mean=1.5
# r_var=1.5
# r_change_prop=0.6
# 
# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
perm_label=0
perm_method="" #c("","b") 
sim_folder="sim_v6"

n_seq=c(50,30,20,10,5)
ncell_seq=c(200,100,50,20)

sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
t_meta=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

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
    selected_index=sample.int(max(n_seq),n)
    total_cell_index=matrix(ncol=1,nrow=0)
    for(i_s in c(selected_index,(max(n_seq)+selected_index))){
      cell_index=(max(ncell_seq)*i_s-ncell+1):(max(ncell_seq)*i_s)
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
    
    
    if(perm_label>0){
      #count cases and controls
      diag_info=paste0(meta$individual,":",meta$phenotype)
      diag_kind=unique(diag_info)
      diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
      
      #permute
      if(perm_method=="b"){
        
        n_exchange=n/2
        n_exchange=floor(n_exchange)+(n_exchange-floor(n_exchange))*2*rbinom(1,1,0.5) #the number changed to other side
        
        i_exchange=sample.int(n,n_exchange)
        perm_order=1:(2*n)
        temp=perm_order[i_exchange]
        perm_order[i_exchange]=perm_order[(i_exchange+n)]
        perm_order[(i_exchange+n)]=temp
      }
      if(perm_method!="b"){
        perm_order=sample.int(2*n)
        
      }
      
      diag_kind[,2]=diag_kind[perm_order,2]
      #match back to each individuals
      ind_index=match(meta$individual,diag_kind[,1])
      diagnosis=as.factor(diag_kind[ind_index,2])
    }
    
    diagnosis2=matrix("Control",ncol=1,nrow=length(diagnosis))
    diagnosis2[which(diagnosis == "1")] = "Case"
    diagnosis= as.factor(diagnosis2)
    
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
    
    print(paste0("print system details, ncell= ",ncell,", n= ",n, ", before b0"))
    print(date())
    print(gc())
    b0 = tryCatch(MAST::zlm(formula = ~ diagnosis, sca = sca, parallel = TRUE), error = function(e) {NA} )
    print(paste0("print system details, ncell= ",ncell,", n= ",n, ", after b0"))
    print(date())
    print(gc())
    b1 = tryCatch(MAST::zlm(formula = ~ diagnosis + ( 1 | ind ), sca = sca, method = 'glmer', ebayes = FALSE, parallel = TRUE), error = function(e) {NA} )
    print(paste0("print system details, ncell= ",ncell,", n= ",n, ", after b1"))
    print(date())
    print(gc())
    lrt0 = tryCatch(MAST::lrTest(b0, "diagnosis"), error = function(e) {NA} )
    print(paste0("print system details, ncell= ",ncell,", n= ",n, ", after lrTest b0"))
    print(date())
    print(gc())
    lrt1 = tryCatch(MAST::lrTest(b1, "diagnosis"), error = function(e) {NA} )
    print(paste0("print system details, ncell= ",ncell,", n= ",n, ", after lrTest b1"))
    print(date())
    print(gc())
    MAST_pval0 = tryCatch(apply(lrt0, 1, function(x){x[3,3]}), error = function(e) {NA} )
    MAST_pval1 = tryCatch(apply(lrt1, 1, function(x){x[3,3]}), error = function(e) {NA} )
    print(paste0("print system details, ncell= ",ncell,", n= ",n, ", after all"))
    print(date())
    print(gc())
    tryCatch(saveRDS(MAST_pval0,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_pval/p",perm_label,perm_method,"_MAST_pval0_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds")), error = function(e) {NA} )
    tryCatch(saveRDS(MAST_pval1,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_pval/p",perm_label,perm_method,"_MAST_pval1_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds")), error = function(e) {NA} )
    print(c(n,ncell))
  }
}



sessionInfo()
q(save="no")
