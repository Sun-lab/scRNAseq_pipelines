#this code is from the 9.3.MAST_rowcount.R
#this code applied to dataset GSE129788
#this code calculate the MAST on the simulated count data from DCA/scvi output

#cluster_tag=1
#file_tag="1k"

covariate_flag=NA #c(NA, "quantile99")
perm_label=1 #perm_label =0 means calculate the observed data other wise, permutated data

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Volumes/SpecialPass/fh_data/GSE129788/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

###########input###############
#input phenotype
tmeta=readRDS(paste0("../GSE129788/meta",file_tag,".rds"))

total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)


#input counts
logcount_data=readRDS(paste0("../GSE129788/exprM",file_tag,".rds"))
logcount_data=as.matrix(logcount_data[,tmeta$cluster==cur_cluster,drop=FALSE])

dim(logcount_data)
dim(meta)

CDR=readRDS("../GSE129788/CDR1k_per_cluster.rds")
CDR_ind=readRDS("../GSE129788/ind_CDR1k_per_cluster.rds")

CDR=CDR[tmeta$cluster==cur_cluster,cur_cluster]
CDR_ind=CDR_ind[,cur_cluster]
dim(CDR)
dim(CDR_ind)
######################Other method comparison: MAST #################################
print("start MAST calculation")

cell_id=as.character(meta$NAME)
gene_id=as.character(dimnames(logcount_data)[[1]])
rownames(logcount_data)=gene_id
colnames(logcount_data)=cell_id


library("MAST")
library("lme4")

fData=data.frame(primerid=gene_id)
rownames(fData) = fData$primerid
cData=data.frame(wellKey=cell_id)
rownames(cData) = cData$wellKey
colnames(meta)
dim(fData)
dim(cData)

for(perm_label in 0:10){
  sca=FromMatrix(logcount_data, cData, fData)
  colData(sca)$cngeneson = as.numeric(CDR) 
  colData(sca)$phenotype =as.factor(meta$phenotype)
  colData(sca)$ind = as.factor(meta$individual)
  
  if(perm_label>0){
    #count cases and controls
    diag_info=paste0(colData(sca)$ind,":",colData(sca)$phenotype)
    diag_kind=unique(diag_info)
    diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
    
    #permute
    perm_order=readRDS(paste0("../GSE129788/ind_perm_order.rds"))
    perm_order=as.numeric(perm_order[,perm_label])
    total_individual_ref=total_individual[perm_order]
    perm_order=match(total_individual_ref,cur_individual)
    perm_order=perm_order[!is.na(perm_order)]
    diag_kind[,2]=diag_kind[perm_order,2]
    
    #match back to each individuals
    ind_index=match(colData(sca)$ind,diag_kind[,1])
    colData(sca)$phenotype=as.factor(diag_kind[ind_index,2])
  }
  
  colData(sca)
  
  
  
  print(paste0("print system details, before b0"))
  date()
  gc()
  b0 = zlm(formula = ~ phenotype, sca = sca, parallel = TRUE)
  print(paste0("print system details, before b1"))
  date()
  gc()
  b1 = zlm(formula = ~ phenotype + ( 1 | ind ), sca = sca, method = 'glmer', 
           ebayes = FALSE, parallel = TRUE)
  print(paste0("print system details, after b1"))
  date()
  gc()
  
  b0
  b1
  print(paste0("print system details, before lrt0"))
  date()
  gc()
  lrt0 = lrTest(b0, "phenotype")
  print(paste0("print system details, before lrt1"))
  date()
  gc()
  lrt1 = lrTest(b1, "phenotype")
  print(paste0("print system details, after lrt1"))
  date()
  gc()
  dim(lrt1)
  lrt1[1,,]
  
  MAST_pval0 = apply(lrt0, 1, function(x){x[3,3]})
  length(MAST_pval0)
  MAST_pval0[1:4]
  
  MAST_pval1 = apply(lrt1, 1, function(x){x[3,3]})
  length(MAST_pval1)
  MAST_pval1[1:4]
  
  saveRDS(MAST_pval0,paste0("../GSE129788/MAST_pval/p",perm_label,"_MAST_pval0_rawcount_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(MAST_pval1,paste0("../GSE129788/MAST_pval/p",perm_label,"_MAST_pval1_rawcount_",cluster_tag,"_",file_tag,".rds"))
  
}



sessionInfo()
q(save="no")
