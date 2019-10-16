#this code try to reproduce the MAST analysi among Velmeshev_2019_autism 
#perform the MAST DE analysis and calculate the permutated version of it.
#Data_PRJNA434002

#something from the header file
#cluster_tag=1 #this tag indicate the clusters it can be choose in 1 to 17
#perm_tag=0    #this tage indicate the permutation tags, 0 means no permutation, otherwise, permutation id

setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

file_name="rawMnorm3k10"

meta=readRDS("../Data_PRJNA434002/meta10.rds")
inputM=readRDS(paste0("../Data_PRJNA434002/",file_name,".rds"))
cur_cluster=as.character(unique(meta$cluster)[cluster_tag])
inputM=as.matrix(inputM[,meta$cluster==cur_cluster])
meta=meta[meta$cluster==cur_cluster,]


dim(inputM)
inputM[1:10,1:10]
cell_list=colnames(inputM)
gene_list=rownames(inputM)


library("MAST")
library("lme4")

fData=data.frame(primerid=gene_list)

cData=data.frame(wellKey=cell_list)
colnames(meta)
length(fData)
length(cData)

sca=FromMatrix(inputM, cData, fData)
colData(sca)$cngeneson = as.numeric((colSums(inputM > 0))/nrow(inputM)) #from Chong and Paul
colData(sca)$diagnosis = meta$diagnosis
colData(sca)$ind = as.factor(meta$individual)
colData(sca)$age = as.numeric(meta$age)
colData(sca)$sex = as.factor(meta$sex)
colData(sca)$RIN = as.numeric(meta$RNA.Integrity.Number)
colData(sca)$PMI = as.numeric(meta$post.mortem.interval..hours.)
colData(sca)$region = as.factor(meta$region)
colData(sca)$Capbatch = as.factor(meta$Capbatch)
colData(sca)$Seqbatch = as.factor(meta$Seqbatch)
colData(sca)$ribo_perc = as.numeric(meta$RNA.ribosomal.percent)

colData(sca)

#do a permutation
if(perm_tag>0){
  #count cases and controls
  diag_info=paste0(colData(sca)$ind,":",colData(sca)$diagnosis)
  diag_kind=unique(diag_info)
  diag_kind=t(apply(as.matrix(diag_kind),1,function(x){return(unlist(strsplit(x,":")))}))
  
  #permute
  diag_kind[,2]=diag_kind[sample.int(nrow(diag_kind),nrow(diag_kind),replace=F),2]
  
  #match back to each individuals
  ind_index=match(colData(sca)$ind,diag_kind[,1])
  colData(sca)$diagnosis=as.factor(diag_kind[ind_index,2])
}

b=zlm(formula=~diagnosis + ( 1 | ind ) + cngeneson + age + sex + RIN + PMI + region + Capbatch + Seqbatch 
      + ribo_perc, sca=sca, method = "glmer", ebayes = F, silent=T)
saveRDS(b, paste0("../Data_PRJNA434002/zlm_",file_name,"_",cluster_tag,"_",perm_tag,".rds"))
bs=summary(b,logFC=TRUE,doLRT = paste0("diagnosis","Control"), level = 0.95)
saveRDS(bs, paste0("../Data_PRJNA434002/zlms_",file_name,"_",cluster_tag,"_",perm_tag,".rds"))

bs$datatable

sessionInfo()
q(save="no")


