#this code is used fr some basic processing of dataset GSE129788
#This dataset is single nucleus of paper "Single-cell transcriptomic profiling of the aging mouse brain (house mouse)". It contains the brain cells of 8 young mouse and 8 old mouse. of 14699 genes and 37070 cells. The data are log transformed.

#This data are pre log transformed with y=log(x+1), and then scaled to 10K transcripts per cell

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

exprM=read.table("../GSE129788/expression_Aging_mouse_brain_portal_data_updated.txt",header=TRUE,row.names=1)

meta=read.table("../GSE129788/meta_Aging_mouse_brain_portal_data.txt",header = TRUE, sep = "\t")


##########count processing#############
dim(exprM)
exprM[1:10,1:10]
quantile(exprM[exprM>0])
#       0%       25%       50%       75%      100%
#0.2957098 1.2729451 1.5769655 1.9853241 8.8877984

#detect max account percentage from the original count
max_count_percentage=apply(exp(exprM)-1,2,function(x){return(max(x)/sum(x))})
quantile(max_count_percentage,prob=c(0:20)/20)
# 0%          5%         10%         15%         20%         25%
#   0.007982522 0.024222012 0.029780691 0.033646277 0.036347304 0.038786362
# 30%         35%         40%         45%         50%         55%
#   0.041118421 0.043518434 0.046133508 0.048793657 0.051737757 0.055243746
# 60%         65%         70%         75%         80%         85%
#   0.059052007 0.063735664 0.069059104 0.075579187 0.083704498 0.094849502
# 90%         95%        100%
# 0.111728923 0.143087785 0.724205516
sort(max_count_percentage,decreasing = TRUE)[1:10]


#generate another 10% percentage cells with the 5000 lowest 0 rate genes, for tSNE plots
exprM_zero_rate=apply(exprM==0,1,function(x){return(sum(x)/length(x))})
o_zero_rate=order(exprM_zero_rate)
mt_index=grep("mt-",rownames(exprM))

exprM1k=exprM[o_zero_rate[1:1000],]
exprM1k[1:10,1:10]
dim(exprM1k)



############meta processing #add individual######################
meta=meta[-1,]
meta$nGene=as.numeric(meta$nGene)
meta$nUMI=as.numeric(meta$nUMI)
table(meta$all_cells_by_age)
table(meta$cell_type)
individual=sapply(as.character(meta$NAME),function(x){return(unlist(strsplit(x,"_"))[6])})
cluster=paste0(meta$cell_class,"_",meta$cell_type)
phenotype=meta$all_cells_by_age
phenotype=((phenotype!="2-3mo")+0) #old==1 young==0

meta=cbind(meta,individual,cluster,phenotype)
meta[1:10,]

table(cluster)

#ind level info
cur_individual=unique(individual)
cur_cluster=unique(cluster)
CDR=matrix(0,ncol=length(cur_cluster),nrow=nrow(meta))
rownames(CDR)=meta$NAME
colnames(CDR)=cur_cluster
CDR_ind=matrix(0,ncol=length(cur_cluster),nrow=length(cur_individual))
rownames(CDR_ind)=cur_individual
colnames(CDR_ind)=cur_cluster

for(i_cluster in 1:length(cur_cluster)){
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]
    cur_c=cur_cluster[i_cluster]
    cur_ind_m=NULL
    #fit org
    cur_ind_m=exprM1k[,meta$individual==cur_ind & meta$cluster==cur_c]
    if(length(cur_ind_m)>0){
      cur_CDR=sum(cur_ind_m>0,na.rm = TRUE)/sum(cur_ind_m>-1,na.rm = TRUE)
      CDR[meta$individual==cur_ind,i_cluster]=cur_CDR
      CDR_ind[i_ind,i_cluster]=cur_CDR
    }
  }
}

saveRDS(CDR,"../GSE129788/CDR1k_per_cluster.rds")
saveRDS(CDR_ind,"../GSE129788/ind_CDR1k_per_cluster.rds")



saveRDS(exprM1k,"../GSE129788/exprM1k.rds")
saveRDS(meta,"../GSE129788/meta1k.rds")
sessionInfo()
q(save="no")
