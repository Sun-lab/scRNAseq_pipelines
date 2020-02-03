#This code generate permutation labels for individuals.

B=2000
n=16

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

perm_order=matrix(ncol=B,nrow=n)

set.seed(77)
#set perm
for(ib in 1:B){
  perm_order[,ib]=sample.int(n,n)
}

saveRDS(perm_order,paste0("../GSE129788/ind_perm_order.rds"))