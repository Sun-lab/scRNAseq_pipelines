#This code generate permutation labels for individuals.

B=2000
n=31

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

perm_order=matrix(ncol=B,nrow=n)

set.seed(7)
#set perm
for(ib in 1:B){
  perm_order[,ib]=sample.int(n,n)
}

saveRDS(perm_order,paste0("../Data_PRJNA434002/7.Result/ind_perm_order.rds"))