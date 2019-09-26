#this code do some analysis on the DCA output and comparied it with the input. To see if the result is reliable.

#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
file_name="3k10"


input_m=as.matrix(read.table(paste0("../Data_PRJNA434002/rawM",file_name,".csv")
                             ,stringsAsFactors = FALSE,header=TRUE,row.names = 1,sep=","))

#1.check the output from DCA
output_mean=as.matrix(read.table(paste0("../Data_PRJNA434002/res_raw_dca",file_name,"/mean.tsv")))
output_latent=as.matrix(read.table(paste0("../Data_PRJNA434002/res_raw_dca",file_name,"/latent.tsv")))
output_dropout=as.matrix(read.table(paste0("../Data_PRJNA434002/res_raw_dca",file_name,"/dropout.tsv")))
output_dispersion=as.matrix(read.table(paste0("../Data_PRJNA434002/res_raw_dca",file_name,"/dispersion.tsv")))


hit_index=match(rownames(output_mean),rownames(input_m))

input_m=input_m[hit_index,]#!all 0

input_m[1:30,1:10]
output_mean[1:30,1:10]

pdf(paste0("DCA_input_output_compare_",file_name,".pdf"),height=8,width=8)

hist(apply(output_mean,2,mean),col=rgb(0,0,1,0.3),breaks=50,ylim=c(0,5000), main="mean expression per cell, DCA input vs output")
hist(apply(input_m,2,mean),col=rgb(1,0,0,0.3),add=T,breaks=50)
legend("topright",legend=c("input","output"),col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch=15,bty="n")

hist(apply(output_mean,1,mean),col=rgb(0,0,1,0.3),breaks=50,ylim=c(0,2000), main="mean expression per gene, DCA input vs output")
hist(apply(input_m,1,mean),col=rgb(1,0,0,0.3),add=T,breaks=50)
legend("topright",legend=c("input","output"),col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch=15,bty="n")

hist(apply(output_mean,2,median),col=rgb(0,0,1,0.3),breaks=50,ylim=c(0,5000), main="median expression per cell, DCA input vs output")
hist(apply(input_m,2,median),col=rgb(1,0,0,0.3),add=T,breaks=50)
legend("topright",legend=c("input","output"),col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch=15,bty="n")

hist(apply(output_mean,1,median),col=rgb(0,0,1,0.3),breaks=50,ylim=c(0,2000), main="median expression per gene, DCA input vs output")
hist(apply(input_m,1,median),col=rgb(1,0,0,0.3),add=T,breaks=50)
legend("topright",legend=c("input","output"),col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch=15,bty="n")

cor_gene=matrix(ncol=1,nrow=nrow(input_m))
names(cor_gene)=rownames(input_m)
for(ig in 1:nrow(input_m)){
  cor_gene[ig]=cor(input_m[ig,],output_mean[ig,])
}
hist(cor_gene,main="correlation of mean expression, input vs output, per gene",breaks=50)

cor_gene[1:10]


plot(apply(output_mean,1,mean), apply(output_mean,1,var), pch=20, col=rgb(0,0,1,0.3), 
     xlab="mean", ylab="var",main="mean vs variance")
points(apply(input_m,1,mean), apply(input_m,1,var), pch=20, col=rgb(1,0,0,0.3), 
     xlab="mean", ylab="var")
legend("topleft",legend=c("input","output"),col=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3)),pch=15,bty="n")

dev.off()


sessionInfo()
q(save="no")


