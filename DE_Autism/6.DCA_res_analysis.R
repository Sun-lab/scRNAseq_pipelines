#this code do some analysis on the DCA output and comparied it with the input. To see if the result is reliable.
library("ggplot2")

#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
setwd("/Users/grasscliff/Desktop/fh/1.Testing_scRNAseq/")

file_name="3k10"

if(is.na(unlist(strsplit(file_name,"k"))[2])){
  meta=readRDS(paste0("../Data_PRJNA434002/meta.rds"))
}
if(!is.na(unlist(strsplit(file_name,"k"))[2])){
  meta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_name,"k"))[2],".rds"))
}



input_m=as.matrix(read.table(paste0("../Data_PRJNA434002/rawM",file_name,".csv")
                             ,stringsAsFactors = FALSE,header=TRUE,row.names = 1,sep=","))

#1.check the output from DCA
output_mean=as.matrix(read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_name,"/mean.tsv")))
output_latent=as.matrix(read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_name,"/latent.tsv")))
output_dropout=as.matrix(read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_name,"/dropout.tsv")))
output_dispersion=as.matrix(read.table(paste0("../Data_PRJNA434002/res_dca_rawM",file_name,"/dispersion.tsv")))


hit_index=match(rownames(output_mean),rownames(input_m))

input_m=input_m[hit_index,]#!all 0

input_m[1:30,1:10]
output_mean[1:30,1:10]
cur_individual=unique(meta$individual)

pdf(paste0("DCA_input_output_compare_",file_name,".pdf"),height=8,width=8)

#1. plot histogram of a few genes of all cells.

set.seed(94)
#ramdomly selected 10 genes
select_flag=sample.int(100,10,replace = FALSE)
selected_output_mean=output_mean[select_flag,]
selected_input=input_m[select_flag,]
#1. plot histogram of all 10 genes of all cells
hist(log(selected_input),
     xlim=c(log(min(selected_input[selected_input>0],selected_output_mean[selected_output_mean>0])),
            log(max(selected_input[selected_input>0],selected_output_mean[selected_output_mean>0]))),
     col=rgb(0,0,1,0.3),main=paste0("expression hist rawCount vs DCACount, 10 genes for all cells"))
hist(log(selected_output_mean),add=T,col=rgb(1,0,0,0.3))
legend("topright",c("raw_count","dca_simulated count"),col=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3)),pch=15)

#2. plot histogram of each gene of all cells

for(i_g in 1:length(select_flag)){
  cur_selected_input=selected_input[i_g,]
  cur_selected_output_mean=selected_output_mean[i_g,]
  
  hist(log(cur_selected_input),
       xlim=c(log(min(cur_selected_input[cur_selected_input>0],cur_selected_output_mean[cur_selected_output_mean>0])),
              log(max(cur_selected_input[cur_selected_input>0],cur_selected_output_mean[cur_selected_output_mean>0]))),
       col=rgb(0,0,1,0.3),main=paste0("expression hist rawCount vs DCACount, gene ",select_flag[i_g]," for all cells"))
  hist(log(cur_selected_output_mean),add=T,col=rgb(1,0,0,0.3))
  legend("topright",c("raw_count","dca_simulated count"),col=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3)),pch=15)
}

#3 calculate the proportion of 0’s and then draw a scatter plot the proportion of 0’s across all the genes(3c). 
zero_rate_gene_input=apply(selected_input>0,1,function(x){return(1-sum(x,na.rm=TRUE)/length(x))})
zero_rate_gene_output=apply(selected_output_mean>0,1,function(x){return(1-sum(x,na.rm=TRUE)/length(x))})
plot(zero_rate_gene_input,zero_rate_gene_output,type="p",main=paste0("0 rate rawCount vs DCACount, all gene all cells"))
hist(zero_rate_gene_input,breaks=20,main=paste0("hist of rawCount, all gene all cells"))




for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  
  
  selected_output_mean=output_mean[select_flag,meta$individual==cur_ind]
  selected_input=input_m[select_flag,meta$individual==cur_ind]
  
  #4. plot histogram of all 10 genes of each individual
  hist(log(selected_input),
       xlim=c(log(min(selected_input[selected_input>0],selected_output_mean[selected_output_mean>0])),
              log(max(selected_input[selected_input>0],selected_output_mean[selected_output_mean>0]))),
       col=rgb(0,0,1,0.3),main=paste0("expression hist rawCount vs DCACount, 10 genes, cells",i_ind))
  hist(log(selected_output_mean),add=T,col=rgb(1,0,0,0.3))
  legend("topright",c("raw_count","dca_simulated count"),col=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3)),pch=15)
  
  #5. plot histogram of each gene of each individual
  
  for(i_g in 1:nrow(selected_output_mean)){
    cur_selected_input=selected_input[i_g,]
    cur_selected_output_mean=selected_output_mean[i_g,]
    
    hist(log(cur_selected_input),
         xlim=c(log(min(cur_selected_input[cur_selected_input>0],cur_selected_output_mean[cur_selected_output_mean>0])),
                log(max(cur_selected_input[cur_selected_input>0],cur_selected_output_mean[cur_selected_output_mean>0]))),
         col=rgb(0,0,1,0.3),main=paste0("expression hist rawCount vs DCACount, gene ",select_flag[i_g],", cells",i_ind))
    hist(log(cur_selected_output_mean),add=T,col=rgb(1,0,0,0.3))
    legend("topright",c("raw_count","dca_simulated count"),col=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3)),pch=15)
  }
  
  #6 calculate the proportion of 0’s and then draw a scatter plot the proportion of 0’s across all the genes(3c) for each individual
  zero_rate_gene_input=apply(selected_input>0,1,function(x){return(1-sum(x,na.rm=TRUE)/length(x))})
  zero_rate_gene_output=apply(selected_output_mean>0,1,function(x){return(1-sum(x,na.rm=TRUE)/length(x))})
  plot(zero_rate_gene_input,zero_rate_gene_output,type="p",main=paste0("0 rate rawCount vs DCACount, all gene, cells",i_ind))
  hist(zero_rate_gene_input,breaks=20,main=paste0("hist of rawCount, all gene, cells",i_ind))
}




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


