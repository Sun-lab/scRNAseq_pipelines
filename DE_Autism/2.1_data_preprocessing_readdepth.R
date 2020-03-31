library("moments")

file_tag="PFC5k"
dataset_folder="Data_PRJNA434002"  #Data_PRJNA434002   MS

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Volumes/SpecialPass/fh_data/Data_PRJNA434002/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########input###############
#input diagnosis
#input phenotype
if(length(grep("PFC",file_tag))>0){
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    meta=readRDS(paste0("../",dataset_folder,"/meta_PFC.rds"))
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    meta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],"_PFC.rds"))
  }
}
if(length(grep("PFC",file_tag))==0){
  if(is.na(unlist(strsplit(file_tag,"k"))[2])){
    meta=read.table(paste0("../",dataset_folder,"/meta.tsv"),header = TRUE, sep = "\t")
  }
  if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
    meta=readRDS(paste0("../",dataset_folder,"/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
  }
}
#name match for MS samples
colnames(meta)[grep("sample",names(meta))]="individual"

cur_individual=unique(meta$individual)


#input counts
rawM=readRDS(paste0("../",dataset_folder,"/rawM",file_tag,".rds"))



############### Section III some basic calculations ###########
#require: results from section II



#calculate features per gene
rawM_nozero_median=matrix(ncol=1,nrow=nrow(rawM))
rawM_nozero_mean=matrix(ncol=1,nrow=nrow(rawM))
rawM_nozero_sd=matrix(ncol=1,nrow=nrow(rawM))
rawM_nozero_skewness=matrix(ncol=1,nrow=nrow(rawM))
rawM_zero_rate=matrix(ncol=1,nrow=nrow(rawM))
cell_4_each_gene=matrix(ncol=1,nrow=nrow(rawM))

for(i in 1:nrow(rawM)){
  cur_count=as.numeric(rawM[i,])
  cur_count=cur_count[cur_count>0]
  rawM_nozero_median[i]=median(cur_count)
  rawM_nozero_mean[i]=mean(cur_count)
  rawM_nozero_sd[i]=sd(cur_count)
  rawM_nozero_skewness[i]=skewness(cur_count)
  cell_4_each_gene[i]=length(cur_count)
  rawM_zero_rate[i]=(ncol(rawM)-cell_4_each_gene[i])/ncol(rawM)
}


#calculate features per cell
percell_gene_expressed_num=matrix(ncol=1,nrow=ncol(rawM))
percell_CDR=matrix(ncol=1,nrow=ncol(rawM))
percell_gene_expressed_count_total=matrix(ncol=1,nrow=ncol(rawM))
percell_gene_expressed_90perc=matrix(ncol=1,nrow=ncol(rawM))
percell_gene_expressed_99perc=matrix(ncol=1,nrow=ncol(rawM))

for(i in 1:floor(ncol(rawM)/1000)){
  cur_count=rawM[,((i-1)*1000+1):(i*1000)]
  cur_count=apply(cur_count,2,as.numeric)
  percell_gene_expressed_num[((i-1)*1000+1):(i*1000)]=apply(cur_count>0,2,sum)
  percell_CDR[((i-1)*1000+1):(i*1000)]=percell_gene_expressed_num[((i-1)*1000+1):(i*1000)]/nrow(cur_count) #CDR: cell expressed num rate, see MAST paper.
  percell_gene_expressed_count_total[((i-1)*1000+1):(i*1000)]=apply(cur_count,2,sum)
  percell_gene_expressed_90perc[((i-1)*1000+1):(i*1000)]=apply(cur_count,2,function(x){return(quantile(x,0.9))})
  percell_gene_expressed_99perc[((i-1)*1000+1):(i*1000)]=apply(cur_count,2,function(x){return(quantile(x,0.99))})
  print(i)
}
cur_count=rawM[,(i*1000+1):ncol(rawM)]
cur_count=apply(cur_count,2,as.numeric)
percell_gene_expressed_num[(i*1000+1):ncol(rawM)]=apply(cur_count>0,2,sum)
percell_CDR[(i*1000+1):ncol(rawM)]=percell_gene_expressed_num[(i*1000+1):ncol(rawM)]/nrow(cur_count) #CDR: cell expressed num rate, see MAST paper.
percell_gene_expressed_count_total[(i*1000+1):ncol(rawM)]=apply(cur_count,2,sum)
percell_gene_expressed_90perc[(i*1000+1):ncol(rawM)]=apply(cur_count,2,function(x){return(quantile(x,0.9))})
percell_gene_expressed_99perc[(i*1000+1):ncol(rawM)]=apply(cur_count,2,function(x){return(quantile(x,0.99))})


cell_features=cbind(percell_gene_expressed_num,percell_CDR,percell_gene_expressed_count_total,
                    percell_gene_expressed_90perc,percell_gene_expressed_99perc)
rownames(cell_features)=colnames(rawM)
colnames(cell_features)=c("percell_gene_expressed_num","percell_CDR","percell_gene_expressed_count_total",
                          "percell_gene_expressed_90perc","percell_gene_expressed_99perc")

saveRDS(cell_features,paste0("../Data_PRJNA434002/rawM",file_tag,"_percell_features.rds"))
cell_features=readRDS(paste0("../Data_PRJNA434002/rawM",file_tag,"_percell_features.rds"))

#load cell cluster for 31 individuals

cust_col=c("#FEE391","#67001F","#F7FCB9","#88419D","#CC4C02","#A6BDDB","#EF6548","#BFD3E6",
           "#7F2704","#DF65B0","#F16913","#74A9CF","#1D91C0","#8C6BB1","#9E9AC8","#FDAE6B",
           "#238B45","#CE1256","#000000","#810F7C","#D94801","#FE9929","#F768A1","#7FCDBB",
           "#41AB5D","#66C2A4","#7F0000","#E0ECF4","#FB6A4A","#662506","#00441B")

pdf(paste0("../Data_PRJNA434002/2.histogram_rawMatrix_cell_features_stat_",file_tag,".pdf"))
pairs(cell_features[,2:5], main = "Cell features comparison",
      pch = 21, bg = cust_col[match(meta$individual,unique(meta$individual))])
for(ic in 1:ncol(cell_features)){
  hist(cell_features[,ic],main=colnames(cell_features)[ic])
}
dev.off()

min_cell_4_each_gene=min(cell_4_each_gene)
min_gene_4_each_cell=min(gene_4_each_cell)

min_cell_4_each_gene
min_gene_4_each_cell

pdf(paste0("../Data_PRJNA434002/2.histogram_rawMatrix_stat_",file_tag,".pdf"))
hist(rawM[rawM>0],breaks=50)
hist(log10(cell_4_each_gene),breaks=50)
hist(log10(gene_4_each_cell),breaks=50)
hist(rawM_nozero_median,breaks=50)
hist(rawM_nozero_mean,breaks=50)
hist(rawM_nozero_sd,breaks=50)
hist(rawM_nozero_skewness,breaks=50)
hist(rawM_zero_rate,breaks=50,
     sub=paste0("0rate==1: ",
                round(sum(rawM_zero_rate==1)/length(rawM_zero_rate),3),
                ", 0rate>0.99: ",   
                round(sum(rawM_zero_rate>0.99)/length(rawM_zero_rate),3),
                ", 0rate>0.98: ", 
                round(sum(rawM_zero_rate>0.98)/length(rawM_zero_rate),3)))


dev.off()


#For read-depth per cell, calculate four measurements, 
#total read count, CDR, 90 percentile and 99 percentile. 
#Then draw the scatter plots of the four variables, using function pairs. Also draw histogram for each measurement. 


total_count_per_cell=apply(rawM,2,sum)




#############################Section IIII calculate the individual based info ####################
#calculate each individual read depth.each cell num,and the zero proportions
#require: results from section II


cur_individual=unique(meta$individual)
cell_num=matrix(ncol=1,nrow=length(cur_individual))
rownames(cell_num)=cur_individual
colnames(cell_num)="cell_num"
read_depth=matrix(ncol=1,nrow=length(cur_individual))
rownames(read_depth)=cur_individual
colnames(read_depth)="read_depth"
zero_rate_ind=matrix(nrow=nrow(rawM),ncol=length(cur_individual))
rownames(zero_rate_ind)=rownames(rawM)
colnames(zero_rate_ind)=cur_individual

for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=rawM[,meta$individual==cur_ind]
  cell_num[i_ind]=ncol(cur_ind_m)
  read_depth[i_ind]=sum(cur_ind_m,na.rm = TRUE)/cell_num[i_ind]*1000
  zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(x,na.rm = TRUE))})/cell_num[i_ind]
}

saveRDS(cell_num,paste0("../Data_PRJNA434002/rawM",file_tag,"_cell_num_per_ind.rds"))
saveRDS(read_depth,paste0("../Data_PRJNA434002/rawM",file_tag,"_read_depth_per_1Kcell_ind.rds"))
saveRDS(zero_rate_ind,paste0("../Data_PRJNA434002/rawM",file_tag,"_zero_rate_per_ind.rds"))

#For each gene, check the correlation (rank-based, spearman correlation) between 
#its average zero-inflation proportion across all the cells of on individual and 
# the  total read depth of this gene of each individual?).  
# Summarize such correlations across genes by a histogram

pdf(paste0("../Data_PRJNA434002/2.histogram_rawM",file_tag,"_correlation_zero_rate_vs_read_depth_per_ind.pdf"))
cor_zero_rate_vs_read_depth_ind=apply(zero_rate_ind,1,function(x){return(cor(x,read_depth,method="spearman"))})
hist(cor_zero_rate_vs_read_depth_ind,sub=paste0("non-NA correlations: ",sum(is.na(cor_zero_rate_vs_read_depth_ind))," out of ",length(cor_zero_rate_vs_read_depth_ind)))
saveRDS(cor_zero_rate_vs_read_depth_ind,paste0("../Data_PRJNA434002/rawM",file_tag,"_cor_zero_rate_vs_read_depth_ind.rds"))
dev.off()

#conclusion: the NA correlation is due to the 100% percent of zeros.