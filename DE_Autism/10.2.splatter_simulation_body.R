#this code used for the

# Simulate scRNA-seq data using splatter
# Using their default simulator for now. We can consider simulate 5000 genes in 40 individuals, with 20 cases vs. 20 controls. We can set DE for 500 genes in terms of mean expression difference, and another 500 genes DE in terms of variance difference. In the simulation, try to set the total read-depth across samples to be the same. 
# As an initial analysis, do not add any covariates. 
# For each gene, calculate density, and then JSD across samples, and then use PERMANOVA to calculate p-value for each gene. 
# Also calculate p-value for each gene using MAST-mixed effect model (less priority for now). 
# Collapse gene expression across cells per individual, then run DESeq2 for differential expression testing. 
# Type I error is the proportion of those 4000 equivalently expressed genes with p-value smaller than 0.05. 
# Power1, the proportion of the first 500 genest with p-value < 0.05
# Power2, the proportion of the next 500 genest with p-value < 0.05


# ##### a. Simulation Plan (method 4 was choisen)#############

# ####### method1: 1 batch, 40 groups, de.prob=rep(0.5,0.5),  "de.facScale", .25,  "de.facLoc", 1
# eg: DEseq2 methods:from https://bioconductor.github.io/BiocWorkshops/rna-seq-data-analysis-with-deseq2.html
#   library("splatter")
# params <- newSplatParams()
# params <- setParam(params, "de.facLoc", 1) 
# params <- setParam(params, "de.facScale", .25)
# params <- setParam(params, "dropout.type", "experiment")
# params <- setParam(params, "dropout.mid", 3)
# set.seed(1)
# sim <- splatSimulate(params, group.prob=c(.5,.5), method="groups")
# 
# 
# pros:
# easy to operated,have examples
# 
# cons:
# 1. no individual differences--batch effects
# 2. hard to characterize the differences.
# 
# 
# ####### method2: times simulation,each is 1 batch, 20 batch have group.prob=1(DE or DV), 20 batch have group.prob=0
# 
# pros:
# easy to operated
# 
# cons:
# can't fix the genes who's differential expressed. Need some further steps to putting them together...
# 
# 
# ####### method3: 40 batch, each has 3 group de.prob=c(0.33,0.33,0.33), 
#                       group1 change mean "de.facScale", .25,  "de.facLoc", 1
#                       group2 change variance  
#                       group3 Set default
#                       
# Then doing the DE analysis with 
# (1).mean comparison: first 20 group1 cell vs latter 20 group3 cell
# (2).mean comparison: first 20 group2 cell vs latter 20 group3 cell
# 
# 
# Based on all of those, I plan to take method 3 for the simulation.
# Something remaining:
# (1)make sure where is the influence place for "de.facScale",  "de.facLoc", what it influence on the gamma distribution mean.rate and mean.shape.
# (2)make sure the DE genes are consistence between batches in your simulation.
# 
#  
#  
# 
# ####### method4: generate 6 sets, each with 20 batch
#After carefully implement of method 3(chosen from the 1,2,3 method),we found that the parameter de.prob/de.fracLoc/de.fracScale are all applied on the place of outlier location, which means they are not able for the adjustment of the variance without non-influence on the mean.

#personally, in this situation, I would prefer to use Lun2 method, which is more straight forward.



#based on that, we plan to simulate 6 parts separately
#               case         |       Control
#MeanDE Genes                |
#VarDE  Genes                |
#NochangeGenes               |

#according to http://www.math.wm.edu/~leemis/chart/UDR/PDFs/Gammapoisson.pdf. The gamma-poission distribution will have the 
#mean= alpha*beta
#variance=alpha*beta+alpha^2*beta
#Thus we will adjust the alphaand beta to make our changes.

##### Start simulation from here #############

setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#functions##########
source("./Command/7.0_ZINB_fit_functions.R")


library("splatter")
library("scater")
library("ggplot2")

nGeneMean=500
nGeneVar=500
nGeneBlank=4000
ncase=20
nctrl=20
ncell=100

file_tag=1
covariate_flag=NA

##############functions#################
####this function returns the parameters for changing mean and variance of the gamma-poission distribution
# meanCase: a'b'=r_m*ab and (1+a')a'b'=(1+a)ab
#(1+a')r_m*ab=(1+a)ab
#(1+a')r_m=(1+a)

#a'=(1+a)/r_m-1
#b'=r_m*ab/a'=r_m*ab/((1+a)/r_m-1)

# var case:  (1+a')a'b'=r_v* (1+a)ab and a'b'=ab
# Then:
# 1+a'=r_v*(1+a)
# So:
#a'=r_v*(1+a)-a
#b'=ab/a'=ab/(r_v*(1+a)-a)

#complex case: a'b'=r_m*ab and (1+a')a'b'=r_v* (1+a)ab
#then
#a'=(a+1)*r_v/r_m - 1
#b'=a*b*r_m*r_m/((a+1)*r_v - r_m)

calc_gapois_param=function(a,b,r_m=1,r_v=1){
  #alpha=shape
  #beta=rate
  alpha=(a+1)*r_v/r_m - 1
  beta=a*b*r_m*r_m/((a+1)*r_v - r_m)
  return(c(alpha,beta))
}

############
#1.meanCtrl
params=newSplatParams(batchCells=rep(ncell,nctrl), nGenes = nGeneMean)
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
sim_meanCtrl=splatSimulate(params)


#2.varCtrl
params=newSplatParams(batchCells=rep(ncell,nctrl), nGenes = nGeneVar)
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
sim_varCtrl=splatSimulate(params)


#3.blankCtrl
params=newSplatParams(batchCells=rep(ncell,nctrl), nGenes = nGeneBlank)
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
sim_blankCtrl=splatSimulate(params)

#1.meanCase 
r_mean=1.5

params=newSplatParams()
ref=getParams(params, c("mean.shape", "mean.rate"))
cur=calc_gapois_param(ref$mean.shape,ref$mean.rate,r_m=r_mean)
params=newSplatParams(batchCells=rep(ncell,ncase), nGenes = nGeneMean)
params=setParams(params, update = list(mean.shape = cur[1],
                                       mean.rate = cur[2]))
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
sim_meanCase=splatSimulate(params)

#3.varCase 
r_var=4

params=newSplatParams()
ref=getParams(params, c("mean.shape", "mean.rate"))
cur=calc_gapois_param(ref$mean.shape,ref$mean.rate,r_v=r_var)
params=newSplatParams(batchCells=rep(ncell,ncase), nGenes = nGeneVar)
params=setParams(params, update = list(mean.shape = cur[1],
                                       mean.rate = cur[2]))
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
sim_varCase=splatSimulate(params)

#6.blankCase
params=newSplatParams(batchCells=rep(ncell,ncase), nGenes = nGeneBlank)
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
sim_blankCase=splatSimulate(params)

op <- par(mfrow = c(4, 3), pty = "s")   
hist(log(apply(counts(sim_blankCtrl),1,mean)))
hist(log(apply(counts(sim_meanCtrl),1,mean)))
hist(log(apply(counts(sim_varCtrl),1,mean)))

hist(log(apply(counts(sim_blankCase),1,mean)),col=rgb(1,0,0,0.3))
hist(log(apply(counts(sim_meanCase),1,mean)),col=rgb(1,0,0,0.3))
hist(log(apply(counts(sim_varCase),1,mean)),col=rgb(1,0,0,0.3))

hist(log(apply(counts(sim_blankCtrl),1,var)))
hist(log(apply(counts(sim_meanCtrl),1,var)))
hist(log(apply(counts(sim_varCtrl),1,var)))

hist(log(apply(counts(sim_blankCase),1,var)),col=rgb(1,0,0,0.3))
hist(log(apply(counts(sim_meanCase),1,var)),col=rgb(1,0,0,0.3))
hist(log(apply(counts(sim_varCase),1,var)),col=rgb(1,0,0,0.3))

par(op)


counts(sim_meanCase)[1:10,1:10]

head(rowData(sim_meanCase))
head(colData(sim_meanCase))
names(assays(sim_meanCase))

individual=c(as.numeric(gsub("Batch", "", colData(sim_meanCase)[,2])),as.numeric(gsub("Batch", "", colData(sim_meanCase)[,2]))+ncase)
phenotype=c(rep(1,ncell*ncase),rep(0,ncell*nctrl))

cell_id=paste0("cell",1:(ncell*ncase+ncell*nctrl))
gene_id=paste0("gene",1:(nGeneMean+nGeneVar+nGeneBlank))

de_mean_flag=c(rep(1,nGeneMean),rep(0,nGeneVar),rep(0,nGeneBlank))
de_var_flag=c(rep(0,nGeneMean),rep(1,nGeneVar),rep(0,nGeneBlank))
de_blank_flag=c(rep(0,nGeneMean),rep(1,nGeneVar),rep(0,nGeneBlank))

sim_matrix=as.matrix(rbind(cbind(counts(sim_meanCase),counts(sim_meanCtrl)),
          cbind(counts(sim_varCase),counts(sim_varCtrl)),
          cbind(counts(sim_blankCase),counts(sim_blankCtrl))))

rownames(sim_matrix)=gene_id
colnames(sim_matrix)=cell_id

cellsum=apply(sim_matrix,2,sum)
genesum=apply(sim_matrix,1,sum)
CDR=apply(sim_matrix>0,2,sum)/nrow(sim_matrix)
meta=data.frame(cbind(cell_id,individual,phenotype,cellsum,CDR))

for(i in 501:510){
  hist(sim_matrix[i,1:100])
  hist(sim_matrix[i,2001:2100])
}


#start to calculation:


cur_individual=unique(meta$individual)
cell_num=matrix(ncol=1,nrow=length(cur_individual))
rownames(cell_num)=cur_individual
colnames(cell_num)="cell_num"
read_depth=matrix(ncol=1,nrow=length(cur_individual))
rownames(read_depth)=cur_individual
colnames(read_depth)="read_depth"
zero_rate_ind=matrix(nrow=nrow(sim_matrix),ncol=length(cur_individual))
rownames(zero_rate_ind)=rownames(sim_matrix)
colnames(zero_rate_ind)=cur_individual

for(i_ind in 1:length(cur_individual)){
  cur_ind=cur_individual[i_ind]
  #fit org
  cur_ind_m=sim_matrix[,meta$individual==cur_ind]
  cell_num[i_ind]=ncol(cur_ind_m)
  read_depth[i_ind]=sum(cur_ind_m,na.rm = TRUE)/cell_num[i_ind]*1000
  zero_rate_ind[,i_ind]=apply(cur_ind_m==0,1,function(x){return(sum(x,na.rm = TRUE))})/cell_num[i_ind]
}

hist(read_depth)
#almost no need to adjust library size.




#fit sim data by individual
if(!is.na(covariate_flag)){
  #logsum_count=log(apply(sim_matrix,2,sum))
  quantile99=log(apply(sim_matrix,2,function(x)return(quantile(x,0.99)+1)))
  covariate=as.matrix(quantile99)
  pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind_raw_",covariate_flag,"_",file_tag,".pdf"),height = 4,width = 6)
}
if(is.na(covariate_flag)){
  pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind_raw_",file_tag,".pdf"),height = 4,width = 6)
}
fit_ind_org=array(dim=c(nrow(sim_matrix),length(cur_individual),3),
                  dimnames = list(rownames(sim_matrix),cur_individual,c("logmean","dispersion","dropout_rate")))

for(i_g in 1:nrow(sim_matrix)){
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]
    #fit org
    cur_org_ind=sim_matrix[i_g,meta$individual==cur_ind]
    
    if(!is.na(covariate_flag)){
      cur_covariate=covariate[meta$individual==cur_ind,]
      fit_ind_org[i_g,i_ind,]=fit_nbzinb(cur_org_ind,cur_covariate)
    }
    if(is.na(covariate_flag)){
      fit_ind_org[i_g,i_ind,]=fit_nbzinb(cur_org_ind)
    }
    
    if(i_g<=10 & i_ind<=5){
      cur_org_ind=data.frame(cur_org_ind)
      colnames(cur_org_ind)="raw_count"
      ggplot(cur_org_ind, aes(x=raw_count),stat="count") + geom_histogram(fill="lightblue")+
        labs(title=paste0("Histogram of rawcount, ",rownames(sim_matrix)[i_g]," of ",cur_ind,", ",cur_cluster),x="Count", y = "Frequency")
      #+theme_classic()
    }
    
  }
  print(i_g)
}
dev.off()











