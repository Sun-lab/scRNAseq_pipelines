#this code used for the simulation

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
#setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")

sim_folder="sim_v3"
library("emdbook")

nGeneMean=300
nGeneVar=300
nGeneDisp=300
nGeneMult=300
nGeneBlank=1800
nGeneTotal=nGeneMean+nGeneVar+nGeneMult+nGeneDisp+nGeneBlank
ncase=100   #!!!!!!!!Different from Wei codes, this is inside the body.
nctrl=100   #!!!!!!!!Different from Wei codes, this is inside the body.
ncell=200
nall=ncase+nctrl

HET=0 #0/1 label: if HET==0, generate case cells from completely case condition
#           if HET==1, generate 1/2 case cells from completely case condition

i_mean=1:nGeneMean
i_var=(nGeneMean+1):(nGeneMean+nGeneVar)
i_disp=(nGeneMean+nGeneVar+1):(nGeneMean+nGeneVar+nGeneDisp)
i_mult=(nGeneMean+nGeneVar+nGeneDisp+1):(nGeneMean+nGeneVar+nGeneDisp+nGeneMult)
i_blank=(nGeneMean+nGeneVar+nGeneDisp+nGeneMult+1):nGeneTotal

#first control, then case
i_ctrl=1:(nctrl*ncell)
i_case=(nctrl*ncell+1):((ncase+nctrl)*ncell)

case_modify_flag=rbinom(nGeneTotal,1,0.5)
i_case_modify=which(case_modify_flag==1) #randomlized settings, switching cases and control samples to make sure the library sizes the same

##############functions#################



####this function returns the parameters for changing mean and variance of the gamma-poission distribution
#require:  r_m/r_v <1+mu/theta
calc_nb_param=function(mu,theta,r_m=1,r_v=1){
  #mu=mean
  #theta=overdispersion
  mu2=r_m*mu
  theta2=theta*mu*r_m/(mu*r_v*r_m+(r_v*r_m-1)*theta)
  return(c(mu2,theta2))
}


#require:  r_m/r_v <1+mu/theta
calc_zinb_param=function(mu,theta,drop=0,r_m=1,r_v=1){
  #mu=mean
  #theta=overdispersion
  mu2=r_m*mu
  theta2=theta*mu*r_m/(mu*r_v*r_m+(r_v*r_m-1)*theta + (r_v-1)*r_m*mu*drop*theta)
  if(theta2 < 0){ stop("negative theta2\n") }
  return(c(mu2,theta2))
}

# the function calc_nb_param_var returns the parameters for changing 
# variance of the negative binomial distribution. 
# to make sure theta is larger than 0, r_v should be > 1. 

calc_nb_param_var = function(mu, theta, r_v) {
  theta2 = theta  * mu / (mu * r_v + (r_v - 1) * theta)
  if(theta2 < 0){ stop("negative theta2\n") }
  theta2
}


# This function calcuate the parameters of 3rd NB/ZINB distribution according to a mixtured 50%-50% 
# Suppose Z have equaliy probability to falls into two NB distribution ZINB(mu1,size1,drop1), ZINB(mu2,size2,drop2). Calculate Z's ZINB(mu3,size3,drop3)

#forward
#input:
# mu1,mu2,size1,size2,drop2,drop1 
# or 
# mean1,mean2,size1,size2,drop2,drop1
#output: rbind(c(mu1,size1,drop1),c(mu2,size2,drop2),c(mu3,size3,drop3))

cal_nbzinb_param_multimodality_forward=function(mu1=NA,size1=4,drop1=0,mu2=NA,size2=4,drop2=0,mean1=NA,mean2=NA){
  drop3=(drop1+drop2)/2
  if(is.na(mu1)||is.na(mu2)){
    mean3=(mean1+mean2)/2
    
    mu1=mean1/(1-drop1)
    mu2=mean2/(1-drop2)
    mu3=mean3/(1-drop3)
  }
  if(is.na(mean1)||is.na(mean2)){
    mu3=(mu1+mu2)/2
    
    mean1=(1-drop1)*mu1
    mean2=(1-drop2)*mu2
    mean3=(1-drop3)*mu3
  }
  
  size3=1/(((1/4*(mean1-mean2)^2+1/2*(mean1*(1+mu1*(drop1+1/size1))+mean2*(1+mu2*(drop2+1/size2))) )/mean3-1 )/mu3-drop3)
  res=rbind(c(mu1,size1,drop1),c(mu2,size2,drop2),c(mu3,size3,drop3))
  return(res)
}

#rerward
#input:
# mu3,size3,drop3,change_proportion
# or 
# mean3,size3,drop3,change_proportion

#Here we have to have some support
#output: size #size==theta in the notation
cal_nbzinb_param_multimodality_reward=function(mu=NA,size,drop,change_proportion=0.8, mean3=NA){
  
  if(is.na(mu)){
    mu=mean3/(1-drop)
  }
  if(is.na(mean3)){
    mean3=(1-drop)*mu
  }
  mean1=(1+change_proportion)*mean3
  mean2=(1+change_proportion)*mean3
  mu1=(1+change_proportion)*mu
  mu2=(1+change_proportion)*mu
  
  size2=(mean1*mu1+mean2*mu2)/( 2*mean3*mu*(drop+1/size)-1/2*(mean1-mean2)^2  - drop*(mean1*mu1+mean2*mu2))
  if(size2 < 0){ stop("negative size2\n") }
  return(size2)
}

#enlarge: using the input as the smaller distribution
#input:
# mu2,size2,drop2,change_proportion
# or 
# mean2,size2,drop2,change_proportion

#Here we have to have some support
#output: size #size==theta in the notation
cal_nbzinb_param_multimodality_enlarge=function(mu=NA,size,drop,change_proportion=0.5, mean2=NA){
  if(is.na(mu)){
    mu=mean2/(1-drop)
  }
  if(is.na(mean2)){
    mean2=(1-drop)*mu
  }
  mean3=mean2/(1-change_proportion)
  mean1=(1+change_proportion)*mean3
  mu3=mu/(1-change_proportion)
  mu1=(1+change_proportion)*mu3
  t=mean3-mean2
  m=mean3
  size3=m^2/(m^2/size + t^2/size + t^2)
  return(size3)
}
# # ############ Simulation data with Method 5 (Pre data analysis)###############################
# 
# #generate the data by ZINB, based on given parameters...
# #input_file_tag="dca_PFC_all"
# t_mean = read.table(paste0("../Data_PRJNA434002/dca_PFC_all/mean.tsv"), 
#                     sep = "\t", header = TRUE)
# 
# t_dispersion = read.table(paste0("../Data_PRJNA434002/dca_PFC_all/dispersion.tsv"),
#                     sep = "\t", header = TRUE)
# 
# t_dropout = read.table(paste0("../Data_PRJNA434002/dca_PFC_all/pi.tsv"),
#                     sep = "\t", header = TRUE)
# 
# 
# dim(t_mean)
# t_mean[1:2,1:5]
# 
# dim(t_dispersion)
# t_dispersion[1:2,1:5]
# 
# dim(t_dropout)
# t_dropout[1:2,1:5]
# 
# rownames(t_mean)=t_mean[, 1]
# t_mean=t_mean[,-1]
# t_mean=apply(t_mean,2,as.numeric)
# t_mean=as.matrix(t_mean)
# 
# rownames(t_dispersion)=t_dispersion[, 1]
# t_dispersion=t_dispersion[,-1]
# t_dispersion=apply(t_dispersion,2,as.numeric)
# t_dispersion=as.matrix(t_dispersion)
# 
# rownames(t_dropout)=t_dropout[, 1]
# t_dropout=t_dropout[,-1]
# t_dropout=apply(t_dropout,2,as.numeric)
# t_dropout=as.matrix(t_dropout)
# 
# 
# saveRDS(t_mean,paste0("../Data_PRJNA434002/dca_PFC_all/mean.rds"))
# saveRDS(t_dispersion,paste0("../Data_PRJNA434002/dca_PFC_all/dispersion.rds"))
# saveRDS(t_dropout,paste0("../Data_PRJNA434002/dca_PFC_all/dropout.rds"))

# dim(t_mean)
# t_mean[1:2,1:5]
# 
# dim(t_dispersion)
# t_dispersion[1:2,1:5]
# 
# dim(t_dropout)
# t_dropout[1:2,1:5]
# 
# 
# #first threshold dropout
# mid_dropout=apply(t_dropout,1,median)
# mid_dropout=as.numeric(mid_dropout)
# quantile(mid_dropout)
# dropout_flag=(mid_dropout<0.2)
# sum(dropout_flag)
# 
# t_mean=t_mean[dropout_flag==1,]
# t_dispersion=t_dispersion[dropout_flag==1,]
# t_dropout=t_dropout[dropout_flag==1,]
# mid_dropout=mid_dropout[dropout_flag==1]
# 
# 
# #then sort mean
# mid_mean=apply(t_mean,1,median)
# mid_mean=as.numeric(mid_mean)
# quantile(mid_mean)
# order_mid_mean=order(mid_mean,decreasing = TRUE)
# 
# #refine to 1k, 3k, 5k
# sub_mean=t_mean[order_mid_mean[1:5000],]
# sub_dispersion=t_dispersion[order_mid_mean[1:5000],]
# sub_dropout=t_dropout[order_mid_mean[1:5000],]
# 
# mid_mean=mid_mean_total[order_mid_mean[1:5000]]
# mid_dispersion=mid_dispersion_total[order_mid_mean[1:5000]]
# mid_dropout=mid_dropout_total[order_mid_mean[1:5000]]
# 
# #mid_mean=apply(sub_mean,1,median)
# #mid_dispersion=apply(sub_dispersion,1,median)
# #mid_dropout=apply(sub_dropout,1,median)
# 
# mid_mean=as.numeric(mid_mean)
# mid_dispersion=as.numeric(mid_dispersion)
# mid_dropout=as.numeric(mid_dropout)
# 
# quantile(mid_mean)
# quantile(mid_dispersion)
# quantile(mid_dropout)
# 
# log_mean_sample=(log(as.numeric(mid_mean)))
# log_disp_sample=(log(as.numeric(mid_dispersion)))
# logit_drop_sample=(log(as.numeric(mid_dropout)/(1-as.numeric(mid_dropout))))
# 
# # Pram_log_mean=rnorm(nGeneMean+nGeneVar+nGeneDisp+nGeneMult+nGeneBlank,mean = (log_mean_sample),sd=sd(log_mean_sample))
# # Pram_log_disp=rnorm(nGeneMean+nGeneVar+nGeneDisp+nGeneMult+nGeneBlank,mean = (log_disp_sample),sd=sd(log_disp_sample))
# # Pram_logit_drop=rnorm(nGeneMean+nGeneVar+nGeneDisp+nGeneMult+nGeneBlank,mean = (logit_drop_sample),sd=sd(logit_drop_sample))
# 
# sample_data=cbind(log_mean_sample,log_disp_sample,logit_drop_sample)
# sample_mean=apply(sample_data,2,mean)
# sample_var=apply(sample_data,2,var)
# 
# cov_matrix=matrix(NA,3,3)
# for(i in 1:3){
#   for(j in 1:3){
#     cov_matrix[i,j]=cov(sample_data[,i],sample_data[,j])
#   }
# }


# saveRDS(cov_matrix,paste0("../Data_PRJNA434002/dca_PFC_all/5k_cov_matrix.rds"))
# saveRDS(sub_mean,paste0("../Data_PRJNA434002/dca_PFC_all/5k_mean.rds"))
# saveRDS(sub_dispersion,paste0("../Data_PRJNA434002/dca_PFC_all/5k_dispersion.rds"))
# saveRDS(sub_dropout,paste0("../Data_PRJNA434002/dca_PFC_all/5k_dropout.rds"))

############# Simulation data with Method 5 (Pre data analysis)###############################

#generate the data by ZINB, based on given parameters...
input_file_tag="3k" 

if(!file.exists(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_data_",sim_folder,".rds"))){
  cov_matrix=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_cov_matrix.rds"))
  t_mean=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_mean.rds"))
  t_dispersion=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_dispersion.rds"))
  t_dropout=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_dropout.rds"))
  
  
  dim(t_mean)
  t_mean[1:2,1:5]
  
  dim(t_dispersion)
  t_dispersion[1:2,1:5]
  
  dim(t_dropout)
  t_dropout[1:2,1:5]
  
  mid_mean=apply(t_mean,1,median)
  mid_dispersion=apply(t_dispersion,1,median)
  mid_dropout=apply(t_dropout,1,median)
  log_mean_sample=(log(as.numeric(mid_mean)))
  log_disp_sample=(log(as.numeric(mid_dispersion)))
  logit_drop_sample=(log(as.numeric(mid_dropout)/(1-as.numeric(mid_dropout))))
  
  #Pram_log_mean=rnorm(nGeneMean+nGeneVar+nGeneBlank,mean = (log_mean_sample),sd=sd(log_mean_sample))
  #Pram_log_disp=rnorm(nGeneMean+nGeneVar+nGeneBlank,mean = (log_disp_sample),sd=sd(log_disp_sample))
  #Pram_logit_drop=rnorm(nGeneMean+nGeneVar+nGeneBlank,mean = (logit_drop_sample),sd=sd(logit_drop_sample))
  
  sample_data=cbind(log_mean_sample,log_disp_sample,logit_drop_sample)
  
  saveRDS(sample_data,paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_data_",sim_folder,".rds"))
}

sample_data=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_data_",sim_folder,".rds"))

sample_mean=apply(sample_data,2,mean)
sample_var=apply(sample_data,2,var)

cov_matrix=matrix(NA,3,3)
for(i in 1:3){
  for(j in 1:3){
    cov_matrix[i,j]=cov(sample_data[,i],sample_data[,j])
  }
}

require(MASS)
gpar_ctrl=exp(mvrnorm(nGeneMean+nGeneVar+nGeneDisp+nGeneMult+nGeneBlank, mu = sample_mean, Sigma = cov_matrix,empirical = TRUE))
gpar_ctrl[,3]=gpar_ctrl[,3]/(1+gpar_ctrl[,3])

colnames(gpar_ctrl)=c("mean","dispersion","dropout")
gpar_case=gpar_ctrl

special_index=sample.int(nGeneTotal,(nGeneMean+nGeneVar+nGeneDisp+nGeneMult))
mean_index=as.numeric(special_index[i_mean])
var_index=as.numeric(special_index[i_var])
disp_index=as.numeric(special_index[i_disp])
mult_index=as.numeric(special_index[i_mult])
# label and save the DE index information.
de.mean = rep(0, nGeneTotal)
de.var  = rep(0, nGeneTotal)
de.disp  = rep(0, nGeneTotal)
de.mult  = rep(0, nGeneTotal)
de.mean[mean_index] = 1
de.var[var_index]   = 1
de.disp[disp_index]   = 1
de.mult[mult_index]   = 1

r_mean2=r_mean
r_var2=r_var  
if(r_mean<1){r_mean2=1/r_mean}
if(r_var<1){r_var2=1/r_var}

gpar_case[mean_index,1:2]=t(apply(gpar_case[mean_index,,drop=FALSE],1,function(x){return(calc_zinb_param(mu=x[1],theta=x[2],drop=x[3],r_m=r_mean2,r_v=1))})) 
gpar_case[var_index,1:2]=t(apply(gpar_case[var_index,,drop=FALSE],1,function(x){return(calc_zinb_param(mu=x[1],theta=x[2],drop=x[3],r_m=1,r_v=r_var2))})) 
gpar_case[disp_index,2]=gpar_case[mean_index,2,drop=FALSE]*r_disp 
gpar_case[mult_index,2]=t(apply(gpar_case[var_index,,drop=FALSE],1,function(x){return(cal_nbzinb_param_multimodality_enlarge(mu=x[1],size=x[2],drop=x[3],change_proportion=r_change_prop))}))


sim_case=matrix(nrow=nGeneTotal,ncol=ncase*ncell)
sim_ctrl=matrix(nrow=nGeneTotal,ncol=nctrl*ncell)

for(ig in 1:nGeneTotal){
  sim_case[ig,]=emdbook::rzinbinom(ncase*ncell,gpar_case[ig,1], gpar_case[ig,2], gpar_case[ig,3])
  sim_ctrl[ig,]=emdbook::rzinbinom(ncase*ncell,gpar_ctrl[ig,1], gpar_ctrl[ig,2], gpar_ctrl[ig,3])
}


# scatter plot 
pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/check_simulation/check_simulation_scatter_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".pdf"), 
    width = 9, height = 12)
op=par(mfrow = c(6, 5))

hist(log(apply(sim_ctrl[de.mean+de.var+de.disp+de.mult==0,],1,mean)))
hist(log(apply(sim_ctrl[de.mean==1,],1,mean)))
hist(log(apply(sim_ctrl[de.var==1,],1,mean)))
hist(log(apply(sim_ctrl[de.disp==1,],1,mean)))
hist(log(apply(sim_ctrl[de.mult==1,],1,mean)))

hist(log(apply(sim_case[de.mean+de.var+de.disp+de.mult==0,],1,mean)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.mean==1,],1,mean)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.var==1,],1,mean)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.disp==1,],1,mean)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.mult==1,],1,mean)),col=rgb(1,0,0,0.3))

hist(log(apply(sim_ctrl[de.mean+de.var+de.disp+de.mult==0,],1,var)))
hist(log(apply(sim_ctrl[de.mean==1,],1,var)))
hist(log(apply(sim_ctrl[de.var==1,],1,var)))
hist(log(apply(sim_ctrl[de.disp==1,],1,var)))
hist(log(apply(sim_ctrl[de.mult==1,],1,var)))

hist(log(apply(sim_case[de.mean+de.var+de.disp+de.mult==0,],1,var)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.mean==1,],1,var)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.var==1,],1,var)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.disp==1,],1,var)),col=rgb(1,0,0,0.3))
hist(log(apply(sim_case[de.mult==1,],1,var)),col=rgb(1,0,0,0.3))

plot(sort(log(apply(sim_ctrl[de.mean+de.var+de.disp+de.mult==0,],1,mean))),sort(log(apply(sim_case[de.mean+de.var+de.disp+de.mult==0,],1,mean))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.mean==1,],1,mean))),sort(log(apply(sim_case[de.mean==1,],1,mean))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.var==1,],1,mean))),sort(log(apply(sim_case[de.var==1,],1,mean))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.disp==1,],1,mean))),sort(log(apply(sim_case[de.disp==1,],1,mean))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.mult==1,],1,mean))),sort(log(apply(sim_case[de.mult==1,],1,mean))))
lines(c(-10,10),c(-10,10),col="red")

plot(sort(log(apply(sim_ctrl[de.mean+de.var+de.disp+de.mult==0,],1,var))),sort(log(apply(sim_case[de.mean+de.var++de.disp+de.mult==0,],1,var))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.mean==1,],1,var))),sort(log(apply(sim_case[de.mean==1,],1,var))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.var==1,],1,var))),sort(log(apply(sim_case[de.var==1,],1,var))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.disp==1,],1,var))),sort(log(apply(sim_case[de.disp==1,],1,var))))
lines(c(-10,10),c(-10,10),col="red")
plot(sort(log(apply(sim_ctrl[de.mult==1,],1,var))),sort(log(apply(sim_case[de.mult==1,],1,var))))
lines(c(-10,10),c(-10,10),col="red")
par(op)
dev.off()

#switch
temp=sim_case[i_case_modify,]
sim_case[i_case_modify,]=sim_ctrl[i_case_modify,]
sim_ctrl[i_case_modify,]=temp

temp=gpar_case[i_case_modify,]
gpar_case[i_case_modify,]=gpar_ctrl[i_case_modify,]
gpar_ctrl[i_case_modify,]=temp

sim_matrix=cbind(sim_case,sim_ctrl)

dim(sim_matrix)
sim_matrix[1:4,1:4]

table(c(sim_matrix) == 0)
table(c(sim_matrix) == 0)/(nrow(sim_matrix)*ncol(sim_matrix))

####################### Meta information collection #################

#the phenotype and individual information of simulated samples.
phenotype = c(rep(0, nctrl * ncell), rep(1, ncase * ncell))
individual = paste0("ind", c(rep(1:nall, each = ncell)))

#Count info for matrix
cell_id = paste0("cell", 1:ncol(sim_matrix))
gene_id = paste0("gene", 1:nrow(sim_matrix))

rownames(sim_matrix) = gene_id
colnames(sim_matrix) = cell_id

sim_param=array(dim=c(length(gene_id),3,2),dimnames=list(gene_id,c("mean","overdisp","dropout"),c("case","ctrl")))
sim_param[,,1]=gpar_case
sim_param[,,2]=gpar_ctrl

#Cell info for meta
cellsum = apply(sim_matrix, 2, sum)
genesum = apply(sim_matrix, 1, sum)
CDR  = apply(sim_matrix > 0, 2, sum) / nrow(sim_matrix)
meta = data.frame(cell_id, individual, phenotype, cellsum, CDR, 
                  stringsAsFactors=FALSE)
dim(meta)
meta[1:2,]

pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/check_covariates/check_covariates_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".pdf"), 
    width=6, height=3)
op=par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
boxplot(meta$cellsum ~ meta$phenotype, xlab="group", ylab="read-depth")
boxplot(meta$CDR ~ meta$phenotype, xlab="group", ylab="CDR")
par(op)
dev.off()

######### Some basic stat & Preparation for Bulk RNAseq analysis ############

# individual level info
cur_individual = unique(meta$individual)
cell_num = matrix(ncol = 1, nrow = length(cur_individual))
rownames(cell_num) = cur_individual
colnames(cell_num) = "cell_num"

read_depth = matrix(ncol = 1, nrow = length(cur_individual))
rownames(read_depth) = cur_individual
colnames(read_depth) = "read_depth"

phenotype_ind = matrix(ncol = 1, nrow = length(cur_individual))
rownames(phenotype_ind) = cur_individual
colnames(phenotype_ind) = "phenotype"

zero_rate_ind = matrix(nrow = nrow(sim_matrix),
                       ncol = length(cur_individual))
rownames(zero_rate_ind) = rownames(sim_matrix)
colnames(zero_rate_ind) = cur_individual

sim_matrix_bulk = matrix(nrow = nrow(sim_matrix),
                         ncol = length(cur_individual))
rownames(sim_matrix_bulk) = rownames(sim_matrix)
colnames(sim_matrix_bulk) = cur_individual

for (i_ind in 1:length(cur_individual)) {
  cur_ind = cur_individual[i_ind]
  #fit org
  cur_ind_m = sim_matrix[, meta$individual == cur_ind]
  cell_num[i_ind]   = ncol(cur_ind_m)
  read_depth[i_ind] = sum(cur_ind_m, na.rm = TRUE) / cell_num[i_ind] * 1000
  phenotype_ind[i_ind] = meta$phenotype[meta$individual == cur_ind][1]
  
  zero_rate_ind[, i_ind] = rowSums(cur_ind_m == 0, na.rm = TRUE)/cell_num[i_ind]
  sim_matrix_bulk[, i_ind] = rowSums(cur_ind_m, na.rm = TRUE)
}

tapply(read_depth, phenotype_ind, summary)



#almost no need to adjust library size.
saveRDS(de.mean,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(de.var,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.var_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(de.disp,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.disp_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(de.mult,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

saveRDS(read_depth,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_ind/sim_ind_readdepth_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(phenotype_ind,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_ind/sim_ind_phenotye_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(zero_rate_ind,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_ind/sim_ind_zero_rate_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

saveRDS(sim_param,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_param_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(sim_matrix,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(meta,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(sim_matrix_bulk,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))



