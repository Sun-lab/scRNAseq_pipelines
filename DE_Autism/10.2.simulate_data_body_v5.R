#this code used for the simulation
#compared with version v1, it extend the number of individuals by estimating the mean and overdisp from the individual level mean and overdisp distribution.

# Simulate scRNA-seq data using splatter
# Using their default simulator for now. We can consider simulate 3000 genes in 40 individuals, with 20 cases vs. 20 controls. We can set DE for 300 genes in terms of mean expression difference, and another 300 genes DE in terms of variance difference. In the simulation, try to set the total read-depth across samples to be the same. 
# As an initial analysis, do not add any covariates. 
# For each gene, calculate density, and then JSD across samples, and then use PERMANOVA to calculate p-value for each gene. 
# Also calculate p-value for each gene using MAST-mixed effect model (less priority for now). 
# Collapse gene expression across cells per individual, then run DESeq2 for differential expression testing. 
# Type I error is the proportion of those 2400 equivalently expressed genes with p-value smaller than 0.05. 
# Power1, the proportion of the first 300 genest with p-value < 0.05
# Power2, the proportion of the next 300 genest with p-value < 0.05


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



#However, this method doesn't work since the parameter mean.rate influence nothing in the package, so finally we plan to simulate 40 batches by the splatter, then manually modify the mean and variance different on the output count data.


# ####### method5: Using ZINB method directly ##########

##### Start simulation from here #############
#something from the head files please uncomment the following lines to run it separately.
# file_tag=1
# sim_method="zinb.naive" #splat.mean or splat.var--method 5, separate the mean and variance using splat
#                         #splat.org--method 4, change the mean.shape and mean.rate originally

# #r_mean/r_var should < 1+mean.shape
# r_mean=1.2  #1.2,1.5,2,4,6
# r_var=1.2 #1.2,1.5,2,4,6
# r_disp=1.2 #1.2,1.5,2,4,6
# r_change_prop=0.2 #c(0.1,0.2,0.3,0.4)
#  dp_minor_prop=0.2 #c(0.1,0.2,0.3,0.4)
  
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
#setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")

library("emdbook")
sim_folder="sim_v5"

nGeneMean=300
nGeneVar=300
nGeneMult=300
nGeneDP=300
nGeneBlank=1800
nGeneTotal=nGeneMean+nGeneVar+nGeneMult+nGeneDP+nGeneBlank
ncase=50   #!!!!!!!!Different from Wei codes, this is inside the body.
nctrl=50   #!!!!!!!!Different from Wei codes, this is inside the body.
ncell=200
nall=ncase+nctrl

HET=0 #0/1 label: if HET==0, generate case cells from completely case condition
      #           if HET==1, generate 1/2 case cells from completely case condition

i_mean=1:nGeneMean
i_var=(nGeneMean+1):(nGeneMean+nGeneVar)
i_mult=(nGeneMean+nGeneVar+1):(nGeneMean+nGeneVar+nGeneMult)
i_dp=(nGeneMean+nGeneVar+nGeneMult+1):(nGeneMean+nGeneVar+nGeneMult+nGeneDP)
i_blank=(nGeneMean+nGeneVar+nGeneMult+nGeneDP+1):nGeneTotal

case_modify_flag=rbinom(nGeneTotal,1,0.5)
i_case_modify=which(case_modify_flag==1) #randomlized settings, switching cases and control samples to make sure the library sizes the same


#first case, then ctrl
i_case=1:(ncase*ncell)
i_ctrl=(ncase*ncell+1):(nall*ncell)

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
  return(c(size3))
}
# ############ Simulation Data Preparation ###############################


#generate the data by ZINB, based on given parameters...
input_file_tag="3k" 

if(!file.exists(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_mean.rds"))){
  #cov_matrix=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_cov_matrix.rds"))
  t_mean=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_mean.rds"))
  t_disp=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_dispersion.rds"))
  t_drop=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_dropout.rds"))
  
  dim(t_mean)
  t_mean[1:2,1:5]
  
  dim(t_disp)
  t_disp[1:2,1:5]
  
  dim(t_drop)
  t_drop[1:2,1:5]
  
  
  
  ############### collect sample information ###################
  
  col_info   = strsplit(colnames(t_mean), split="_")
  sample_ids = sapply(col_info, function(x)x[2])
  table(sample_ids)
  
  tapply_mean <- function(x){tapply(x, sample_ids, function(x)return(mean(x,na.rm = TRUE)))}
  tapply_median <- function(x){tapply(x, sample_ids, function(x)return(median(x,na.rm = TRUE)))}
  tapply_sd <- function(x){tapply(x, sample_ids, function(x)return(sd(x,na.rm = TRUE)))}
  
  mid_mean = t(apply(t_mean, 1, tapply_median))
  mid_disp = t(apply(t_disp, 1, tapply_median))
  mid_drop = t(apply(t_drop, 1, tapply_median))
  sample_log_mean_sd = t(apply(t_mean, 1, function(x)tapply_sd(log(as.numeric(x)))))
  
  sample_log_mean=log(mid_mean)
  sample_log_disp=log(mid_disp)
  sample_logit_drop=log((mid_drop)/(1-mid_drop))
  
  
  
  dim(sample_log_mean)
  dim(sample_log_disp)
  dim(sample_logit_drop)
  dim(sample_log_mean_sd)
  
  sample_log_mean[1:2,1:5]
  sample_log_mean_sd[1:2,1:5]
  
  # the sd within an individual, across cells is large, 
  # probably becaue the cells are from different cell types
  # so we reduce them by a factor of 10 here. 
  summary(apply(sample_log_mean, 1, sd))
  summary(apply(sample_log_mean_sd,1,mean))
  
  #sample_log_mean_sd = sample_log_mean_sd/10
  
  
  #simulate individual level parameters from correlated mutinormal distribution
  # sample_data=cbind(c(sample_log_mean),c(sample_log_disp),c(sample_logit_drop))
  # sample_data_mean=apply(sample_data,2,mean)
  # 
  # cov_matrix=matrix(NA,3,3)
  # for(i in 1:3){
  #   for(j in 1:3){
  #     cov_matrix[i,j]=cov(sample_data[,i],sample_data[,j])
  #   }
  # }
  
  dim(sample_log_mean)
  dim(sample_log_disp)
  dim(sample_logit_drop)
  dim(sample_log_mean_sd)
  
  saveRDS(sample_log_mean,paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_mean.rds"))
  saveRDS(sample_log_disp,paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_disp.rds"))
  saveRDS(sample_logit_drop,paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_logit_drop.rds"))
  saveRDS(sample_log_mean_sd,paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_mean_sd.rds"))
}
###########################

sample_log_mean=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_mean.rds"))
sample_log_disp=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_disp.rds"))
sample_logit_drop=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_logit_drop.rds"))
sample_log_mean_sd=readRDS(paste0("../Data_PRJNA434002/dca_PFC_all/",input_file_tag,"_sample_log_mean_sd.rds"))



require(MASS)
sample_ctrl=array(dim=c(nGeneTotal,nall,4),dimnames=list(paste0("gene", 1:nGeneTotal),paste0("ind", 1:nall),c("mean","dispersion","dropout","mean_sd")))

ind_strength=0.5 #from 0 to 1, use this to adjust individual mean expression strength in the simulation.
for(ig in 1:nGeneTotal){
  sample_data=cbind(c(sample_log_mean[ig,]),c(sample_log_disp[ig,]),c(sample_logit_drop[ig,]),c(sample_log_mean_sd[ig,]))
  
  sample_data_reg=apply(sample_data,2,mean)
  sample_data2=sample_data*ind_strength+matrix(rep(sample_data_reg,each=nrow(sample_data)),nrow=nrow(sample_data))*(1-ind_strength)
  
  sample_data_mean=apply(sample_data2,2,function(x)mean(x,na.rm = TRUE))
  
  cov_matrix=matrix(NA,4,4)
  for(i in 1:4){
    for(j in 1:4){
      cov_matrix[i,j]=cov(sample_data[,i],sample_data[,j])
    }
  }
  sample_ctrl[ig,,]=exp(mvrnorm(nall, mu = sample_data_mean, Sigma = cov_matrix,empirical = TRUE))
}

sample_ctrl[,,3]=sample_ctrl[,,3]/(1+sample_ctrl[,,3])




###################### simulation based on real data ######################

random_idx_gene = sample.int(nGeneTotal)
random_idx_sam  = sample.int(nall)

sample_ctrl=sample_ctrl[random_idx_gene, random_idx_sam,]

#################  calculate parameters for the case data #################

# sample gene index for genes differential expressed by mean or variance.
special_index = sample.int(nGeneTotal, (nGeneMean + nGeneVar +nGeneMult+nGeneDP))
mean_index    = as.numeric(special_index[i_mean])
var_index     = as.numeric(special_index[i_var])
mult_index     = as.numeric(special_index[i_mult])
dp_index     = as.numeric(special_index[i_dp])
# label and save the DE index information.
de.mean = rep(0, nGeneTotal)
de.var  = rep(0, nGeneTotal)
de.mult  = rep(0, nGeneTotal)
de.dp  = rep(0, nGeneTotal)
de.mean[mean_index] = 1
de.var[var_index]   = 1
de.mult[mult_index]   = 1
de.dp[dp_index]   = 1
# To make sure all parameters are non-negative, we do some transformation
r_mean2 = r_mean
r_var2  = r_var

if (r_mean < 1) {
  r_mean2 = 1 / r_mean
}

if (r_var < 1) {
  r_var2 = 1 / r_var
}



# modify parameters for cases mean and var
# now r_m (r_mean2) is smaller than 1. We need to modify gene expression in 
# x proportion of genes by fold r_m and (1-x) proportion of genes by 
# fold of 1/r_m, so that x(1 - r_m) = (1-x)(1/r_m - 1) 
# x (1 - r_m + 1/r_m -1) = (1/r_m -1) => x = 1/(1 + r_m) 


sample_param_case=sample_ctrl[,1:ncase,] #gene x ind x (mean,overdispersion, droupout)
sample_param_ctrl=sample_ctrl[,(ncase+1):nall,]
temp=NA

for(i in mean_index){
  for(j in 1:ncase){
    x=sample_param_case[i,j,]
    sample_param_case[i,j,1:2]=calc_zinb_param(mu=x[1],theta=x[2],drop=x[3],r_m=r_mean2,r_v=1)
  }
}
for(i in var_index){
  for(j in 1:ncase){
    x=sample_param_case[i,j,]
    sample_param_case[i,j,1:2]=calc_zinb_param(mu=x[1],theta=x[2],drop=x[3],r_m=1,r_v=r_var2)
  }
}

for(i in mult_index){
  for(j in 1:ncase){
    sample_param_ctrl[i,j,2]=cal_nbzinb_param_multimodality_enlarge(mu=sample_param_ctrl[i,j,1]*(1-r_change_prop),
                                                                    size=sample_param_ctrl[i,j,2],
                                                                    drop=sample_param_ctrl[i,j,3],
                                                                    change_proportion=r_change_prop)
    sample_param_case[i,j,1]=sample_param_case[i,j,1]+sample_param_case[i,j,1]*2*(0.5-rbinom(1,1,0.5))*(r_change_prop)
  }
}

for(i in dp_index){
  for(j in 1:ncase){
    cur_dp=rbinom(1,1,dp_minor_prop)
    if(cur_dp==0){
      x=sample_param_case[i,j,]
      sample_param_case[i,j,1:2]=calc_zinb_param(mu=x[1],theta=x[2],drop=x[3],r_m=r_mean2,r_v=1)
    }
    if(cur_dp==1){
      x= sample_param_ctrl[i,j,]
      sample_param_ctrl[i,j,1:2]=calc_zinb_param(mu=x[1],theta=x[2],drop=x[3],r_m=r_mean2,r_v=1)
    }
  }
}

temp=sample_param_case[i_case_modify,,]
sample_param_case[i_case_modify,,]=sample_param_ctrl[i_case_modify,,]
sample_param_ctrl[i_case_modify,,]=temp


#################  check the parameters  ###################

# check the total read depth
#summary(colSums(sample_mean))
#summary(colSums(sample_mean_cases))
#ratio.adjust = median(colSums(sample_param_case[,,1]))/median(colSums(sample_mean))
#sample_mean_cases = sample_mean_cases/ratio.adjust

quantile(colSums(sample_param_case[mean_index,,1]))
quantile(colSums(sample_param_ctrl[mean_index,,1]))
quantile(colSums(sample_param_case[var_index,,1]))
quantile(colSums(sample_param_ctrl[var_index,,1]))
quantile(colSums(sample_param_case[mult_index,,1]))
quantile(colSums(sample_param_ctrl[mult_index,,1]))
quantile(colSums(sample_param_case[dp_index,,1]))
quantile(colSums(sample_param_ctrl[dp_index,,1]))

# scatter plot 
pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/check_simulation/check_simulation_scatter_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".pdf"), 
    width = 8, height = 6)
par(mfrow = c(3, 4), pty = "s")

plot(log10(apply(sample_param_ctrl[de.mean + de.var + de.mult + de.dp == 0,,1], 1, mean)),
     log10(apply(sample_param_case[de.mean + de.var + de.mult + de.dp == 0,,1], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, non-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mean== 1,,1], 1, mean)),
     log10(apply(sample_param_case[de.mean== 1,,1], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Mean-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.var== 1,,1], 1, mean)),
     log10(apply(sample_param_case[de.var== 1,,1], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Var-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mult== 1,,1], 1, mean)),
     log10(apply(sample_param_case[de.mult== 1,,1], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, MULT-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.dp== 1,,1], 1, mean)),
     log10(apply(sample_param_case[de.dp== 1,,1], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, DP-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mean + de.var + de.mult + de.dp == 0,,2], 1, mean)),
     log10(apply(sample_param_case[de.mean + de.var + de.mult + de.dp == 0,,2], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 disp, non-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mean== 1,,2], 1, mean)),
     log10(apply(sample_param_case[de.mean== 1,,2], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 disp, Mean-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.var== 1,,2], 1, mean)),
     log10(apply(sample_param_case[de.var== 1,,2], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Var-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mult== 1,,2], 1, mean)),
     log10(apply(sample_param_case[de.mult== 1,,2], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, MULT-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.dp== 1,,2], 1, mean)),
     log10(apply(sample_param_case[de.dp== 1,,2], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, DP-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mean + de.var + de.mult + de.dp == 0,,1], 1, var)),
     log10(apply(sample_param_case[de.mean + de.var + de.mult + de.dp == 0,,1], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, non-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mean== 1,,1], 1, var)),
     log10(apply(sample_param_case[de.mean== 1,,1], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, Mean-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.var== 1,,1], 1, var)),
     log10(apply(sample_param_case[de.var== 1,,1], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, Var-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.mult== 1,,1], 1, var)),
     log10(apply(sample_param_case[de.mult== 1,,1], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, MULT-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_param_ctrl[de.dp== 1,,1], 1, var)),
     log10(apply(sample_param_case[de.dp== 1,,1], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, DP-DE genes")
abline(0, 1, col = "red")

dev.off()


# simulate scRNAseq based on zinb parameters of cases and controls:
sim_matrix = matrix(nrow = nGeneTotal, ncol = nall * ncell)
sim_param = array(dim=c(nGeneTotal,nall * ncell,3))
date()
for(i in 1:nall){
  idx_i = ((i-1)*ncell+1):(i*ncell)
  for(k in 1:ncell){
    if(HET){
      if(i > ncase && k > ncell/2 ){
        mean_i = sample_param_ctrl[,(i-ncase),1]
        disp_i = sample_param_ctrl[,(i-ncase),2]
        drop_i = sample_param_ctrl[,(i-ncase),3]
        sample_mean_sd_i = sample_param_ctrl[,(i-ncase),4]
      }else{
        mean_i = sample_param_case[,i,1]
        disp_i = sample_param_case[,i,2]
        drop_i = sample_param_case[,i,3]
        sample_mean_sd_i = sample_param_case[,i,4]
      }
    }
    if(!HET){
      if(i > nctrl){
        mean_i = sample_param_ctrl[,(i-ncase),1]
        disp_i = sample_param_ctrl[,(i-ncase),2]
        drop_i = sample_param_ctrl[,(i-ncase),3]
        sample_mean_sd_i = sample_param_ctrl[,(i-ncase),4]
      }else{
        mean_i = sample_param_case[,i,1]
        disp_i = sample_param_case[,i,2]
        drop_i = sample_param_case[,i,3]
        sample_mean_sd_i = sample_param_case[,i,4]
      }
    }
    
    sample_mean_k = exp(rnorm(nGeneTotal, log(mean_i), log(sample_mean_sd_i)))
    for (ig in 1:nGeneTotal) {
      sim_matrix[ig,idx_i[k]] = emdbook::rzinbinom(1, sample_mean_k[ig], 
                                                   disp_i[ig], drop_i[ig])
      sim_param[ig,idx_i[k],]=c(sample_mean_k[ig], disp_i[ig], drop_i[ig])
    }
  }
  print(i)
}
date()

dim(sim_matrix)
sim_matrix[1:4,1:4]

table(c(sim_matrix) == 0)
table(c(sim_matrix) == 0)/(nrow(sim_matrix)*ncol(sim_matrix))

####################### Meta information collection #################

#the phenotype and individual information of simulated samples.
phenotype = c(rep(0, nctrl * ncell), rep(1, ncase * ncell))
individual = paste0("ind", c(rep(1:nall, each = ncell)))

#Count info for matrix
cell_id = paste0("cell", 1:(nall*ncell))
gene_id = paste0("gene", 1:nGeneTotal)

rownames(sim_matrix) = gene_id
colnames(sim_matrix) = cell_id

dimnames(sim_param)=list(gene_id,cell_id,c("mean","overdisp","dropout"))

#Cell info for meta
cellsum = apply(sim_matrix, 2, sum)
genesum = apply(sim_matrix, 1, sum)
CDR  = apply(sim_matrix > 0, 2, sum) / nrow(sim_matrix)
meta = data.frame(cell_id, individual, phenotype, cellsum, CDR, 
                  stringsAsFactors=FALSE)
dim(meta)
meta[1:2,]

pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/check_covariates/check_covariates_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".pdf"), 
    width=6, height=3)
par(mfrow=c(1,2), mar=c(5,4,1,1), bty="n")
boxplot(meta$cellsum ~ meta$phenotype, xlab="group", ylab="read-depth")
boxplot(meta$CDR ~ meta$phenotype, xlab="group", ylab="CDR")
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
saveRDS(de.mean,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(de.var,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.var_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(de.mult,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(de.dp,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.dp_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

saveRDS(read_depth,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_ind/sim_ind_readdepth_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(phenotype_ind,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_ind/sim_ind_phenotye_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(zero_rate_ind,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_ind/sim_ind_zero_rate_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

saveRDS(sim_param,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_param_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(sim_matrix,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(meta,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
saveRDS(sim_matrix_bulk,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

write.csv(sim_matrix,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_matrix_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".csv"))








