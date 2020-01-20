#this code used for the simulation

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
# #something from the head files please uncomment the following lines to run it separately.
# file_tag=1
# sim_method="zinb.naive" #splat.mean or splat.var--method 5, separate the mean and variance using splat
#                         #splat.org--method 4, change the mean.shape and mean.rate originally
# 
# #r_mean/r_var should < 1+mean.shape
#r_mean=1.5  #1.2,1.5,2,4,6
#r_var=1.5 #1.2,1.5,2,4,6
#r_disp=1.5 #1.2,1.5,2,4,6
#r_change_prop=0.75 #c(0.1,0.2,0.4,0.6,0.8)
  
  
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

library("emdbook")

nGeneMean=300
nGeneVar=300
nGeneDisp=300
nGeneMult=300
nGeneBlank=1800
nGeneTotal=nGeneMean+nGeneVar+nGeneMult+nGeneDisp+nGeneBlank
ncase=20   #!!!!!!!!Different from Wei codes, this is inside the body.
nctrl=20   #!!!!!!!!Different from Wei codes, this is inside the body.
ncell=100
nall=ncase+nctrl

HET=0 #0/1 label: if HET==0, generate case cells from completely case condition
      #           if HET==1, generate 1/2 case cells from completely case condition

i_mean=1:nGeneMean
i_var=(nGeneMean+1):(nGeneMean+nGeneVar)
i_disp=(nGeneMean+nGeneVar+1):(nGeneMean+nGeneVar+nGeneDisp)
i_mult=(nGeneMean+nGeneVar+nGeneDisp+1):(nGeneMean+nGeneVar++nGeneDisp+nGeneMult)
i_blank=(nGeneMean+nGeneVar+nGeneDisp+nGeneMult+1):nGeneTotal

case_modify_flag=rbinom(nGeneTotal,1,0.5)
i_case_modify=which(case_modify_flag==1) #randomlized settings, switching cases and control samples to make sure the library sizes the same


#first control, then case
i_ctrl=1:(nctrl*ncell)
i_case=(nctrl*ncell+1):((ncase+nctrl)*ncell)

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
# ############ Simulation data with Method 5 ###############################

#generate the data by ZINB, based on given parameters...
input_file_tag="3k10"
t_mean = read.table(paste0("../Data_PRJNA434002/res_dca_rawM",input_file_tag,"/mean_signif4.tsv.gz"), 
                    sep = "\t", header = TRUE)

t_disp = read.table(paste0("../Data_PRJNA434002/res_dca_rawM",input_file_tag,"/dispersion_signif4.tsv.gz"),
                    sep = "\t", header = TRUE)

t_drop = read.table(paste0("../Data_PRJNA434002/res_dca_rawM",input_file_tag,"/dropout_signif4.tsv.gz"),
                    sep = "\t", header = TRUE)

dim(t_mean)
t_mean[1:2,1:5]

dim(t_disp)
t_disp[1:2,1:5]

dim(t_drop)
t_drop[1:2,1:5]

############### collect sample information ###################

col_info   = strsplit(names(t_mean), split="_")
sample_ids = sapply(col_info, function(x){paste(x[-1], collapse="_")})
table(sample_ids)

tapply_mean <- function(x){tapply(x, sample_ids, mean)}
tapply_sd <- function(x){tapply(x, sample_ids, sd)}

sample_log_mean = t(apply(log(t_mean), 1, tapply_mean))
sample_disp = t(apply(t_disp, 1, tapply_mean))
sample_drop = t(apply(t_drop, 1, tapply_mean))

sample_log_mean_sd = t(apply(log(t_mean), 1, tapply_sd))

dim(sample_log_mean)
dim(sample_disp)
dim(sample_drop)
dim(sample_log_mean_sd)

sample_log_mean[1:2,1:5]
sample_log_mean_sd[1:2,1:5]

# the sd within an individual, across cells is large, 
# probably becaue the cells are from different cell types
# so we reduce them by a factor of 10 here. 
summary(apply(sample_log_mean, 1, sd))
summary(apply(sample_log_mean_sd,1,mean))

sample_log_mean_sd = sample_log_mean_sd/10

###################### simulation based on real data ######################

# extract information from the data
# We use the gene specific medians of each parameter.
# and simulate parameters assuming log-normal or logic-normal distribution.

# randomlize the orders of genes and samples, and take 40 samples
# set.seed(2019)

random_idx_gene = sample.int(nrow(sample_log_mean))
random_idx_sam  = sample.int(ncol(sample_log_mean))[1:nall]

sample_log_mean = sample_log_mean[random_idx_gene, random_idx_sam]
sample_drop = sample_drop[random_idx_gene, random_idx_sam]
sample_disp = sample_disp[random_idx_gene, random_idx_sam]

sample_log_mean_sd = sample_log_mean_sd[random_idx_gene, random_idx_sam]

dim(sample_log_mean)
dim(sample_drop)
dim(sample_disp)

dim(sample_log_mean_sd)

#################  calculate parameters for the case data #################

# sample gene index for genes differential expressed by mean or variance.
special_index = sample.int(nGeneTotal, (nGeneMean + nGeneVar +nGeneDisp+nGeneMult))
mean_index    = as.numeric(special_index[i_mean])
var_index     = as.numeric(special_index[i_var])
disp_index     = as.numeric(special_index[i_disp])
mult_index     = as.numeric(special_index[i_mult])
# label and save the DE index information.
de.mean = rep(0, nGeneTotal)
de.var  = rep(0, nGeneTotal)
de.disp  = rep(0, nGeneTotal)
de.mult  = rep(0, nGeneTotal)
de.mean[mean_index] = 1
de.var[var_index]   = 1
de.disp[disp_index]   = 1
de.mult[mult_index]   = 1
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


sample_mean = exp(sample_log_mean)

sample_mean_ctrls = sample_mean[,1:nctrl]
sample_disp_ctrls = sample_disp[,1:nctrl]
sample_drop_ctrls = sample_drop[,1:nctrl]
sample_mean_cases = sample_mean[,(nctrl+1):nall]
sample_disp_cases = sample_disp[,(nctrl+1):nall]
sample_drop_cases = sample_drop[,(nctrl+1):nall]


for(j in mean_index){
  cur_param=c(NA,NA) #initializaiton
  if(case_modify_flag[j]==1){
    for(i in 1:ncol(sample_disp_cases)){
      mu_ji    = sample_mean_cases[j,i]
      theta_ji = sample_disp_cases[j,i]
      drop_ji  = sample_drop_cases[j,i]
      cur_param = calc_zinb_param(mu=mu_ji,theta=theta_ji,drop=drop_ji,r_m=r_mean2,r_v=1)
      sample_mean_cases[j,i]=cur_param[1]
      sample_disp_cases[j,i]=cur_param[2]
    }
  }
  if(case_modify_flag[j]==0){
    for(i in 1:ncol(sample_disp_ctrls)){
      mu_ji    = sample_mean_ctrls[j,i]
      theta_ji = sample_disp_ctrls[j,i]
      drop_ji  = sample_drop_ctrls[j,i]
      cur_param = calc_zinb_param(mu=mu_ji,theta=theta_ji,drop=drop_ji,r_m=r_mean2,r_v=1)
      sample_mean_ctrls[j,i]=cur_param[1]
      sample_disp_ctrls[j,i]=cur_param[2]
    }
  }
}

for(j in var_index){
  cur_param=c(NA,NA) #initializaiton
  if(case_modify_flag[j]==1){
    for(i in 1:ncol(sample_disp_cases)){
      mu_ji    = sample_mean_cases[j,i]
      theta_ji = sample_disp_cases[j,i]
      drop_ji  = sample_drop_cases[j,i]
      cur_param = calc_zinb_param(mu=mu_ji,theta=theta_ji,drop=drop_ji,r_m=1,r_v=r_var2)
      sample_mean_cases[j,i]=cur_param[1]
      sample_disp_cases[j,i]=cur_param[2]
    }
  }
  if(case_modify_flag[j]==0){
    for(i in 1:ncol(sample_disp_ctrls)){
      mu_ji    = sample_mean_ctrls[j,i]
      theta_ji = sample_disp_ctrls[j,i]
      drop_ji  = sample_drop_ctrls[j,i]
      cur_param = calc_zinb_param(mu=mu_ji,theta=theta_ji,drop=drop_ji,r_m=1, r_v=r_var2)
      sample_mean_ctrls[j,i]=cur_param[1]
      sample_disp_ctrls[j,i]=cur_param[2]
    }
  }
}


for(j in disp_index){
  cur_param=c(NA,NA) #initializaiton
  if(case_modify_flag[j]==1){
    for(i in 1:ncol(sample_disp_cases)){
      theta_ji = sample_disp_cases[j,i]
      sample_disp_cases[j,i]=theta_ji*r_disp
    }
  }
  if(case_modify_flag[j]==0){
    for(i in 1:ncol(sample_disp_ctrls)){
      theta_ji = sample_disp_ctrls[j,i]
      sample_disp_ctrls[j,i]=theta_ji*r_disp
    }
  }
}

for(j in mult_index){
  cur_param=c(NA,NA) #initializaiton
  if(case_modify_flag[j]==1){
    for(i in 1:ncol(sample_disp_cases)){
      mu_ji    = sample_mean_cases[j,i]
      theta_ji = sample_disp_cases[j,i]
      drop_ji  = sample_drop_cases[j,i]
      sample_disp_cases[j,i]=cal_nbzinb_param_multimodality_enlarge(mu=mu_ji,size=theta_ji,drop=drop_ji,
                                                                   change_proportion=r_change_prop)
    }
  }
  if(case_modify_flag[j]==0){
    for(i in 1:ncol(sample_disp_ctrls)){
      mu_ji    = sample_mean_ctrls[j,i]
      theta_ji = sample_disp_ctrls[j,i]
      drop_ji  = sample_drop_ctrls[j,i]
      sample_disp_ctrls[j,i]=cal_nbzinb_param_multimodality_enlarge(mu=mu_ji,size=theta_ji,drop=drop_ji,
                                                                   change_proportion=r_change_prop)
    }
  }
}
#################  check the parameters  ###################

# check the total read depth
summary(colSums(sample_mean))
summary(colSums(sample_mean_cases))

ratio.adjust = median(colSums(sample_mean_cases))/median(colSums(sample_mean))
sample_mean_cases = sample_mean_cases/ratio.adjust

# scatter plot 
pdf(paste0("../Data_PRJNA434002/10.Result/check_simulation_scatter_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".pdf"), 
    width = 9, height = 6)
par(mfrow = c(3, 5), pty = "s")

plot(log10(apply(sample_mean_ctrls[de.mean + de.var + de.disp + de.mult == 0,], 1, mean)),
     log10(apply(sample_mean_cases[de.mean + de.var + de.disp + de.mult == 0,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, non-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.mean== 1,], 1, mean)),
     log10(apply(sample_mean_cases[de.mean== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Mean-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.var== 1,], 1, mean)),
     log10(apply(sample_mean_cases[de.var== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Var-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.disp== 1,], 1, mean)),
     log10(apply(sample_mean_cases[de.disp== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, DISP-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.mult== 1,], 1, mean)),
     log10(apply(sample_mean_cases[de.mult== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, MULT-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_disp_ctrls[de.mean + de.var + de.disp + de.mult == 0,], 1, mean)),
     log10(apply(sample_disp_cases[de.mean + de.var + de.disp + de.mult == 0,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 disp, non-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_disp_ctrls[de.mean== 1,], 1, mean)),
     log10(apply(sample_disp_cases[de.mean== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Mean-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_disp_ctrls[de.var== 1,], 1, mean)),
     log10(apply(sample_disp_cases[de.var== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, Var-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_disp_ctrls[de.disp== 1,], 1, mean)),
     log10(apply(sample_disp_cases[de.disp== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, DISP-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_disp_ctrls[de.mult== 1,], 1, mean)),
     log10(apply(sample_disp_cases[de.mult== 1,], 1, mean)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 mean, MULT-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.mean + de.var + de.disp + de.mult == 0,], 1, var)),
     log10(apply(sample_mean_cases[de.mean + de.var + de.disp + de.mult == 0,], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, non-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.mean== 1,], 1, var)),
     log10(apply(sample_mean_cases[de.mean== 1,], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, Mean-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.var== 1,], 1, var)),
     log10(apply(sample_mean_cases[de.var== 1,], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, Var-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.disp== 1,], 1, var)),
     log10(apply(sample_mean_cases[de.disp== 1,], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, DISP-DE genes")
abline(0, 1, col = "red")

plot(log10(apply(sample_mean_ctrls[de.mult== 1,], 1, var)),
     log10(apply(sample_mean_cases[de.mult== 1,], 1, var)),
     cex = .2, xlab = "control cells", ylab = "case cells",
     main = "log10 var, MULT-DE genes")
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
      if(i > nctrl && k > ncell/2 ){
        mean_i = sample_mean_cases[,i-nctrl]
        disp_i = sample_disp_cases[,i-nctrl]
        drop_i = sample_drop_cases[,i-nctrl]
        mean_i[de.mult & case_modify_flag]=mean_i[de.mult & case_modify_flag]+
          mean_i[de.mult & case_modify_flag]*((r_change_prop)/(1-r_change_prop))*rbinom(sum(de.mult & case_modify_flag),1,0.5)
      }else{
        mean_i = sample_mean_ctrls[,i]
        disp_i = sample_disp_ctrls[,i]
        drop_i = sample_drop_ctrls[,i]
        mean_i[de.mult & !case_modify_flag]=mean_i[de.mult & !case_modify_flag]+
          mean_i[de.mult & !case_modify_flag]*((r_change_prop)/(1-r_change_prop))*rbinom(sum(de.mult & !case_modify_flag),1,0.5)
      }
    }
    else{
      if(i > nctrl){
        mean_i = sample_mean_cases[,i-nctrl]
        disp_i = sample_disp_cases[,i-nctrl]
        drop_i = sample_drop_cases[,i-nctrl]
        mean_i[de.mult & case_modify_flag]=mean_i[de.mult & case_modify_flag]+
          mean_i[de.mult & case_modify_flag]*((r_change_prop)/(1-r_change_prop))*rbinom(sum(de.mult & case_modify_flag),1,0.5)
      }else{
        mean_i = sample_mean_ctrls[,i]
        disp_i = sample_disp_ctrls[,i]
        drop_i = sample_drop_ctrls[,i]
        mean_i[de.mult & !case_modify_flag]=mean_i[de.mult & !case_modify_flag]+
          mean_i[de.mult & !case_modify_flag]*((r_change_prop)/(1-r_change_prop))*rbinom(sum(de.mult & !case_modify_flag),1,0.5)
      }
    }
    
    sample_mean_k = exp(rnorm(nGeneTotal, log(mean_i), sample_log_mean_sd[,i]))
    
    for (ig in 1:nGeneTotal) {
      sim_matrix[ig,idx_i[k]] = emdbook::rzinbinom(1, sample_mean_k[ig], 
                                                   disp_i[ig], drop_i[ig])
      sim_param[ig,idx_i[k],]=c(sample_mean_k[ig], disp_i[ig], drop_i[ig])
    }
  }
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
cell_id = paste0("cell", 1:ncol(sim_matrix))
gene_id = paste0("gene", 1:nrow(sim_matrix))

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

pdf(paste0("../Data_PRJNA434002/10.Result/check_covariates_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".pdf"), 
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
saveRDS(de.mean,paste0("../Data_PRJNA434002/10.Result/sim_de.mean_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(de.var,paste0("../Data_PRJNA434002/10.Result/sim_de.var_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(de.disp,paste0("../Data_PRJNA434002/10.Result/sim_de.disp_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(de.mult,paste0("../Data_PRJNA434002/10.Result/sim_de.mult_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

saveRDS(read_depth,paste0("../Data_PRJNA434002/10.Result/sim_ind_readdepth_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(phenotype_ind,paste0("../Data_PRJNA434002/10.Result/sim_ind_phenotye_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(zero_rate_ind,paste0("../Data_PRJNA434002/10.Result/sim_ind_zero_rate_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

saveRDS(sim_matrix,paste0("../Data_PRJNA434002/10.Result/sim_param_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(sim_matrix,paste0("../Data_PRJNA434002/10.Result/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(meta,paste0("../Data_PRJNA434002/10.Result/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
saveRDS(sim_matrix_bulk,paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))




