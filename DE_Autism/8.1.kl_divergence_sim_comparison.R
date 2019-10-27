#this code construct some distance measures between 2 vectors (no integer needed, no same length needed)
# 
# We will do simulations on the following situation to see if our method is robust.
# (1)Same variants, different mean
# Expectation: longer distance between means, larger kl divergence
# 
# #conclusion: distance longer, no distinguish 
# 
# (2)fix mean, different variant
# Expectation: longer distance between variance, larger kl divergence
# 
# #conclusion: It works to detect the difference of variances
# #It distinguish some complex situations when mu and sd are all differences as general concepts of human being.
# 
# (3)Discrete vs continuous data patterns.
# Expectation: robust to it.
# 
# #Generally, Thus this influence can be regarded as robust.
# 
# (4)Sample size differences
# Expectation: robust to it.
# 
# #conclusion: bigger sample size, better estimation,no huge influence in relative size differences.
# #generally robust
# 
# 
# (5)Bin lengths
# Expectation: robust to it.
# 
# #generally it can be controlled when the bin number are small. 
# #the bin number are recommended no bigger than 10% of the sample number
# 
# (6)Outliers and alt proportions
# Expectation: robust to it.
# #conclusion: unfortunately, this method is not very sensitive to the mixture situations.
# #even if 50% of the points from another distinct distribution, the distance value is not very huge.
# #its reasonable since the mixture situation means at least half of the points from the same distribution.



#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

############functions #####################################
source("./Command/8.0_kl_divergence_functions.R")

##############Simulation for testing#############################
# We will do simulations on the following situation to see if our method is robust.
# (1). Same variants, different mean
# Expectation: longer distance between means, larger kl divergence
pdf("8.kl_divergence_sim_comparison_1_for_mean.pdf",height=8,width = 6)
op=par(mfrow = c(2, 1))
n=1000
a1=rnorm(n,0,1)
distance1=matrix(ncol=1,nrow=length(0:10))
distance2=matrix(ncol=1,nrow=length(0:10))
for(mu in 0:10){
  a2=rnorm(n,mu,1)
  distance1[mu]=mean_KL_dens(a1,a2,empirical.n=n,alter="mean",fit_model="empirical")
  distance2[mu]=mean_KL_dens(a1,a2,empirical.n=n,alter="JSD",fit_model="empirical")
}
plot(0:10,distance1,pch=4,xlab="mu",ylab="Mean distance",main="Same sd, different mean",sub=paste0("norm KLdis,n=",n,",sd=1"))
plot(0:10,distance2,pch=4,xlab="mu",ylab="JSD distance",main="Same sd, different mean",sub=paste0("norm KLdis,n=",n,",sd=1"))
par(op)
dev.off()
#conclusion: distance longer, no distinguish 

# (2). fix mean, different variant
# Expectation: longer distance between variance, larger kl divergence
pdf("8.kl_divergence_sim_comparison_2_for_sd.pdf",height=8,width = 6)
op=par(mfrow = c(2, 1))
n=1000
a1=rnorm(n,0,1)
distance1=matrix(ncol=1,nrow=length(0:50))
distance2=matrix(ncol=1,nrow=length(0:50))
for(mu in c(0,1,2,3,5,10)){
  for(theta in 0:50){
    a2=rnorm(n,mu,theta/10)
    distance1[theta]=mean_KL_dens(a1,a2,empirical.n=n,alter="mean",fit_model="empirical")
    distance2[theta]=mean_KL_dens(a1,a2,empirical.n=n,alter="JSD",fit_model="empirical")
  }
  plot((0:50)/10,distance1,pch=4,xlab="sd",ylab="Mean distance",main="fix mean, different sd",sub=paste0("norm KLdis,n=",n,",ref_sd=1,mu=0 vs ",mu))
  plot((0:50)/10,distance2,pch=4,xlab="sd",ylab="JSD distance",main="fix mean, different sd",sub=paste0("norm KLdis,n=",n,",ref_sd=1,mu=0 vs ",mu))
}
par(op)
dev.off()
#conclusion: It works to detect the difference of variances
#It distinguish some complex situations when mu and sd are all differences as general concepts of human being.


# (3). Discrete vs continuous data patterns.
# Expectation: robust to it.
pdf("8.kl_divergence_sim_comparison_3_robust_discrete.pdf",height=8,width = 6)
op=par(mfrow = c(2, 1))
n=1000
a1=rnorm(n,0,1)

d=c(1,2,5,10,20,50,100)
distance1=matrix(ncol=1,nrow=length(d))
distance2=matrix(ncol=1,nrow=length(d))
for(i_d in 1:length(d)){
  a2=rep(a1[1:round(n/d[i_d])],d[i_d])
  distance1[i_d]=mean_KL_dens(a1,a2,empirical.n=n,alter="mean",fit_model="empirical")
  distance2[i_d]=mean_KL_dens(a1,a2,empirical.n=n,alter="JSD",fit_model="empirical")
}
plot(d,distance1,pch=4,xlab="divisor",ylab="Mean distance",main="Discrete vs continuous data patterns",sub=paste0("norm KLdis,n=",n,",mu=0,sd=1"))
plot(d,distance2,pch=4,xlab="divisor",ylab="JSD distance",main="Discrete vs continuous data patterns",sub=paste0("norm KLdis,n=",n,",mu=0,sd=1"))
par(op)
dev.off()

#Generally, Thus this influence can be regarded as robust.

# (4). Sample size differences
# Expectation: robust to it.

pdf("8.kl_divergence_sim_comparison_4_robust_set_size_diff.pdf",height=8,width = 6)
op=par(mfrow = c(2, 1))
n_seq=c(10,20,50,100,500,1000)
a1=rnorm(n,3,1)
distance1=matrix(ncol=1,nrow=100)
for(n in n_seq){
  for(n_sub in 1:100){
    a2=rnorm(n_sub*10,0,1)
    distance1[n_sub]=mean_KL_dens(a1,a2,empirical.n=n,alter="mean",fit_model="empirical")
    distance2[n_sub]=mean_KL_dens(a1,a2,empirical.n=n,alter="JSD",fit_model="empirical")
  }
  plot((1:100)*10,distance1,pch=4,xlab=paste0("n ",n," vs *"),ylab="Mean distance",main="Sample size differences",sub=paste0("norm KLdis,n_ref=",n,",mu=",mu," vs 0,sd=1"))
  plot((1:100)*10,distance2,pch=4,xlab=paste0("n ",n," vs *"),ylab="JSD distance",main="Sample size differences",sub=paste0("norm KLdis,n_ref=",n,",mu=",mu," vs 0,sd=1"))
}
par(op)
dev.off()
#conclusion: bigger sample size, better estimation,no huge influence in relative size differences.
#generally robust
#conclusion: distance quantile provide more distinguish resuls.


# (5). Bin lengths
# Expectation: robust to it.
pdf("8.kl_divergence_sim_comparison_5_robust_bin_size_diff.pdf",height=8,width = 6)
op=par(mfrow = c(2, 1))
n_seq=c(100,200,500,1000)
n_bin_seq=c(5,10,20,50,100)
a1=rnorm(n,0,1)

distance1=matrix(ncol=1,nrow=length(n_bin_seq))
distance2=matrix(ncol=1,nrow=length(n_bin_seq))
for(n in n_seq){
  for(mu in c(0,1,3,5,10)){
    a2=rnorm(n,mu,1)
    for(i_nbin in 1:length(n_bin_seq)){
      n_bin=n_bin_seq[i_nbin]
      distance1[i_nbin]=mean_KL_dens(a1,a2,empirical.n=n_bin*100,alter="mean",fit_model="empirical")
      distance2[i_nbin]=mean_KL_dens(a1,a2,empirical.n=n_bin*100,alter="JSD",fit_model="empirical")
    }
    plot(n_bin_seq,distance1,pch=4,xlab="dens pointx100",ylab="Mean distance",main="density sample number effect",sub=paste0("norm KLdis,n=",n,",mu=0 vs ",mu,",sd=1"))
    plot(n_bin_seq,distance2,pch=4,xlab="dens pointx100",ylab="JSD distance",main="density sample number effect",sub=paste0("norm KLdis,n=",n,",mu=0 vs ",mu,",sd=1"))
  }
}
par(op)
dev.off()
#generally it can be controlled when the bin number are small. 
#the bin number are recommended no bigger than 10% of the sample number

# (6). Outliers and alt proportions
# Expectation: robust to it.
pdf("8.kl_divergence_sim_comparison_6_robust_outlier_altproportion.pdf",height=8,width = 12)
op=par(mfrow = c(2, 2))
n=1000
a1=rnorm(n,0,1)
alt_num_seq=c(1,2,5,10,100,200,500,800,1000)

distance1=matrix(ncol=1,nrow=length(alt_num_seq))
distance2=matrix(ncol=1,nrow=length(alt_num_seq))
distance3=matrix(ncol=1,nrow=length(alt_num_seq))
distance4=matrix(ncol=1,nrow=length(alt_num_seq))
for(mu in c(1,3,5,10,100)){
  for(i_alt in 1:length(alt_num_seq)){
    alt_num=alt_num_seq[i_alt]
    a2=rnorm(n,0,1)
    a2[1:alt_num]=rnorm(alt_num,mu,1)
    distance1[i_alt]=mean_KL_dens(a1,a2,empirical.n=n,alter="mean",fit_model="empirical")
    distance2[i_alt]=mean_KL_dens(a1,a2,empirical.n=n,alter="JSD",fit_model="empirical")
    a2[1:max(round(alt_num/2),1)]=rnorm(max(round(alt_num/2),1),-mu,1)
    distance3[i_alt]=mean_KL_dens(a1,a2,empirical.n=n,alter="mean",fit_model="empirical")
    distance4[i_alt]=mean_KL_dens(a1,a2,empirical.n=n,alter="JSD",fit_model="empirical")
  }
  plot(alt_num_seq,distance1,pch=4,xlab="alt number",ylab="Mean distance",main="Outliers and alt proportions, single-side",sub=paste0("norm KLdis,n=",n,",mu=0 vs ",mu,",sd=1"))
  plot(alt_num_seq,distance2,pch=4,xlab="alt number",ylab="JSD distance",main="Outliers and alt proportions, single-side",sub=paste0("norm KLdis,n=",n,",mu=0 vs ",mu,",sd=1"))
  plot(alt_num_seq,distance3,pch=4,xlab="alt number",ylab="Mean distance",main="Outliers and alt proportions, two-side",sub=paste0("norm KLdis,n=",n,",mu=0 vs ",mu,",sd=1"))
  plot(alt_num_seq,distance4,pch=4,xlab="alt number",ylab="JSD distance",main="Outliers and alt proportions, two-side",sub=paste0("norm KLdis,n=",n,",mu=0 vs ",mu,",sd=1"))
}
par(op)
dev.off()

#conclusion: unfortunately, this method is not very sensitive to the mixture situations.
#even if 50% of the points from another distinct distribution, the distance value is not very huge.
#its reasonable since the mixture situation means at least half of the points from the same distribution.
#the double side alternative situation performs similar to the single side alternative

#sessionInfo()
#q(save="no")


