#this check is especially designed for the dist array.

#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#here we did it for
#seedtag=5574587 #Jz*_1.2_1.2_1.5_0.4_1_10_100*.R
#seedtag=7638446 #Je*1.2_1.2_1.5_0.4_1_40_60*.R
#seedtag=9917471 #Jz_1.2_1.2_1.5_0.4_1_40_80.R


perm_label_seq=1:10
perm_label=10
perm_method=""
param_tag="disp"

r_mean=1.2
r_var=1.2
r_disp=1.5
r_mult=0.4
file_tag=1
n_ind=40
n_cell=80

ncell=n_cell
n=n_ind/2
r_change_prop=r_mult

mean_index=readRDS(paste0("../Data_PRJNA434002/10.Result/de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))
var_index=readRDS(paste0("../Data_PRJNA434002/10.Result/de_label/sim_de.var_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))
disp_index=readRDS(paste0("../Data_PRJNA434002/10.Result/de_label/sim_de.disp_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))
mult_index=readRDS(paste0("../Data_PRJNA434002/10.Result/de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))

dist_method="JSD"
fit_method="zinb"




for(cur_n_cell in (n_cell-30):n_cell){
  cur_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_pval/p10",perm_method,"_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",cur_n_cell,".rds"))
  cur_pval=cur_pval[,perm_label]

  op=par(mfrow = c(3, 2))
  tryCatch({hist(cur_pval[mean_index==1,1],main=paste0("pval of mean-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 1st perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(cur_pval[var_index==1,1],main=paste0("pval of var-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 1st perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(cur_pval[disp_index==1,1],main=paste0("pval of disp-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 1st perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(cur_pval[mult_index==1,1],main=paste0("pval of mult-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 1st perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  tryCatch({hist(cur_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,1],main=paste0("pval of non-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 1st perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  par(op)

  # op=par(mfrow = c(3, 2))
  # tryCatch({hist(cur_pval[mean_index==1,],main=paste0("pval of mean-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 10 perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  # tryCatch({hist(cur_pval[var_index==1,],main=paste0("pval of var-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 10 perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  # tryCatch({hist(cur_pval[disp_index==1,],main=paste0("pval of disp-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 10 perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  # tryCatch({hist(cur_pval[mult_index==1,],main=paste0("pval of mult-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 10 perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  # tryCatch({hist(cur_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,],main=paste0("pval of non-DE genes ",dist_method,"_",fit_method," cell ",cur_n_cell," 10 perm"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  # par(op)
}
  























t_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_sim_param=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_param_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
read_depth=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_ind/sim_ind_readdepth_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
phenotype_ind=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_ind/sim_ind_phenotye_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
zero_rate_ind=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_ind/sim_ind_zero_rate_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

perm_label=1
cur_n_cell=79
cur_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_pval/p10",perm_method,"_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",cur_n_cell,".rds"))
cur_pval=cur_pval[,1]
dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",cur_n_cell,".rds"))
sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

selected_index=1:20
#selected_index=c(10,15,11,16,8) #jz disp 1.2 1.2 1.1 0.4 2 10 60
#selected_index=c(17,13,5,18,8,6,3,11,15,4) #jz mean 1.1 1.2 1.2 0.4 1 20 80
#selected_index=c(1,11,14,12,17,5,9,4,7,3) # je mult 1.2 1.2 1.2 0.2 1 20 100

total_cell_index=matrix(ncol=1,nrow=0)
sub_cell_index_case=matrix(ncol=1,nrow=0)
sub_cell_index_ctrl=matrix(ncol=1,nrow=0)
for(i_s in c(selected_index,(20+selected_index))){
  cell_index=(100*i_s-cur_n_cell+1):(100*i_s)
  #cell_index=(100*i_s-ncell+1):(100*i_s-ncell+20)
  total_cell_index=c(total_cell_index,cell_index)
  cell_index_case=(100*i_s-cur_n_cell+1)
  cell_index_ctrl=(100*i_s-cur_n_cell+2):(100*i_s)
  sub_cell_index_case=c(sub_cell_index_case,cell_index_case)
  sub_cell_index_ctrl=c(sub_cell_index_ctrl,cell_index_ctrl)
}

#calculation
sim_matrix=t_sim_matrix[,total_cell_index]
sim_param=t_sim_param[,total_cell_index,]
meta=t_meta[total_cell_index,]

dist_array=dist_array[order(cur_pval),,]
sim_param=sim_param[order(cur_pval),,]
sim_matrix=sim_matrix[order(cur_pval),]
cur_pval2=cur_pval[order(cur_pval)]

sim_matrix_case=t_sim_matrix[,sub_cell_index_case]
sim_param_case=t_sim_param[,sub_cell_index_case,]
meta_case=t_meta[sub_cell_index_case,]
sim_param_case=sim_param_case[order(cur_pval),,]
sim_matrix_case=sim_matrix_case[order(cur_pval),]

sim_matrix_ctrl=t_sim_matrix[,sub_cell_index_ctrl]
sim_param_ctrl=t_sim_param[,sub_cell_index_ctrl,]
meta_ctrl=t_meta[sub_cell_index_ctrl,]
sim_param_ctrl=sim_param_ctrl[order(cur_pval),,]
sim_matrix_ctrl=sim_matrix_ctrl[order(cur_pval),]


#check param
op=par(mfrow = c(2, 3))
for(i in 1:10){
  j=i+2000
  hist(sim_param[i,,1],breaks=10)
  hist(sim_param_case[i,,1],breaks=10)
  hist(sim_param_ctrl[i,,1],breaks=10)
  hist(sim_param[j,,1],breaks=10)
  hist(sim_param_case[j,,1],breaks=10)
  hist(sim_param_ctrl[j,,1],breaks=10)
  hist(sim_matrix[i,],breaks=10)
  hist(sim_matrix_case[i,],breaks=10)
  hist(sim_matrix_ctrl[i,],breaks=10)
  hist(sim_matrix[j,],breaks=10)
  hist(sim_matrix_case[j,],breaks=10)
  hist(sim_matrix_ctrl[j,],breaks=10)
} 
par(op)

#check homogeneity














res=matrix(ncol=3,nrow=0)
res_case=matrix(ncol=3,nrow=0)
res_ctrl=matrix(ncol=3,nrow=0)
res_out=matrix(ncol=3,nrow=0)
zero_num=matrix(ncol=1,nrow=3000)
dist_diff_vector=zero_num
for(ig in 1:3000){
  dist_matrix=dist_array[ig,,]
  for(i in 1:(nrow(dist_matrix)-1)){
    for(j in (i+1):nrow(dist_matrix)){
      dist_range=dist_matrix[i,-c(i,j)]-dist_matrix[j,-c(i,j)]
      dist_diff=abs(max(dist_range,na.rm = TRUE)-min(dist_range,na.rm = TRUE))
       if(dist_diff<0.00001){
          res=rbind(res,c(ig,i,j))
      #   print(c(ig,i,j))
      #   print(quantile(dist_range))
       }
    }
  }
  zero_num[ig]=sum(is.na(dist_matrix))
  dist_diff_vector[ig]=dist_diff
}

tie_table=table(res[,1])

op=par(mfrow = c(2, 3))
hist(cur_pval2)
hist(cur_pval2[as.numeric(names(tie_table[tie_table>2]))])
hist(cur_pval2[-as.numeric(names(tie_table[tie_table>2]))])

hist(cur_pval2)
hist(cur_pval2[which(zero_num>0)])
hist(cur_pval2[-which(zero_num>0)])

hist(cur_pval2)
hist(cur_pval2[unique(as.numeric(names(tie_table)),which(zero_num>0))])
hist(cur_pval2[-unique(as.numeric(names(tie_table)),which(zero_num>0))])

hist(cur_pval2[-unique(as.numeric(names(tie_table)))])
hist(cur_pval2[-unique(as.numeric(names(tie_table[tie_table>2])),which(zero_num>0))])
hist(cur_pval2[-unique(as.numeric(names(tie_table[tie_table>4])),which(zero_num>0))])
par(op)










#redo all the method

#levene test from Anderson 2005
library("ggplot2")
library("emdbook")
library("Rcpp")
library("Barycenter")

source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/7.0_ZINB_fit_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/8.0_kl_divergence_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/9.0_Fstat_functions.R")

library("vegan")


#dist_matrix=nxn symmetric matrix with diag =0
#label nx1 vector about categories.


cal_betadisper_pval=function(dist_matrix,label,disper_type=c("median","centroid"),div_method=c("anova","TukeyHSD")){
  flag=which(apply(dist_matrix,1,function(x){sum(is.na(x))==length(x)}))
  if(length(flag)>0){
    label=as.factor(label[-flag])
    dist_matrix=dist_matrix[-flag,-flag]
    dist_matrix[is.na(dist_matrix)]=median(dist_matrix)
  }
  
  dis=as.dist(dist_matrix)
  mod=betadisper(dis, label,type=disper_type)

  if(div_method=="anova"){
    pval=anova(mod)$"Pr(>F)"[1]
  }
  if(div_method=="TukeyHSD"){
    pval=TukeyHSD(mod)$group[4]
  }
  return(pval)
}




pdf(paste0("~/Desktop/11.check_simulation_",dist_method,"_",fit_method,"_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".pdf"),height = 25,width = 25)

for(i_s in 1:20){
  set.seed(i_s)
  # part I homogenies
  cur_phenotype_ind=phenotype_ind[sample.int(length(phenotype_ind),length(phenotype_ind))]
  dist_res=cal_permanova_pval2(dist_array,cur_phenotype_ind,perm_num.min =1000)
  dist_pval=dist_res$pval
  
  
  
  #part II ties
  res=matrix(ncol=3,nrow=0)
  for(ig in 1:3000){
    dist_matrix=dist_array[ig,,]
    for(i in 1:(nrow(dist_matrix)-1)){
      for(j in (i+1):nrow(dist_matrix)){
        dist_range=dist_matrix[i,-c(i,j)]-dist_matrix[j,-c(i,j)]
        dist_diff=abs(max(dist_range,na.rm = TRUE)-min(dist_range,na.rm = TRUE))
        if(dist_diff<0.00001){
          res=rbind(res,c(ig,i,j))
        }
      }
    }
  }
  tie_table=table(res[,1])

  
  #if(sum(dist_pval<0.1)>500){

    print(i)
    
    # pval=matrix(ncol=1,nrow=3000)
    # for(i in 1:3000){
    #   print(i)
    #   pval[i]=cal_betadisper_pval(dist_array[i,,],label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD")
    # }
    tukey_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
    anova_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )})
    tukey_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
    anova_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="anova"), error = function(e) {NA} )}) 
    
    
    op=par(mfrow = c(6, 5))
    tryCatch({hist(dist_pval[mean_index==1],main="pval of mean-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[var_index==1],main="pval of var-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[disp_index==1],main="pval of disp-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[mult_index==1],main="pval of mult-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    tryCatch({hist(tukey_centroid_pval[mean_index==1],main="pval of mean-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_centroid_pval[var_index==1],main="pval of var-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_centroid_pval[disp_index==1],main="pval of disp-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_centroid_pval[mult_index==1],main="pval of mult-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,tukey_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

    tryCatch({hist(anova_centroid_pval[mean_index==1],main="pval of mean-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_centroid_pval[var_index==1],main="pval of var-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_centroid_pval[disp_index==1],main="pval of disp-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_centroid_pval[mult_index==1],main="pval of mult-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_centroid",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

    tryCatch({hist(tukey_median_pval[mean_index==1],main="pval of mean-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_median_pval[var_index==1],main="pval of var-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_median_pval[disp_index==1],main="pval of disp-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_median_pval[mult_index==1],main="pval of mult-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(tukey_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,tukey_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

    tryCatch({hist(anova_median_pval[mean_index==1],main="pval of mean-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_median_pval[var_index==1],main="pval of var-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_median_pval[disp_index==1],main="pval of disp-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_median_pval[mult_index==1],main="pval of mult-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(anova_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_median",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    
    
    
    
    tryCatch({plot(dist_pval[mean_index==1],tukey_centroid_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],tukey_centroid_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[var_index==1],tukey_centroid_pval[var_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,var-DE genes, cor ",round(cor(dist_pval[var_index==1],tukey_centroid_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[disp_index==1],tukey_centroid_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],tukey_centroid_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mult_index==1],tukey_centroid_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs tukey_centroid,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],tukey_centroid_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs tukey_centroid,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
    
    tryCatch({plot(dist_pval[mean_index==1],anova_centroid_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],anova_centroid_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[var_index==1],anova_centroid_pval[var_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,var-DE genes, cor ",round(cor(dist_pval[var_index==1],anova_centroid_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[disp_index==1],anova_centroid_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],anova_centroid_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mult_index==1],anova_centroid_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs anova_centroid,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],anova_centroid_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs anova_centroid,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_centroid_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
    
    
    tryCatch({plot(dist_pval[mean_index==1],tukey_median_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],tukey_median_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[var_index==1],tukey_median_pval[var_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,var-DE genes, cor ",round(cor(dist_pval[var_index==1],tukey_median_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[disp_index==1],tukey_median_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],tukey_median_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mult_index==1],tukey_median_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs tukey_median,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],tukey_median_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs tukey_median,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],tukey_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
    
    tryCatch({plot(dist_pval[mean_index==1],anova_median_pval[mean_index==1],cex=.2, main=paste0("dist_pval vs anova_median,mean-DE genes, cor ",round(cor(dist_pval[mean_index==1],anova_median_pval[mean_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[var_index==1],anova_median_pval[var_index==1],cex=.2, main=paste0("dist_pval vs anova_median,var-DE genes, cor ",round(cor(dist_pval[var_index==1],anova_median_pval[var_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[disp_index==1],anova_median_pval[disp_index==1],cex=.2, main=paste0("dist_pval vs anova_median,disp-DE genes, cor ",round(cor(dist_pval[disp_index==1],anova_median_pval[disp_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mult_index==1],anova_median_pval[mult_index==1],cex=.2, main=paste0("dist_pval vs anova_median,mult-DE genes, cor ",round(cor(dist_pval[mult_index==1],anova_median_pval[mult_index==1],use="complete.obs"),3)))}, error = function(e) {NA} )
    tryCatch({plot(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],cex=.2, main=paste0("dist_pval vs anova_median,non-DE genes, cor ",round(cor(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],anova_median_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],use="complete.obs"),3)))}, error = function(e) {NA} )
    par(op)
    
    #Part excusion
    op=par(mfrow = c(5, 5))
    tryCatch({hist(dist_pval[mean_index==1],main="pval of mean-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[var_index==1],main="pval of var-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[disp_index==1],main="pval of disp-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[mult_index==1],main="pval of mult-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,dist method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    tryCatch({hist(dist_pval[anova_median_pval<=0.05 & mean_index==1],main="pval of mean-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval<=0.05 & var_index==1],main="pval of var-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval<=0.05 & disp_index==1],main="pval of disp-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval<=0.05 & mult_index==1],main="pval of mult-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval<=0.05 & mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_median_pval<=0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    
    tryCatch({hist(dist_pval[anova_median_pval>0.05 & mean_index==1],main="pval of mean-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval>0.05 & var_index==1],main="pval of var-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval>0.05 & disp_index==1],main="pval of disp-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval>0.05 & mult_index==1],main="pval of mult-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[anova_median_pval>0.05 & mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,anova_median_pval>0.05",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    
    
    
    tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(mean_index==1)))],main="pval of mean-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(var_index==1)))],main="pval of var-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(disp_index==1)))],main="pval of disp-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which(mult_index==1)))],main="pval of mult-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[unique(c(as.numeric(names(tie_table[tie_table>2])),which( mean_index==0 & var_index==0 & disp_index==0 & mult_index==0)))],main="pval of non-DE genes,tie>2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(mean_index==0)))],main="pval of mean-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(var_index==0)))],main="pval of var-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(disp_index==0)))],main="pval of disp-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(mult_index==0)))],main="pval of mult-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval[-unique(c(as.numeric(names(tie_table[tie_table>2])),which(mean_index==1 | var_index==1 | disp_index==1 | mult_index==1)))],main="pval of non-DE genes,tie<=2",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    
    par(op)
  #}
}
dev.off()


op=par(mfrow = c(2, 3))
hist(dist_pval)
hist(dist_pval[as.numeric(names(tie_table[tie_table>2]))])
hist(dist_pval[-as.numeric(names(tie_table[tie_table>2]))])
hist(dist_pval[-which(anova_median_pval<0.05) ])
hist(dist_pval[which(anova_median_pval<0.05) ])

