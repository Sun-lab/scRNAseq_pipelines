
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
library("vegan")

sim_folder="sim_v1"
perm_num=500
covariate_flag=NA #c(NA, "quantile99")
tol=0.2
perm_label_seq=0:10



##############functions#################
library("ggplot2")
library("emdbook")

source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/7.0_ZINB_fit_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/8.0_kl_divergence_functions.R")
source("~/Desktop/github/scRNAseq_pipelines/DE_Autism/9.0_Fstat_functions.R")


#########functions###############
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
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




perm_label=1
perm_method=""

r_mean=1.2
r_var=1.2
r_disp=1.5
r_mult=0.4
file_tag=1
n_ind=40
n_cell=80
n_cell_seq=(n_cell-20):n_cell

ncell=n_cell
n=n_ind/2
r_change_prop=r_mult

dist_method="JSD"
fit_method="zinb"

mean_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))
var_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.var_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))
disp_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.disp_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))
mult_index=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".rds"))

#pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/fig_check_pval/pval_hist_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,".pdf"),width=25,height = 15)

for(ncell in n_cell_seq){
  dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",n_ind,"_",ncell,".rds"))
  hist(apply(dist_array,2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),xlim=c(0,1),sub=ncell)
  
}


tail_modification=function(m,prop_thres=0.8,fold_change=0.1,alter=c("greater","less")){
  thres_value=quantile(m,prop_thres,na.rm = TRUE)
  if(alter=="greater"){
    res=m-pmax((m-thres_value),0)*(1-fold_change)
  }
  if(alter=="less"){
    res=m-pmin((m-thres_value),0)*(1-fold_change)
  }
  res
}




for(ncell in n_cell_seq){
  
  phenotype=c(rep(1,n),rep(0,n))
  
  tryCatch({cur_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/",dist_method,"_",fit_method,"_pval/p10_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",ncell,".rds"))}, error = function(e) {NA} )
  tryCatch({het_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/tukey_median_pval/p10_",dist_method,"_",fit_method,"_perm_tukey_median_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",ncell,".rds"))}, error = function(e) {NA} )

  
  dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",n_ind,"_",ncell,".rds"))
  perm_order=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/perm_order/p10_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",ncell,".rds"))
  
  dist_array2=array(dim=dim(dist_array),dimnames=dimnames(dist_array))
  dist_array3=array(dim=dim(dist_array),dimnames=dimnames(dist_array))
  dist_array4=array(dim=dim(dist_array),dimnames=dimnames(dist_array))
  dist_array5=array(dim=dim(dist_array),dimnames=dimnames(dist_array))
  
  for(ig in 1:dim(dist_array)[1]){
    dist_array2[ig,,]=tail_modification(dist_array[ig,,],prop_thres = 0.8,fold_change = 0.1,alter="greater")
    dist_array3[ig,,]=tail_modification(dist_array[ig,,],prop_thres = 0.6,fold_change = 0.05,alter="greater")
    dist_array4[ig,,]=tail_modification(dist_array[ig,,],prop_thres = 0.8,fold_change = 3,alter="greater")
    dist_array5[ig,,]=tail_modification(dist_array[ig,,],prop_thres = 0.6,fold_change = 3,alter="greater")
  }
  
  
  
  
  
  op=par(mfrow = c(5, 5))
  tryCatch({hist(apply(dist_array[mean_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, ",ncell," dist_array"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array[var_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, ",ncell," dist_array"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )        
  tryCatch({hist(apply(dist_array[disp_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, ",ncell," dist_array"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array[mult_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, ",ncell," dist_array"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, ",ncell," dist_array"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  
  
  tryCatch({hist(apply(dist_array2[mean_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, ",ncell," dist_array2"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array2[var_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, ",ncell," dist_array2"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )        
  tryCatch({hist(apply(dist_array2[disp_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, ",ncell," dist_array2"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array2[mult_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, ",ncell," dist_array2"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array2[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, ",ncell," dist_array2"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  
  
  tryCatch({hist(apply(dist_array3[mean_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, ",ncell," dist_array3"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array3[var_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, ",ncell," dist_array3"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )        
  tryCatch({hist(apply(dist_array3[disp_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, ",ncell," dist_array3"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array3[mult_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, ",ncell," dist_array3"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array3[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, ",ncell," dist_array3"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  
  tryCatch({hist(apply(dist_array4[mean_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, ",ncell," dist_array4"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array4[var_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, ",ncell," dist_array4"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )        
  tryCatch({hist(apply(dist_array4[disp_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, ",ncell," dist_array4"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array4[mult_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, ",ncell," dist_array4"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array4[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, ",ncell," dist_array4"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  
  tryCatch({hist(apply(dist_array5[mean_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, ",ncell," dist_array5"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array5[var_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, ",ncell," dist_array5"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )        
  tryCatch({hist(apply(dist_array5[disp_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, ",ncell," dist_array5"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array5[mult_index==1,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, ",ncell," dist_array5"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  tryCatch({hist(apply(dist_array5[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,,],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, ",ncell," dist_array5"),xlab="dist",breaks = 20,xlim=c(0,1.2),ylim=c(0,9))}, error = function(e) {NA} )
  
  for(ip in 1:5){
    cur_phenotype=phenotype[perm_order[,ip]]
    dist_pval1=cal_permanova_pval2(dist_array,cur_phenotype,perm_num.min = perm_num)$pval
    dist_pval2=cal_permanova_pval2(dist_array2,cur_phenotype,perm_num.min = perm_num)$pval
    dist_pval3=cal_permanova_pval2(dist_array3,cur_phenotype,perm_num.min = perm_num)$pval
    dist_pval4=cal_permanova_pval2(dist_array4,cur_phenotype,perm_num.min = perm_num)$pval
    dist_pval5=cal_permanova_pval2(dist_array5,cur_phenotype,perm_num.min = perm_num)$pval
    dist_pval6=cal_permanova_pval2(dist_array,1-cur_phenotype,perm_num.min = perm_num)$pval
    
    tryCatch({hist(dist_pval1[mean_index==1],main=paste0("mean-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval1",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval1[var_index==1],main=paste0("var-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval1",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(dist_pval1[disp_index==1],main=paste0("disp-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval1",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval1[mult_index==1],main=paste0("mult-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval1",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval1[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main=paste0("non-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval1",breaks = 20)}, error = function(e) {NA} )
    # 
    tryCatch({hist(dist_pval6[mean_index==1],main=paste0("mean-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval6",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval6[var_index==1],main=paste0("var-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval6",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(dist_pval6[disp_index==1],main=paste0("disp-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval6",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval6[mult_index==1],main=paste0("mult-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval6",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval6[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main=paste0("non-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval6",breaks = 20)}, error = function(e) {NA} )
    
    tryCatch({hist(dist_pval2[mean_index==1],main=paste0("mean-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval2",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval2[var_index==1],main=paste0("var-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval2",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(dist_pval2[disp_index==1],main=paste0("disp-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval2",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval2[mult_index==1],main=paste0("mult-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval2",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval2[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main=paste0("non-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval2",breaks = 20)}, error = function(e) {NA} )
    # 
    
    tryCatch({hist(dist_pval3[mean_index==1],main=paste0("mean-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval3",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval3[var_index==1],main=paste0("var-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval3",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(dist_pval3[disp_index==1],main=paste0("disp-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval3",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval3[mult_index==1],main=paste0("mult-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval3",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval3[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main=paste0("non-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval3",breaks = 20)}, error = function(e) {NA} )
    
    tryCatch({hist(dist_pval4[mean_index==1],main=paste0("mean-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval4",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval4[var_index==1],main=paste0("var-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval4",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(dist_pval4[disp_index==1],main=paste0("disp-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval4",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval4[mult_index==1],main=paste0("mult-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval4",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval4[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main=paste0("non-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval4",breaks = 20)}, error = function(e) {NA} )
    
    tryCatch({hist(dist_pval5[mean_index==1],main=paste0("mean-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval5",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval5[var_index==1],main=paste0("var-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval5",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(dist_pval5[disp_index==1],main=paste0("disp-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval5",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval5[mult_index==1],main=paste0("mult-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval5",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(dist_pval5[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main=paste0("non-DE, ",ncell," ", dist_method,"_",fit_method),xlab="dist_pval5",breaks = 20)}, error = function(e) {NA} )
    

    # 
  }
  par(op)
 
  
  
  
  ########heatmap#####################

  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)

  
  for(ip in 1:5){
    cur_pheno=phenotype[perm_order[,ip]]
    o_pheno=order(cur_pheno)
    o_pval=order(cur_pval[,ip])
    for(io in c(1:10,10001:1010)){
      cur_o=o_pval[io]
      heatmap.2(dist_array[cur_o,o_pheno,o_pheno],
                main = paste0("Distance of pval ",round(cur_pval[cur_o],5)," p ",ip,", o ",io), # heat map title
                notecol="black",      # change font color of cell labels to black
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(12,9),     # widens margins around plot
                col=my_palette,       # use on color palette defined earlier
                #breaks=col_breaks,    # enable color transition at specified limits
                dendrogram='none', 
                Rowv=FALSE, 
                Colv=FALSE)
      hist(apply(dist_array[cur_o,,],1,sum))
    }
  }

  # ############Distance #########
  # mean_dist_perc_plot=function(m,label,cur_range=NA){
  # 
  #   mean_1=apply(m[label==1,],2,function(x){return(sum(x,na.rm = TRUE))}) #the distance to cases
  #   mean_0=apply(m[label==0,],2,function(x){return(sum(x,na.rm = TRUE))}) #the distance to ctrls
  #   mean1=mean_1/(sum(label==1)-label)
  #   mean0=mean_0/(sum(label==0)+label-1)
  #   if(sum(!is.na(cur_range))==0){
  #     cur_range=c(0,max(c(mean0,mean1)))
  #   }
  #   
  #   plot(mean1[label==1],mean0[label==1],col="red",pch=16,xlim=cur_range,ylim=cur_range,xlab="mean dist to in-network",ylab="mean dist to out-network")
  #   points(mean0[label==0],mean1[label==0],pch=16,col="blue")
  #   lines(cur_range,cur_range)
  #   
  #   res=matrix(ncol=1,nrow=length(mean1))
  #   res[label==1]=mean1[label==1]/mean0[label==1]
  #   res[label==0]=mean0[label==0]/mean1[label==0]
  #   res
  # }
  # 
  # 
  # # to check if remove the ones with median extreme pval in cases distance or controls distance.... It does remove signals.
  # cal_permanova_pval_rm=function(m,label,case_remove_num=2,zm=NA,Fstat_method="p",perm_num.min=500){
  #     mean_1=apply(m[label==1,],2,function(x){return(sum(x,na.rm = TRUE))}) #the distance to cases
  #     mean_0=apply(m[label==0,],2,function(x){return(sum(x,na.rm = TRUE))}) #the distance to ctrls
  #     mean1=mean_1/(sum(label==1)-label)
  #     mean0=mean_0/(sum(label==0)+label-1)
  #     
  #     mean_max=pmax(mean1,mean0)
  #     mean_max_case_temp=mean_max
  #     mean_max_case_temp[label==0]=-1
  #     mean_max_ctrl_temp=mean_max
  #     mean_max_ctrl_temp[label==1]=-1
  #     order_case=order(c(mean_max_case_temp,decreasing = TRUE))[1:case_remove_num]
  #     order_ctrl=order(c(mean_max_ctrl_temp,decreasing = TRUE))[1:case_remove_num]
  #     order_rm=unique(c(order_case,order_ctrl))
  #     m=m[-order_rm,-order_rm]
  #     label=label[-order_rm]
  #     res=cal_permanova_pval(m,label,zm=zm,Fstat_method=Fstat_method,perm_num.min=perm_num.min)
  #     res
  # }
  
  
  for(ip in 1:5){
    cur_pheno=phenotype[perm_order[,ip]]
    tryCatch({cur_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/",dist_method,"_",fit_method,"_pval/p10_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",ncell,".rds"))}, error = function(e) {NA} )
    dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",n_ind,"_",ncell,".rds"))
    
    #het_pval[,ip]=apply(dist_array,1,function(x)return(cal_betadisper_pval(x,cur_pheno,disper_type="median",div_method="anova")))
    hist(cur_pval[,ip])
    cur_pheno=phenotype[perm_order[,ip]]
    o_pheno=order(cur_pheno)
    o_pval=order(cur_pval[,ip])
    
    # op=par(mfrow = c(4, 5))
    # for(ig in order(cur_pval[,ip])[1:20]){
    #   prop=mean_dist_perc_plot(dist_array[ig,,],cur_pheno,cur_range=c(0,0.05))
    # }
    # for(ig in order(cur_pval[,ip])[2001:2020]){
    #   mean_dist_perc_plot(dist_array[ig,,],cur_pheno,cur_range=c(0,0.05))
    # }
    # par(op)
    
    dist_pval=matrix(ncol=1,nrow=dim(dist_array)[1])
    dist_pval_500=dist_pval

    for(ig in 1:dim(dist_array)[1]){
      dist_pval[ig]=cal_permanova_pval(dist_array[ig,,],cur_pheno,perm_num.min = 500)$pval
      dist_pval_500[ig]=cal_permanova_pval(dist_array[ig,,],cur_pheno[sample.int(length(cur_pheno),length(cur_pheno))],perm_num.min = 500)$pval

    }
    
    op=par(mfrow = c(1, 2))
    hist(dist_pval_500,breaks=20)
    hist(dist_pval,breaks=20)
    par(op)
    
    dist_array_flag=apply(dist_array,1,function(x)sum(is.na(x)))
    dist_array_sub=dist_array[dist_array_flag==0,,]
    dist_pval_eigen_vector=apply(dist_array_sub,1,function(x)eigen(x)$values)
  }

  
  ###########################
  
  op=par(mfrow = c(4, 5))
  for(ip in 1:5){
    cur_pheno=phenotype[perm_order[,ip]]
    #het_pval[,ip]=apply(dist_array,1,function(x)return(cal_betadisper_pval(x,cur_pheno,disper_type="median",div_method="anova")))
    
    tryCatch({hist(cur_pval[mean_index==1,ip],main=paste0("mean-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(cur_pval[var_index==1,ip],main=paste0("var-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
    tryCatch({hist(cur_pval[disp_index==1,ip],main=paste0("disp-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(cur_pval[mult_index==1,ip],main=paste0("mult-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    tryCatch({hist(cur_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip],main=paste0("non-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    # 
    # tryCatch({hist(het_pval[mean_index==1,ip],main=paste0("mean-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20,xlim=c(0,0.8))}, error = function(e) {NA} )
    # tryCatch({hist(het_pval[var_index==1,ip],main=paste0("var-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20,xlim=c(0,0.8))}, error = function(e) {NA} )        
    # tryCatch({hist(het_pval[disp_index==1,ip],main=paste0("disp-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20,xlim=c(0,0.8))}, error = function(e) {NA} )
    # tryCatch({hist(het_pval[mult_index==1,ip],main=paste0("mult-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20,xlim=c(0,0.8))}, error = function(e) {NA} )
    # tryCatch({hist(het_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip],main=paste0("non-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20,xlim=c(0,0.8))}, error = function(e) {NA} )
    
    
    tryCatch({hist(apply(dist_array[mean_index==1,cur_pheno==1,cur_pheno==1],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, p",ip," ",ncell," dist_array_case"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    tryCatch({hist(apply(dist_array[var_index==1,cur_pheno==1,cur_pheno==1],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, p",ip," ",ncell," dist_array_case"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )        
    tryCatch({hist(apply(dist_array[disp_index==1,cur_pheno==1,cur_pheno==1],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, p",ip," ",ncell," dist_array_case"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    tryCatch({hist(apply(dist_array[mult_index==1,cur_pheno==1,cur_pheno==1],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, p",ip," ",ncell," dist_array_case"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    tryCatch({hist(apply(dist_array[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,cur_pheno==1,cur_pheno==1],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, p",ip," ",ncell," dist_array_case"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    
    tryCatch({hist(apply(dist_array[mean_index==1,cur_pheno==0,cur_pheno==0],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mean-DE, p",ip," ",ncell," dist_array_ctrl"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    tryCatch({hist(apply(dist_array[var_index==1,cur_pheno==0,cur_pheno==0],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("var-DE, p",ip," ",ncell," dist_array_ctrl"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )        
    tryCatch({hist(apply(dist_array[disp_index==1,cur_pheno==0,cur_pheno==0],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("disp-DE, p",ip," ",ncell," dist_array_ctrl"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    tryCatch({hist(apply(dist_array[mult_index==1,cur_pheno==0,cur_pheno==0],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("mult-DE, p",ip," ",ncell," dist_array_ctrl"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    tryCatch({hist(apply(dist_array[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,cur_pheno==0,cur_pheno==0],2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),main=paste0("non-DE, p",ip," ",ncell," dist_array_ctrl"),xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    
    tryCatch({hist(c(apply(dist_array[mean_index==1,cur_pheno==0,cur_pheno==1],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),
                     apply(dist_array[mean_index==1,cur_pheno==1,cur_pheno==0],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1])),
      main=paste0("mean-DE, p",ip," ",ncell," dist_array_out"),
      xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )

    tryCatch({hist(c(apply(dist_array[var_index==1,cur_pheno==0,cur_pheno==1],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),
                     apply(dist_array[var_index==1,cur_pheno==1,cur_pheno==0],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1])),
      main=paste0("var-DE, p",ip," ",ncell," dist_array_out"),
      xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )  

    tryCatch({hist(c(apply(dist_array[disp_index==1,cur_pheno==0,cur_pheno==1],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),
                     apply(dist_array[disp_index==1,cur_pheno==1,cur_pheno==0],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1])),
      main=paste0("disp-DE, p",ip," ",ncell," dist_array_out"),
      xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )

    tryCatch({hist(c(apply(dist_array[mult_index==1,cur_pheno==0,cur_pheno==1],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),
                     apply(dist_array[mult_index==1,cur_pheno==1,cur_pheno==0],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1])),
      main=paste0("mult-DE, p",ip," ",ncell," dist_array_out"),
      xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )

    tryCatch({hist(c(apply(dist_array[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,cur_pheno==0,cur_pheno==1],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1]),
                     apply(dist_array[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,cur_pheno==1,cur_pheno==0],
                           2,function(x)sum(x,na.rm = TRUE)/dim(x)[1])),
      main=paste0("non-DE, p",ip," ",ncell," dist_array_out"),
      xlab="dist",breaks = 20,xlim=c(0,0.8),ylim=c(0,9))}, error = function(e) {NA} )
    
  }
  par(op)
  
  # op=par(mfrow = c(5, 5))
  # for(ip in 1:5){
  #   tryCatch({hist(cur_pval[mean_index==1,ip],main=paste0("mean-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  #   tryCatch({hist(cur_pval[var_index==1,ip],main=paste0("var-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
  #   tryCatch({hist(cur_pval[disp_index==1,ip],main=paste0("disp-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  #   tryCatch({hist(cur_pval[mult_index==1,ip],main=paste0("mult-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  #   tryCatch({hist(cur_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip],main=paste0("non-DE, p",ip," ",ncell," ", dist_method,"_",fit_method),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
  #   
    # tryCatch({hist(het_pval[mean_index==1,ip],main=paste0("mean-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    # tryCatch({hist(het_pval[var_index==1,ip],main=paste0("var-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
    # tryCatch({hist(het_pval[disp_index==1,ip],main=paste0("disp-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    # tryCatch({hist(het_pval[mult_index==1,ip],main=paste0("mult-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    # tryCatch({hist(het_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip],main=paste0("non-DE, p",ip," ",ncell," anova_median"),xlab="p-values",breaks = 20)}, error = function(e) {NA} )
    # 
    # tryCatch({plot(-log10(cur_pval[mean_index==1,ip]),-log10(het_pval[mean_index==1,ip]),xlab="permanova pval",ylab="anova_median_pval",cex=.1,main=paste0("mean-DE, p",ip," ",ncell," permanova vs anova_median"),sub=round(cor(-log10(cur_pval[mean_index==1,ip]),-log10(het_pval[mean_index==1,ip]),use="complete.obs"),3))}, error = function(e) {NA} )
    # tryCatch({plot(-log10(cur_pval[var_index==1,ip]),-log10(het_pval[var_index==1,ip]),xlab="permanova pval",ylab="anova_median_pval",cex=.1,main=paste0("var-DE, p",ip," ",ncell," permanova vs anova_median"),sub=round(cor(-log10(cur_pval[var_index==1,ip]),-log10(het_pval[var_index==1,ip]),use="complete.obs"),3))}, error = function(e) {NA} )       
    # tryCatch({plot(-log10(cur_pval[disp_index==1,ip]),-log10(het_pval[disp_index==1,ip]),xlab="permanova pval",ylab="anova_median_pval",cex=.1,main=paste0("disp-DE, p",ip," ",ncell," permanova vs anova_median"),sub=round(cor(-log10(cur_pval[disp_index==1,ip]),-log10(het_pval[disp_index==1,ip]),use="complete.obs"),3))}, error = function(e) {NA} )
    # tryCatch({plot(-log10(cur_pval[mult_index==1,ip]),-log10(het_pval[mult_index==1,ip]),xlab="permanova pval",ylab="anova_median_pval",cex=.1,main=paste0("mult-DE, p",ip," ",ncell," permanova vs anova_median"),sub=round(cor(-log10(cur_pval[mult_index==1,ip]),-log10(het_pval[mult_index==1,ip]),use="complete.obs"),3))}, error = function(e) {NA} )
    # tryCatch({plot(-log10(cur_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip]),-log10(het_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip]),xlab="permanova pval",ylab="anova_median_pval",cex=.1,main=paste0("non-DE, p",ip," ",ncell," permanova vs anova_median"),sub=round(cor(-log10(cur_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip]),-log10(het_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0,ip]),use="complete.obs"),3))}, error = function(e) {NA} )
  #   
  # }
  # par(op)
  print(ncell)
}

#dev.off()


