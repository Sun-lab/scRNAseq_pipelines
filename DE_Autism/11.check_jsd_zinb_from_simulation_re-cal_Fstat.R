
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


perm_label=1
perm_method=""
param_tag="var"

r_mean=1.2
r_var=1.2
r_disp=1.2
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

jsd_zinb_pval=NA
jsd_empirical_pval=NA
jsd_direct_pval=NA
klmean_zinb_pval=NA
klmean_empirical_pval=NA
klmean_direct_pval=NA
deseq2_pval=NA
MAST_pval=NA

if(perm_label>0){
  tryCatch({jsd_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_zinb_pval/p",perm_label,perm_method,"_JSD_zinb_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({jsd_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_empirical_pval/p",perm_label,perm_method,"_JSD_empirical_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({jsd_direct_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_direct_pval/p",perm_label,perm_method,"_JSD_direct_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_zinb_pval/p",perm_label,perm_method,"_mean_zinb_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_empirical_pval/p",perm_label,perm_method,"_mean_empirical_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_direct_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_direct_pval/p",perm_label,perm_method,"_mean_direct_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({deseq2_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,".rds"))}, error = function(e) {NA} )
  
  #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
  tryCatch({MAST_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/MAST_pval/p",perm_label,perm_method,"_MAST_pval1_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
}
if(perm_label==0){
  tryCatch({jsd_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_zinb_pval/JSD_zinb_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({jsd_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_empirical_pval/JSD_empirical_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({jsd_direct_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/JSD_direct_pval/p",perm_label,perm_method,"_JSD_direct_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_zinb_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_zinb_pval/mean_zinb_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_empirical_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_empirical_pval/mean_empirical_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({klmean_direct_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/mean_direct_pval/p",perm_label,perm_method,"_mean_direct_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
  tryCatch({deseq2_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/DESeq2_pval/DESeq2_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,".rds"))}, error = function(e) {NA} )
  
  #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
  tryCatch({MAST_pval=readRDS(paste0("../Data_PRJNA434002/10.Result/MAST_pval/MAST_pval1_",r_mean,"_",r_var,"_",r_disp,"_",r_mult,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
}



op=par(mfrow = c(3, 2))

# tryCatch({hist(klmean_zinb_pval[mean_index==1],main="pval of mean-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# tryCatch({hist(klmean_zinb_pval[var_index==1],main="pval of var-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
# tryCatch({hist(klmean_zinb_pval[disp_index==1],main="pval of disp-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# tryCatch({hist(klmean_zinb_pval[mult_index==1],main="pval of mult-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# tryCatch({hist(klmean_zinb_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# 
# tryCatch({hist(jsd_zinb_pval[mean_index==1],main="pval of mean-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# tryCatch({hist(jsd_zinb_pval[var_index==1],main="pval of var-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
# tryCatch({hist(jsd_zinb_pval[disp_index==1],main="pval of disp-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# tryCatch({hist(jsd_zinb_pval[mult_index==1],main="pval of mult-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
# tryCatch({hist(jsd_zinb_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

tryCatch({hist(klmean_empirical_pval[mean_index==1],main="pval of mean-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_empirical_pval[var_index==1],main="pval of var-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
tryCatch({hist(klmean_empirical_pval[disp_index==1],main="pval of disp-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_empirical_pval[mult_index==1],main="pval of mult-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_empirical_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

tryCatch({hist(jsd_empirical_pval[mean_index==1],main="pval of mean-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_empirical_pval[var_index==1],main="pval of var-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
tryCatch({hist(jsd_empirical_pval[disp_index==1],main="pval of disp-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_empirical_pval[mult_index==1],main="pval of mult-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_empirical_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

tryCatch({hist(jsd_direct_pval[mean_index==1],main="pval of mean-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_direct_pval[var_index==1],main="pval of var-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
tryCatch({hist(jsd_direct_pval[disp_index==1],main="pval of disp-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_direct_pval[mult_index==1],main="pval of mult-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(jsd_direct_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

tryCatch({hist(klmean_direct_pval[mean_index==1],main="pval of mean-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_direct_pval[var_index==1],main="pval of var-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
tryCatch({hist(klmean_direct_pval[disp_index==1],main="pval of disp-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_direct_pval[mult_index==1],main="pval of mult-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(klmean_direct_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

tryCatch({hist(deseq2_pval[mean_index==1],main="pval of mean-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(deseq2_pval[var_index==1],main="pval of var-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )  
tryCatch({hist(deseq2_pval[disp_index==1],main="pval of disp-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(deseq2_pval[mult_index==1],main="pval of mult-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(deseq2_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

tryCatch({hist(MAST_pval[mean_index==1],main="pval of mean-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(MAST_pval[var_index==1],main="pval of var-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )  
tryCatch({hist(MAST_pval[disp_index==1],main="pval of disp-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(MAST_pval[mult_index==1],main="pval of mult-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
tryCatch({hist(MAST_pval[mean_index==0 & var_index==0 & disp_index==0 & mult_index==0],main="pval of non-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )

par(op)









t_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_sim_param=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_param_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))


dist_method="jsd"
fit_method="empirical"
cur_pval=jsd_empirical_pval
dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))



read_depth=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_ind/sim_ind_readdepth_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
phenotype_ind=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_ind/sim_ind_phenotye_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
zero_rate_ind=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_ind/sim_ind_zero_rate_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))


sim_matrix_bulk=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_bulk_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))

selected_index=1:20
#selected_index=c(10,15,11,16,8) #jz disp 1.2 1.2 1.1 0.4 2 10 60
#selected_index=c(17,13,5,18,8,6,3,11,15,4) #jz mean 1.1 1.2 1.2 0.4 1 20 80
#selected_index=c(1,11,14,12,17,5,9,4,7,3) # je mult 1.2 1.2 1.2 0.2 1 20 100

total_cell_index=matrix(ncol=1,nrow=0)
for(i_s in c(selected_index,(20+selected_index))){
  cell_index=(100*i_s-ncell+1):(100*i_s)
  #cell_index=(100*i_s-ncell+1):(100*i_s-ncell+20)
  total_cell_index=c(total_cell_index,cell_index)
}

#calculation
sim_matrix=t_sim_matrix[,total_cell_index]
sim_param=t_sim_param[,total_cell_index,]
meta=t_meta[total_cell_index,]

dist_array=dist_array[order(cur_pval),,]
sim_param=sim_param[order(cur_pval),,]
sim_matrix=sim_matrix[order(cur_pval),]

op=par(mfrow = c(2, 2))
for(i in 1:10){
  j=i+2000
  plot(sim_param[i,,1],sim_param[j,,1],main=paste0("param mean",i))
  lines(c(0,10),c(0,10),col="red")
  plot(sim_param[i,,2],sim_param[j,,2],main=paste0("param disp",i))
  lines(c(0,10),c(0,10),col="red")
  plot(sim_param[i,,3],sim_param[j,,3],main=paste0("param drop",i))
  lines(c(0,10),c(0,10),col="red")
  plot(sim_matrix[i,],sim_matrix[j,],main=paste0("count",i))
  lines(c(0,10),c(0,10),col="red")
}
par(op)



param12_quantile=t(apply(sim_param[,,2]/sim_param[,,1],1,function(x){quantile(x,na.rm = TRUE)}))
param12_var=t(apply(sim_param[,,2]/sim_param[,,1],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param12_mean=t(apply(sim_param[,,2]/sim_param[,,1],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

param1_quantile=t(apply(sim_param[,,1],1,function(x){quantile(x,na.rm = TRUE)}))
param1_var=t(apply(sim_param[,,1],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param1_mean=t(apply(sim_param[,,1],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

param2_quantile=t(apply(sim_param[,,2],1,function(x){quantile(x,na.rm = TRUE)}))
param2_var=t(apply(sim_param[,,2],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param2_mean=t(apply(sim_param[,,2],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

param3_quantile=t(apply(sim_param[,,3],1,function(x){quantile(x,na.rm = TRUE)}))
param3_var=t(apply(sim_param[,,3],1,function(x){var(as.numeric(x),na.rm = TRUE)}))
param3_mean=t(apply(sim_param[,,3],1,function(x){mean(as.numeric(x),na.rm = TRUE)}))

op=par(mfrow = c(4, 2))
for(i in c(1,3,5)){
  plot(1:100,param1_quantile[1:100,i],ylim=c(0,5))
  plot(2001:2100,param1_quantile[2001:2100,i],ylim=c(0,5))
}
plot(1:100,param1_var[1:100],ylim=c(0,5))
plot(2001:2100,param1_var[2001:2100],ylim=c(0,5))

for(i in c(1,3,5)){
  plot(1:100,param2_quantile[1:100,i],ylim=c(0,40))
  plot(2001:2100,param2_quantile[2001:2100,i],ylim=c(0,40))
}
plot(1:100,param2_var[1:100],ylim=c(0,80))
plot(2001:2100,param2_var[2001:2100],ylim=c(0,80))

#plot(1:100,param2_quantile[1:100,5]-param2_quantile[1:100,1],ylim=c(0,60))
#plot(2001:2100,param2_quantile[2001:2100,5]-param2_quantile[2001:2100,1],ylim=c(0,60))

for(i in c(1,3,5)){
  plot(1:100,param3_quantile[1:100,i],ylim=c(0,1))
  plot(2001:2100,param3_quantile[2001:2100,i],ylim=c(0,1))
}
plot(1:100,param3_var[1:100],ylim=c(0,1))
plot(2001:2100,param3_var[2001:2100],ylim=c(0,1))

for(i in c(1,3,5)){
  plot(1:100,param12_quantile[1:100,i],ylim=c(0,300))
  plot(2001:2100,param12_quantile[2001:2100,i],ylim=c(0,300))
}
plot(1:100,param12_var[1:100],ylim=c(0,300))
plot(2001:2100,param12_var[2001:2100],ylim=c(0,300))

for(i in c(1,3,5)){
  plot(1:100,param2_quantile[1:100,5]-param2_quantile[1:100,1],ylim=c(0,10))
  plot(2001:2100,param2_quantile[2001:2100,5]-param2_quantile[2001:2100,1],ylim=c(0,10))
}
plot(1:100,param2_var[1:100],ylim=c(0,50))
plot(2001:2100,param2_var[2001:2100],ylim=c(0,50))
par(op)


plot(sort(param2_quantile[1:100,5]),sort(param2_quantile[2001:2100,5]))
lines(c(0,30),c(0,30))


op=par(mfrow = c(4, 2))
for(i in (order(param2_quantile[1:100,5],decreasing = TRUE)[1:8])){
  j=i+2000
  hist(sim_param[i,,2])
  hist(sim_param[j,,2])
}
par(op)

op=par(mfrow = c(2, 2))
for(i in 1:10){
  j=i+1500
  heatmap(dist_array[i,,])
  print(dist_array[i,,])
  print(quantile(dist_array[i,,]))
  hist(dist_array[i,,])
  hist(dist_array[j,,])
  plot(dist_array[i,,],dist_array[j,,])
  lines(c(0,10),c(0,10),col="red")
}
par(op)

max_count=apply(sim_matrix,1,max)
head(dist_array)
dist_quantile=t(apply(dist_array,1,function(x){quantile(x,na.rm = TRUE)}))
dist_var=t(apply(dist_array,1,function(x){var(as.numeric(x),na.rm = TRUE)}))
nlog10_dist_var=-log10(dist_var)
View(dist_quantile)

op=par(mfrow = c(4, 2))
for(i in c(1,3,5)){
  plot(1:100,dist_quantile[1:100,i],ylim=c(0,0.5))
  plot(2001:2100,dist_quantile[2001:2100,i],ylim=c(0,0.5))
}
plot(1:100,nlog10_dist_var[1:100],ylim=c(0,6))
plot(2001:2100,nlog10_dist_var[2001:2100],ylim=c(0,6))
par(op)


plot(sort(dist_quantile[1:100,5]),sort(dist_quantile[2001:2100,5]))
lines(c(0,30),c(0,30))


op=par(mfrow = c(2, 2))
#param1 max
plot(dist_quantile[1:100,3],param1_quantile[1:100,3],sub=paste0("cor=", round(cor(dist_quantile[1:100,3],param1_quantile[1:100,3]),3)))
plot(dist_quantile[1:1000,3],param1_quantile[1:1000,3],sub=paste0("cor=", round(cor(dist_quantile[1:1000,3],param1_quantile[1:1000,3]),3)))
plot(dist_quantile[1501:1600,3],param1_quantile[1501:1600,3],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,3],param1_quantile[1501:1600,3]),3)))
plot(dist_quantile[,3],param1_quantile[,3],sub=paste0("cor=", round(cor(dist_quantile[,3],param1_quantile[,3]),3)))

#param 1 median
plot(dist_quantile[1:100,5],param1_quantile[1:100,5],sub=paste0("cor=", round(cor(dist_quantile[1:100,5],param1_quantile[1:100,5]),3)))
plot(dist_quantile[1:1000,5],param1_quantile[1:1000,5],sub=paste0("cor=", round(cor(dist_quantile[1:1000,5],param1_quantile[1:1000,5]),3)))
plot(dist_quantile[1501:1600,5],param1_quantile[1501:1600,5],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,5],param1_quantile[1501:1600,5]),3)))
plot(dist_quantile[,5],param1_quantile[,5],sub=paste0("cor=", round(cor(dist_quantile[,5],param1_quantile[,5]),3)))


#param1 max
plot(dist_quantile[1:100,3],param2_quantile[1:100,3],sub=paste0("cor=", round(cor(dist_quantile[1:100,3],param2_quantile[1:100,3]),3)))
plot(dist_quantile[1:1000,3],param2_quantile[1:1000,3],sub=paste0("cor=", round(cor(dist_quantile[1:1000,3],param2_quantile[1:1000,3]),3)))
plot(dist_quantile[1501:1600,3],param2_quantile[1501:1600,3],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,3],param2_quantile[1501:1600,3]),3)))
plot(dist_quantile[2001:3000,3],param2_quantile[2001:3000,3],sub=paste0("cor=", round(cor(dist_quantile[,3],param2_quantile[,3]),3)))

#param 1 median
plot(dist_quantile[1:100,5],param2_quantile[1:100,5],sub=paste0("cor=", round(cor(dist_quantile[1:100,5],param2_quantile[1:100,5]),3)))
plot(dist_quantile[1:1000,5],param2_quantile[1:1000,5],sub=paste0("cor=", round(cor(dist_quantile[1:1000,5],param2_quantile[1:1000,5]),3)))
plot(dist_quantile[1501:1600,5],param2_quantile[1501:1600,5],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,5],param2_quantile[1501:1600,5]),3)))
plot(dist_quantile[2001:3000,5],param2_quantile[2001:3000,5],sub=paste0("cor=", round(cor(dist_quantile[,5],param2_quantile[,5]),3)))

#param1 max
plot(dist_quantile[1:100,3],param12_quantile[1:100,3],sub=paste0("cor=", round(cor(dist_quantile[1:100,3],param12_quantile[1:100,3]),3)))
plot(dist_quantile[1:1000,3],param12_quantile[1:1000,3],sub=paste0("cor=", round(cor(dist_quantile[1:1000,3],param12_quantile[1:1000,3]),3)))
plot(dist_quantile[1501:1600,3],param12_quantile[1501:1600,3],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,3],param12_quantile[1501:1600,3]),3)))
plot(dist_quantile[2001:3000,3],param12_quantile[2001:3000,3],sub=paste0("cor=", round(cor(dist_quantile[,3],param12_quantile[,3]),3)))

#param 1 median
plot(dist_quantile[1:100,5],param12_quantile[1:100,5],sub=paste0("cor=", round(cor(dist_quantile[1:100,5],param12_quantile[1:100,5]),3)))
plot(dist_quantile[1:1000,5],param12_quantile[1:1000,5],sub=paste0("cor=", round(cor(dist_quantile[1:1000,5],param12_quantile[1:1000,5]),3)))
plot(dist_quantile[1501:1600,5],param12_quantile[1501:1600,5],sub=paste0("cor=", round(cor(dist_quantile[1501:1600,5],param12_quantile[1501:1600,5]),3)))
plot(dist_quantile[2001:3000,5],param12_quantile[2001:3000,5],sub=paste0("cor=", round(cor(dist_quantile[,5],param12_quantile[,5]),3)))
par(op)