#this code care with the results of simulation data, and calculate the KL

# dist_method="JSD"        #c("mean","JSD")
# fit_method="empirical"   #c("empirical","zinb")
# file_tag=1
# r_mean=1.5
# r_var=1.5
# r_disp=1.5
# r_change_prop=0.75
# n=5      #c(20,15,10,5)
# ncell=20 #c(100,80,60,40,20)

# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
# setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")
perm_label=1

perm_num=500
covariate_flag=NA #c(NA, "quantile99")
tol=0.2


t_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))


##############functions#################
library("ggplot2")
library("emdbook")

source("./Command/7.0_ZINB_fit_functions.R")
source("./Command/8.0_kl_divergence_functions.R")
source("./Command/9.0_Fstat_functions.R")


#set labels
selected_index=sample.int(20,n)
total_cell_index=matrix(ncol=1,nrow=0)
for(i_s in c(selected_index,(20+selected_index))){
  cell_index=(100*i_s-ncell+1):(100*i_s)
  total_cell_index=c(total_cell_index,cell_index)
}

#calculation
sim_matrix=t_sim_matrix[,total_cell_index]
meta=t_meta[total_cell_index,]





##########################fit sim data by individual: Refer section 7 ############################

cur_individual=unique(meta$individual)
phenotype=meta$phenotype[match(cur_individual,meta$individual)]

if(file.exists(paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))){
  dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
}

if(!file.exists(paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))){
  if(fit_method!="empirical"){
    if(!is.na(covariate_flag)){
      #logsum_count=log(apply(sim_matrix,2,sum))
      quantile99=log(apply(sim_matrix,2,function(x)return(quantile(x,0.99)+1)))
      covariate=as.matrix(quantile99)
      pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind_raw_",covariate_flag,"_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".pdf"),height = 4,width = 6)
    }
    if(is.na(covariate_flag)){
      pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind_raw_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".pdf"),height = 4,width = 6)
    }
    sim_fit=array(dim=c(nrow(sim_matrix),length(cur_individual),3),
                  dimnames = list(rownames(sim_matrix),cur_individual,c("logmean","dispersion","dropout_rate")))
    
    for(i_g in 1:nrow(sim_matrix)){
      for(i_ind in 1:length(cur_individual)){
        cur_ind=cur_individual[i_ind]
        #fit org
        cur_org_ind=sim_matrix[i_g,meta$individual==cur_ind]
        
        if(!is.na(covariate_flag)){
          cur_covariate=covariate[meta$individual==cur_ind,]
          sim_fit[i_g,i_ind,]=fit_nbzinb(cur_org_ind,cur_covariate)
        }
        if(is.na(covariate_flag)){
          sim_fit[i_g,i_ind,]=fit_nbzinb(cur_org_ind)
        }
        
        if(i_g<=10 & i_ind<=5){
          cur_org_ind=data.frame(cur_org_ind)
          colnames(cur_org_ind)="raw_count"
          ggplot(cur_org_ind, aes(x=raw_count),stat="count") + geom_histogram(fill="lightblue")+
            labs(title=paste0("Histogram of rawcount, ",rownames(sim_matrix)[i_g]," of ",cur_ind),x="Count", y = "Frequency")
          #+theme_classic()
        }
        
      }
      if(i_g%%100==0){
        print(c("ind fit",i_g))
      }
    }
    dev.off()
    
    sim_fit[,,1]=exp(sim_fit[,,1]) #change the log mean to mean!!!
    saveRDS(sim_fit,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_sim_fit_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  }
  
  
  ###################JSD and KL calculation: Refer section 8 #######################
  ###################calculation t#################################
  print(paste0("start calculation: ",fit_method,"_",dist_method))
  
  if(fit_method=="empirical"){
    dist_array=array(dim=c(nrow(sim_matrix),length(cur_individual),length(cur_individual)),
                     dimnames = list(rownames(sim_matrix),cur_individual,cur_individual))
    for(i_g in 1:nrow(sim_matrix)){
      cur_sim=sim_matrix[i_g,]
      for(i_ind_a in 1:length(cur_individual)){
        for(i_ind_b in 1:length(cur_individual)){
          cur_ind_a=cur_individual[i_ind_a]
          cur_ind_b=cur_individual[i_ind_b]
          #fit sim
          cur_sim_ind_a=as.numeric(cur_sim[meta$individual==cur_ind_a])
          cur_sim_ind_b=as.numeric(cur_sim[meta$individual==cur_ind_b])
          
          dist_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_sim_ind_a,cur_sim_ind_b,alter=dist_method,fit_model=fit_method), error = function(e) {NA} )
        }
      }
      if(i_g%%100==0){
        print(i_g)
      }
      
    }
  }
  if(fit_method!="empirical"){
    dist_array=array(dim=c(nrow(sim_fit),length(cur_individual),length(cur_individual)),
                     dimnames = list(rownames(sim_fit),cur_individual,cur_individual))
    for(i_g in 1:nrow(sim_fit)){
      cur_fit=sim_fit[i_g,,]
      for(i_ind_a in 1:length(cur_individual)){
        for(i_ind_b in 1:length(cur_individual)){
          cur_a=cur_fit[i_ind_a,]
          cur_b=cur_fit[i_ind_b,]
          #kl and jsd
          dist_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_a,cur_b,alter=dist_method,zinb.quantile=0.975,fit_model=fit_method), error = function(e) {NA} )
        }
      }
      if(i_g%%100==0){
        print(i_g)
      }
    }
  }
  
  
  saveRDS(dist_array,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
}


###################Fstat calculation: Refer section 9 #######################
if(perm_label>0){
  phenotype=phenotype[sample.int(length(phenotype))]
}

print("start manova calculation")

#real data test
#dist_array=round(dist_array,5)
dist_pval=cal_permanova_pval2(dist_array,phenotype,perm_num.min = perm_num)
thres=tol/perm_num
second_index=which(dist_pval<thres)
if(length(second_index)>0){
  sub_dist_array=dist_array[second_index,,,drop=FALSE]
  sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*10)
  dist_pval[second_index]=sub_dist_pval
  thres=tol/(perm_num*10)
  second_index=which(dist_pval<thres)
  if(length(second_index)>0){
    sub_dist_array=dist_array[second_index,,,drop=FALSE]
    sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*100)
    dist_pval[second_index]=sub_dist_pval
    thres=tol/(perm_num*100)
    second_index=which(dist_pval<thres)
    if(length(second_index)>0){
      sub_dist_array=dist_array[second_index,,,drop=FALSE]
      sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*1000)
      dist_pval[second_index]=sub_dist_pval
    }
  }
}

#dist_pval=apply(dist_array,1,function(x){tryCatch(return(cal_permanova_pval(x,phenotype)), error = function(e) {NA} )})
saveRDS(dist_pval,paste0("../Data_PRJNA434002/10.Result/p",perm_label,"_",dist_method,"_",fit_method,"_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))


sessionInfo()
q(save="no")

