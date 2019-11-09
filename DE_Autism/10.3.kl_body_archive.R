#this code care with the results of simulation data, and calculate the KL


# file_tag=1
# sim_method="zinb.naive" #splat.mean or splat.var--method 3, separate the mean and variance using splat
#                         #splat.org--method 4, change the mean.shape and mean.rate originally
# #zinb.naive--method 5, using naive zinb models to do so.
# 
# dist_method="JSD"        #c("mean","JSD")
# fit_method="empirical"   #c("empirical","zinb")


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#r_mean=1.5  #r_mean/r_var should < 1+mean.shape
#r_var=4

perm_num=500
covariate_flag=NA #c(NA, "quantile99")

sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_meta_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))

cur_info=meta[,c("individual","phenotype")]
cur_info=unique(cur_info)
rownames(cur_info)=cur_info$individual
phenotype=as.numeric(cur_info[,2])-1
#sim_matrix=sim_matrix[1:10,]
##############functions#################
library("ggplot2")
library("emdbook")

source("./Command/7.0_ZINB_fit_functions.R")
source("./Command/8.0_kl_divergence_functions.R")
source("./Command/9.0_Fstat_functions.R")

##########################fit sim data by individual: Refer section 7 ############################
cur_individual=unique(meta$individual)

if(fit_method!="empirical"){
  if(!is.na(covariate_flag)){
    #logsum_count=log(apply(sim_matrix,2,sum))
    quantile99=log(apply(sim_matrix,2,function(x)return(quantile(x,0.99)+1)))
    covariate=as.matrix(quantile99)
    pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind_raw_",covariate_flag,"_",file_tag,".pdf"),height = 4,width = 6)
  }
  if(is.na(covariate_flag)){
    pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind_raw_",file_tag,".pdf"),height = 4,width = 6)
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
    print(c("ind fit",i_g))
  }
  dev.off()
  
  sim_fit[,,1]=exp(sim_fit[,,1]) #change the log mean to mean!!!
  saveRDS(sim_fit,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_sim_fit_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
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
    print(i_g)
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
    print(i_g)
  }
}


saveRDS(dist_array,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_dist_array_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
###################Fstat calculation: Refer section 9 #######################
print("start permanova calculation")

#real data test
dist_array=round(dist_array,5)
dist_pval=cal_permanova_pval2(dist_array,phenotype)
#dist_pval=apply(dist_array,1,function(x){tryCatch(return(cal_permanova_pval(x,phenotype)), error = function(e) {NA} )})
saveRDS(dist_pval,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_raw_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))











