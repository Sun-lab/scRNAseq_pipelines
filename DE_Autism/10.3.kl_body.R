#this code care with the results of simulation data, and calculate the KL

# dist_method="JSD"        #c("mean","JSD")
# fit_method="empirical"   #c("empirical","zinb")
# file_tag=1
# r_mean=1.5
# r_var=1.5
# r_disp=1.5
# r_change_prop=0.4
# n=5      #c(20,15,10,5)
# ncell=20 #c(100,80,60,40,20) +10,5,3

# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
# setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#perm_label=1
#perm_method="b" #c("","b","bs","s") #b means balanced permutation #s means remove samples with limited cells.
cell_thres=5 #

perm_num=500
covariate_flag=NA #c(NA, "quantile99")
tol=0.2
perm_label_seq=0:10

t_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_matrix_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
t_sim_param=readRDS(paste0("../Data_PRJNA434002/10.Result/sim_param_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,".rds"))
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
sim_param=t_sim_param[,total_cell_index,]
meta=t_meta[total_cell_index,]


##########################fit sim data by individual: Refer section 7 ############################

cur_individual=unique(meta$individual)
phenotype=meta$phenotype[match(cur_individual,meta$individual)]

if(file.exists(paste0("../Data_PRJNA434002/10.Result/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))){
  dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
}

if(!file.exists(paste0("../Data_PRJNA434002/10.Result/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))){
  if(fit_method=="zinb"){
    if(!is.na(covariate_flag)){
      #logsum_count=log(apply(sim_matrix,2,sum))
      quantile99=log(apply(sim_matrix,2,function(x)return(quantile(x,0.99)+1)))
      covariate=as.matrix(quantile99)
      pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind/hist_sim_ind_raw_",covariate_flag,"_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".pdf"),height = 4,width = 6)
    }
    if(is.na(covariate_flag)){
      pdf(paste0("../Data_PRJNA434002/10.Result/hist_sim_ind/hist_sim_ind_raw_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".pdf"),height = 4,width = 6)
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
  if(fit_method=="zinb"){
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
  
  if(fit_method=="direct"){
    dist_array=array(dim=c(nrow(sim_matrix),length(cur_individual),length(cur_individual)),
                     dimnames = list(rownames(sim_matrix),cur_individual,cur_individual))
    
    for(i_g in 1:nrow(sim_matrix)){
      cell_param=sim_param[i_g,,]
      dist_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter=dist_method), error = function(e) {NA} )
      if(i_g%%100==0){
        print(i_g)
      }
    }
  }
  
  saveRDS(dist_array,paste0("../Data_PRJNA434002/10.Result/dist_array/",dist_method,"_",fit_method,"_dist_array_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
}


###################Fstat calculation: Refer section 9 #######################
pval_perm=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])
perm_tukey_centroid_pval=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])
perm_anova_centroid_pval=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])
perm_tukey_median_pval=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])
perm_anova_median_pval=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])

perm_order=matrix(ncol=(length(perm_label_seq)-1),nrow=length(phenotype))

for(perm_label in perm_label_seq){
  #for(perm_method in c("","b")){
   perm_method="" 
    
    dist_res=NA
    dist_pval=NA
    tukey_centroid_pval=NA
    anova_centroid_pval=NA
    tukey_median_pval=NA
    anova_median_pval=NA
    
    cur_phenotype=phenotype
    if(perm_label>0){
      if(length(grep("b",perm_method))==0){ #naive permutation
        cur_order=sample.int(length(cur_phenotype))
        cur_phenotype=cur_phenotype[cur_order]
        perm_order[,perm_label]=cur_order
        
      }
      if(length(grep("b",perm_method))>0){ #balanced permutation
        n_exchange=n/2
        n_exchange=floor(n_exchange)+(n_exchange-floor(n_exchange))*2*rbinom(1,1,0.5) #the number changed to other side
        
        i_exchange=sample.int(n,n_exchange)
        i_exchange=c(i_exchange,n+i_exchange)
        cur_phenotype[i_exchange]=1-cur_phenotype[i_exchange]
        
      }
    }
    
    
    #print("start manova calculation")
    
    #real data test
    #dist_array=round(dist_array,5)
    dist_res=cal_permanova_pval2(dist_array,cur_phenotype,perm_num.min = perm_num)
    dist_pval=dist_res$pval
    F_ob0=dist_res$F_ob
    F_perm0=dist_res$F_perm
    
    #plot 
    
    
    
    # png(paste0("../Data_PRJNA434002/10.Result/fig_Fstat/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,"_Fperm.png"),height = 900,width = 1800,type="cairo")
    # op=par(mfrow = c(2, 4))
    # pval_order=order(dist_pval)[c(1:4,100,200,300,400)]
    # for(iplot in pval_order){
    #   hist(F_perm0[iplot,],xlim=range(min(F_perm0[iplot,],F_ob0[iplot]),max(F_perm0[iplot,],F_ob0[iplot])),main=paste0("Fstat distr of gene ",dimnames(dist_array)[[1]][iplot],", pval ",round(dist_pval[iplot],3),", rank ",iplot),breaks=50)
    #   abline(v=F_ob0[iplot],col="red")
    # }
    # par(op)
    # dev.off()
    # 
    # png(paste0("../Data_PRJNA434002/10.Result/fig_Fstat/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,"_Fob.png"),height = 600,width = 1200,type="cairo")
    # hist(F_ob0,main=paste0("Fstat distr of gene, Fob"),breaks=50)
    # dev.off()
    # 
    
    thres=tol/perm_num
    second_index=which(dist_pval<thres)
    if(length(second_index)>0){
      sub_dist_array=dist_array[second_index,,,drop=FALSE]
      sub_dist_pval=cal_permanova_pval2(sub_dist_array,cur_phenotype,perm_num.min = perm_num*10)$pval
      dist_pval[second_index]=sub_dist_pval
      thres=tol/(perm_num*10)
      second_index=which(dist_pval<thres)
      if(length(second_index)>0){
        sub_dist_array=dist_array[second_index,,,drop=FALSE]
        sub_dist_pval=cal_permanova_pval2(sub_dist_array,cur_phenotype,perm_num.min = perm_num*100)$pval
        dist_pval[second_index]=sub_dist_pval
        thres=tol/(perm_num*100)
        second_index=which(dist_pval<thres)
        if(length(second_index)>0){
          sub_dist_array=dist_array[second_index,,,drop=FALSE]
          sub_dist_pval=cal_permanova_pval2(sub_dist_array,cur_phenotype,perm_num.min = perm_num*1000)$pval
          dist_pval[second_index]=sub_dist_pval
        }
      }
    }
    
    #dist_pval=apply(dist_array,1,function(x){tryCatch(return(cal_permanova_pval(x,cur_phenotype)), error = function(e) {NA} )})
    if(perm_label==0){
      saveRDS(dist_pval,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_raw_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
      
      
      #########start pre analysis calculation #########################
      tukey_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
      anova_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )}) 
      tukey_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
      anova_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="anova"), error = function(e) {NA} )}) 
      saveRDS(tukey_centroid_pval,paste0("../Data_PRJNA434002/10.Result/tukey_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_tukey_centroid_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
      saveRDS(anova_centroid_pval,paste0("../Data_PRJNA434002/10.Result/anova_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_anova_centroid_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
      saveRDS(tukey_median_pval,paste0("../Data_PRJNA434002/10.Result/tukey_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_tukey_median_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
      saveRDS(anova_median_pval,paste0("../Data_PRJNA434002/10.Result/anova_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_anova_median_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
      
    }
    
    if(perm_label>0){
      pval_perm[,perm_label]=dist_pval
      print(perm_label)
      
      perm_tukey_centroid_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
      perm_anova_centroid_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )}) 
      perm_tukey_median_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
      perm_anova_median_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="anova"), error = function(e) {NA} )})
    }
  #}
}



if(perm_label>0){
  saveRDS(pval_perm,paste0("../Data_PRJNA434002/10.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  saveRDS(perm_order,paste0("../Data_PRJNA434002/10.Result/perm_order/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  
  
  saveRDS(perm_tukey_centroid_pval,paste0("../Data_PRJNA434002/10.Result/tukey_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_tukey_centroid_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  saveRDS(perm_anova_centroid_pval,paste0("../Data_PRJNA434002/10.Result/anova_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_anova_centroid_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  saveRDS(perm_tukey_median_pval,paste0("../Data_PRJNA434002/10.Result/tukey_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_tukey_median_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  saveRDS(perm_anova_median_pval,paste0("../Data_PRJNA434002/10.Result/anova_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_anova_median_pval_",r_mean,"_",r_var,"_",r_disp,"_",r_change_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
  
  
  
}





sessionInfo()
q(save="no")

