#this code care with the results of simulation data, and calculate the KL

# file_tag=1
# dist_method="JSD"
# fit_method="empirical"
# r_mean=1.2
# r_var=1.1
# r_change_prop=0.2
# dp_minor_prop=0.2
# n=30
# ncell=20

# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
# setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#perm_label=1
#perm_method="b" #c("","b","bs","s") #b means balanced permutation #s means remove samples with limited cells.

library("abind")

#####################parameters from simulation 
total_ncase=50  
total_ncell=800
###############################Other parameters
cell_thres=5 #

perm_num=500
covariate_flag=NA #c(NA, "quantile99")
tol=10
perm_label_seq=0:10
perm_label_seq=0:1
sub_index_num=30

sim_folder="sim_v6"
pre_tag_seq=c("", "dca")
fit_tag_seq=c("", "nb") 
covariate_flag_seq=c("","readdepth")
resid_flag_seq=c("","resid")
#pre_tag_seq="" #c("" "dca")
#fit_tag_seq="" #c("" "nb") 
#covariate_flag_seq="readdepth" #c("","readdepth","q50","q75","q100")
#resid_flag_seq=""#c("","resid")

#pre_tag_seq="dca"
#noise_flag="0.01_10" #("",0.05_10,0.01_10,0.01_5,0.05_5)
noise_flag=""

for(pre_tag in pre_tag_seq){
  for(fit_tag in fit_tag_seq){
    for(covariate_flag in covariate_flag_seq){
      for(resid_flag in resid_flag_seq){
        print(c(pre_tag,fit_tag,covariate_flag,resid_flag))
        gc()
        meta=NULL
        dist_array=NULL
        
        if(sim_folder=="sim_v3" && fit_method=="direct"){
          q(save="no")
        }
        
        t_meta=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
        if(pre_tag==""){
          t_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
          t_sim_param=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_param_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
          
        }
        
        if(pre_tag=="dca"){
          t_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_ind_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_1.rds")) #actually sim_data, with dim= n_gene, n_cell, sim_n
          t_sim_param=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_param_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_1.rds"))
          for(sub_index in 2:sub_index_num){
            sub_sim_matrix=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_ind_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",sub_index,".rds")) #actually sim_data, with dim= n_gene, n_cell, sim_n
            sub_sim_param=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_param_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",sub_index,".rds"))
            t_sim_matrix=abind(t_sim_matrix,sub_sim_matrix,along=1)
            t_sim_param=abind(t_sim_param,sub_sim_param,along=1)
          }
          dimnames(t_sim_param)=list(rownames(t_sim_matrix),colnames(t_sim_matrix),c("mean","overdisp","dropout"))
          
        }
        
        if(covariate_flag=="readdepth"){
          t_covariate=as.matrix(readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_cell_readdepth_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds")))
          
        }
        if(covariate_flag=="q50"){
          t_covariate=as.matrix(readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_cell_count_quantile_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds")))[3,]
          
        }
        if(covariate_flag=="q75"){
          t_covariate=as.matrix(readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_cell_count_quantile_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds")))[4,]
          
        }
        if(covariate_flag=="q100"){
          t_covariate=as.matrix(readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_cell_count_quantile_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds")))[5,]
          
        }
        
        
        if(noise_flag!=""){
          sim_matrix_extreme_flag=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_extreme_flag_",noise_flag,"_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
          if(length(dim(t_sim_matrix))==2){
            t_sim_matrix=t_sim_matrix*sim_matrix_extreme_flag
          }
          if(length(dim(t_sim_matrix))==3){
            for(i_sim in 1:dim(t_sim_matrix)[[3]]){
              t_sim_matrix[,,i_sim] =t_sim_matrix[,,i_sim]*sim_matrix_extreme_flag
            }
          }
        }
        
        ##############functions#################
        library("ggplot2")
        library("emdbook")
        
        source("./Command/7.0_ZINB_fit_functions.R")
        source("./Command/8.0_kl_divergence_functions.R")
        source("./Command/9.0_Fstat_functions.R")
        
        
        #set labels
        
        selected_index=sample.int(total_ncase,n)
        total_cell_index=matrix(ncol=1,nrow=0)
        for(i_s in c(selected_index,(total_ncase+selected_index))){
          cell_index=(total_ncell*i_s-ncell+1):(total_ncell*i_s)
          total_cell_index=c(total_cell_index,cell_index)
        }
        
        #calculation
        covariate=matrix(1,ncol=1,nrow=length(total_cell_index))
        
        if(pre_tag==""){
          sim_matrix=t_sim_matrix[,total_cell_index]
        }
        if(pre_tag=="dca"){
          sim_matrix=t_sim_matrix[,total_cell_index,]
        }
        if(covariate_flag!=""){
          covariate=t_covariate[total_cell_index,]
        }
        
        meta=t_meta[total_cell_index,]
        
        if(sim_folder!="sim_v3"){
          sim_param=t_sim_param[,total_cell_index,]
        }
        

        
        ########################## resid pre-covariate adjustment ############################
        if(fit_tag!="zinb"){
          if(resid_flag=="resid"){
            sim_matrix_y  = log(sim_matrix+1)
            sim_param_y  = log(sim_param[,,1])
            
            sim_matrix_resid   =sim_matrix_y
            sim_matrix_resid[] =NA
            sim_param_resid    =sim_param_y
            sim_param_resid[]  =NA
            
            
            for(ig in 1:nrow(sim_matrix_y)){
              if(length(dim(sim_matrix_y))==2){
                sim_matrix_resid[ig,] = lm(sim_matrix_y[ig,]~log(covariate))$residuals
              }
              if(length(dim(sim_matrix_y))==3){
                for(i_sim in 1:dim(sim_matrix_y)[[3]]){
                  sim_matrix_resid[ig,,i_sim] = lm(sim_matrix_y[ig,,i_sim]~log(covariate))$residuals
                }
              }
              sim_param_resid[ig,]  = lm(sim_param_y[ig,]~log(covariate))$residuals
            }
            
            sim_matrix=exp(sim_matrix_resid) # only for empirical method, we didn't log transformed back 
            sim_param[,,1]=exp(sim_param_resid)
            
          }
        }
        
        if(resid_flag=="adj"){
          if(covariate_flag!=""){
            for(ig in 1:nrow(sim_matrix)){
              if(length(dim(sim_matrix))==2){
                sim_matrix[ig,] = sim_matrix[ig,]/covariate*10000
              }
              if(length(dim(sim_matrix))==3){
                for(i_sim in 1:dim(sim_matrix_y)[[3]]){
                  sim_matrix[ig,,i_sim] = sim_matrix[ig,,i_sim]/covariate*10000
                }
              }
              sim_param[ig,,1]  = sim_param[ig,,1]/covariate*10000
            }
          }
        }
        
        ##########################fit sim data by individual: Refer section 7 ############################
        
        cur_individual=unique(meta$individual)
        phenotype=meta$phenotype[match(cur_individual,meta$individual)]
        
        
        if(file.exists(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))){
          dist_array=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
        }
        
        if(!file.exists(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))){
          
          #cell level covariate
          if(pre_tag==""){
            
            if(fit_method=="zinb"){
              pdf(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/hist_sim_ind/hist_sim_ind_raw_",resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".pdf"),height = 4,width = 6)
              
              sim_fit=array(dim=c(nrow(sim_matrix),length(cur_individual),3),
                            dimnames = list(rownames(sim_matrix),cur_individual,c("logmean","dispersion","dropout_rate")))
              
              for(i_g in 1:nrow(sim_matrix)){
                for(i_ind in 1:length(cur_individual)){
                  cur_ind=cur_individual[i_ind]
                  #fit org
                  cur_org_ind=sim_matrix[i_g,meta$individual==cur_ind]
                  
                  if(covariate_flag!=""){
                    cur_covariate=covariate[meta$individual==cur_ind]
                    sim_fit[i_g,i_ind,]=fit_nbzinb(cur_org_ind,cur_covariate)
                  }
                  if(covariate_flag==""){
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
              saveRDS(sim_fit,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/",dist_method,"_",fit_method,"_sim_fit_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
            }
            
            ###################calculation t#################################
            print(paste0("start calculation: ",fit_method,"_",dist_method))
            
            dist_array=array(dim=c(nrow(sim_matrix),length(cur_individual),length(cur_individual)), dimnames = list(rownames(sim_matrix),cur_individual,cur_individual))
            
            if(fit_method=="empirical"){
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
              for(i_g in 1:nrow(sim_matrix)){
                cell_param=sim_param[i_g,,]
                dist_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter=dist_method), error = function(e) {NA} )
                if(i_g%%100==0){
                  print(i_g)
                }
              }
            }
          }
          if(pre_tag=="dca"){
            if(fit_method=="zinb"){
              sim_fit=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/fit_ind_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_1.rds"))
              for(sub_index in 2:sub_index_num){
                sub_sim_fit=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/fit_ind_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",sub_index,".rds"))
                sim_fit=abind(sim_fit,sub_sim_fit,along=1)
              }
              dimnames(sim_fit)=list(rownames(t_sim_matrix),dimnames(sub_sim_fit)[[2]],c("logmean","dispersion","dropout_rate"))
            }
            
            ###################JSD and KL calculation: Refer section 8 #######################
            print(paste0("start calculation: ",fit_method,"_",dist_method))
            dist_array=array(dim=c(nrow(sim_matrix),length(cur_individual),length(cur_individual)),  dimnames = list(rownames(sim_matrix),cur_individual,cur_individual))
            
            if(fit_method=="empirical"){
              for(i_g in 1:nrow(sim_matrix)){
                cur_sim=sim_matrix[i_g,,]
                for(i_ind_a in 1:length(cur_individual)){
                  for(i_ind_b in 1:length(cur_individual)){
                    cur_ind_a=cur_individual[i_ind_a]
                    cur_ind_b=cur_individual[i_ind_b]
                    #fit sim
                    cur_sim_ind_a=as.numeric(cur_sim[meta$individual==cur_ind_a,])
                    cur_sim_ind_b=as.numeric(cur_sim[meta$individual==cur_ind_b,])
                    
                    dist_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_sim_ind_a,cur_sim_ind_b,alter=dist_method,fit_model=fit_method), error = function(e) {NA} )
                  }
                }
                if(i_g%%100==0){
                  print(i_g)
                }
              }
            }
            if(fit_method=="zinb"){
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
              for(i_g in 1:nrow(sim_matrix)){
                cell_param=sim_param[i_g,,]
                dist_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter=dist_method), error = function(e) {NA} )
                if(i_g%%100==0){
                  print(i_g)
                }
              }
            }
            
            
            saveRDS(dist_array,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dist_array/",dist_method,"_",fit_method,"_dist_array_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          }
          
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
              cur_order=sample.int(length(cur_phenotype),length(cur_phenotype))
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
          
          
          
          # png(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/fig_Fstat/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,"_Fperm.png"),height = 900,width = 1800,type="cairo")
          # op=par(mfrow = c(2, 4))
          # pval_order=order(dist_pval)[c(1:4,100,200,300,400)]
          # for(iplot in pval_order){
          #   hist(F_perm0[iplot,],xlim=range(min(F_perm0[iplot,],F_ob0[iplot]),max(F_perm0[iplot,],F_ob0[iplot])),main=paste0("Fstat distr of gene ",dimnames(dist_array)[[1]][iplot],", pval ",round(dist_pval[iplot],3),", rank ",iplot),breaks=50)
          #   abline(v=F_ob0[iplot],col="red")
          # }
          # par(op)
          # dev.off()
          # 
          # png(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/fig_Fstat/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,"_Fob.png"),height = 600,width = 1200,type="cairo")
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
            saveRDS(dist_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_raw_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
            
            
            # #########start pre analysis calculation #########################
            # tukey_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
            # anova_centroid_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )}) 
            # tukey_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
            # anova_median_pval=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="anova"), error = function(e) {NA} )}) 
            # saveRDS(tukey_centroid_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/tukey_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_tukey_centroid_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
            # saveRDS(anova_centroid_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/anova_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_anova_centroid_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
            # saveRDS(tukey_median_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/tukey_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_tukey_median_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
            # saveRDS(anova_median_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/anova_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_anova_median_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
            
          }
          
          if(perm_label>0){
            pval_perm[,perm_label]=dist_pval
            print(perm_label)
            # 
            # perm_tukey_centroid_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="TukeyHSD"), error = function(e) {NA} )})
            # perm_anova_centroid_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )}) 
            # perm_tukey_median_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="TukeyHSD"), error = function(e) {NA} )})
            # perm_anova_median_pval[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "median",div_method="anova"), error = function(e) {NA} )})
          }
          #}
        }
        
        
        
        if(perm_label>0){
          saveRDS(pval_perm,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          saveRDS(perm_order,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/perm_order/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          
          # 
          # saveRDS(perm_tukey_centroid_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/tukey_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_tukey_centroid_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          # saveRDS(perm_anova_centroid_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/anova_centroid_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_anova_centroid_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          # saveRDS(perm_tukey_median_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/tukey_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_tukey_median_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          # saveRDS(perm_anova_median_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/anova_median_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_perm_anova_median_pval_",noise_flag,resid_flag,covariate_flag,pre_tag,fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"_",(2*n),"_",ncell,".rds"))
          
          
        }
        
      }
    }
  }
}






sessionInfo()
q(save="no")

