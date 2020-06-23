#this code is modify from code series 7.2,8.2 and 9.1
#this code is the original code for the data analysis

cluster_tag=4
file_tag="PFC3k"
pre_tag="dca"

setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

fit_tag="" # ""(zinb) or "nb"
fit_tag2="" # "" or "raw"(directely use dca res)
resid_flag="logresid" #c("", "logresid","adj")
sim_n=10
covariate_flag="readdepth" #c(NA, "quantile99","quantile99_readdepth","readdepth")
dataset_folder="Data_PRJNA434002"  #Data_PRJNA434002   MS

perm_label_seq=0:10
dist_method_seq=c("jsd","klmean")
fit_method_seq=c("nbzinb","empirical","direct")
F_method_seq=c("p","ps")
ind_covariate_flag_seq=c("","ind") 
perm_method="" 
perm_num=500
tol=1

library("ggplot2")
library("emdbook")



###########functions#############
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/7.0_ZINB_fit_functions.R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/8.0_kl_divergence_functions.R")
source("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/9.0_Fstat_functions.R")
###########input###############

#dca input
t_mean=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
t_dispersion=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/dispersion.tsv"),stringsAsFactors = FALSE,row.names = 1)
t_dropout=read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/dropout.tsv"),stringsAsFactors = FALSE,row.names = 1)

# test
# t_mean=t_mean[1:10,]
# t_dispersion=t_dispersion[1:10,]
# t_dropout=t_dropout[1:10,]

#meta input
tmeta=readRDS(paste0("/fh/fast/sun_w/mengqi/Data_PRJNA434002/meta_PFC.rds"))

##########data processing############
#1. restrict the cluster
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
cur_cluster

sub_mean=as.matrix(t_mean[,tmeta$cluster==cur_cluster])

sub_dispersion=as.matrix(t_dispersion[,tmeta$cluster==cur_cluster])
sub_dropout=as.matrix(t_dropout[,tmeta$cluster==cur_cluster])

meta=tmeta[tmeta$cluster==cur_cluster,]

#2. generate idividual label
cur_individual=unique(meta$individual)

#########################################Section 7 fit models ################################

fit_ind_sub_sim=array(dim=c(nrow(sub_mean),length(cur_individual),3),
                      dimnames = list(rownames(sub_mean),cur_individual,c("logmean","dispersion","dropout_rate")))

sim_ind=array(dim=c(nrow(sub_mean),ncol(sub_mean),sim_n),
                      dimnames = list(rownames(sub_mean),colnames(sub_mean),1:sim_n))
sim_ind_adj=sim_ind #adj for covariate
sim_ind_logresid=sim_ind #adj for covariate

if("readdepth" %in% covariate_flag){  #only 1 covaraite allowed
  covariate=apply(sub_mean,2,function(x)return(sum(x,na.rm = TRUE)))
  log_covariate=log(covariate)
  covariate_ratio=covariate/median(covariate)
}
  
pdf(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/hist_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".pdf"),height = 4,width = 6)
for(i_g in 1:nrow(sub_mean)){
  cur_sim=matrix(ncol=sim_n,nrow=ncol(sub_mean))
  for(i_s in 1:ncol(sub_mean)){
    if(fit_tag==""){
      cur_sim[i_s,]=emdbook::rzinbinom(sim_n,sub_mean[i_g,i_s], sub_dispersion[i_g,i_s], sub_dropout[i_g,i_s])
    }
    if(fit_tag=="nb"){
      cur_sim[i_s,]=rnbinom(n=sim_n,mu=sub_mean[i_g,i_s], size=sub_dispersion[i_g,i_s])
    }
    
  }
  sim_ind[i_g,,]=cur_sim
  if(covariate_flag!=""){
    cur_sim_adj=cur_sim
    cur_sim_logresid=cur_sim
    #resid
    for(i_sim in 1:sim_n){
      cur_sim_logresid[,i_sim]=lm(log(cur_sim_logresid[,i_sim]+1)~log_covariate)$residuals
    }
    #adj
    for(i_s in 1:ncol(sub_mean)){
      cur_sim_adj[i_s,]=cur_sim_adj[i_s,]/covariate_ratio[i_s] 
    }
    sim_ind_logresid[i_g,,]=cur_sim_logresid
    sim_ind_adj[i_g,,]=cur_sim_adj
  }
  
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]

    #fit sim
    cur_sim_ind=as.numeric(cur_sim[meta$individual==cur_ind,])
    
    if(covariate_flag!=""){
      cur_covariate=rep(log_covariate[meta$individual==cur_ind],sim_n)
      if(fit_tag==""){
        fit_ind_sub_sim[i_g,i_ind,]=fit_nbzinb(cur_sim_ind,cur_covariate)
      }
      if(fit_tag=="nb"){
        fit_ind_sub_sim[i_g,i_ind,]=fit_nb(cur_sim_ind,cur_covariate)
      }

    }
    if(covariate_flag==""){
      if(fit_tag==""){
        fit_ind_sub_sim[i_g,i_ind,]=fit_nbzinb(cur_sim_ind)
      }
      if(fit_tag=="nb"){
        fit_ind_sub_sim[i_g,i_ind,]=fit_nb(cur_sim_ind)
      }
    }
    
    #plot
    if(i_g<=10 & i_ind<=5){
      cur_sim_ind=data.frame(cur_sim_ind)
      colnames(cur_sim_ind)="simulated_count"
      ggplot(cur_sim_ind, aes(x=simulated_count),stat="count") + geom_histogram(fill="lightblue")+
        labs(title=paste0("Histogram of ",pre_tag,", ",rownames(sub_mean)[i_g]," of ",cur_ind,", ",cur_cluster),x="Simulated Count", y = "Frequency")
      #+theme_classic()
    }
  }
  print(i_g)
}
dev.off()

saveRDS(fit_ind_sub_sim,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/fit_ind_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(sim_ind,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/sim_ind_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(sim_ind_adj,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/sim_ind_adj",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds")) #if covariate flag works, then adj can be defined from here
saveRDS(sim_ind_logresid,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/sim_ind_logresid",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(covariate,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/covariate_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(covariate_ratio,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/covariate_ratio_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))



###############################Section 8: distance calculation ##############################

#input counts
sim_fit=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/fit_ind_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

if(fit_tag2=="raw"){
  sim_data=as.matrix(read.table(paste0("/fh/scratch/delete90/sun_w/mengqi/Data_PRJNA434002/res_dca_rawM",file_tag,"/mean.tsv"),stringsAsFactors = FALSE))
  if(resid_flag=="adj"){
    covariate_ratio=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/covariate_ratio_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    sim_data_adj=matrix(as.numeric(sim_data)*rep(covariate_ratio,each=nrow(sim_data)),nrow=nrow(sim_data),ncol=ncol(sim_data))
    rownames(sim_data_adj)=rownames(sim_data)
    colnames(sim_data_adj)=colnames(sim_data)
    sim_data=sim_data_adj
  }
  if(resid_flag=="resid"){
    covariate=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/covariate_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    sim_data_adj=sim_data
    sim_data_adj[]=NA
    for(i_g in 1:nrow(sim_data)){
      sim_data_adj[ig,]=lm(log(sim_data[ig,]+0.1)~log(covariate))$residuals
    }
    sim_data=sim_data_adj
  }
  
  sim_data=array(sim_data,dim=c(dim(sim_data),1))
}
if(fit_tag2!="raw"){
  if(resid_flag==""){
    sim_data=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/sim_ind_",covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  }
  if(resid_flag!=""){
    sim_data=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/sim_ind_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  }
}


print("start calculation: Part I: Empirical KLmean and JSD")

klmean_empirical_array=array(dim=c(nrow(sim_data),length(cur_individual),length(cur_individual)),
                             dimnames = list(rownames(sim_data),cur_individual,cur_individual))
jsd_empirical_array=array(dim=c(nrow(sim_data),length(cur_individual),length(cur_individual)),
                          dimnames = list(rownames(sim_data),cur_individual,cur_individual))

for(i_g in 1:nrow(sim_data)){
  cur_sim=sim_data[i_g,,,drop=FALSE]
  for(i_ind_a in 1:length(cur_individual)){
    for(i_ind_b in 1:length(cur_individual)){
      cur_ind_a=cur_individual[i_ind_a]
      cur_ind_b=cur_individual[i_ind_b]
      #fit sim
      cur_sim_ind_a=as.numeric(cur_sim[,meta$individual==cur_ind_a,])
      cur_sim_ind_b=as.numeric(cur_sim[,meta$individual==cur_ind_b,])
      klmean_empirical_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_sim_ind_a,cur_sim_ind_b,alter="mean",fit_model="empirical"), error = function(e) {NA} )
      jsd_empirical_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_sim_ind_a,cur_sim_ind_b,alter="JSD",fit_model="empirical"), error = function(e) {NA} )
    }
  }
  print(i_g)
}

print("start calculation: Part II: model fitted KLmean and JSD")

klmean_nbzinb_array=array(dim=c(nrow(sim_fit),length(cur_individual),length(cur_individual)),
                          dimnames = list(rownames(sim_fit),cur_individual,cur_individual))
jsd_nbzinb_array=array(dim=c(nrow(sim_fit),length(cur_individual),length(cur_individual)),
                       dimnames = list(rownames(sim_fit),cur_individual,cur_individual))

sim_fit[,,1]=exp(sim_fit[,,1]) #change the log mean to mean!!!

for(i_g in 1:nrow(sim_fit)){
  cur_fit=sim_fit[i_g,,]
  for(i_ind_a in 1:length(cur_individual)){
    for(i_ind_b in 1:length(cur_individual)){
      cur_a=cur_fit[i_ind_a,]
      cur_b=cur_fit[i_ind_b,]
      #kl and jsd
      klmean_nbzinb_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_a,cur_b,alter="mean",zinb.quantile=0.975,fit_model="zinb"), error = function(e) {NA} )
      jsd_nbzinb_array[i_g,i_ind_a,i_ind_b]=tryCatch(mean_KL_dens(cur_a,cur_b,alter="JSD",zinb.quantile=0.975,fit_model="zinb"), error = function(e) {NA} )
    }
  }
  print(i_g)
}

print("calculation end")

if(fit_tag2=="raw"){
  fit_tag=fit_tag2
}



print("start calculation: Part III: direct KLmean and JSD")
if(covariate_flag!=""){
  if(resid_flag=="adj"){
    covariate_ratio=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/covariate_ratio_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    sub_mean_adj=matrix(as.numeric(sub_mean)*rep(covariate_ratio,each=nrow(sub_mean)),nrow=nrow(sub_mean),ncol=ncol(sub_mean))
    rownames(sub_mean_adj)=rownames(sub_mean)
    colnames(sub_mean_adj)=colnames(sub_mean)
    sub_mean=sub_mean_adj
  }
  if(resid_flag=="resid"){
    covariate=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/covariate_",covariate_flag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
    sub_mean_adj=sub_mean
    sub_mean_adj[]=NA
    for(i_g in 1:nrow(sub_mean)){
      sub_mean_adj[ig,]=lm(log(sub_mean[ig,]+0.1)~log(covariate))$residuals
    }
    sub_mean=sub_mean_adj
  }
}

klmean_direct_array=array(dim=c(nrow(sub_mean),length(cur_individual),length(cur_individual)),
                          dimnames = list(rownames(sub_mean),cur_individual,cur_individual))
jsd_direct_array=array(dim=c(nrow(sub_mean),length(cur_individual),length(cur_individual)),
                       dimnames = list(rownames(sub_mean),cur_individual,cur_individual))

for(i_g in 1:nrow(sub_mean)){
  if(fit_tag==""){
    cell_param=cbind(sub_mean[i_g,], sub_dispersion[i_g,], sub_dropout[i_g,])
  }
  if(fit_tag=="nb"){
    cell_param=cbind(sub_mean[i_g,], sub_dispersion[i_g,], 0)
  }
  klmean_direct_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter="mean"), error = function(e) {NA} )
  jsd_direct_array[i_g,,]=tryCatch(mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter="JSD"), error = function(e) {NA} )
  print(i_g)
}

###################calculation end, output#################################
print("calculation end")

saveRDS(klmean_empirical_array,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/klmean_empirical_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(jsd_empirical_array,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/jsd_empirical_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(klmean_nbzinb_array,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/klmean_nbzinb_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(jsd_nbzinb_array,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/jsd_nbzinb_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

saveRDS(klmean_direct_array,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/klmean_direct_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
saveRDS(jsd_direct_array,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/jsd_direct_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

###############################Section 9: Fstat ############################################


for(dist_method in dist_method_seq){
  for(fit_method in fit_method_seq){
    for(F_method in F_method_seq){
      for(ind_covariate_flag in ind_covariate_flag_seq){
        if(file.exists(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/",dist_method,"_",fit_method,"_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))){
          dist_array=NA
          dist_pval=NA
          
          dist_array=readRDS(paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/",dist_method,"_",fit_method,"_array_",resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          
          phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
          phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0
          
          pval_perm=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])
          rownames(pval_perm)=dimnames(dist_array)[[1]]
          anova_centroid_pval_perm=pval_perm
          
          for(perm_label in perm_label_seq){
            
            #set perm
            if(perm_label>0){
              phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
              phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0
              
              total_individual=unique(tmeta$individual)
              perm_order=readRDS(paste0("../",dataset_folder,"/7.Result/ind_perm_order.rds"))
              perm_order=as.numeric(perm_order[,perm_label])
              total_individual_ref=total_individual[perm_order]
              perm_order=match(total_individual_ref,cur_individual)
              perm_order=perm_order[!is.na(perm_order)]
              phenotype=phenotype[perm_order,drop=FALSE]
            }
            
            #set covariate
            if(ind_covariate_flag==""){
              covariate_model_matrix=NA
            }
            if(ind_covariate_flag=="ind"){
              #
              cur_covariate=meta[match(cur_individual,meta$individual),c("age","sex","RNA.Integrity.Number","Seqbatch")]
              #read depth
              if(length(grep("rdpadj",file_tag))==0){ #if not justed before
                read_depth=readRDS(paste0("../",dataset_folder,"/rawM",file_tag,"_read_depth_per_1Kcell_ind.rds"))
                read_depth=read_depth[match(cur_individual,rownames(read_depth)),]
                cur_covariate=cbind(cur_covariate, read_depth)
              }
              rownames(cur_covariate)=cur_individual
              covariate_model_matrix=model.matrix(~.,cur_covariate)
            }
            
            dist_res=cal_permanova_pval2(dist_array,phenotype,Fstat_method=F_method,perm_num.min = perm_num,zm=covariate_model_matrix)
            dist_pval=dist_res$pval
            F_ob0=dist_res$F_ob
            F_perm0=dist_res$F_perm
            
            
            #print("0 level complete")
            thres=tol/perm_num
            second_index=which(dist_pval<thres)
            if(length(second_index)>0){
              #print("1st level")
              sub_dist_array=dist_array[second_index,,,drop=FALSE]
              sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*10,Fstat_method=F_method,zm=covariate_model_matrix)$pval
              dist_pval[second_index]=sub_dist_pval
              thres=tol/(perm_num*10)
              second_index=which(dist_pval<thres)
              if(length(second_index)>0){
                #print("2nd level")
                sub_dist_array=dist_array[second_index,,,drop=FALSE]
                sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*100,Fstat_method=F_method,zm=covariate_model_matrix)$pval
                dist_pval[second_index]=sub_dist_pval
                thres=tol/(perm_num*100)
                second_index=which(dist_pval<thres)
                if(length(second_index)>0){
                  print("3rd level")
                  sub_dist_array=dist_array[second_index,,,drop=FALSE]
                  sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*1000,Fstat_method=F_method,zm=covariate_model_matrix)$pval
                  dist_pval[second_index]=sub_dist_pval
                }
              }
            }
            ##here we give the NAs a second chance by removing missing samples/ missing distances, 
            # or fix the missing samples with median distances.
            
            second_index=which(is.na(dist_pval))
            tol_missing_dist=0 # tolerate missing p-values number, if missing numer is no bigger than it, we will make it up with median of distance.
            tol_missing_sample=dim(dist_array)[2]/2 #if effective sample less than this number, stop calculation
            
            for(i2 in second_index){
              print(i2)
              x=dist_array[i2,,]
              #calculate zeros
              zero_sum=apply(is.na(x),2,sum)
              
              #first thres: to remove all zero inds
              flag=(zero_sum<nrow(x))
              if(sum(flag)>=tol_missing_sample){
                x=x[flag,flag]
                cur_pheno=phenotype[flag]
                
                #second thres:to remove inds with more than tolerate missing values
                zero_sum=apply(is.na(x),2,sum)
                flag=(zero_sum<=tol_missing_dist)
                if(sum(flag)>=tol_missing_sample){ 
                  #third thres:to 
                  x=x[flag,flag]
                  cur_pheno=cur_pheno[flag]
                  #add missing values:
                  fill_index=which(!complete.cases(x))
                  if(length(fill_index)>0){
                    for(i_f in fill_index){
                      for(j_f in fill_index){
                        if(j_f>i_f){
                          x[i_f,j_f]=median(c(x[,i_f],x[,j_f]),na.rm = TRUE) #here is a little recurrence, but that's OK...
                          x[j_f,i_f]=x[i_f,j_f]
                        }
                      }
                    }
                  }
                  dist_pval[i2]=tryCatch(cal_permanova_pval(x,cur_pheno,Fstat_method=F_method,perm_num.min = perm_num,zm=covariate_model_matrix)$pval, error = function(e) {NA} )
                  thres=tol/perm_num
                  if(!is.na(dist_pval[i2]) && dist_pval[i2]<thres){
                    dist_pval[i2]=tryCatch(cal_permanova_pval(x,cur_pheno,Fstat_method=F_method,perm_num.min = (perm_num*10),zm=covariate_model_matrix)$pval, error = function(e) {NA} )
                    thres=tol/(perm_num*10)
                    if(!is.na(dist_pval[i2]) && dist_pval[i2]<thres){
                      dist_pval[i2]=tryCatch(cal_permanova_pval(x,cur_pheno,Fstat_method=F_method,perm_num.min = (perm_num*100),zm=covariate_model_matrix)$pval, error = function(e) {NA} )
                      thres=tol/(perm_num*100)
                      if(!is.na(dist_pval[i2]) && dist_pval[i2]<thres){
                        dist_pval[i2]=tryCatch(cal_permanova_pval(x,cur_pheno,Fstat_method=F_method,perm_num.min = (perm_num*1000),zm=covariate_model_matrix)$pval, error = function(e) {NA} )
                      }
                    }
                  }
                }
              }
            }
            
            if(perm_label==0){
              pval_ob=dist_pval
            }
            if(perm_label>0){
              pval_perm[,perm_label]=dist_pval
              print(perm_label)
            }
          }
          
          rownames(pval_ob)=dimnames(dist_array)[[1]]
          rownames(pval_perm)=dimnames(dist_array)[[1]]
          
          saveRDS(pval_ob,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          saveRDS(pval_perm,paste0("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/Command/11.1.mengqi_code_check/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_perm_pval_",ind_covariate_flag,resid_flag,covariate_flag,fit_tag,pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
          gc()
          print(paste0(dist_method,"_",fit_method,"_",F_method))
        }
      }
    }
  }
}





sessionInfo()
#q(save="no")
