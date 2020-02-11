#this code calculate the permutation based manova or permanovaS
#The input are the dist_array from step 8.

# cluster_tag=1
# file_tag="3k10"
# pre_tag="dca"
# dist_method="jsd"
# fit_method="nbzinb"
# F_method="p"

perm_label_seq=0:10
perm_method="" # "" or "s"(set cell threshold, any individuals with less cells would be removed)
cell_threshold=4

ind_covariate_flag="ind" 
perm_num=500
tol=0.2
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


###########functions#############
source("./Command/9.0_Fstat_functions.R")
###########input###############
#set input
dist_array=readRDS(paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_array_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))

#input phenotype
if(is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=meta=read.table("../Data_PRJNA434002/meta.tsv",header = TRUE, sep = "\t")
}
if(!is.na(unlist(strsplit(file_tag,"k"))[2])){
  tmeta=readRDS(paste0("../Data_PRJNA434002/meta",unlist(strsplit(file_tag,"k"))[2],".rds"))
}
total_individual=unique(tmeta$individual)
cur_cluster=as.character(unique(tmeta$cluster)[cluster_tag])
meta=tmeta[tmeta$cluster==cur_cluster,]
cur_individual=unique(meta$individual)

if(perm_method=="s"){
  cell_num=table(meta$individual)
  cell_index=which(cell_num>=cell_threshold)
  cur_individual=cur_individual[cell_index]
  dist_array=dist_array[,cell_index,cell_index,drop=FALSE]
}


phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0



###################calculation t#################################
print("start calculation: Part I: Empirical KLmean and JSD")

pval_perm=matrix(ncol=(length(perm_label_seq)-1),nrow=dim(dist_array)[1])
rownames(pval_perm)=dimnames(dist_array)[[1]]
anova_centroid_pval_perm=pval_perm


for(perm_label in perm_label_seq){
  
  #set perm
  if(perm_label>0){
    phenotype=matrix(1,ncol=1,nrow=length(cur_individual))
    phenotype[which(meta$diagnosis[match(cur_individual,meta$individual)]=="Control")]=0
    
    perm_order=readRDS(paste0("../Data_PRJNA434002/7.Result/ind_perm_order.rds"))
    perm_order=as.numeric(perm_order[,perm_label])
    total_individual_ref=total_individual[perm_order]
    perm_order=match(total_individual_ref,cur_individual)
    perm_order=perm_order[!is.na(perm_order)]
    phenotype=phenotype[perm_order,drop=FALSE]
  }
  
  #set covariate
  if(is.na(ind_covariate_flag)){
    covariate_model_matrix=NA
  }
  if(ind_covariate_flag=="ind"){
    #
    cur_covariate=meta[match(cur_individual,meta$individual),c("age","sex","Capbatch","Seqbatch")]
    #read depth
    if(length(grep("rdpadj",file_tag))==0){ #if not justed before
      read_depth=readRDS(paste0("../Data_PRJNA434002/rawM",file_tag,"_read_depth_per_1Kcell_ind.rds"))
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
  
  #plot 
  png(paste0("../Data_PRJNA434002/8.Result/fig_Fstat/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".png"),height = 900,width = 1800,type="cairo")
  op=par(mfrow = c(2, 4))
  pval_order=order(dist_pval)[c(1:4,100,200,300,400)]
  for(iplot in pval_order){
    hist(F_perm0[iplot,],xlim=range(min(F_perm0[iplot,],F_ob0[iplot]),max(F_perm0[iplot,],F_ob0[iplot])),main=paste0("Fstat distr of gene ",dimnames(dist_array)[[1]][iplot],", pval ",round(dist_pval[iplot],3),", rank ",iplot),breaks=50)
    abline(v=F_ob0[iplot],col="red")
  }
  par(op)
  dev.off()
  
  png(paste0("../Data_PRJNA434002/8.Result/fig_Fstat/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,"_Fob.png"),height = 600,width = 1200,type="cairo")
  hist(F_ob0,main=paste0("Fstat distr of gene, Fob"),breaks=50)
  dev.off()
  
  print("0 level complete")
  print(Sys.time())
  print(gc())
  thres=tol/perm_num
  second_index=which(dist_pval<thres)
  if(length(second_index)>0){
    print("1st level")
    print(Sys.time())
    print(gc())
    sub_dist_array=dist_array[second_index,,,drop=FALSE]
    sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*10,Fstat_method=F_method,zm=covariate_model_matrix)$pval
    dist_pval[second_index]=sub_dist_pval
    thres=tol/(perm_num*10)
    second_index=which(dist_pval<thres)
    if(length(second_index)>0){
      print("2nd level")
      print(Sys.time())
      print(gc())
      sub_dist_array=dist_array[second_index,,,drop=FALSE]
      sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*100,Fstat_method=F_method,zm=covariate_model_matrix)$pval
      dist_pval[second_index]=sub_dist_pval
      thres=tol/(perm_num*100)
      second_index=which(dist_pval<thres)
      if(length(second_index)>0){
        print("3rd level")
        print(Sys.time())
        print(gc())
        sub_dist_array=dist_array[second_index,,,drop=FALSE]
        sub_dist_pval=cal_permanova_pval2(sub_dist_array,phenotype,perm_num.min = perm_num*1000,Fstat_method=F_method,zm=covariate_model_matrix)$pval
        dist_pval[second_index]=sub_dist_pval
      }
    }
  }
  
  
  # if(!is.na(ind_covariate_flag)){
  #   saveRDS(dist_pval,paste0("../Data_PRJNA434002/8.Result/"dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  #   
  # }
  # if(is.na(ind_covariate_flag)){
  #   saveRDS(dist_pval,paste0("../Data_PRJNA434002/8.Result/"dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  # }
  
  
  ####Second Chance, to see if we can get more from our method (OPTIONAL)##############################
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
    anova_centroid_pval_ob=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )}) 
    
  }
  if(perm_label>0){
    pval_perm[,perm_label]=dist_pval
    anova_centroid_pval_perm[,perm_label]=apply(dist_array,1,function(x){tryCatch(cal_betadisper_pval(x,label = cur_phenotype_ind,disper_type = "centroid",div_method="anova"), error = function(e) {NA} )}) 
    print(perm_label)
  }
  
}

if(!is.na(ind_covariate_flag)){
  saveRDS(pval_ob,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(pval_perm,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_perm_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  
  saveRDS(anova_centroid_pval_ob,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_anova_centroid_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(anova_centroid_pval_perm,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_perm_anova_centroid_pval_",ind_covariate_flag,"_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}
if(is.na(ind_covariate_flag)){
  saveRDS(pval_ob,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(pval_perm,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_perm_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  
  saveRDS(anova_centroid_pval_ob,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_anova_centroid_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
  saveRDS(anova_centroid_pval_perm,paste0("../Data_PRJNA434002/8.Result/",dist_method,"_",fit_method,"_pval/p",perm_label,perm_method,"_",dist_method,"_",fit_method,"_",F_method,"_perm_anova_centroid_pval_",pre_tag,"_sim_",cluster_tag,"_",file_tag,".rds"))
}

sessionInfo()
q(save="no")
  
