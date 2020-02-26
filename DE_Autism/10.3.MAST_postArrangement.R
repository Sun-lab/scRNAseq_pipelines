#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

perm_num=500
sim_folder="sim_v1"
sim_method_seq=c("splat.org","zinb.naive") #"splat.mean","splat.var"
file_tag_seq=1:5
perm_tag_seq=0:500
#round 1:
r_mean_seq=c(1.2,1.5,2,4,6)
r_var_seq=2

#round 2:
r_mean_seq=1.5
r_var_seq=c(1.2,1.5,2,4,6)


for(r_mean in r_mean_seq){
  for(r_var in r_var_seq){
    for(sim_method in sim_method_seq){
      for(file_tag in file_tag_seq){
        
        MAST_pval_ob=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_0.rds"))
        
        MAST_pval_perm=matrix(ncol=perm_num,nrow=length(MAST_pval_ob))
        for(i in 1:perm_num){
          if(file.exists(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",i,".rds"))){
            MAST_pval_perm[,i]=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",i,".rds"))
          }
        }
        
        MAST_pval_perm_flag=apply(MAST_pval_perm,2,function(x){return(x-MAST_pval_ob)})
        
        MAST_pval=apply(MAST_pval_perm_flag>=0,1,function(x){return(sum(x,na.rm = TRUE))})
        MAST_pval_length=apply(MAST_pval_perm_flag,1,function(x){return(sum(!is.na(x)))})
        MAST_pval=MAST_pval/MAST_pval_length
        
        saveRDS(MAST_pval,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
        saveRDS(MAST_pval_ob,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_pval_ob_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
        saveRDS(MAST_pval_perm,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_pval_perm_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,".rds"))
        print(paste0(sim_method,"_",r_mean,"_",r_var,"_",file_tag))
      }
    }
  }
}

print("strat to remove files")

#file remove

for(r_mean in r_mean_seq){
  for(r_var in r_var_seq){
    for(sim_method in sim_method_seq){
      for(file_tag in file_tag_seq){
        for(i in 0:perm_num){
          if(file.exists(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",i,".rds"))){
            file.remove(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/MAST_org_pval_",sim_method,"_",r_mean,"_",r_var,"_",file_tag,"_",i,".rds"))
          }
        }
      }
    }
  }
}




