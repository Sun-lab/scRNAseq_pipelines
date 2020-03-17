#this code convert the simulated and dca translated results to sim_fit and sim_data

#input: dca parameters
#output sim_fit and sim_data

# file_tag=4
# r_mean=1.1
# r_var=1.1
# r_change_prop=0.1
# dp_minor_prop=0.3

# setwd("~/Desktop/fh/1.Testing_scRNAseq/")
# setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#perm_label=1
#perm_method="b" #c("","b","bs","s") #b means balanced permutation #s means remove samples with limited cells.
cell_thres=5 #

perm_num=500
covariate_flag=NA #c(NA, "quantile99")
tol=10
perm_label_seq=0:10
sim_folder="sim_v6"
fit_tag="nb" # ""(zinb) or "nb"
sim_n=3

##############functions#################
library("ggplot2")
library("emdbook")
library("abind")
source("./Command/7.0_ZINB_fit_functions.R")
source("./Command/8.0_kl_divergence_functions.R")
source("./Command/9.0_Fstat_functions.R")


######construct param####################
meta=readRDS(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_meta_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

sim_mean=read.table(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/dca_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"/mean.tsv"),stringsAsFactors = FALSE)
sim_dispersion=read.table(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/dca_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"/dispersion.tsv"),stringsAsFactors = FALSE,row.names = 1)
sim_dropout=read.table(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/dca_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,"/dropout.tsv"),stringsAsFactors = FALSE,row.names = 1)

sim_param=abind(sim_mean,sim_dispersion,sim_dropout,along=3)
dimnames(sim_param)[[2]]=colnames(sim_mean)
dimnames(sim_param)[[3]]=c("mean","overdisp","dropout")

#####construct reconstructed data #############

cur_individual=unique(meta$individual)

fit_ind=array(dim=c(nrow(sim_mean),length(cur_individual),3),
                      dimnames = list(rownames(sim_mean),cur_individual,c("logmean","dispersion","dropout_rate")))
sim_ind=array(dim=c(nrow(sim_mean),ncol(sim_mean),sim_n),
              dimnames = list(rownames(sim_mean),colnames(sim_mean),1:sim_n))

for(i_g in 1:nrow(sim_mean)){
  cur_sim=matrix(ncol=sim_n,nrow=ncol(sim_mean))
  for(i_s in 1:ncol(sim_mean)){
    if(fit_tag==""){
      cur_sim[i_s,]=emdbook::rzinbinom(sim_n,sim_mean[i_g,i_s], sim_dispersion[i_g,i_s], sim_dropout[i_g,i_s])
    }
    if(fit_tag=="nb"){
      cur_sim[i_s,]=rnbinom(n=sim_n,mu=sim_mean[i_g,i_s], size=sim_dispersion[i_g,i_s])
    }
  }
  sim_ind[i_g,,]=cur_sim
  for(i_ind in 1:length(cur_individual)){
    cur_ind=cur_individual[i_ind]
    
    #fit sim
    cur_sim_ind=as.numeric(cur_sim[meta$individual==cur_ind,])
    
    if(!is.na(covariate_flag)){
      cur_covariate=rep(covariate[meta$individual==cur_ind,],sim_n)
      if(fit_tag==""){
        fit_ind[i_g,i_ind,]=fit_nbzinb(cur_sim_ind,cur_covariate)
      }
      if(fit_tag=="nb"){
        fit_ind[i_g,i_ind,]=fit_nb(cur_sim_ind,cur_covariate)
      }
      
    }
    if(is.na(covariate_flag)){
      if(fit_tag==""){
        fit_ind[i_g,i_ind,]=fit_nbzinb(cur_sim_ind)
      }
      if(fit_tag=="nb"){
        fit_ind[i_g,i_ind,]=fit_nb(cur_sim_ind)
      }
    }

  }
  print(i_g)
}

saveRDS(fit_ind,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/fit_ind_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

saveRDS(sim_ind,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_ind_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))


saveRDS(sim_param,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/dca_data/sim_param_",fit_tag,r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))

sessionInfo()
#q(save="no")

