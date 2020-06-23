setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

nGeneMean=300
nGeneVar=300
nGeneMult=300
nGeneDP=300
nGeneBlank=1800
nGeneTotal=nGeneMean+nGeneVar+nGeneMult+nGeneDP+nGeneBlank
ncase=50   #!!!!!!!!Different from Wei codes, this is inside the body.
nctrl=50   #!!!!!!!!Different from Wei codes, this is inside the body.
ncell=800
nall=ncase+nctrl

sim_folder="sim_v6"
file_tag_seq=1
r_mean_seq=c(1.1,1.2,1.5,2,4)
r_var_seq=c(1.1,1.2,1.5,2,4)
r_mult_seq=5:9/10
r_dp_seq=1:4/10
  
#r_mean_seq=1.2
#r_var_seq=1.2
#r_mult_seq=0.6
#r_dp_seq=0.2

noise_perc=0.01
extreme_fold_change=10 #the particular cell and genes who have extreme values than expected.

for(r_mean in r_mean_seq){
  for(r_var in r_var_seq){
    for(r_dp in r_dp_seq){
      for(r_mult in r_mult_seq){
        for(file_tag in file_tag_seq){
          if(file.exists(paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))){
            sim_matrix_extreme_flag=rbinom(nGeneTotal*nall*ncell,1,0.05)*(extreme_fold_change-1)+1
            dim(sim_matrix_extreme_flag) = c(nGeneTotal, nall * ncell)
            
            saveRDS(sim_matrix_extreme_flag,paste0("../Data_PRJNA434002/10.Result/",sim_folder,"/sim_data/sim_matrix_extreme_flag_",noise_perc,"_",extreme_fold_change,"_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds")) 
            print(paste(r_mean,r_var,r_dp,r_mult))
          }
        }
      }
    }
  }
}


