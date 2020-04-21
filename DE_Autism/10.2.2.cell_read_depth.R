setwd("/fh/fast/sun_w/mengqi/Data_PRJNA434002/10.Result/")

file_tag_seq=1:4
r_mean_seq=c(1.1,1.2,1.5,2,4)
r_var_seq=c(1.1,1.2,1.5,2,4)
r_change_prop_seq=c(1:10)/10
dp_minor_prop_seq=1:4/10

#r_mean_seq=1.2
#r_var_seq=1.2
#r_change_prop_seq=0.6
#dp_minor_prop=0.2
sim_folder="sim_v6"

for(r_mean in r_mean_seq){
  for(r_var in r_var_seq){
    #for(r_disp in r_disp_seq){
    for(r_change_prop in r_change_prop_seq){
      for(file_tag in file_tag_seq){
        for(dp_minor_prop in dp_minor_prop_seq){
          if(file.exists(paste0(sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))){
            sim_matrix=readRDS(paste0(sim_folder,"/sim_data/sim_matrix_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
            cell_readdepth=apply(sim_matrix,2,sum)
            cell_count_quantile=apply(sim_matrix,2,quantile)
            saveRDS(cell_readdepth,paste0(sim_folder,"/sim_data/sim_cell_readdepth_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
            saveRDS(cell_count_quantile,paste0(sim_folder,"/sim_data/sim_cell_count_quantile_",r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag,".rds"))
            print(paste0(r_mean,"_",r_var,"_",r_change_prop,"_",dp_minor_prop,"_",file_tag))
          }
        }
      }
    }
  }
}


