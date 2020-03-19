
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/10.Result/sim_v6/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

perm_label=0
perm_method=""
param_tag=1 #c(1,2,3,4)

perm_label_seq=0
param_tag_seq=1:4
perm_method_seq="" #c("","b")


dist_method_seq=c("mean","JSD")
fit_method_seq=c("empirical","zinb","direct")


for(pre_tag in c("","dca")){
  if(pre_tag==""){
    fit_tag=""
  }
  if(pre_tag=="dca"){
    fit_tag="nb"
  }
  
  for(perm_method in perm_method_seq){
    for(perm_label in perm_label_seq){
      for(param_tag in param_tag_seq){ 
        file_tag_seq=1
        
        r_mean_seq=1.2
        r_var_seq=1.2
        r_dp_seq=0.2
        r_mult_seq=0.2
        
        ind_seq=c(10,20,40,60)
        cell_seq=c(20,50,100)
        
        if(param_tag==1){
          r_mean_seq=c(1.1,1.2,1.5,2,4)
          param_tag="mean"
        }
        if(param_tag==2){
          r_var_seq=c(1.1,1.2,1.5,2,4)
          param_tag="var"
        }
        if(param_tag==3){
          r_mult_seq=1:5/10
          param_tag="mult"
        }
        if(param_tag==4){
          r_dp_seq=1:4/10
          param_tag="dp"
        }

        #cell_seq=100
        #ind_seq=40
        
        # #test
        # perm_num=500
        # r_mean=1.5  #r_mean/r_var should < 1+mean.shape
        # r_var=1.5
        # file_tag=1
        
        
        
        ###################functions###################
        
        #power calculation
        cal_power=function(x,threshold){
          return(sum(x<=threshold,na.rm = TRUE)/length(x))
        }
        
        cal_range=function(x,threshold1=0,threshold2=1){
          return(sum(x<=threshold2 & x>=threshold1,na.rm = TRUE)/sum(x>-1,na.rm=TRUE))
        }
        ###############################################
        power_array=array(dim=c(
          length(file_tag_seq),
          length(r_mean_seq),
          length(r_var_seq),
          length(r_mult_seq),
          length(r_dp_seq),
          length(ind_seq),
          length(cell_seq),
          8,5),
          dimnames = list(
            file_tag_seq,
            r_mean_seq,
            r_var_seq,
            r_mult_seq,
            r_dp_seq,
            ind_seq,
            cell_seq,
            c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct"),
            c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")))
        
        range01_array=power_array
        range09_array=power_array
        range46_array=power_array
        
        count=1
        max_count=length(file_tag_seq)*length(r_mean_seq)*length(r_var_seq)*length(r_dp_seq)*length(r_mult_seq)*length(cell_seq)*length(ind_seq)
        zeros=matrix(ncol=8,nrow=max_count)
        rownames_zeros=matrix(ncol=1,nrow=max_count)
        colnames(zeros)=c("jsd_zinb_pval","jsd_empirical_pval","jsd_direct_pval","klmean_zinb_pval","klmean_empirical_pval","klmean_direct_pval","MAST_pval","deseq2_pval")
        for(i_file in 1:length(file_tag_seq)){
          for(i_mean in 1:length(r_mean_seq)){
            for(i_var in 1:length(r_var_seq)){
              for(i_dp in 1:length(r_dp_seq)){
                for(i_mult in 1:length(r_mult_seq)){
                  for(i_cell in 1:length(cell_seq)){
                    for(i_ind in 1:length(ind_seq)){
                      file_tag=file_tag_seq[i_file]
                      r_mean=r_mean_seq[i_mean]
                      r_var=r_var_seq[i_var]
                      r_mult=r_mult_seq[i_mult]
                      r_dp=r_dp_seq[i_dp]
                      n_ind=ind_seq[i_ind]
                      n_cell=cell_seq[i_cell]
                      
                      # file_tag=1
                      # r_mean=1.1
                      # r_var=1.2
                      # r_mult=0.2
                      # r_dp=0.3
                      # n_ind=20
                      # n_cell=50
                      
                      fit_tag=""
                      pre_tag=""
                      
                      mean_index=readRDS(paste0("./de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      var_index=readRDS(paste0("./de_label/sim_de.var_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      dp_index=readRDS(paste0("./de_label/sim_de.dp_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      mult_index=readRDS(paste0("./de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      
                      #tryCatch({     }, error = function(e) {NA} )
                      
                      jsd_zinb_pval=NA
                      jsd_empirical_pval=NA
                      jsd_direct_pval=NA
                      klmean_zinb_pval=NA
                      klmean_empirical_pval=NA
                      klmean_direct_pval=NA
                      deseq2_pval=NA
                      MAST_pval=NA
                      
                      if(perm_label>0){
                        tryCatch({jsd_zinb_pval=readRDS(paste0("./JSD_zinb_pval/p10",perm_method,"_JSD_zinb_perm_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
                        tryCatch({jsd_empirical_pval=readRDS(paste0("./JSD_empirical_pval/p10",perm_method,"_JSD_empirical_perm_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
                        tryCatch({jsd_direct_pval=readRDS(paste0("./JSD_direct_pval/p10",perm_method,"_JSD_direct_perm_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
                        tryCatch({klmean_zinb_pval=readRDS(paste0("./mean_zinb_pval/p10",perm_method,"_mean_zinb_perm_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
                        tryCatch({klmean_empirical_pval=readRDS(paste0("./mean_empirical_pval/p10",perm_method,"_mean_empirical_perm_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
                        tryCatch({klmean_direct_pval=readRDS(paste0("./mean_direct_pval/p10",perm_method,"_mean_direct_perm_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))[,perm_label]}, error = function(e) {NA} )
                        tryCatch({deseq2_pval=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,".rds"))[,perm_label]}, error = function(e) {NA} )
                        
                        #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
                        tryCatch({MAST_pval=readRDS(paste0("./MAST_pval/p",perm_label,perm_method,"_MAST_pval1_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                      }
                      if(perm_label==0){
                        tryCatch({jsd_zinb_pval=readRDS(paste0("./JSD_zinb_pval/p",perm_label,perm_method,"_JSD_zinb_raw_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                        tryCatch({jsd_empirical_pval=readRDS(paste0("./JSD_empirical_pval/p",perm_label,perm_method,"_JSD_empirical_raw_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                        tryCatch({jsd_direct_pval=readRDS(paste0("./JSD_direct_pval/p",perm_label,perm_method,"_JSD_direct_raw_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                        tryCatch({klmean_zinb_pval=readRDS(paste0("./mean_zinb_pval/p",perm_label,perm_method,"_mean_zinb_raw_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                        tryCatch({klmean_empirical_pval=readRDS(paste0("./mean_empirical_pval/p",perm_label,perm_method,"_mean_empirical_raw_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                        tryCatch({klmean_direct_pval=readRDS(paste0("./mean_direct_pval/p",perm_label,perm_method,"_mean_direct_raw_pval_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                        tryCatch({deseq2_pval=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,".rds"))}, error = function(e) {NA} )
                        
                        #note! please be sure to use 10.3.MAST_postArrangment.R when all permutation results are ready.
                        tryCatch({MAST_pval=readRDS(paste0("./MAST_pval/p0_MAST_pval1_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".rds"))}, error = function(e) {NA} )
                      }
                      
                      zeros[count,]=c(sum(is.na(jsd_zinb_pval)),sum(is.na(jsd_empirical_pval)),sum(is.na(jsd_direct_pval)),sum(is.na(klmean_zinb_pval)),sum(is.na(klmean_empirical_pval)),sum(is.na(klmean_direct_pval)),sum(is.na(MAST_pval)),sum(is.na(deseq2_pval)))
                      rownames_zeros[count]=paste0(r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell)
                      count=count+1  
                      
                      #histogram
                      if((!is.na(jsd_direct_pval)) || (!is.na(jsd_empirical_pval)) || (!is.na(jsd_zinb_pval)) ){
                        png(paste0("./fig_pval_hist/p",perm_label,perm_method,"_pval_hist_",param_tag,"_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 4800,width = 3000)
                        op=par(mfrow = c(8, 5))
                        
                        tryCatch({hist(klmean_zinb_pval[mean_index==1],main="pval of mean-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_zinb_pval[var_index==1],main="pval of var-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
                        tryCatch({hist(klmean_zinb_pval[dp_index==1],main="pval of dp-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_zinb_pval[mult_index==1],main="pval of mult-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,klmean_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(jsd_zinb_pval[mean_index==1],main="pval of mean-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_zinb_pval[var_index==1],main="pval of var-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
                        tryCatch({hist(jsd_zinb_pval[dp_index==1],main="pval of dp-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_zinb_pval[mult_index==1],main="pval of mult-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,jsd_zinb method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(klmean_empirical_pval[mean_index==1],main="pval of mean-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_empirical_pval[var_index==1],main="pval of var-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
                        tryCatch({hist(klmean_empirical_pval[dp_index==1],main="pval of dp-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_empirical_pval[mult_index==1],main="pval of mult-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,klmean_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(jsd_empirical_pval[mean_index==1],main="pval of mean-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_empirical_pval[var_index==1],main="pval of var-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
                        tryCatch({hist(jsd_empirical_pval[dp_index==1],main="pval of dp-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_empirical_pval[mult_index==1],main="pval of mult-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,jsd_empirical method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(jsd_direct_pval[mean_index==1],main="pval of mean-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_direct_pval[var_index==1],main="pval of var-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
                        tryCatch({hist(jsd_direct_pval[dp_index==1],main="pval of dp-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_direct_pval[mult_index==1],main="pval of mult-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,jsd_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(klmean_direct_pval[mean_index==1],main="pval of mean-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_direct_pval[var_index==1],main="pval of var-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )        
                        tryCatch({hist(klmean_direct_pval[dp_index==1],main="pval of dp-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_direct_pval[mult_index==1],main="pval of mult-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,klmean_direct method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(deseq2_pval[mean_index==1],main="pval of mean-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval[var_index==1],main="pval of var-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )  
                        tryCatch({hist(deseq2_pval[dp_index==1],main="pval of dp-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval[mult_index==1],main="pval of mult-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(MAST_pval[mean_index==1],main="pval of mean-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(MAST_pval[var_index==1],main="pval of var-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )  
                        tryCatch({hist(MAST_pval[dp_index==1],main="pval of dp-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(MAST_pval[mult_index==1],main="pval of mult-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],main="pval of non-DE genes,MAST method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        par(op)
                        dev.off()
                      }
                      
                      
                      
                      power_matrix=matrix(nrow=8,ncol=5)
                      
                      rownames(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
                      colnames(power_matrix)=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
                      
                      tryCatch({power_matrix[1,]=c(cal_range(deseq2_pval[mean_index==1],0,0.1),
                                                   cal_range(deseq2_pval[var_index==1],0,0.1),
                                                   cal_range(deseq2_pval[dp_index==1],0,0.1),
                                                   cal_range(deseq2_pval[mult_index==1],0,0.1),
                                                   cal_range(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[2,]=c(cal_range(MAST_pval[mean_index==1],0,0.1),
                                                   cal_range(MAST_pval[var_index==1],0,0.1),
                                                   cal_range(MAST_pval[dp_index==1],0,0.1),
                                                   cal_range(MAST_pval[mult_index==1],0,0.1),
                                                   cal_range(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[3,]=c(cal_range(jsd_empirical_pval[mean_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[var_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[dp_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[mult_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[4,]=c(cal_range(klmean_empirical_pval[mean_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[var_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[dp_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[mult_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[5,]=c(cal_range(jsd_zinb_pval[mean_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[var_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[dp_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[mult_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[6,]=c(cal_range(klmean_zinb_pval[mean_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[var_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[dp_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[mult_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[7,]=c(cal_range(jsd_direct_pval[mean_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[var_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[dp_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[mult_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[8,]=c(cal_range(klmean_direct_pval[mean_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[var_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[dp_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[mult_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      power_matrix
                      
                      # #barplot
                      # png(paste0("./fig_barplot/p",perm_label,perm_method,"_barplot_",param_tag,"_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 2400,width = 1200)
                      # op=par(mfrow = c(4, 2), pty = "s")
                      # barplot(power_matrix[1,],ylab="power",main=rownames(power_matrix)[1],ylim=c(0,1))
                      # barplot(power_matrix[2,],ylab="power",main=rownames(power_matrix)[2],ylim=c(0,1))
                      # barplot(power_matrix[3,],ylab="power",main=rownames(power_matrix)[3],ylim=c(0,1))
                      # barplot(power_matrix[4,],ylab="power",main=rownames(power_matrix)[4],ylim=c(0,1))
                      # barplot(power_matrix[5,],ylab="power",main=rownames(power_matrix)[5],ylim=c(0,1))
                      # barplot(power_matrix[6,],ylab="power",main=rownames(power_matrix)[6],ylim=c(0,1))
                      # barplot(power_matrix[7,],ylab="power",main=rownames(power_matrix)[7],ylim=c(0,1))
                      # barplot(power_matrix[8,],ylab="power",main=rownames(power_matrix)[8],ylim=c(0,1))
                      # par(op)
                      # dev.off()
                      if((!is.na(jsd_direct_pval)) || (!is.na(jsd_empirical_pval)) || (!is.na(jsd_zinb_pval)) ){
                        png(paste0("./fig_power_point/p",perm_label,perm_method,"_power_point_",param_tag,"_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 800,width = 800)
                        plot(power_matrix[1,5],power_matrix[1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
                        points(power_matrix[1,5],power_matrix[1,2],col="red",pch=4,cex=3)
                        points(power_matrix[1,5],power_matrix[1,3],col="red",pch=5,cex=3)
                        points(power_matrix[1,5],power_matrix[1,4],col="red",pch=6,cex=3)
                        
                        points(power_matrix[2,5],power_matrix[2,1],col="blue",pch=3,cex=3)
                        points(power_matrix[2,5],power_matrix[2,2],col="blue",pch=4,cex=3)
                        points(power_matrix[2,5],power_matrix[2,3],col="blue",pch=5,cex=3)
                        points(power_matrix[2,5],power_matrix[2,4],col="blue",pch=6,cex=3)
                        
                        points(power_matrix[3,5],power_matrix[3,1],col="pink",pch=3,cex=3)
                        points(power_matrix[3,5],power_matrix[3,2],col="pink",pch=4,cex=3)
                        points(power_matrix[3,5],power_matrix[3,3],col="pink",pch=5,cex=3)
                        points(power_matrix[3,5],power_matrix[3,4],col="pink",pch=6,cex=3)
                        
                        points(power_matrix[4,5],power_matrix[4,1],col="brown",pch=3,cex=3)
                        points(power_matrix[4,5],power_matrix[4,2],col="brown",pch=4,cex=3)
                        points(power_matrix[4,5],power_matrix[4,3],col="brown",pch=5,cex=3)
                        points(power_matrix[4,5],power_matrix[4,4],col="brown",pch=6,cex=3)
                        
                        points(power_matrix[5,5],power_matrix[5,1],col="orange",pch=3,cex=3)
                        points(power_matrix[5,5],power_matrix[5,2],col="orange",pch=4,cex=3)
                        points(power_matrix[5,5],power_matrix[5,3],col="orange",pch=5,cex=3)
                        points(power_matrix[5,5],power_matrix[5,4],col="orange",pch=6,cex=3)
                        
                        points(power_matrix[6,5],power_matrix[6,1],col="green",pch=3,cex=3)
                        points(power_matrix[6,5],power_matrix[6,2],col="green",pch=4,cex=3)
                        points(power_matrix[6,5],power_matrix[6,3],col="green",pch=5,cex=3)
                        points(power_matrix[6,5],power_matrix[6,4],col="green",pch=6,cex=3)
                        
                        points(power_matrix[7,5],power_matrix[7,1],col="goldenrod",pch=3,cex=3)
                        points(power_matrix[7,5],power_matrix[7,2],col="goldenrod",pch=4,cex=3)
                        points(power_matrix[7,5],power_matrix[7,3],col="goldenrod",pch=5,cex=3)
                        points(power_matrix[7,5],power_matrix[7,4],col="goldenrod",pch=6,cex=3)
                        
                        points(power_matrix[8,5],power_matrix[8,1],col="deeppink",pch=3,cex=3)
                        points(power_matrix[8,5],power_matrix[8,2],col="deeppink",pch=4,cex=3)
                        points(power_matrix[8,5],power_matrix[8,3],col="deeppink",pch=5,cex=3)
                        points(power_matrix[8,5],power_matrix[8,4],col="deeppink",pch=6,cex=3)
                        
                        legend("topright",c(rownames(power_matrix),"mean diff","var diff","dp_diff","mult_diff"),pch=c(rep(15,8),3:8),cex=1,col=c("red","blue","pink","brown","orange","green","goldenrod","deeppink","black","black","black","black"))
                        
                        dev.off()
                      }
                      
                      
                      power_array[i_file,i_mean,i_var,i_mult,i_dp,i_ind,i_cell,,]=power_matrix
                      
                      
                      
                      
                      
                      power_matrix=matrix(nrow=8,ncol=5)
                      
                      rownames(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
                      colnames(power_matrix)=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
                      
                      tryCatch({power_matrix[1,]=c(cal_range(deseq2_pval[mean_index==1],0,0.1),
                                                   cal_range(deseq2_pval[var_index==1],0,0.1),
                                                   cal_range(deseq2_pval[dp_index==1],0,0.1),
                                                   cal_range(deseq2_pval[mult_index==1],0,0.1),
                                                   cal_range(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[2,]=c(cal_range(MAST_pval[mean_index==1],0,0.1),
                                                   cal_range(MAST_pval[var_index==1],0,0.1),
                                                   cal_range(MAST_pval[dp_index==1],0,0.1),
                                                   cal_range(MAST_pval[mult_index==1],0,0.1),
                                                   cal_range(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[3,]=c(cal_range(jsd_empirical_pval[mean_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[var_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[dp_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[mult_index==1],0,0.1),
                                                   cal_range(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[4,]=c(cal_range(klmean_empirical_pval[mean_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[var_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[dp_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[mult_index==1],0,0.1),
                                                   cal_range(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[5,]=c(cal_range(jsd_zinb_pval[mean_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[var_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[dp_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[mult_index==1],0,0.1),
                                                   cal_range(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[6,]=c(cal_range(klmean_zinb_pval[mean_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[var_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[dp_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[mult_index==1],0,0.1),
                                                   cal_range(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[7,]=c(cal_range(jsd_direct_pval[mean_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[var_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[dp_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[mult_index==1],0,0.1),
                                                   cal_range(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[8,]=c(cal_range(klmean_direct_pval[mean_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[var_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[dp_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[mult_index==1],0,0.1),
                                                   cal_range(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0,0.1))}, error = function(e) {NA} )
                      power_matrix
                      range01_array[i_file,i_mean,i_var,i_mult,i_dp,i_ind,i_cell,,]=power_matrix           
                      
                      
                      
                      power_matrix=matrix(nrow=8,ncol=5)
                      
                      rownames(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
                      colnames(power_matrix)=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
                      
                      tryCatch({power_matrix[1,]=c(cal_range(deseq2_pval[mean_index==1],0.9,1),
                                                   cal_range(deseq2_pval[var_index==1],0.9,1),
                                                   cal_range(deseq2_pval[dp_index==1],0.9,1),
                                                   cal_range(deseq2_pval[mult_index==1],0.9,1),
                                                   cal_range(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[2,]=c(cal_range(MAST_pval[mean_index==1],0.9,1),
                                                   cal_range(MAST_pval[var_index==1],0.9,1),
                                                   cal_range(MAST_pval[dp_index==1],0.9,1),
                                                   cal_range(MAST_pval[mult_index==1],0.9,1),
                                                   cal_range(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[3,]=c(cal_range(jsd_empirical_pval[mean_index==1],0.9,1),
                                                   cal_range(jsd_empirical_pval[var_index==1],0.9,1),
                                                   cal_range(jsd_empirical_pval[dp_index==1],0.9,1),
                                                   cal_range(jsd_empirical_pval[mult_index==1],0.9,1),
                                                   cal_range(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[4,]=c(cal_range(klmean_empirical_pval[mean_index==1],0.9,1),
                                                   cal_range(klmean_empirical_pval[var_index==1],0.9,1),
                                                   cal_range(klmean_empirical_pval[dp_index==1],0.9,1),
                                                   cal_range(klmean_empirical_pval[mult_index==1],0.9,1),
                                                   cal_range(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[5,]=c(cal_range(jsd_zinb_pval[mean_index==1],0.9,1),
                                                   cal_range(jsd_zinb_pval[var_index==1],0.9,1),
                                                   cal_range(jsd_zinb_pval[dp_index==1],0.9,1),
                                                   cal_range(jsd_zinb_pval[mult_index==1],0.9,1),
                                                   cal_range(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[6,]=c(cal_range(klmean_zinb_pval[mean_index==1],0.9,1),
                                                   cal_range(klmean_zinb_pval[var_index==1],0.9,1),
                                                   cal_range(klmean_zinb_pval[dp_index==1],0.9,1),
                                                   cal_range(klmean_zinb_pval[mult_index==1],0.9,1),
                                                   cal_range(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[7,]=c(cal_range(jsd_direct_pval[mean_index==1],0.9,1),
                                                   cal_range(jsd_direct_pval[var_index==1],0.9,1),
                                                   cal_range(jsd_direct_pval[dp_index==1],0.9,1),
                                                   cal_range(jsd_direct_pval[mult_index==1],0.9,1),
                                                   cal_range(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      tryCatch({power_matrix[8,]=c(cal_range(klmean_direct_pval[mean_index==1],0.9,1),
                                                   cal_range(klmean_direct_pval[var_index==1],0.9,1),
                                                   cal_range(klmean_direct_pval[dp_index==1],0.9,1),
                                                   cal_range(klmean_direct_pval[mult_index==1],0.9,1),
                                                   cal_range(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.9,1))}, error = function(e) {NA} )
                      power_matrix
                      range09_array[i_file,i_mean,i_var,i_mult,i_dp,i_ind,i_cell,,]=power_matrix
                      
                      
                      
                      
                      power_matrix=matrix(nrow=8,ncol=5)
                      
                      rownames(power_matrix)=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
                      colnames(power_matrix)=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
                      
                      tryCatch({power_matrix[1,]=c(cal_range(deseq2_pval[mean_index==1],0.4,0.6),
                                                   cal_range(deseq2_pval[var_index==1],0.4,0.6),
                                                   cal_range(deseq2_pval[dp_index==1],0.4,0.6),
                                                   cal_range(deseq2_pval[mult_index==1],0.4,0.6),
                                                   cal_range(deseq2_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[2,]=c(cal_range(MAST_pval[mean_index==1],0.4,0.6),
                                                   cal_range(MAST_pval[var_index==1],0.4,0.6),
                                                   cal_range(MAST_pval[dp_index==1],0.4,0.6),
                                                   cal_range(MAST_pval[mult_index==1],0.4,0.6),
                                                   cal_range(MAST_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[3,]=c(cal_range(jsd_empirical_pval[mean_index==1],0.4,0.6),
                                                   cal_range(jsd_empirical_pval[var_index==1],0.4,0.6),
                                                   cal_range(jsd_empirical_pval[dp_index==1],0.4,0.6),
                                                   cal_range(jsd_empirical_pval[mult_index==1],0.4,0.6),
                                                   cal_range(jsd_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[4,]=c(cal_range(klmean_empirical_pval[mean_index==1],0.4,0.6),
                                                   cal_range(klmean_empirical_pval[var_index==1],0.4,0.6),
                                                   cal_range(klmean_empirical_pval[dp_index==1],0.4,0.6),
                                                   cal_range(klmean_empirical_pval[mult_index==1],0.4,0.6),
                                                   cal_range(klmean_empirical_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[5,]=c(cal_range(jsd_zinb_pval[mean_index==1],0.4,0.6),
                                                   cal_range(jsd_zinb_pval[var_index==1],0.4,0.6),
                                                   cal_range(jsd_zinb_pval[dp_index==1],0.4,0.6),
                                                   cal_range(jsd_zinb_pval[mult_index==1],0.4,0.6),
                                                   cal_range(jsd_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[6,]=c(cal_range(klmean_zinb_pval[mean_index==1],0.4,0.6),
                                                   cal_range(klmean_zinb_pval[var_index==1],0.4,0.6),
                                                   cal_range(klmean_zinb_pval[dp_index==1],0.4,0.6),
                                                   cal_range(klmean_zinb_pval[mult_index==1],0.4,0.6),
                                                   cal_range(klmean_zinb_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[7,]=c(cal_range(jsd_direct_pval[mean_index==1],0.4,0.6),
                                                   cal_range(jsd_direct_pval[var_index==1],0.4,0.6),
                                                   cal_range(jsd_direct_pval[dp_index==1],0.4,0.6),
                                                   cal_range(jsd_direct_pval[mult_index==1],0.4,0.6),
                                                   cal_range(jsd_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      tryCatch({power_matrix[8,]=c(cal_range(klmean_direct_pval[mean_index==1],0.4,0.6),
                                                   cal_range(klmean_direct_pval[var_index==1],0.4,0.6),
                                                   cal_range(klmean_direct_pval[dp_index==1],0.4,0.6),
                                                   cal_range(klmean_direct_pval[mult_index==1],0.4,0.6),
                                                   cal_range(klmean_direct_pval[mean_index==0 & var_index==0 & dp_index==0 & mult_index==0],0.4,0.6))}, error = function(e) {NA} )
                      power_matrix
                      range46_array[i_file,i_mean,i_var,i_mult,i_dp,i_ind,i_cell,,]=power_matrix
                      
                      
                      
                      
                      print(paste0(pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell))
                    }
                  }
                }
              }
            }
          }
        }
        
        saveRDS(power_array,paste0("./p",perm_label,perm_method,pre_tag,fit_tag,"_final_power_array_",param_tag,".rds"))
        saveRDS(range01_array,paste0("./p",perm_label,perm_method,pre_tag,fit_tag,"_final_range01_array_",param_tag,".rds"))
        saveRDS(range09_array,paste0("./p",perm_label,perm_method,pre_tag,fit_tag,"_final_range09_array_",param_tag,".rds"))
        saveRDS(range46_array,paste0("./p",perm_label,perm_method,pre_tag,fit_tag,"_final_range46_array_",param_tag,".rds"))
        
        rownames(zeros)=rownames_zeros
        View(zeros)
        saveRDS(power_array,paste0("./p",perm_label,perm_method,pre_tag,fit_tag,"_pval_NAs_",param_tag,".rds"))
        
        #more plot
        
        for(i_file in 1:length(file_tag_seq)){
          for(i_cell in 1:length(cell_seq)){
            for(i_ind in 1:length(ind_seq)){
              file_tag=file_tag_seq[i_file]
              n_ind=ind_seq[i_ind]
              n_cell=cell_seq[i_cell]
              
              cur_power_array=power_array[i_file,,,,,i_ind,i_cell,,]
              
              if((!is.na(jsd_direct_pval)) || (!is.na(jsd_empirical_pval)) || (!is.na(jsd_zinb_pval)) ){
                png(paste0("./fig_final_power/p",perm_label,perm_method,pre_tag,fit_tag,"_final_power_",param_tag,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 800,width = 800)
                plot(cur_power_array[,1,5],cur_power_array[,1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3,main=paste0("power scatter: ",file_tag))
                points(cur_power_array[,1,5],cur_power_array[,1,2],col="red",pch=4,cex=3)
                points(cur_power_array[,1,5],cur_power_array[,1,3],col="red",pch=5,cex=3)
                points(cur_power_array[,1,5],cur_power_array[,1,4],col="red",pch=6,cex=3)
                
                points(cur_power_array[,2,5],cur_power_array[,2,1],col="blue",pch=3,cex=3)
                points(cur_power_array[,2,5],cur_power_array[,2,2],col="blue",pch=4,cex=3)
                points(cur_power_array[,2,5],cur_power_array[,2,3],col="blue",pch=5,cex=3)
                points(cur_power_array[,2,5],cur_power_array[,2,4],col="blue",pch=6,cex=3)
                
                points(cur_power_array[,3,5],cur_power_array[,3,1],col="pink",pch=3,cex=3)
                points(cur_power_array[,3,5],cur_power_array[,3,2],col="pink",pch=4,cex=3)
                points(cur_power_array[,3,5],cur_power_array[,3,3],col="pink",pch=5,cex=3)
                points(cur_power_array[,3,5],cur_power_array[,3,4],col="pink",pch=6,cex=3)
                
                points(cur_power_array[,4,5],cur_power_array[,4,1],col="brown",pch=3,cex=3)
                points(cur_power_array[,4,5],cur_power_array[,4,2],col="brown",pch=4,cex=3)
                points(cur_power_array[,4,5],cur_power_array[,4,3],col="brown",pch=5,cex=3)
                points(cur_power_array[,4,5],cur_power_array[,4,4],col="brown",pch=6,cex=3)
                
                points(cur_power_array[,5,5],cur_power_array[,5,1],col="orange",pch=3,cex=3)
                points(cur_power_array[,5,5],cur_power_array[,5,2],col="orange",pch=4,cex=3)
                points(cur_power_array[,5,5],cur_power_array[,5,3],col="orange",pch=5,cex=3)
                points(cur_power_array[,5,5],cur_power_array[,5,4],col="orange",pch=6,cex=3)
                
                points(cur_power_array[,6,5],cur_power_array[,6,1],col="green",pch=3,cex=3)
                points(cur_power_array[,6,5],cur_power_array[,6,2],col="green",pch=4,cex=3)
                points(cur_power_array[,6,5],cur_power_array[,6,3],col="green",pch=5,cex=3)
                points(cur_power_array[,6,5],cur_power_array[,6,4],col="green",pch=6,cex=3)
                
                points(cur_power_array[,7,5],cur_power_array[,7,1],col="goldenrod",pch=3,cex=3)
                points(cur_power_array[,7,5],cur_power_array[,7,2],col="goldenrod",pch=4,cex=3)
                points(cur_power_array[,7,5],cur_power_array[,7,3],col="goldenrod",pch=5,cex=3)
                points(cur_power_array[,7,5],cur_power_array[,7,4],col="goldenrod",pch=6,cex=3)
                
                points(cur_power_array[,8,5],cur_power_array[,8,1],col="deeppink",pch=3,cex=3)
                points(cur_power_array[,8,5],cur_power_array[,8,2],col="deeppink",pch=4,cex=3)
                points(cur_power_array[,8,5],cur_power_array[,8,3],col="deeppink",pch=5,cex=3)
                points(cur_power_array[,8,5],cur_power_array[,8,4],col="deeppink",pch=6,cex=3)
                
                legend("topright",c(rownames(power_matrix),"mean diff","var diff","dp_diff","mult_diff"),pch=c(rep(15,8),3:8),cex=1,col=c("red","blue","pink","brown","orange","green","goldenrod","deeppink","black","black","black","black"))
                dev.off()
              }
            }
          }
        }
        
      }
    }
  }
  
}
