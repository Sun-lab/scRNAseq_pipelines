#this code follows 10.4.1.*R and use the output from it.


setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/10.Result/sim_v6/")
#sim_method_seq=c("splat.org","zinb.naive") #"splat.mean","splat.var"



perm_label_seq=0
param_tag_seq=1:4
perm_method_seq="" #c("","b")
pre_tag_seq=c("", "dca")
fit_tag_seq=c("", "nb") 
covariate_flag_seq=c("","readdepth")
resid_flag_seq=c("","resid")
dist_method_seq=c("mean","JSD")
fit_method_seq=c("empirical","zinb","direct")
method_tag_seq=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
de_tag_seq=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
file_tag=1 
perm_label=0
perm_method=""
param_tag=1

factor_len=c(5,5,4,5)
cell_seq=c(20,50,100,200,400)
cell_seq=50
ind_seq=20
file_tag=1
#SO
i_file=1
i_cell=1
i_ind=1
i_perm=1

library("RColorBrewer")
#power_plot is used for plotting the power curves
#note: y_matrix max column=29 since max_col_num=29, max pch_num=14

power_plot1=function(y_matrix,cur_main="",cur_xlab="",cur_ylab="Power",cur_sub="",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA){
  if(sum(!is.na(x_seq))==0){
    x_seq=as.numeric(rownames(y_matrix))
  }
  if(sum(!is.na(cur_xlim))==0){
    cur_xlim=c((min(x_seq)-0.1),(max(x_seq)+0.1))
  }
  if(sum(!is.na(cur_ylim))==0){
    cur_ylim=c(0,1)
  }
  if(sum(!is.na(cur_legend))==0){
    cur_legend=colnames(y_matrix)
  }
  if(sum(!is.na(cur_col))==0){
    total_col=c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))
    cur_col=total_col[(1:ncol(y_matrix))%%29+1]
  }
  plot(x_seq,y_matrix[,1],type="l",main=cur_main,sub=cur_sub,xlim=cur_xlim,ylim=cur_ylim,xlab=cur_xlab,ylab=cur_ylab,col=cur_col[1],lwd=3)
  for(i_c in 1:ncol(y_matrix)){
    lines(x_seq,y_matrix[,i_c],col=cur_col[i_c],lwd=3)
    points(x_seq,y_matrix[,i_c],col=cur_col[i_c],pch=(i_c%%14+1),cex=2)
  }
  
  legend("bottomright",cur_legend,pch=((1:ncol(y_matrix))%%14+1),col=cur_col,cex=.5)
}



#if(pre_tag==""){fit_tag=""}
#if(pre_tag=="dca"){fit_tag="nb"}



param_tag_seq=c("mean","var","dp","mult")

namelist=matrix(ncol=1,nrow=0)
for(i_pre in 1:length(pre_tag_seq)){
  for(i_fit in 1:length(fit_tag_seq)){
    for(i_covariate in 1:length(covariate_flag_seq)){
      for(i_resid in 1:length(resid_flag_seq)){
        pre_tag=pre_tag_seq[i_pre]
        fit_tag=fit_tag_seq[i_fit]
        covariate_flag=covariate_flag_seq[i_covariate]
        resid_flag=resid_flag_seq[i_resid]
        namelist=c(namelist,paste0(pre_tag,fit_tag,covariate_flag,resid_flag))
      }
    }
  }
}
namelist[1]="raw"

power_array=list()
power_array2=list()

for(i_param in 1:length(param_tag_seq)){ 
  param_tag=param_tag_seq[i_param]
  #power_array[[param_tag]]=array(dim=c(length(pre_tag_seq),length(fit_tag_seq),length(covariate_flag_seq),length(resid_flag_seq),length(de_tag_seq),factor_len[i_param],length(method_tag_seq)), dimnames=list(pre_tag_seq,fit_tag_seq,covariate_flag_seq,resid_flag_seq,de_tag_seq,1:factor_len[i_param],method_tag_seq))
  
  power_array[[param_tag]]=array(dim=c(length(pre_tag_seq)*length(fit_tag_seq)*length(covariate_flag_seq)*length(resid_flag_seq),length(de_tag_seq),factor_len[i_param],length(method_tag_seq)), dimnames=list(namelist,de_tag_seq,1:factor_len[i_param],method_tag_seq))
  for(i_pre in 1:length(pre_tag_seq)){
    for(i_fit in 1:length(fit_tag_seq)){
      for(i_covariate in 1:length(covariate_flag_seq)){
        for(i_resid in 1:length(resid_flag_seq)){
          pre_tag=pre_tag_seq[i_pre]
          fit_tag=fit_tag_seq[i_fit]
          covariate_flag=covariate_flag_seq[i_covariate]
          resid_flag=resid_flag_seq[i_resid]
          cur_array=readRDS(paste0("./p",perm_label,perm_method,resid_flag,covariate_flag,pre_tag,fit_tag,"_final_range005_array_",param_tag,".rds"))
          for(i_de in 1:length(de_tag_seq)){
            cur_sub_array=cur_array[i_file,,,,,i_ind,i_cell,,i_de]
            #power_array[[param_tag]][i_pre,i_fit,i_covariate,i_resid,i_de,,]=cur_sub_array
            #dimnames(power_array[[param_tag]])[[6]]=rownames(cur_sub_array)
            power_array[[param_tag]][((i_pre-1)*length(fit_tag_seq)*
                                       length(covariate_flag_seq)*length(resid_flag_seq)+
                                     (i_fit-1)*length(covariate_flag_seq)*length(resid_flag_seq)+
                                     (i_covariate-1)*length(resid_flag_seq)+
                                     i_resid),
                                     i_de,
                                     ,]=cur_sub_array
            dimnames(power_array[[param_tag]])[[3]]=rownames(cur_sub_array)
          }
        }
      }
    }
  }
}



# for(i_param in 1:length(param_tag_seq)){ 
#   param_tag=param_tag_seq[i_param]
#   power_array2[[param_tag]]=apply(power_array[[param_tag]],c(5,6,7),c)
#   dimnames(power_array2[[param_tag]])[[1]]=namelist
# } 




###ploting method comparison
namelist
##a. nb vs zinb resid
model_seq=c(2,6,10,14) #resid
pdf(paste0("./fig_power_curve_method/p",perm_label,"_nb_vs_zinb_resid_",file_tag,"_",ind_seq,"_",cell_seq,"_",param_tag,".pdf"),height = 12,width = 36)

##b. nb vs zinb raw
#model_seq=c(1,5,9,13) #raw
#pdf(paste0("./fig_power_curve_method/p",perm_label,"_nb_vs_zinb_raw_",file_tag,"_",ind_seq,"_",cell_seq,"_",param_tag,".pdf"),height = 12,width = 36)

##c. resid method compare dca-zinb
#model_seq=c(9:12) #raw
#pdf(paste0("./fig_power_curve_method/p",perm_label,"_residcompare_dca_",file_tag,"_",ind_seq,"_",cell_seq,"_",param_tag,".pdf"),height = 12,width = 36)

op=par(mfrow = c(2, 5),cex=1.5)
for(i_method in 1:length(method_tag_seq)){
  for(i_param in 1:length(param_tag_seq)){
    method_tag=method_tag_seq[i_method]
    param_tag=param_tag_seq[i_param]
    power_cur=t(power_array[[param_tag]][model_seq,i_param,,i_method])
    power_cur[power_cur==0]=NA
    print(power_cur)
    power_plot1(power_cur,
                cur_main=paste0("power of ",method_tag,", n ",ind_seq,", ncell ",cell_seq),
                cur_xlab=paste0("DE ",param_tag," change rate"),cur_ylab="Power",
                cur_sub=paste0(pre_tag," ",fit_tag," ",covariate_flag," ",resid_flag),
                x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
  }
  i_param=i_param+1
  power_cur=t(power_array[[param_tag]][,i_param,,i_method])
  print(power_cur)
  ##barplot
  barplot(colMeans(power_cur, na.rm = TRUE, dims = 1),main=paste0("power of ",method_tag,", n ",ind_seq,", ncell ",cell_seq),xlab="",las=2,ylab="Type I error",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA,ylim=c(0,0.2),col=c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"))[(1:ncol(power_cur))%%29+1])
  abline(h=0.05,col="red")
}

par(op)
dev.off()


