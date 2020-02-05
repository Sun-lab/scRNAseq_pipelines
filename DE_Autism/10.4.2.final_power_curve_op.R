#this code follows 10.4.1.*R and use the output from it.


setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")

perm_label=0
perm_method=""

file_tag_seq=1:2

sim_method="zinb.naive"

cell_seq=1:5*20
ind_seq=1:4*10

r_mean_seq=c(1.1,1.2,1.5,2,4)
r_var_seq=c(1.1,1.2,1.5,2,4)
r_disp_seq=c(1.1,1.2,1.5,2,4)
r_mult_seq=c(0.2,0.4,0.6,0.8)

param_tag_seq=c("mean","var","disp","mult")

power_array=list()
for(param_tag in param_tag_seq){
  power_array[[param_tag]]=readRDS(paste0("../Data_PRJNA434002/10.Result/op/p",perm_label,perm_method,"_final_power_array_",param_tag,".rds"))
  power_array[[param_tag]][power_array[[param_tag]]==0]=NA
  
}



library("RColorBrewer")
#power_plot is used for plotting the power curves
#note: y_matrix max column=12

power_plot1=function(y_matrix,cur_main="",cur_xlab="",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA){
  if(sum(!is.na(x_seq))==0){
    x_seq=as.numeric(rownames(y_matrix))
  }
  if(sum(!is.na(cur_xlim))==0){
    cur_xlim=c((min(x_seq)-0.2),(max(x_seq)+0.5))
  }
  if(sum(!is.na(cur_ylim))==0){
    cur_ylim=c(0,1)
  }
  if(sum(!is.na(cur_legend))==0){
    cur_legend=colnames(y_matrix)
  }
  if(sum(!is.na(cur_col))==0){
    cur_col=brewer.pal(ncol(y_matrix),"Paired")
  }
  plot(x_seq,y_matrix[,1],type="l",main=cur_main,xlim=cur_xlim,ylim=cur_ylim,xlab=cur_xlab,ylab=cur_ylab,col=cur_col[1],lwd=3)
  for(i_c in 1:ncol(y_matrix)){
    lines(x_seq,y_matrix[,i_c],col=cur_col[i_c],lwd=3)
    points(x_seq,y_matrix[,i_c],col=cur_col[i_c],pch=i_c,cex=2)
  }
  
  legend("bottomright",cur_legend,pch=1:ncol(y_matrix),col=cur_col,cex=.5)
}


pdf(paste0("../Data_PRJNA434002/10.Result/op/power_curve_final.pdf"),height = 16,width = 8)
op=par(mfrow = c(4, 2), pty = "s")
for(i_file in 1:length(file_tag_seq)){
  for(i_cell in 1:length(cell_seq)){
    for(i_ind in 1:length(ind_seq)){
      file_tag=file_tag_seq[i_file]
      n_ind=ind_seq[i_ind]
      n_cell=cell_seq[i_cell]

      for(i_param in 1:length(param_tag_seq)){
        param_tag=param_tag_seq[i_param]
        power_cur=power_array[[param_tag]][i_file,,,,,i_ind,i_cell,,]
        
        i_cur=i_param
        #for(i_cur in 1:4){
          power_plot1(power_cur[,,i_cur],cur_main=paste0("power of ",dimnames(power_cur)[[3]][i_cur],", for ",sim_method,", rep ",file_tag,", n ",n_ind,", ncell ",n_cell),cur_xlab="DE mean change rate",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
          power_plot1(power_cur[,,5],cur_main=paste0("power of ",dimnames(power_cur)[[3]][5],", for ",sim_method,", rep ",file_tag,", n ",n_ind,", ncell ",n_cell),cur_xlab="DE mean change rate",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
        #}
      }       
    }
  }
}
par(op)
dev.off()
