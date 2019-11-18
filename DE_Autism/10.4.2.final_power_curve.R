#this code follows 10.4.1.*R and use the output from it.


setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")

#sim_method_seq=c("splat.org","zinb.naive") #"splat.mean","splat.var"

sim_method_seq="zinb.naive"
file_tag_seq=1:5


power_array1=readRDS(paste0("../Data_PRJNA434002/10.Result/final_power_array_param1.rds"))
power_array2=readRDS(paste0("../Data_PRJNA434002/10.Result/final_power_array_param2.rds"))

power_array1[power_array1==0]=NA
power_array2[power_array2==0]=NA

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


pdf(paste0("../Data_PRJNA434002/10.Result/power_curve_final.pdf"),height = 8,width = 8)
op=par(mfrow = c(2, 2), pty = "s")
for(i_sim in length(sim_method_seq):1){
  for(i_file in 1:length(file_tag_seq)){
    file_tag=file_tag_seq[i_file]
    sim_method=sim_method_seq[i_sim]
    
    power_cur1=power_array1[i_file,i_sim,,,,]
    power_cur2=power_array2[i_file,i_sim,,,,]
    
    
    power_plot1(power_cur1[,,1],cur_main=paste0("power of ",dimnames(power_cur1)[[3]][1],", for ",sim_method,", rep ",file_tag),cur_xlab="DE mean change rate",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
    power_plot1(power_cur1[,,3],cur_main=paste0("power of ",dimnames(power_cur1)[[3]][3],", for ",sim_method,", rep ",file_tag),cur_xlab="DE mean change rate",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
    
    power_plot1(power_cur2[,,2],cur_main=paste0("power of ",dimnames(power_cur2)[[3]][2],", for ",sim_method,", rep ",file_tag),cur_xlab="DE var change rate",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
    power_plot1(power_cur2[,,3],cur_main=paste0("power of ",dimnames(power_cur2)[[3]][3],", for ",sim_method,", rep ",file_tag),cur_xlab="DE mean change rate",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)

  }
}
par(op)
dev.off()