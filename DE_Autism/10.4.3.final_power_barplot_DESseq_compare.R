
#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/10.Result/sim_v6/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

#complete param
perm_label=0
perm_method=""
param_tag=1 #c(1,2,3,4)

perm_label_seq=0
param_tag_seq=1:4
perm_method_seq="" #c("","b")
pre_tag_seq=c("","dca")

dist_method_seq=c("mean","JSD")
fit_method_seq=c("empirical","zinb","direct")

#shrink param
perm_method_seq=""
pre_tag_seq=""

for(pre_tag in pre_tag_seq){
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
        r_mult_seq=0.6
        
        ind_seq=c(10,20,40,60,100)
        cell_seq=c(20,50,100,200)
        
        if(param_tag==1){
          r_mean_seq=c(1.1,1.2,1.5,2,4)
          param_tag="mean"
        }
        if(param_tag==2){
          r_var_seq=c(1.1,1.2,1.5,2,4)
          param_tag="var"
        }
        if(param_tag==3){
          r_mult_seq=1:4/5
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
          2,5),
          dimnames = list(
            file_tag_seq,
            r_mean_seq,
            r_var_seq,
            r_mult_seq,
            r_dp_seq,
            ind_seq,
            cell_seq,
            c("DESeq","DESeq_sep"),
            c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")))
        
        range005_array=power_array
        
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
                      #fit_tag=""
                      #pre_tag=""
                      
                      mean_index=readRDS(paste0("./de_label/sim_de.mean_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      var_index=readRDS(paste0("./de_label/sim_de.var_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      dp_index=readRDS(paste0("./de_label/sim_de.dp_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      mult_index=readRDS(paste0("./de_label/sim_de.mult_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,".rds"))
                      
                      #tryCatch({     }, error = function(e) {NA} )

                      deseq2_pval=NA
                      deseq2_pval_sep=NA
                      if(perm_label==0){
                        tryCatch({deseq2_pval=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,".rds"))}, error = function(e) {NA} )
                        
                        
                        deseq2_pval_sep=matrix(ncol=1,nrow=3000)
                        tryCatch({deseq2_pval_sep[which(mean_index==1)]=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_mean.rds"))}, error = function(e) {NA} )
                        tryCatch({deseq2_pval_sep[which(var_index==1)]=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_var.rds"))}, error = function(e) {NA} )
                        tryCatch({deseq2_pval_sep[which(dp_index==1)]=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_dp.rds"))}, error = function(e) {NA} )
                        tryCatch({deseq2_pval_sep[which(mult_index==1)]=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_mult.rds"))}, error = function(e) {NA} )
                        tryCatch({deseq2_pval_sep[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)]=readRDS(paste0("./DESeq2_pval/p",perm_label,perm_method,"_DESeq2_pval_",r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_none.rds"))}, error = function(e) {NA} )
                      }
                      
                      #histogram
                      if((!is.na(deseq2_pval)) || (!is.na(deseq2_pval_sep))){
                        png(paste0("./fig_DESeq2_compare_pval_hist/p",perm_label,perm_method,"_pval_hist_",param_tag,"_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 1200,width = 3000)
                        op=par(mfrow = c(2, 5),cex=2)
                        tryCatch({hist(deseq2_pval[which(mean_index==1)],main="pval of mean-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval[which(var_index==1)],main="pval of var-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )  
                        tryCatch({hist(deseq2_pval[which(dp_index==1)],main="pval of dp-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval[which(mult_index==1)],main="pval of mult-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)],main="pval of non-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        tryCatch({hist(deseq2_pval_sep[which(mean_index==1)],main="pval_sep of mean-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval_sep[which(var_index==1)],main="pval_sep of var-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )  
                        tryCatch({hist(deseq2_pval_sep[which(dp_index==1)],main="pval_sep of dp-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval_sep[which(mult_index==1)],main="pval_sep of mult-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        tryCatch({hist(deseq2_pval_sep[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)],main="pval_sep of non-DE genes,deseq2 method",xlab="p-values",breaks = 20)}, error = function(e) {NA} )
                        
                        par(op)
                        dev.off()
                      }
                      
                      power_matrix=matrix(nrow=2,ncol=5)
                      rownames(power_matrix)=c("DESeq","DESeq_sep")
                      colnames(power_matrix)=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
                      
                      tryCatch({power_matrix[1,]=c(cal_range(deseq2_pval[which(mean_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(var_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(dp_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(mult_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)],0,0.05))}, error = function(e) {NA} )
                      tryCatch({power_matrix[2,]=c(cal_range(deseq2_pval_sep[which(mean_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(var_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(dp_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(mult_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)],0,0.05))}, error = function(e) {NA} )
                      power_matrix
                      
                      #barplot
                      png(paste0("./fig_DESeq2_compare_barplot/p",perm_label,perm_method,"_barplot_",param_tag,"_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 600,width = 1200)
                      op=par(mfrow = c(1, 2), pty = "s",cex=2)
                      barplot(power_matrix[1,],ylab="power",main=rownames(power_matrix)[1],ylim=c(0,1))
                      abline(h=0.05,col="red")
                      barplot(power_matrix[2,],ylab="power",main=rownames(power_matrix)[2],ylim=c(0,1))
                      abline(h=0.05,col="red")
                      par(op)
                      dev.off()
                      if((!is.na(deseq2_pval)) || (!is.na(deseq2_pval_sep)) ){
                        png(paste0("./fig_DESeq2_compare_power_point/p",perm_label,perm_method,"_power_point_",param_tag,"_",pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 800,width = 800)
                        plot(power_matrix[1,5],power_matrix[1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3)
                        points(power_matrix[1,5],power_matrix[1,2],col="red",pch=4,cex=3)
                        points(power_matrix[1,5],power_matrix[1,3],col="red",pch=5,cex=3)
                        points(power_matrix[1,5],power_matrix[1,4],col="red",pch=6,cex=3)
                        
                        plot(power_matrix[2,5],power_matrix[2,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="blue",pch=3,cex=3)
                        points(power_matrix[2,5],power_matrix[2,2],col="blue",pch=4,cex=3)
                        points(power_matrix[2,5],power_matrix[2,3],col="blue",pch=5,cex=3)
                        points(power_matrix[2,5],power_matrix[2,4],col="blue",pch=6,cex=3)
                        
                        legend("topright",c(rownames(power_matrix),"mean diff","var diff","dp_diff","mult_diff"),pch=c(rep(15,2),3:8),cex=1,col=c("red","blue"))
                        
                        dev.off()
                      }
                      
                      
                      power_array[i_file,i_mean,i_var,i_mult,i_dp,i_ind,i_cell,,]=power_matrix
                      power_matrix=matrix(nrow=2,ncol=5)
                      
                      rownames(power_matrix)=c("DESeq","DESeq_sep")
                      colnames(power_matrix)=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
                      
                      tryCatch({power_matrix[1,]=c(cal_range(deseq2_pval[which(mean_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(var_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(dp_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(mult_index==1)],0,0.05),
                                                   cal_range(deseq2_pval[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)],0,0.05))}, error = function(e) {NA} )
                      tryCatch({power_matrix[2,]=c(cal_range(deseq2_pval_sep[which(mean_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(var_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(dp_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(mult_index==1)],0,0.05),
                                                   cal_range(deseq2_pval_sep[which(mean_index==0 & var_index==0 & dp_index==0 & mult_index==0)],0,0.05))}, error = function(e) {NA} )
                      power_matrix
                      range005_array[i_file,i_mean,i_var,i_mult,i_dp,i_ind,i_cell,,]=power_matrix           
                      print(paste0(pre_tag,fit_tag,r_mean,"_",r_var,"_",r_mult,"_",r_dp,"_",file_tag,"_",n_ind,"_",n_cell))
                    }
                  }
                }
              }
            }
          }
        }
        saveRDS(power_array,paste0("./DESeq2_compare_p",perm_label,perm_method,pre_tag,fit_tag,"_final_power_array_",param_tag,".rds"))
        saveRDS(range005_array,paste0("./DESeq2_compare_p",perm_label,perm_method,pre_tag,fit_tag,"_final_range005_array_",param_tag,".rds"))
        
        #more plot
        for(i_file in 1:length(file_tag_seq)){
          for(i_cell in 1:length(cell_seq)){
            for(i_ind in 1:length(ind_seq)){
              file_tag=file_tag_seq[i_file]
              n_ind=ind_seq[i_ind]
              n_cell=cell_seq[i_cell]
              
              cur_power_array=power_array[i_file,,,,,i_ind,i_cell,,]
              
              if((!is.na(deseq2_pval)) || (!is.na(deseq2_pval_sep)) ){
                png(paste0("./fig_DESeq2_compare_final_power/p",perm_label,perm_method,pre_tag,fit_tag,"_final_power_",param_tag,"_",file_tag,"_",n_ind,"_",n_cell,".png"),height = 800,width = 800)
                plot(cur_power_array[,1,5],cur_power_array[,1,1],xlim=c(0,1),ylim=c(0,1),xlab="False positive rate (FPR)",ylab="True positive rate (TPR)",type="p",col="red",pch=3,cex=3,main=paste0("power scatter: ",file_tag))
                points(cur_power_array[,1,5],cur_power_array[,1,2],col="red",pch=4,cex=3)
                points(cur_power_array[,1,5],cur_power_array[,1,3],col="red",pch=5,cex=3)
                points(cur_power_array[,1,5],cur_power_array[,1,4],col="red",pch=6,cex=3)
                
                points(cur_power_array[,2,5],cur_power_array[,2,1],col="blue",pch=3,cex=3)
                points(cur_power_array[,2,5],cur_power_array[,2,2],col="blue",pch=4,cex=3)
                points(cur_power_array[,2,5],cur_power_array[,2,3],col="blue",pch=5,cex=3)
                points(cur_power_array[,2,5],cur_power_array[,2,4],col="blue",pch=6,cex=3)
                
                legend("topright",c(rownames(power_matrix),"mean diff","var diff","dp_diff","mult_diff"),pch=c(rep(15,2),3:8),cex=1,col=c("red","blue","black","black","black","black"))
                dev.off()
              }
            }
          }
        }
      }
    }
  }
}












#this code follows 10.4.1.*R and use the output from it.

file_tag_seq=1
perm_label=0
perm_method=""

param_tag=1
power_array=list()
param_tag_seq=c("mean","var","dp","mult")

perm_method_seq="" #c("","b")
pre_tag_seq=c("","dca")

pre_tag=""

if(pre_tag==""){fit_tag=""}
if(pre_tag=="dca"){fit_tag="nb"}

for(i_param in 1:length(param_tag_seq)){
  param_tag=param_tag_seq[i_param]
  power_array[[param_tag]]=readRDS(paste0("./DESeq2_compare_p",perm_label,perm_method,pre_tag,fit_tag,"_final_range005_array_",param_tag,".rds"))
  power_array[[param_tag]][power_array[[param_tag]]==0]=NA
}


ind_seq=c(10,20,40,60,100)
cell_seq=c(20,50,100,200)

library("RColorBrewer")
#power_plot is used for plotting the power curves
#note: y_matrix max column=12

power_plot1=function(y_matrix,cur_main="",cur_xlab="",cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA){
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
    cur_col=brewer.pal(ncol(y_matrix),"Paired")
  }
  plot(x_seq,y_matrix[,1],type="l",main=cur_main,xlim=cur_xlim,ylim=cur_ylim,xlab=cur_xlab,ylab=cur_ylab,col=cur_col[1],lwd=3)
  for(i_c in 1:ncol(y_matrix)){
    lines(x_seq,y_matrix[,i_c],col=cur_col[i_c],lwd=3)
    points(x_seq,y_matrix[,i_c],col=cur_col[i_c],pch=i_c,cex=2)
  }
  
  legend("topleft",cur_legend,pch=1:ncol(y_matrix),col=cur_col,cex=.5)
}


pdf(paste0("./fig_DESeq2_compare_power_curve/p",perm_label,perm_method,pre_tag,fit_tag,"_power_curve.pdf"),height = 12,width = 30)

op=par(mfrow = c(2, 5),cex=1.5)
for(i_ind in 1:length(ind_seq)){
  for(i_cell in 1:length(cell_seq)){
    for(i_file in 1:length(file_tag_seq)){
      
      file_tag=file_tag_seq[i_file]
      
      for(i_param in 1:length(param_tag_seq)){
        param_tag=param_tag_seq[i_param]
        power_cur=power_array[[param_tag]][i_file,,,,,i_ind,i_cell,,i_param]
        print(power_cur)
        power_plot1(power_cur,cur_main=paste0("power of ",param_tag,", n ",ind_seq[i_ind],", ncell ",cell_seq[i_cell],", rep ",file_tag),cur_xlab=paste0("DE ",param_tag_seq[i_param]," change rate"),cur_ylab="Power",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
      }
      i_param=i_param+1
      power_cur=power_array[[param_tag]][i_file,,,,,i_ind,i_cell,,i_param]
      print(power_cur)
      ##option1: still curve plot
      #power_plot1(power_cur,cur_main=paste0("False Positive of ",param_tag,", n ",ind_seq[i_ind],", ncell ",cell_seq[i_cell],", rep ",file_tag),cur_xlab=paste0("DE ",param_tag," change rate"),cur_ylab="Type I error",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
      ##option2: barplot
      barplot(colMeans(power_cur, na.rm = TRUE, dims = 1),main=paste0("FDR of ",param_tag,", n ",ind_seq[i_ind],", ncell ",cell_seq[i_cell],", rep ",file_tag),xlab="",las=2,ylab="Type I error",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA,ylim=c(0,0.2),col=brewer.pal(ncol(power_cur),"Paired"))
      abline(h=0.05,col="red")
      
    }
  }
}


par(op)
dev.off()
