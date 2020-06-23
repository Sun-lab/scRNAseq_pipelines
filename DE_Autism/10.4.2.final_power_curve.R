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

file_tag_seq=1
perm_label=0
perm_method=""
param_tag=1

ind_seq=c(10,20,40) #need to be correspond to its in 10.4.1
fit_tag_seq=""
cell_seq=c(100,200,400)#need to be correspond to its in 10.4.1
#i_cell=1
library("RColorBrewer")
#power_plot is used for plotting the power curves
#note: y_matrix max column=12

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
    #cur_col=brewer.pal(ncol(y_matrix),"Paired")
    #cur_col=brewer.pal(ncol(y_matrix),"Pastel1")
    cur_col=brewer.pal(ncol(y_matrix),"Set1")
  }
  plot(x_seq,y_matrix[,1],type="l",main=cur_main,sub=cur_sub,xlim=cur_xlim,ylim=cur_ylim,xlab=cur_xlab,ylab=cur_ylab,col=cur_col[1],lwd=3)
  for(i_c in 1:ncol(y_matrix)){
    lines(x_seq,y_matrix[,i_c],col=cur_col[i_c],lwd=3)
    points(x_seq,y_matrix[,i_c],col=cur_col[i_c],pch=i_c,cex=2)
  }
  
  legend("topleft",cur_legend,pch=1:ncol(y_matrix),col=cur_col,cex=.5)
}



#if(pre_tag==""){fit_tag=""}
#if(pre_tag=="dca"){fit_tag="nb"}

pdf(paste0("./fig_power_curve/p",perm_label,perm_method,"_",ind_seq,"_power_curve.pdf"),height = 25,width = 30)
op=par(mfrow = c(4, 5),cex=1.5,mar=c(8,5,3,3))
for(pre_tag in pre_tag_seq){
  for(fit_tag in fit_tag_seq){
    for(covariate_flag in covariate_flag_seq){
      for(resid_flag in resid_flag_seq){
        
        power_array=list()
        param_tag_seq=c("mean","var","dp","mult")
        
        for(i_param in 1:length(param_tag_seq)){
          param_tag=param_tag_seq[i_param]
          power_array[[param_tag]]=readRDS(paste0("./p",perm_label,perm_method,resid_flag,covariate_flag,pre_tag,fit_tag,"_final_range005_array_",param_tag,".rds"))
          power_array[[param_tag]][power_array[[param_tag]]==0]=NA
        }
        
        #pdf(paste0("./fig_power_curve/p",perm_label,perm_method,resid_flag,covariate_flag,pre_tag,fit_tag,"_power_curve.pdf"),height = 48,width = 30)
        
        #op=par(mfrow = c(8, 5),cex=1.5)
        for(i_ind in 1:length(ind_seq)){
          for(i_cell in 1:length(cell_seq)){
            for(i_file in 1:length(file_tag_seq)){
              
              file_tag=file_tag_seq[i_file]
              
              for(i_param in 1:length(param_tag_seq)){
                param_tag=param_tag_seq[i_param]
                power_cur=power_array[[param_tag]][i_file,,,,,i_ind,i_cell,,i_param]
                power_cur=power_cur[,c(1,2,7)] #no klmean
                print(power_cur)
                colnames(power_cur)=c("DESeq2","MAST","IDEAS")
                power_plot1(power_cur,
                            cur_main=paste0("Power of ",param_tag,"-DE, n ",ind_seq[i_ind],", ncell ",cell_seq[i_cell]),
                            cur_xlab=paste0("DE ",param_tag_seq[i_param]," change rate"),cur_ylab="Power",
                            cur_sub=paste0(pre_tag," ",fit_tag," ",covariate_flag," ",resid_flag),
                            x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
              }
              i_param=i_param+1
              power_cur=power_array[[param_tag]][i_file,,,,,i_ind,i_cell,,i_param]
              power_cur=power_cur[,c(1,2,7)] #no klmean
              colnames(power_cur)=c("DESeq2","MAST","IDEAS")
              print(power_cur)
              ##option1: still curve plot
              #power_plot1(power_cur,cur_main=paste0("False Positive of ",param_tag,", n ",ind_seq[i_ind],", ncell ",cell_seq[i_cell]),cur_xlab=paste0("DE ",param_tag," change rate"),cur_ylab="Type I error",x_seq=NA,cur_xlim=NA,cur_ylim=NA,cur_legend=NA,cur_col=NA)
              ##option2: barplot
              barplot(colMeans(power_cur, na.rm = TRUE, dims = 1),main=paste0("FDR, n ",ind_seq[i_ind],", ncell ",cell_seq[i_cell]),xlab="",las=2,ylab="Type I error",x_seq=NA,cur_xlim=NA,cur_ylim=c(0,1),cur_legend=NA,cur_col=NA,ylim=c(0,0.2),col=brewer.pal(ncol(power_cur),"Set1"))
              abline(h=0.05,col="red")
              
            }
          }
        }
        
        #par(op)
        #dev.off()
        
      }
    }
  }
}
par(op)
dev.off()

