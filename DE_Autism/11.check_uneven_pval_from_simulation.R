#this code check the uneven p-values proportion under centrain senarios of simulation.
#It follows up 10.4.1.final_power_barplot.R 
#Tt use the range01_array,range09_array,range46_array as its input.
#range01_array,range09_array,range46_array records the proportions of p-values falls into c(0,0.05),c(0.95,1) and c(0.4,0.6)

#The output is the element of each matrix shows the proportion of uneven pval,under certain conditions
#such as
#     nrow= categories of uneven pval   c("FP_inflation",               #pval<=0.05 have a unexpected large porpotion
#                                          "Loose Power",               #pval>=0.95 have a unexpected large porpotion
#                                          "strange median",            #pval falls into 0.4-0.6 are too less or too much.
#                                          "any")                       #any of above
#
#     ncol= degree of extreme c( "sign",                                #standard are low for the juedgement of categories of uneven pval, its easy for pvals to falls into the certain uneven cateogries.
#                                "moderate",                            #standard are moderate for the juedgement of categories of uneven pval
#                                "extreme")                             #standard are high for the juedgement of categories of uneven pval, pvals pass it should really extreme shape to fall into uneven categories.
#  
# The ath is the theshold for 
#                 "sign"  "moderate" "extreme"
# FP_inflation:    0.07    0.1        0.2
# Loose Power:     0.07    0.1        0.2
# strange median:  0.05    0.1        0.2  (this is for single side distance between the expected proportion, here 0.6-0.4=0.2)
#These output can be defined for each method("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct"), each datasubset("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/10.Result/sim_v6")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

perm_label=0
perm_method=""
param_tag=2 #c(1,2,3,4)

perm_label_seq=0
param_tag_seq=1:4

#The functions, given the 3 same size input array, it helps calculate the proportion of each array that their number less than the given threshold, and the proportion of the union of those extreme pvalues. 

#input: 3 same shape array, a01,a09, a46
#output, a 4x3 matrix with 
#     nrow= categories of uneven pval   c("FP_inflation","Loose Power","strange median","any") 
#     ncol= degree of extreme c("sign","moderate","extreme")
#the element of each matrix shows the proportion of uneven pval,under certain conditions

#cal_prob_array=function(a005,a095,a46,ath=matrix(c(0.07,0.07,0.05,0.1,0.1,0.1,0.2,0.2,0.2),3,3)){
cal_prob_array=function(a01,a09,a46,ath=matrix(c(0.15,0.15,0.05,0.2,0.2,0.1,0.4,0.4,0.2),3,3)){
  x=matrix(ncol=3,nrow=4)
  rownames(x)=c("FP_inflation","Loose Power","strange median","any")
  colnames(x)=c("sign","moderate","extreme")
  
  x[1,1]=sum(a01>=ath[1,1],na.rm = TRUE)/sum(a01<=2,na.rm = TRUE)
  x[2,1]=sum(a09>=ath[2,1],na.rm = TRUE)/sum(a09<=2,na.rm = TRUE)
  x[3,1]=sum(a46>=(0.2+ath[3,1]) | a46<=(0.2-ath[3,1]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  x[4,1]=sum(a01>=ath[1,1] | a09>=ath[2,1] | a46>=(0.2+ath[3,1]) | a46<=(0.2-ath[3,1]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  
  x[1,2]=sum(a01>=ath[1,2],na.rm = TRUE)/sum(a01<=2,na.rm = TRUE)
  x[2,2]=sum(a09>=ath[2,2],na.rm = TRUE)/sum(a09<=2,na.rm = TRUE)
  x[3,2]=sum(a46>=(0.2+ath[3,2]) | a46<=(0.2-ath[3,2]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  x[4,2]=sum(a01>=ath[1,2] | a09>=ath[2,2] | a46>=(0.2+ath[3,2]) | a46<=(0.2-ath[3,2]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  
  x[1,3]=sum(a01>=ath[1,3],na.rm = TRUE)/sum(a01<=2,na.rm = TRUE)
  x[2,3]=sum(a09>=ath[2,3],na.rm = TRUE)/sum(a09<=2,na.rm = TRUE)
  x[3,3]=sum(a46>=(0.2+ath[3,3]) | a46<=(0.2-ath[3,3]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  x[4,3]=sum(a01>=ath[1,3] | a09>=ath[2,3] | a46>=(0.2+ath[3,3]) | a46<=(0.2-ath[3,3]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  return(x)
}


for(param_tag in param_tag_seq){ 
  for(perm_label in perm_label_seq){
    file_tag_seq=1
    
    r_mean_seq=1.2
    r_var_seq=1.2
    r_dp_seq=0.2
    r_mult_seq=0.2
    
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
      r_mult_seq=1:5/10
      param_tag="mult"
    }
    if(param_tag==4){
      r_dp_seq=1:4/10
      param_tag="dp"
    }
    
    range01_array=readRDS(paste0("./p",perm_label,perm_method,"_final_range01_array_",param_tag,".rds"))
    range09_array=readRDS(paste0("./p",perm_label,perm_method,"_final_range09_array_",param_tag,".rds"))
    range46_array=readRDS(paste0("./p",perm_label,perm_method,"_final_range46_array_",param_tag,".rds"))
    
    
    method_seq=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")
    diff_seq=c("mean_diff","var_diff","dp_diff","mult_diff","control(FDR)")
    
    range_method=array(dim=c(length(method_seq),4,3),dimnames=list(method_seq,c("FP_inflation","Loose Power","strange median","any"),c("sign","moderate","extreme")))
    range_data=array(dim=c(length(diff_seq),4,3),dimnames=list(diff_seq,c("FP_inflation","Loose Power","strange median","any"),c("sign","moderate","extreme")))
    range_method_data=array(dim=c(length(method_seq),length(diff_seq),4,3),dimnames=list(method_seq,diff_seq,c("FP_inflation","Loose Power","strange median","Total"),c("sign","moderate","extreme")))
    
    pdf(paste0("./fig_uneven_pval_prop/p",perm_label,perm_method,"_final_unevenpval_prop_",param_tag,".pdf"),width = 9,height = 15)
    op=par(mfrow = c(5,3))
    for(i_m in 1:8){
      for(i_diff in 1:5){
        if(length(range01_array[,,,,,,,i_m,i_diff])>sum(is.na(range01_array[,,,,,,,i_m,i_diff]))){
          hist(as.numeric(range01_array[,,,,,,,i_m,i_diff]),main=paste0("hist of prop pval<0.05, ",method_seq[i_m],"_",diff_seq[i_diff]),breaks=10,xlim=range(0,1))
          abline(v=0.1,col="red")
          hist(as.numeric(range09_array[,,,,,,,i_m,i_diff]),main=paste0("hist of prop pval>0.95, ",method_seq[i_m],"_",diff_seq[i_diff]),breaks=10,xlim=range(0,1))
          abline(v=0.1,col="red")
          hist(as.numeric(range46_array[,,,,,,,i_m,i_diff]),main=paste0("hist of prop pval 0.4-0.6, ",method_seq[i_m],"_",diff_seq[i_diff]),breaks=10,xlim=range(0,1))
          abline(v=0.3,col="red")
          
          
          
          range_method_data[i_m,i_diff,,]=cal_prob_array(a01=range01_array[,,,,,,,i_m,i_diff],a09=range09_array[,,,,,,,i_m,i_diff],a46=range46_array[,,,,,,,i_m,i_diff],ath=matrix(c(0.15,0.15,0.05,0.2,0.2,0.1,0.4,0.4,0.2),3,3))
          
          range_method[i_m,,]=cal_prob_array(a01=range01_array[,,,,,,,i_m,],a09=range09_array[,,,,,,,i_m,],a46=range46_array[,,,,,,,i_m,],ath=matrix(c(0.15,0.15,0.05,0.2,0.2,0.1,0.4,0.4,0.2),3,3))
          
          range_data[i_diff,,]=cal_prob_array(a01=range01_array[,,,,,,,,i_diff],a09=range09_array[,,,,,,,,i_diff],a46=range46_array[,,,,,,,,i_diff],ath=matrix(c(0.15,0.15,0.05,0.2,0.2,0.1,0.4,0.4,0.2),3,3))
          
        }
      }
    }
    par(op)
    dev.off()
    
    write.csv(range_method_data, paste0("./fig_uneven_pval_prop/p",perm_label,perm_method,"_unevenpval_prop_method_data_",param_tag,".csv"))
    write.csv(range_method, paste0("./fig_uneven_pval_prop/p",perm_label,perm_method,"_unevenpval_prop_method_",param_tag,".csv"))
    write.csv(range_data, paste0("./fig_uneven_pval_prop/p",perm_label,perm_method,"_unevenpval_prop_data_",param_tag,".csv"))
  }
}

View(range_method)
View(range_data[5,,])
View(range_method_data[,5,,])

