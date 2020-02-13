#this code check the uneven p-values proportion under centrain senarios of real data. It is generated from the following ones.
#------------------------------from 11.check_uneven_pval_from_simulation.R-----------
#this code check the uneven p-values proportion under centrain senarios of simulation.
#It follows up 10.4.1.final_power_barplot.R 
#Tt use the range005_array,range095_array,range46_array as its input.
#range005_array,range095_array,range46_array records the proportions of p-values falls into c(0,0.05),c(0.95,1) and c(0.4,0.6)

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
#These output can be defined for each method("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct"), each datasubset("mean_diff","var_diff","disp_diff","mult_diff","control(FDR)")


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
#setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")


#setwd("~/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Users/mzhang24/Desktop/fh/1.Testing_scRNAseq/")
setwd("/Volumes/SpecialData/fh_data/Data_PRJNA434002/")
#setwd("/fh/fast/sun_w/mengqi/1.Testing_scRNAseq/")

cluster_tag_seq=1:17
file_tag_seq=c("3k10","5k")
pre_tag_seq=c("dca","scvi")
dist_method_seq=c("klmean","jsd")
fit_method_seq=c("empirical","nbzinb")
F_method_seq=c("p","ps")


file_tag_seq=c("5k","3k10")
pre_tag_seq=c("scvi","dca")

perm_label_seq=0:5
ind_covariate_flag=NA

perm_method=""

##################functions #################################

#The functions, given the 3 same size input array, it helps calculate the proportion of each array that their number less than the given threshold, and the proportion of the union of those extreme pvalues. 

#input: 3 same shape array, a005,a095, a46
#output, a 4x3 matrix with 
#     nrow= categories of uneven pval   c("FP_inflation","Loose Power","strange median","any") 
#     ncol= degree of extreme c("sign","moderate","extreme")
#the element of each matrix shows the proportion of uneven pval,under certain conditions

cal_prob_array=function(a005,a095,a46,ath=matrix(c(0.07,0.07,0.05,0.1,0.1,0.1,0.2,0.2,0.2),3,3)){
  x=matrix(ncol=3,nrow=4)
  rownames(x)=c("FP_inflation","Loose Power","strange median","any")
  colnames(x)=c("sign","moderate","extreme")
  
  x[1,1]=sum(a005>=ath[1,1],na.rm = TRUE)/sum(a005<=2,na.rm = TRUE)
  x[2,1]=sum(a095>=ath[2,1],na.rm = TRUE)/sum(a095<=2,na.rm = TRUE)
  x[3,1]=sum(a46>=(0.2+ath[3,1]) | a46<=(0.2-ath[3,1]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  x[4,1]=sum(a005>=ath[1,1] | a095>=ath[2,1] | a46>=(0.2+ath[3,1]) | a46<=(0.2-ath[3,1]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  
  x[1,2]=sum(a005>=ath[1,2],na.rm = TRUE)/sum(a005<=2,na.rm = TRUE)
  x[2,2]=sum(a095>=ath[2,2],na.rm = TRUE)/sum(a095<=2,na.rm = TRUE)
  x[3,2]=sum(a46>=(0.2+ath[3,2]) | a46<=(0.2-ath[3,2]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  x[4,2]=sum(a005>=ath[1,2] | a095>=ath[2,2] | a46>=(0.2+ath[3,2]) | a46<=(0.2-ath[3,2]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  
  x[1,3]=sum(a005>=ath[1,3],na.rm = TRUE)/sum(a005<=2,na.rm = TRUE)
  x[2,3]=sum(a095>=ath[2,3],na.rm = TRUE)/sum(a095<=2,na.rm = TRUE)
  x[3,3]=sum(a46>=(0.2+ath[3,3]) | a46<=(0.2-ath[3,3]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  x[4,3]=sum(a005>=ath[1,3] | a095>=ath[2,3] | a46>=(0.2+ath[3,3]) | a46<=(0.2-ath[3,3]) ,na.rm = TRUE)/sum(a46<=2,na.rm = TRUE)
  return(x)
}

method_seq=c("DESeq","MAST","jsd_empirical","klmean_empirical","jsd_zinb","klmean_zinb","jsd_direct","klmean_direct")

range005_array=readRDS("../Data_PRJNA434002/8.Result/final_range005_array.rds")
range095_array=readRDS("../Data_PRJNA434002/8.Result/final_range095_array.rds")
range46_array=readRDS("../Data_PRJNA434002/8.Result/final_range46_array.rds")

range_total=cal_prob_array(range005_array,range095_array,range46_array)

write.csv(range_total, paste0("../Data_PRJNA434002/8.Result/fig_uneven_pval_prop/unevenpval_prop_total.csv"))

#[i_file,i_F,i_pre,i_cluster,i_perm_label,]

#for file
range_file=array(dim=c(length(method_seq),length(file_tag_seq),4,3),dimnames=list(method_seq,file_tag_seq,c("FP_inflation","Loose Power","strange median","any"),c("sign","moderate","extreme")))
for(i_m in 1:length(method_seq)){
  for(i_file in 1:length(file_tag_seq)){
    range_file[i_m,i_file,,]=cal_prob_array(range005_array[i_file,,,,,i_m],range095_array[i_file,,,,,i_m],range46_array[i_file,,,,,i_m])
  }
}
write.csv(range_file, paste0("../Data_PRJNA434002/8.Result/fig_uneven_pval_prop/unevenpval_prop_file.csv"))

#for F method
range_F=array(dim=c(length(method_seq),length(F_method_seq),4,3),dimnames=list(method_seq,F_method_seq,c("FP_inflation","Loose Power","strange median","any"),c("sign","moderate","extreme")))
for(i_m in 1:length(method_seq)){
  for(i_F in 1:length(F_method_seq)){
    range_F[i_m,i_F,,]=cal_prob_array(range005_array[,i_F,,,,i_m],range095_array[,i_F,,,,i_m],range46_array[,i_F,,,,i_m])
  }
}
write.csv(range_F, paste0("../Data_PRJNA434002/8.Result/fig_uneven_pval_prop/unevenpval_prop_F.csv"))

#for pre method
range_pre=array(dim=c(length(method_seq),length(pre_tag_seq),4,3),dimnames=list(method_seq,pre_tag_seq,c("FP_inflation","Loose Power","strange median","any"),c("sign","moderate","extreme")))
for(i_m in 1:length(method_seq)){
  for(i_pre in 1:length(pre_tag_seq)){
    range_pre[i_m,i_pre,,]=cal_prob_array(range005_array[,,i_pre,,,i_m],range095_array[,,i_pre,,,i_m],range46_array[,,i_pre,,,i_m])
  }
}
write.csv(range_pre, paste0("../Data_PRJNA434002/8.Result/fig_uneven_pval_prop/unevenpval_prop_pre.csv"))

#for cluster method
range_cluster=array(dim=c(length(method_seq),length(cluster_tag_seq),4,3),dimnames=list(method_seq,cluster_tag_seq,c("FP_inflation","Loose Power","strange median","any"),c("sign","moderate","extreme")))
for(i_m in 1:length(method_seq)){
  for(i_cluster in 1:length(cluster_tag_seq)){
    range_cluster[i_m,i_cluster,,]=cal_prob_array(range005_array[,,,i_cluster,,i_m],range095_array[,,,i_cluster,,i_m],range46_array[,,,i_cluster,,i_m])
  }
}
write.csv(range_pre, paste0("../Data_PRJNA434002/8.Result/fig_uneven_pval_prop/unevenpval_prop_pre.csv"))




