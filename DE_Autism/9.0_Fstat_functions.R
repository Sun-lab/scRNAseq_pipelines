
########################################Functions############################################
#############pre analysis section ####################################

library("vegan")
#dist_matrix=nxn symmetric matrix with diag =0
#label nx1 vector about categories.


cal_betadisper_pval=function(dist_matrix,label,disper_type=c("median","centroid"),div_method=c("anova","TukeyHSD")){
  flag=which(apply(dist_matrix,1,function(x){sum(is.na(x))==length(x)}))
  if(length(flag)>0){
    label=as.factor(label[-flag])
    dist_matrix=dist_matrix[-flag,-flag]
    dist_matrix[is.na(dist_matrix)]=median(dist_matrix)
  }
  
  dis=as.dist(dist_matrix)
  mod=betadisper(dis, label,type=disper_type)
  
  if(div_method=="anova"){
    pval=anova(mod)$"Pr(>F)"[1]
  }
  if(div_method=="TukeyHSD"){
    pval=TukeyHSD(mod)$group[4]
  }
  return(pval)
}




####MANOVA Section ###############

#dist_matrix=nxn symmetric matrix with diag =0
#label nx1 vector about categories.
calc_F_manova=function(dist_matrix,label){
  label=as.numeric(label)
  a=length(unique(label))
  N=length(label)
  
  epsilon=matrix((rep(label,time=N)==rep(label,each=N)),ncol=N,nrow=N)+0
  d2=dist_matrix*dist_matrix
  sst=sum(d2)/(2*N)
  ssw=0
  for(ik in 1:a){
    cur_index=which(label==(unique(label)[ik]))
    ssw=ssw+sum(d2[cur_index,cur_index])/(2*length(cur_index))
  }
  Fstat=((sst-ssw)*(N-a))/(ssw*(a-1))
  return(Fstat)
}

# calc_F_manova2_previous=function(dist_array,label){
#   label=as.numeric(label)
#   a=length(unique(label))
#   N=length(label)
#   d2=dist_array*dist_array
#   sst=apply(d2/(2*N),1,function(x){sum(x)})
#   ssw=matrix(0,nrow=dim(dist_array)[1],ncol=1)
#   for(ik in 1:a){
#     cur_index=which(label==(unique(label)[ik]))
#     ssw=ssw+apply(d2,1,function(x){sum(x[cur_index,cur_index])/(2*length(cur_index))})
#   }
#   #ssw=apply(d2,1,function(x){sum(x*epsilon,na.rm = TRUE)})
#   Fstat=((sst-ssw)*(N-a))/(ssw*(a-1))
#   return(Fstat)
# }

#d2=array(exp(rnorm(10000)),dim=c(100,10,10))
#d2label=c(1,1,1,2,2,2,3,3,3,3)
calc_F_manova2=function(dist_array,label){
  label=as.numeric(label)
  a=length(unique(label))
  N=length(label)
  Ng=dim(dist_array)[1]
  d2=dist_array*dist_array
  
  dim(d2)=c(Ng,N^2)
  sst=rowSums(d2)/(2*N)
  divid_vector=matrix(0,N,N)   #this is nxn indicator with 1/n_1, 1/n_2,...for genexnxn --> genex (nxn) 
  for(ik in 1:a){
    cur_index=which(label==(unique(label)[ik]))
    divid_vector[cur_index,cur_index]=1/(2*length(cur_index))
  }
  divid_vector2=c(divid_vector)
  ssw=d2 %*% divid_vector2
  #ssw=apply(d2,1,function(x){sum(x*epsilon,na.rm = TRUE)})
  Fstat=((sst-ssw)*(N-a))/(ssw*(a-1))
  return(Fstat)
}

####PERMANOVAS Section ###############
##1.distance matrix realted

#almost no difference in cal speed between cal1 and cal2, for possible better numerical resuls(no evidence, by feeling), we may prefer cal1.
calG_1=function(m){#m is a matrix dim2=dim3=n_individual
  m=m*m
  n=nrow(m)
  preGn=(diag(n)-matrix(1/n,n,n))
  -1/2*preGn%*%m%*%preGn
}

calG_2=function(m){ #m is a matrix dim2=dim3=n_individual
  m=m*m
  n=nrow(m)
  preGn=(diag(n)*n-matrix(1,n,n))
  -1/(2*n*n)*preGn%*%m%*%preGn
}

#more suitable for large matrix, say, larger than 40
calG_3=function(m){#m is a matrix dim2=dim3=n_individual
  m=m*m
  n=nrow(m)
  res=m+matrix(sum(m)/(n*n),n,n)-1/n *(matrix((rep(rowSums(m),times=n)+rep(colSums(m),each=n)),n,n))
  -1/2*res
}

calG_1a=function(m){#m is array dim1=n_genes dim2=dim3=n_individual
  m=m*m
  n=dim(m)[2]
  preGn=(diag(n)-matrix(1/n,n,n))
  res=apply(m,1,function(x){return(-1/2*preGn%*%x%*%preGn)})
  dim(res)=c(n,n,dim(m)[1])
  res
}

calG_2a=function(m){ #m is array dim1=n_genes dim2=dim3=n_individual
  m=m*m
  n=dim(m)[2]
  preGn=(diag(n)*n-matrix(1,n,n))
  res=apply(m,1,function(x){return(-1/(2*n*n)*preGn%*%x%*%preGn)})
  dim(res)=c(n,n,dim(m)[1])
  res
}

#more suitable for large matrix, say, larger than 40
calG_3a=function(m){#m is array dim1=n_genes dim2=dim3=n_individual
  m=m*m
  n=dim(m)[2]
  res=-1/2*apply(m,1,function(x){return(x+matrix(sum(x)/(n*n),n,n)-1/n *(matrix((rep(rowSums(x),times=n)+rep(colSums(x),each=n)),n,n)))})
  dim(res)=c(n,n,dim(m)[1])
  res
}

##2.meta data realted

#this function calculate the residuals
#eg
#cov_z = data.frame(cbind(seq(1:50),rnorm(50),rpois(50,1),runif(50,0,1)))
#zm=model.matrix(~.,cov_z)
#x = cbind(rbinom(50,1,prob=0.5),c(rbinom(25,1,prob=0.2),rbinom(25,1,prob=0.8)) ,sort(rnorm(50,1)))
#res=cal_residual(zm,x)
#returns a list, first element residual, second element estimated x, third is the glm flag(discrete flag)
cal_residual_perm=function(zm=NA,x,perm_num=500){
  x=as.matrix(x)
  if(sum(!is.na(zm))==0){
    zm=rep(1,nrow(x))
  }
  perm_res=array(dim=c(nrow(x),ncol(x),perm_num),dimnames=list(rownames(x),colnames(x),1:perm_num))
  for(ix in 1:ncol(x)){
    cur_x=as.numeric(x[,ix])
    if(length(unique(cur_x))==2){
      #step1: fit
      m1=glm(cur_x~0+zm,family=binomial(link="logit"))  #logit model
      cur_res=residuals(m1,type = "response")  #cal residuals
      cur_fit=fitted(m1) #cal fit
      
      #step2: permutation
      for(ip in 1:perm_num){
        perm_x=rbinom(length(cur_fit),1,prob=cur_fit)
        m2=glm(perm_x~0+zm,family=binomial(link="logit"))
        perm_res[,ix,ip]=residuals(m2,type = "response")  
      }
    }
    else{
      #step1: fit
      m1=lm(cur_x~0+zm)  #linear model
      cur_res=residuals(m1,type = "response")   #cal residuals
      cur_fit=fitted(m1) #cal fit
      
      #step2: permutation
      for(ip in 1:perm_num){
        perm_x=cur_res[sample.int(length(cur_res),length(cur_res))]+cur_fit 
        m2=lm(perm_x~0+zm)  #logit model
        perm_res[,ix,ip]=residuals(m2,type = "response")  #cal residuals
      }
    }
  }
  return(perm_res)
}

cal_residual_ob=function(zm=NA,x){
  x=as.matrix(x)
  if(sum(!is.na(zm))==0){
    zm=rep(1,nrow(x))
  }
  ob_res=array(dim=c(nrow(x),ncol(x),1),dimnames=list(rownames(x),colnames(x),1))
  for(ix in 1:ncol(x)){
    cur_x=as.numeric(x[,ix])
    if(length(unique(cur_x))==2){
      m1=glm(cur_x~0+zm,family=binomial(link="logit"))  #logit model
      ob_res[,ix,1]=residuals(m1,type = "response")   #cal residuals
    }
    else{
      m1=lm(cur_x~0+zm)  #linear model
      ob_res[,ix,1]=residuals(m1,type = "response")   #cal residuals
    }
  }
  return(ob_res)
}

##3.final res realted
calH=function(x){
  x=as.matrix(x)
  x%*% solve(crossprod(x))%*%t(x)
}

library("psych")
calTrace=function(G,H){
  tr(H%*%G%*%H)
}

#dist_matrix=nxn symmetric matrix with diag =0
#covariate_x=nxp covariates to test, it should be adjusted if tehre is covariates
#label is a observed(nxp2x1) or permutated(nxp2xperm_num) array
#label=diagnose

calc_F_permanovaS=function(dist_matrix,label_array,G_method=NA){ #G_method=c(NA,1,2,3)
  if(!is.na(G_method)){
    if(G_method==1){
      G=calG_1(dist_matrix)
    }
    if(G_method==2){
      G=calG_2(dist_matrix)
    }
    if(G_method==3){
      G=calG_3(dist_matrix)
    }
  }
  if(is.na(G_method)){
    if(nrow(dist_matrix)<=40){
      G=calG_1(dist_matrix)
    }
    if(nrow(dist_matrix)>40){
      G=calG_3(dist_matrix)
    }
  }
  Fstat_total=matrix(ncol=1,nrow=dim(label_array)[3])
  for(i_label in 1:dim(label_array)[3]){
    label=label_array[,,i_label,drop=FALSE]
    #label=as.matrix(label)
    label=as.matrix(apply(label,2,as.numeric))
    H=calH(label)
    IH=diag(nrow(H))-H
    Fstat_total[i_label]=calTrace(G,H)/calTrace(G,IH)
  }
  return(Fstat_total)
}

calc_F_permanovaS2=function(dist_array,label_array,G_method=NA){ #G_method=c(NA,1,2,3)
  if(!is.na(G_method)){
    if(G_method==1){
      G=calG_1a(dist_array)
    }
    if(G_method==2){
      G=calG_2a(dist_array)
    }
    if(G_method==3){
      G=calG_3a(dist_array)
    }
  }
  if(is.na(G_method)){
    if(nrow(dist_array)<=40){
      G=calG_1a(dist_array)
    }
    if(nrow(dist_array)>40){
      G=calG_3a(dist_array)
    }
  }
  
  Fstat_total=matrix(nrow=dim(dist_array)[1],ncol=dim(label_array)[3])
  for(i_label in 1:dim(label_array)[3]){
    label=label_array[,,i_label]
    label=as.matrix(label)
    label=as.matrix(apply(label,2,as.numeric))
    H=calH(label)
    IH=diag(nrow(H))-H
    Fstat_total[,i_label]=apply(G,3,function(x){return(calTrace(x,H)/calTrace(x,IH))})
  }
  return(Fstat_total)
}


#do balanced permutation. 
#require case num < control num, and case num are expected to be equal or similar to control num
#input label is an indicator vactor with 1 and 0

#eg:
#a_label=c(rep(1,10),rep(0,10))
#balanced_perm(a_label)

balanced_perm=function(label){
  n_len=length(label)
  
  case_index=which(label==1)
  ctrl_index=which(label==0)
  
  n_case=length(case_index)
  n_ctrl=length(ctrl_index)
  
  n_exchange=n_case/2
  n_exchange=floor(n_exchange)+(n_exchange-floor(n_exchange))*2*rbinom(1,1,0.5) #the number changed to other side
  
  switch_index_case=sample(case_index,n_exchange)
  switch_index_ctrl=sample(ctrl_index,n_exchange)
  total_index=c(switch_index_case,switch_index_ctrl)
  label[total_index]=1-label[total_index]

  label
}



# #example:
# zm=NA
# cov_z=matrix(rnorm(50,5,1),10,5)
# cov_z = data.frame(cbind(seq(1:50),rnorm(50),rpois(50,1),runif(50,0,1)))
# zm = model.matrix(~.,cov_z)

# perm_num.min=100
# 
# tol=1
# diagnose=rbinom(10,1,.5)
# #diagnose2=rnorm(10,1,1)
# 
# dist_array=array(dim=c(5,10,10))
# for(k in 1:5){
#   dist_matrix=matrix(rnorm(100),10,10)
#   for(i in 1:10){
#     dist_matrix[i,i]=0
#     for(j in i:10){
#       dist_matrix[j,i]=dist_matrix[i,j]
#     }
#   }
#   dist_array[k,,]=dist_matrix
# }
# 
# cov_z = data.frame(cbind(seq(1:50),rnorm(50),rpois(50,1),runif(50,0,1)))
# zm=model.matrix(~.,cov_z)
# 
# R_array=cal_residual_perm(zm,diagnose,perm_num = 50)
# R_ob=cal_residual_ob(zm,diagnose)
# calc_F_permanovaS(dist_matrix,label_array=R_ob)
# calc_F_permanovaS2(dist_array,label_array=R_array)


cal_permanova_pval=function(dist_matrix,diagnose,zm=NA,Fstat_method="p",perm_num.min=500,perm_method=""){
  n=length(diagnose)
  if(Fstat_method=="p"){
    F_ob=calc_F_manova(dist_matrix,label=diagnose)
    F_perm=t(matrix(ncol=1,nrow=perm_num.min))
    
    if(length(grep("b",perm_method))>0){ #balanced permutation
      F_perm=apply(F_perm,2,function(x){
        return(
          calc_F_manova(dist_matrix,label=balanced_perm(diagnose))
        )})
    }
    
    if(length(grep("b",perm_method))==0){ #normal permutation
      F_perm=apply(F_perm,2,function(x){
        return(
          calc_F_manova(dist_matrix,label=diagnose[sample.int(length(diagnose))])
        )})
      }
    
  }
  if(Fstat_method=="ps"){
    R_array=cal_residual_perm(zm,diagnose,perm_num = perm_num.min)
    R_ob=cal_residual_ob(zm,diagnose)
    F_ob=as.numeric(calc_F_permanovaS(dist_matrix,label_array=R_ob))
    F_perm=as.numeric(calc_F_permanovaS(dist_matrix,label_array=R_array))
  }
  pval=mean(F_perm-F_ob>=0)
  #pval=sum(F_perm-F_ob>=0,na.rm = TRUE)/sum(!is.na(F_perm))
  res=list()
  res[["pval"]]=pval
  res[["F_ob"]]=F_ob
  res[["F_perm"]]=F_perm
  return(res)
}

cal_permanova_pval2=function(dist_array,diagnose,zm=NA,Fstat_method="p",perm_num.min=500,perm_method=""){
  n=length(diagnose)
  if(Fstat_method=="p"){
    F_ob=calc_F_manova2(dist_array,label=diagnose)
    F_perm=t(matrix(ncol=1,nrow=perm_num.min))
    
    if(length(grep("b",perm_method))>0){ #balanced permutation
      F_perm=apply(F_perm,2,function(x){
        return(
          calc_F_manova2(dist_array,label=balanced_perm(diagnose))
        )})
    }
    
    if(length(grep("b",perm_method))==0){ #normal permutation
      F_perm=apply(F_perm,2,function(x){
        return(
          calc_F_manova2(dist_array,label=diagnose[sample.int(length(diagnose))])
        )})
    }
    
    

  }
  if(Fstat_method=="ps"){
    R_array=cal_residual_perm(zm,diagnose,perm_num = perm_num.min)
    R_ob=cal_residual_ob(zm,diagnose)
    F_ob=calc_F_permanovaS2(dist_array,label_array=R_ob)
    F_perm=calc_F_permanovaS2(dist_array,label_array=R_array)
  }
  F_perm=as.matrix(F_perm)
  F_perm0=apply(F_perm,2,function(x){return(x-F_ob)})
  if(length(F_ob)==1){
    F_perm0=t(F_perm0)
  }
  pval=apply(F_perm0>=0,1,mean)
  #pval=apply(F_perm0>=0,1,function(x){return(sum(x,na.rm = TRUE)/sum(!is.na(x)))})
  res=list()
  res[["pval"]]=pval
  res[["F_ob"]]=F_ob
  res[["F_perm"]]=F_perm
  return(res)
}

