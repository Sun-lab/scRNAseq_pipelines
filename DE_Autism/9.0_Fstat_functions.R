
########################################Functions############################################
#cal F version 1#############
#dist_matrix=nxn symmetric matrix with diag =0
#label nx1 vector about categories.
calc_F_permanova=function(dist_matrix,label,covariate_x=NA){
  label=as.numeric(label)
  a=length(unique(label))
  n=length(label)
  
  epsilon=matrix((rep(label,time=n)==rep(label,each=n)),ncol=n,nrow=n)+0
  d2=dist_matrix*dist_matrix
  sst=sum(d2/(2*n),na.rm = TRUE)
  ssw=sum((d2*epsilon)/(2*n),na.rm = TRUE)
  
  Fstat=((sst-ssw)*(n-a))/(ssw*(a-1))
  return(Fstat)
}

calc_F_permanova2=function(dist_array,label,covariate_x=NA){
  label=as.numeric(label)
  a=length(unique(label))
  n=length(label)
  
  epsilon=matrix((rep(label,time=n)==rep(label,each=n)),ncol=n,nrow=n)+0
  d2=dist_array*dist_array/(2*n)
  sst=apply(d2,1,function(x){sum(x,na.rm = TRUE)})
  ssw=apply(d2,1,function(x){sum(x*epsilon,na.rm = TRUE)})
  Fstat=((sst-ssw)*(n-a))/(ssw*(a-1))
  return(Fstat)
}

#cal F version 2############
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

calH=function(x){
  x=as.matrix(x)
  x%*% solve(crossprod(x))%*%t(x)
}

library("psych")
calTrace=function(H,G){
  tr(H%*%G%*%H)
}


#dist_matrix=nxn symmetric matrix with diag =0
#covariate_x=nxp covariates

calc_F_permanovaS=function(dist_matrix,label,covariate_x=NA,G_method=NA){ #G_method=c(NA,1,2,3)
  if(!is.na(covariate_x)){
    label=cbind(label,covariate_x)
  }
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
  
  label=as.matrix(label)
  label=as.matrix(apply(label,2,as.numeric))
  H=calH(label)
  IH=diag(nrow(H))-H
  Fstat=calTrace(G,H)/calTrace(G,IH)
  return(Fstat)
}



calc_F_permanovaS2=function(dist_array,label,covariate_x=NA,G_method=NA){ #G_method=c(NA,1,2,3)
  if(!is.na(covariate_x)){
    label=cbind(label,covariate_x)
  }
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
  label=as.matrix(label)
  label=as.matrix(apply(label,2,as.numeric))
  H=calH(label)
  IH=diag(nrow(H))-H
  Fstat=apply(G,3,function(x){return(calTrace(x,H)/calTrace(x,IH))})
  return(Fstat)
}



#example:
#cov_x=NA
#F_method=c("p","ps")
#perm_num.min=100
#perm_num.max=10000
#tol=1
#diagnose=rbinom(10,1,.5)
#dist_matrix=matrix(rnorm(100),10,10)

cal_permanova_pval=function(dist_matrix,diagnose,cov_x=NA,F_method="p",perm_num.min=500,perm_num.max=500000,tol=0.2){
  if(F_method=="p"){
    cur_cal_F=calc_F_permanova
  }
  if(F_method=="ps"){
    cur_cal_F=calc_F_permanovaS
  }
  n=length(diagnose)
  F_ob=cur_cal_F(dist_matrix,label=diagnose,covariate_x=cov_x)
  
  B=perm_num.min
  F_perm=t(matrix(ncol=1,nrow=B))
  F_perm=apply(F_perm,2,function(x){
    return(
      cur_cal_F(dist_matrix,covariate_x=cov_x,label=diagnose[sample.int(length(diagnose))])
           )})
  pval=sum(F_perm>F_ob,na.rm = TRUE)/sum(!is.na(F_perm))
  while(B<=perm_num.max){
    if(pval>=1/(tol*B)){
      break
    }
    else{
      B=B*10
      F_perm=t(matrix(ncol=1,nrow=B))
      F_perm=apply(F_perm,2,function(x){
        return(cur_cal_F(dist_matrix,covariate_x=cov_x,label=diagnose[sample.int(length(diagnose))]))})
      pval=sum(F_perm>F_ob,na.rm = TRUE)/sum(!is.na(F_perm))
    }
  }
  return(pval)
}

cal_permanova_pval2=function(dist_array,diagnose,cov_x=NA,F_method="p",perm_num.min=500,perm_num.max=500000,tol=0.2){
  if(F_method=="p"){
    cur_cal_F=calc_F_permanova2
  }
  if(F_method=="ps"){
    cur_cal_F=calc_F_permanovaS2
  }
  n=length(diagnose)
  F_ob=cur_cal_F(dist_array,label=diagnose,covariate_x=cov_x)
  
  B=perm_num.min
  F_perm=t(matrix(ncol=1,nrow=B))
  F_perm=apply(F_perm,2,function(x){
    return(
      cur_cal_F(dist_array,covariate_x=cov_x,label=diagnose[sample.int(length(diagnose))])
    )})
  F_perm=apply(F_perm,2,function(x){return(x-F_ob)})
  pval=apply(F_perm>=0,1,function(x){return(sum(x,na.rm = TRUE)/sum(!is.na(x)))})
  pval0_flag=(pval<1/(tol*B))
  while(B<=perm_num.max){
    if(sum(pval0_flag,na.rm = TRUE)==0){
      break
    }
    else{
      B=B*10
      cur_dist_array=dist_array[which(pval0_flag),,,drop=FALSE]
      cur_F_ob=as.numeric(F_ob[pval0_flag])
      F_perm=matrix(ncol=1,nrow=B)
      F_perm=as.matrix(apply(F_perm,1,function(x){
        return(cur_cal_F(cur_dist_array,covariate_x=cov_x,label=diagnose[sample.int(length(diagnose))]))}))
      F_perm=as.matrix(apply(F_perm,2,function(x){return(x-cur_F_ob)}))
      
      if(length(cur_F_ob)==1){
        cur_pval=sum(F_perm>0,na.rm = TRUE)/sum(!is.na(F_perm))
      }
      if(length(cur_F_ob)>1){
        cur_pval=apply(F_perm>0,1,function(x){return(sum(x,na.rm = TRUE)/sum(!is.na(x)))})
      }
      pval[which(pval0_flag)]=cur_pval
      pval0_flag=(pval<1/(tol*B))
    }
  }
  return(pval)
}


#real data test
#res=matrix(ncol=1,nrow=1000)
#res=apply(t(res),2,function(x){return(cal_permanova_pval(matrix(rnorm(n*n),n,n),rbinom(n,1,.5)))})
