
library("emdbook")
library("MASS")

############functions #####################################
calc_KL = function(px,qx){
  -sum(px * (log(qx/px)))
}

calc_JSD=function(px,qx){
  mx=(px+qx)/2
  -(sum(px * (log(mx/px)))+sum(qx * (log(mx/qx))))/2
}

#version 1: bin calculation
makeBins=function(xx,numbin=10,breakmethod="equal"){ #breakmethod=c("equal","quantile")
  min_xx=min(xx[!is.infinite(xx)& !is.na(xx)])
  max_xx=max(xx[!is.infinite(xx)& !is.na(xx)])
  if(breakmethod=="equal"){
    rangepoint=c(min_xx,(1:(numbin-1))/numbin *max_xx,max_xx)
  }
  if(breakmethod=="quantile"){
    rangepoint=c(min_xx,quantile(xx,probs=c(1:(numbin-1))/numbin),max_xx)
  }
  return(rangepoint)
}

calcBinCounts = function(x,y,nbin=10,rangepoint=NA,bmethod="equal"){
  if(is.na(rangepoint)){
    rangepoint=makeBins(c(x,y),numbin=nbin,breakmethod = bmethod)
  }
  xhist=hist(x,breaks=rangepoint,plot=F)
  yhist=hist(y,breaks=rangepoint,plot=F)
  res=list(count_x=xhist$count,count_y=yhist$count)
  return(res)
}

mean_KL_hist=function(x,y,nbins=10,rangep=NA,bmeth="equal"){ 
  res=calcBinCounts(x,y,nbin=nbins,rangepoint=rangep,bmethod=bmeth)
  px=res$count_x+1
  qx=res$count_y+1
  px=px/sum(px)
  qx=qx/sum(qx)
  kl_count=(calc_KL(px,qx)+calc_KL(qx,px))/2
  return(kl_count)
}
#version 2:


calDensity=function(x,y,n=1024){
  cur_range=range(c(x,y),na.rm=TRUE,finite = TRUE)
  x_dens=density(x,from=cur_range[1],to=cur_range[2],n=n,na.rm = TRUE)
  y_dens=density(y,from=cur_range[1],to=cur_range[2],n=n,na.rm = TRUE)
  res=list(dens_x=x_dens$y,dens_y=y_dens$y)
  return(res)
  
}

#p_triple=c("mu", "size", "zprob")
#q_triple=c("mu", "size", "zprob")
#range_max: integer for density estimation max number.
#quantile: the minimium quantile of probability for characterizing the distribution. ignored if range_max>0
calDens_nbzinb=function(x_triple,y_triple,range_max=0,quantile=0.975){
  if(range_max==0){
    range_max=max(qnbinom(p=max(0,(quantile-x_triple[3])), size=x_triple[2], mu=x_triple[1]),
                  qnbinom(p=max(0,(quantile-y_triple[3])), size=y_triple[2], mu=y_triple[1]),
                  na.rm = TRUE)
  }
  x=emdbook::dzinbinom(0:range_max,mu=x_triple[1],size=x_triple[2],zprob=x_triple[3],log=FALSE)
  y=emdbook::dzinbinom(0:range_max,mu=y_triple[1],size=y_triple[2],zprob=y_triple[3],log=FALSE)
  res=list(dens_x=x,dens_y=y)
  return(res)
  
}


mean_KL_dens=function(x,y,empirical.n=1024,zinb.range_max=0,zinb.quantile=0.975,alter=c("mean","JSD"),fit_model=c("empirical","zinb")){ 
  if(fit_model=="zinb"){
    res=calDens_nbzinb(x,y,quantile=zinb.quantile,range_max = zinb.range_max)
    px=res$dens_x
    qx=res$dens_y
  }
  if(fit_model=="empirical"){
    res=calDensity(x,y,n=empirical.n)
    px=res$dens_x+1/empirical.n
    qx=res$dens_y+1/empirical.n
  }
  px=px/sum(px)
  qx=qx/sum(qx)
  if(alter=="mean"){
    kl_dens=(calc_KL(px,qx)+calc_KL(qx,px))/2
  }
  if(alter=="JSD"){
    kl_dens=calc_JSD(px,qx)
  }
  return(kl_dens)
}



##################Versions for direct #################################

#calculate density of each cell based on nbzinb
#vector_triple 
#nrow= eah cell for a given gene
#ncol= mean, dispersion dropout

calDens_nbzinb_cell=function(vector_triple,range_max=0,quantile=0.975){
  if(range_max==0){
    range_candidate=apply(vector_triple,1,function(x){return(qnbinom(p=max(0,(quantile-x[3])), size=x[2], mu=x[1]))})
    range_max=max(range_candidate,na.rm = TRUE)
  }
  dens_vector=apply(vector_triple,1,function(x){return(emdbook::dzinbinom(0:range_max,mu=x[1],size=x[2],zprob=x[3],log=FALSE))})
  return(dens_vector)
}


#combine the density info from each cell for each individual
#
#
calDens_nbzinb_ind=function(dens_vector,cell_ind=meta$individual){
  total_ind=unique(cell_ind)
  dens_vector_avg=matrix(ncol=length(total_ind),nrow=nrow(dens_vector))
  for(i_ind in 1:length(total_ind)){
    cur_ind=total_ind[i_ind]
    dens_vector_avg[,i_ind]=apply(dens_vector[,cell_ind==cur_ind,drop=FALSE],1,mean)
  }
  return(dens_vector_avg)
}

#calculate the dist array based on the zinb parameter of a group of cell

#usage
#dist_matrix=mean_KL_dens2(vector_triple=cell_param,cell_ind_label=meta$individual,alter="JSD")

mean_KL_dens2=function(vector_triple,cell_ind_label,cur_quantile=0.975,alter=c("mean","JSD")){ 
  cur_dens_vector=calDens_nbzinb_cell(vector_triple,range_max=0,quantile=cur_quantile)
  cur_dens_ind=calDens_nbzinb_ind(dens_vector=cur_dens_vector,cell_ind=cell_ind_label)
  cur_dens_ind=apply(cur_dens_ind,2,function(x){return(x/sum(x,na.rm = TRUE))})
  cur_dens_ind=cur_dens_ind+min(cur_dens_ind[cur_dens_ind>0])
  total_ind=unique(cell_ind_label)
  n_ind=length(total_ind)
  kl_dens=matrix(0,ncol=n_ind,nrow=n_ind)
  rownames(kl_dens)=total_ind
  colnames(kl_dens)=total_ind
  if(alter=="mean"){
    for(ia in 1:(n_ind-1)){
      for(ib in (ia+1):n_ind){
        px=cur_dens_ind[,ia]
        qx=cur_dens_ind[,ib]
        kl_dens[ia,ib]=(calc_KL(px,qx)+calc_KL(qx,px))/2
        kl_dens[ib,ia]=kl_dens[ia,ib]
      }
    }
    
  }
  if(alter=="JSD"){
    for(ia in 1:(n_ind-1)){
      for(ib in (ia+1):n_ind){
        px=cur_dens_ind[,ia]
        qx=cur_dens_ind[,ib]
        kl_dens[ia,ib]=calc_JSD(px,qx)
        kl_dens[ib,ia]=kl_dens[ia,ib]
      }
    }
  }
  return(kl_dens)
}


