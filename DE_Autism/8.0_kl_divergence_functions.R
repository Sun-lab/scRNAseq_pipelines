############functions #####################################


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

calc_KL = function(px,qx){
  -sum(px * (log(qx)-log(px)))
}


mean_KL=function(x,y,nbins=10,rangep=NA,bmeth="equal"){ 
  res=calcBinCounts(x,y,nbin=nbins,rangepoint=rangep,bmethod=bmeth)
  px=res$count_x+1
  qx=res$count_y+1
  px=px/sum(px)
  qx=qx/sum(qx)
  kl_count=(calc_KL(px,qx)+calc_KL(qx,px))/2
  return(kl_count)
}
