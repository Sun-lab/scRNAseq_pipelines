

##Part I TWO NB same Mean and Var, 0.5 prob mixture, compared with 3rd NB #################
# Instructions about calculation of mixtured NB
# 
# Our question:
# Suppose Z have equaliy probability to falls into two NB distribution NB(mu1,size1), NB(mu2,size2). Calculate Z's mean and var.
# 
# Our answer:
# Let
# mu1=mu3+t
# mu2=mu3-t
# size2=size1
# size3=size1/(1+(t/mu3)^2*(1+size1))
# 
# a1=rnbinom(n,mu=mu1,size=size1)
# a2=rnbinom(n,mu=mu2,size=size2)
# a3=rnbinom(2*n,mu=mu3,size=size3)
# mean(a1)
# mean(a2)
# mean(c(a1,a2))
# mean(a3)
# var(a1)
# var(a2)
# var(c(a1,a2))
# var(a3)
# 
# 
# (mean(a1))^2+(mean(a2))^2
# (mean(c(a1,a2)))^2
# 
# var(a1)+var(a2)
#   mean(a1^2)+mean(a2^2)-mean(a1)^2-mean(a1)^2
# 
# var(c(a1,a2))
#   mean(c(a1^2,a2^2))-mean(c(a1,a2))^2
#   mean(a1^2)/2+mean(a2^2)/2-0.25*(mean(a1)+mean(a2))^2
#   mean(a1^2)/2+mean(a2^2)/2-0.25*mean(a1)^2-0.25*mean(a2)^2-0.5*mean(a1)*mean(a2)
#   0.5*(var(a1)+var(a2))+0.25*mean(a1)^2+0.25*mean(a2)^2-0.5*mean(a1)*mean(a2)
#   0.5*(var(a1)+var(a2))+0.25*(mean(a1)-mean(a2))^2
#   
#   
#   0.5*(mu1+mu1*mu1/size1+mu2+mu2*mu2/size2)+0.25*mu1*mu1+0.25*mu2*mu2-0.5*mu1*mu2
#   let mu1=mu3+t, mu2=mu3-t
#   mean: mu3=(mu1+mu2)/2
#  variance:
#   mu3+0.5*(mu3+t)*(mu3+t)/size1+0.5*(mu3-t)*(mu3-t)/size2+t^2
#   mu3+0.5*mu3*mu3/size1+mu3*t/size1+0.5*t*t/size1 +0.5*mu3*mu3/size2-mu3*t/size2+0.5*t*t/size2+t^2
# 
#   mu3+0.5*mu3^2*(1/size1+1/size2)+mu3*t*(1/size1-1/size2)+t^2(1+0.5/size1+0.5/size2)
#   so mu3*mu3/size3=0.5*mu3^2*(1/size1+1/size2)+mu3*t*(1/size1-1/size2)+t^2(1+0.5/size1+0.5/size2)
#   so size3=1/(0.5*(1/size1+1/size2)+t*(1/size1-1/size2)/mu3+t^2(1+0.5/size1+0.5/size2)/mu3^2)
#   
#   let size1=size2 
#   so size3=1/(1/size1+(t/mu3)^2*(1+1/size1))
#      size3=size1/(1+(t/mu3)^2*(1+size1))

n=100000
size1=1
size2=1
#size2=size1

mu3_seq=c(2,3,4,5)
t_seq=c(0.1,0.2,0.5,0.8)

lm3=length(mu3_seq)
lt=length(t_seq)

res_mean_1=matrix(,lm3,lt)
res_mean_2=matrix(,lm3,lt)
res_var_1=matrix(,lm3,lt)
res_var_2=matrix(,lm3,lt)

op=par(mfrow = c(2, 1))
for(i_mu3 in 1:lm3){
  for(i_t in 1:lt){
    mu3=mu3_seq[i_mu3]
    t=t_seq[i_t]*mu3
    mu1=mu3+t
    mu2=mu3-t
    if(size1==size2){
      size3=size1/(1+(t/mu3)^2*(1+size1))
    }
    if(size1!=size2){
      size3=2/(2*(t/mu3)^2 + (1+t/mu3)^2/size1 + (1-t/mu3)^2/size2)
    }
    
    
    a1=rnbinom(n,mu=mu1,size=size1)
    a2=rnbinom(n,mu=mu2,size=size2)
    a3=rnbinom(2*n,mu=mu3,size=size3)
    mean(a1)
    mean(a2)
    res_mean_1[i_mu3,i_t]=mean(c(a1,a2))
    res_mean_2[i_mu3,i_t]=mean(a3)
    var(a1)
    var(a2)
    res_var_1[i_mu3,i_t]=var(c(a1,a2))
    res_var_2[i_mu3,i_t]=var(a3)
    hist(a3,breaks=50)
    hist(c(a1,a2),col=rgb(1,0,0,0.3),breaks=50)
   
  }
}
par(op)

plot(res_mean_1,res_mean_2)
lines(c(0,20),c(0,20))
plot(res_var_1,res_var_2)
lines(c(0,20),c(0,20))

##Part II TWO ZINB same Mean and Var, 0.5 prob mixture, compared with 3rd ZINB ##################
#for the ZINB situation, we have:
a1~ZINB(mu1,size1,drop1) n samples
a1~ZINB(mu2,size2,drop2) n samples
a3~ZINB(mu3, size3,drop3) 2n samples

suppose 
drop1=drop2=drop3

mean1=mean3+t
mean2=mean3-t

library("emdbook")
#param_nbzinb returns the parameter set for nb/zinb to make sure that the parameters has the same mean and variance as the given 2 nb/zinb distributions

cal_nbzinb_param_multimodality=function(mu1=NA,size1=4,drop1=0,mu2=NA,size2=4,drop2=0,mean1=NA,mean2=NA){
  drop3=(drop1+drop2)/2
  if(is.na(mu1)||is.na(mu2)){
    mean3=(mean1+mean2)/2
    
    mu1=mean1/(1-drop1)
    mu2=mean2/(1-drop2)
    mu3=mean3/(1-drop3)
    
    
  }
  if(is.na(mean1)||is.na(mean2)){
    mu3=(mu1+mu2)/2
    
    mean1=(1-drop1)*mu1
    mean2=(1-drop2)*mu2
    mean3=(1-drop3)*mu3
  }
  
  size3=1/(((1/4*(mean1-mean2)^2+1/2*(mean1*(1+mu1*(drop1+1/size1))+mean2*(1+mu2*(drop2+1/size2))) )/mean3-1 )/mu3-drop3)
  res=list(c(mu1,size1,drop1),c(mu2,size2,drop2),c(mu3,size3,drop3))
  return(res)
}

n=100000
mean3_seq=c(2,3,4,5)
t_seq=c(0.1,0.2,0.5,0.8)
drop1=0.1
drop2=0.2
size1=10
size2=4

lm3=length(mean3_seq)
lt=length(t_seq)

res_mean_a1=matrix(,lm3,lt)
res_mean_a2=matrix(,lm3,lt)
res_var_a1=matrix(,lm3,lt)
res_var_a2=matrix(,lm3,lt)

res_mean_b1=matrix(,lm3,lt)
res_mean_b2=matrix(,lm3,lt)
res_var_b1=matrix(,lm3,lt)
res_var_b2=matrix(,lm3,lt)

op=par(mfrow = c(2, 1))
for(i_mean3 in 1:lm3){
  for(i_t in 1:lt){
    
    mean3=mean3_seq[i_mean3]
    t=t_seq[i_t]*mean3
    
    mean1=mean3+t
    mean2=mean3-t
    
    drop3=(drop1+drop2)/2
    
    mu1=mean1/(1-drop1)
    mu2=mean2/(1-drop2)
    mu3=mean3/(1-drop3)
    
    #var(a3)=t^2+1/2*(mean1*(1+mu1*(drop1+1/size1)) + mean2*(1+mu2*(drop2+1/size2)))
    size3=1/((( t*t+1/2*(mean1*(1+mu1*(drop1+1/size1)) + mean2*(1+mu2*(drop2+1/size2))) )/mean3-1 )/mu3-drop3)
    
    
    a1=emdbook::rzinbinom(n,mu=mu1, size=size1,zprob=drop1)
    a2=emdbook::rzinbinom(n,mu=mu2, size=size2,zprob=drop2)
    a3=emdbook::rzinbinom(2*n,mu=mu3, size=size3,zprob=drop3)
    
    mean(a1)
    mean(a2)
    mean(a3)
    res_mean_a1[i_mean3,i_t]=mean(c(a1,a2))
    res_mean_a2[i_mean3,i_t]=mean(a3)
    var(a1)
    var(a2)
    var(a3)
    res_var_a1[i_mean3,i_t]=var(c(a1,a2))
    res_var_a2[i_mean3,i_t]=var(a3)
    hist(a3,breaks=50,range=c(0,max(a1,a2,a3)))
    hist(c(a1,a2),col=rgb(1,0,0,0.3),breaks=50,range=c(0,max(a1,a2,a3)))
    
    
    print(c(mu1,mu2,mu3,size1,size2,size3,drop1,drop2,drop3))
    param=cal_nbzinb_param_multimodality(mean1=mean1,mean2=mean2,size1=size1,size2=size2,drop1=drop1,drop2=drop2)
    print(param)
    
    b1=emdbook::rzinbinom(n,mu=param[[1]][1], size=param[[1]][2],zprob=param[[1]][3])
    b2=emdbook::rzinbinom(n,mu=param[[2]][1], size=param[[2]][2],zprob=param[[2]][3])
    b3=emdbook::rzinbinom(2*n,mu=param[[3]][1], size=param[[3]][2],zprob=param[[3]][3])
    
    
    
    
    mean(b1)
    mean(b2)
    mean(b3)
    res_mean_b1[i_mean3,i_t]=mean(c(b1,b2))
    res_mean_b2[i_mean3,i_t]=mean(b3)
    var(b1)
    var(b2)
    var(b3)
    res_var_b1[i_mean3,i_t]=var(c(b1,b2))
    res_var_b2[i_mean3,i_t]=var(b3)
    hist(b3,breaks=50,range=c(0,max(b1,b2,b3)))
    hist(c(b1,b2),col=rgb(1,0,0,0.3),breaks=50,range=c(0,max(b1,b2,b3)))
    
  }
}
par(op)

plot(res_mean_a1,res_mean_a2)
lines(c(0,20),c(0,20))
plot(res_var_a1,res_var_a2)
lines(c(0,20),c(0,20))

plot(res_mean_b1,res_mean_b2)
lines(c(0,20),c(0,20))
plot(res_var_b1,res_var_b2)
lines(c(0,20),c(0,20))








