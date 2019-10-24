
#this code provide zinb or nb fit functions.
#this code is required by 7 series for fitting.

library("pscl")
library("MASS")

###########functions#############
#fit_zinb fits the zinb model(if 0 included in data) and nb model (if 0 is not included)

#input:
#input_numeric: the numerical vectors as observation
#input_x: the covariates.If only need to estimate the observation distribution, input_x=1(default)

#output:
#fit_zinb returns 3 numbers, the fitted log transformed mean, fitted dispersion(theta)

#package requirement:
#fit_zinb needs package pscl::zeroinfl


fit_nbzinb=function(input_numeric,input_x=NA){
  if((max(input_numeric,na.rm = TRUE)-min(input_numeric,na.rm = TRUE)==0)){
    warning("all observation shows the same, returns the log-transformed values themselves or 0 as mean")
    if(max(input_numeric,na.rm = TRUE)==0){
      fit_total=c(0,NA,1)
    }
    if(max(input_numeric,na.rm = TRUE)>0){
      fit_total=c(log(max(input_numeric,na.rm = TRUE)),NA,0)
    }
  }
  else{
    if(min(input_numeric)==0){
      fm_zinb=NA
      if(sum(!is.na(input_x))==0){
        fm_zinb=tryCatch(pscl::zeroinfl(as.numeric(input_numeric) ~  1, dist = "negbin"), error = function(e) {NA} )
        #print("zinb0")
      }
      if(sum(!is.na(input_x))>0){
        fm_zinb=tryCatch(pscl::zeroinfl(as.numeric(input_numeric) ~  input_x, dist = "negbin"), error = function(e) {print(e);NA} )
        #print("zinb1")
      }
      if(sum(!is.na(fm_zinb))>0){
        fit_logitdropout=as.numeric(fm_zinb$coefficients$zero[1])
        fit_dropout=exp(fit_logitdropout)
        fit_dropout=fit_dropout/(1+fit_dropout)
        fit_logmean=as.numeric(fm_zinb$coefficients$count[1])
        fit_dispersion=as.numeric(fm_zinb$theta)
        fit_total=c(fit_logmean,fit_dispersion,fit_dropout)
      }
      else{
        fit_total=c(NA,NA,NA)
      }
    }
    else{
      fm_nb=NA
      if(sum(!is.na(input_x))==0){
        fm_nb=tryCatch(MASS::glm.nb(as.numeric(input_numeric) ~  1), error = function(e) {print(e);NA} )
        #print("nb0")
      }
      if(sum(!is.na(input_x))>0){
        fm_nb=tryCatch(MASS::glm.nb(as.numeric(input_numeric) ~  input_x), error = function(e) {print(e);NA} )
        #print("nb1")
      }
      if(sum(!is.na(fm_nb))>0){
        fit_dropout=0
        fit_logmean=as.numeric(fm_nb$coefficients[1])
        fit_dispersion=as.numeric(fm_nb$theta)
        fit_total=c(fit_logmean,fit_dispersion,fit_dropout)
      }
      else{
        fit_total=c(NA,NA,NA)
      }
      
    }
  }
  names(fit_total)=c("logmean","dispersion","dropout_rate")
  return(fit_total)
}

