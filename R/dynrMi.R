##' Multiple Imputation of dynrModel objects
##' 
##' @param model dynrModel object
##' @param aux.variable names of auxiliary variables used in imputation
##' @param m number of multiple imputations
##' @param imp.obs logical. whether to impute the observed variables
##' @param imp.exo logical. whether to impute the exogenous variables
##' @param imp.method names of imputation methods
##' @param lag numeric. the number of lags to use
##' @param lag.variable names of variables to create lags on
##' @param leads logical. whether to use lags or leads
##' @param cook.save logical. whether to save the dynr.cook objects
##' 
##' @details
##' This function is in alpha-testing form.  Please do not use or rely on it for now. A full implementation is in progress.
dynr.mi <- function(model, aux.variable, m=5, imp.obs=FALSE, imp.exo=FALSE, imp.method = "mice", 
                    lag, lag.variable, leads = FALSE, cook.save = FALSE){    #multiple lag; #factor  #get variable names
  
  
  data=model@data$original.data
  k=length(model@param.names)             #number of parameters estimated
  
  
  ynames=model@data$observed.names
  xnames=model@data$covariate.names
  y=data[,colnames(data)==ynames]
  x=data[,colnames(data)==xnames]
  ID=model@data$id
  n=length(unique(ID))
  time=model@data$time
  
  
  au=data[,colnames(data)==aux.variable]  
  
  datanolag=cbind(ID,y,x)
  dataforlag = subset(datanolag,select = c("ID", lag.variable))
  
  if (lag < 1){
    warning("No lags/leads introduced.")
    datalag = dataforlag
  }
  
  datalag = NULL
  for(i in 1:n){
    tmp = dataforlag[dataforlag$ID == i,]
    tmp = as.matrix(tmp[,-1])
    P = ncol(tmp)
    nt = nrow(tmp)
    if (lag > nt-1){
      lag = nt-1
      warning("The number of lags/leads should be smaller than the number of measurements.")}
    datalag1 = NULL
    for(t in 1:lag){    #number of lags recommendation
      tmp1 = matrix(NA,nrow = nt, ncol = P)
      if (leads == TRUE){
        tmp[1:(nt-t),] = tmp1[(t+1):nt,]
        colnames(tmp1)=paste0(lag.variable,"lead",t)
      }
      else{
        tmp1[(t+1):nt,] = tmp[1:(nt-t),]
        colnames(tmp1)=paste0(lag.variable,"lag",t)
      }
      datalag1 = cbind(datalag1, tmp1)
    }
    datalag = rbind(datalag, datalag1)
  }
  datalag = as.data.frame(datalag)
  
  
  
  
  dataformice=cbind(datanolag[,-1],datalag,au)
  dataformice=data.frame(dataformice)
  
  imp=mice::mice(dataformice,m) #only imputate the initial data
  #fit = with(imp)
  
  pmcarqhat=matrix(NA, nrow=m,ncol=k) #parameter estimates from each imputation
  pmcaru=array(NA,dim=c(k,k,m)) #vcov of par estimates from each imputation
  
  for (j in 1:m){
    
    completedata=mice::complete(imp,action=j) #obtain the jth imputation
    #colnames(completedata)=c(ynames,xnames,
    #                       paste("lag",ynames,sep=''),
    #                       paste("lag",xnames,sep=''),
    #                       colnames(au))
    
    if (imp.obs==TRUE){
      imp.data.obs= completedata[,colnames(completedata)==ynames]  
    }else{
      imp.data.obs=y
    }
    
    if (imp.exo==TRUE){
      imp.data.exo= completedata[,colnames(completedata)==xnames]  
    }else{
      imp.data.exo=x
    }
    
    newdata=cbind(ID,time,imp.data.obs,imp.data.exo)
    
    save(newdata,file="test.rdata")
    
    colnames(newdata)=c("ID","Time",ynames,xnames)
    
    data <- dynr.data(newdata, id="ID", time="Time", 
                      observed=ynames,covariates=xnames)
    
    modelnew=model
    modelnew@data=data
    
    
    trial <- dynr.cook(modelnew)  #names(trial) get names of the params
    #summary(trial)
    if (cook.save == TRUE)
      save(trial, file=paste0("cookresult",j,".rdata"))
    
    #getting parameter estimates
    pmcarqhat[j,]=coef(trial)[1:k]
    pmcaru[, ,j]= vcov(trial)[c(1:k),c(1:k)]
  }
  
  pqbarmcarimpute <- apply(pmcarqhat, 2, mean) 
  pubarmcarimpute <- apply(pmcaru, 1:2, mean)
  #ubar <- apply(u, c(2, 3), mean)
  pe.mcarimpute <- pmcarqhat - matrix(pqbarmcarimpute, nrow = m, ncol = k, byrow = TRUE)
  pb.mcarimpute <- (t(pe.mcarimpute) %*% pe.mcarimpute)/(m - 1)
  pvcovmcarimpute <- pubarmcarimpute + (1 + 1/m) * pb.mcarimpute #vcov for estimates
  psemcarimpute=sqrt(diag(pvcovmcarimpute))
  
  t=pqbarmcarimpute/psemcarimpute
  ci.upper=pqbarmcarimpute+2*psemcarimpute
  ci.lower=pqbarmcarimpute-2*psemcarimpute
  p <- pt( abs(t), df=nobs(model) - k, lower.tail=FALSE)
  
  result=cbind(pqbarmcarimpute,psemcarimpute,t,ci.lower,ci.upper,p)
  
  colnames(result) = c("Estimate", "Std. Error", "t value", "ci.lower", "ci.upper", #"",
                       "Pr(>|t|)")
  row.names(result) <- names(coef(model))
  
  
  return(result)
}
