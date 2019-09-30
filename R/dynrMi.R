##' Multiple Imputation of dynrModel objects
##' 
##' @param model dynrModel object
##' @param m number of multiple imputations
##' @param aux.variable names of auxiliary variables used in imputation
##' @param imp.obs logical. whether to impute the observed variables
##' @param imp.exo logical. whether to impute the exogenous variables
##' @param lag numeric. the number of lags to use
##' 
##' @details
##' This function is in alpha-testing form.  Please do not use or rely on it for now. A full implementation is in progress.
dynr.mi <- function(model, m=5, aux.variable, imp.obs=FALSE, imp.exo=FALSE, lag){    #multiple lag; #factor  #get variable names
	
	
	data=model@data$original.data
	k=length(model@param.names) #- length(model@initial@paramnames)             #number of parameters estimated
	nolag=TRUE
	
	
	ynames=model@data$observed.names
	xnames=model@data$covariate.names
	y=data[,colnames(data)==ynames]
	x=data[,colnames(data)==xnames]
	ID=model@data$id
	time=model@data$time
	
	au=data[,colnames(data)==aux.variable]  
	
	dataforlag=cbind(ID,y,x)
	
	#dataforlag=ts(dataforlag)
	
	#head(stats::lag(as.xts(dataforlag),2))
	#head(lag(dataforlag,2))
	
	datalag <- dataforlag
	#  dataforlag %>%
	#  dplyr::group_by(ID) %>%
	#  dplyr::mutate_all(lag) #Linying, what does this "lag" do?
	
	dataformice=cbind(dataforlag[,-1],datalag[,-1],au)
	dataformice=data.frame(dataformice)

	colnames(dataformice)=c()
	m=m
	imp=mice::mice(dataformice,m=m)
	
	
	pmcarqhat=matrix(NA, nrow=m,ncol=k) #parameter estimates from each imputation
	pmcaru=array(NA,dim=c(k,k,m)) #vcov of par estimates from each imputation
	
	for (j in 1:m){
	  
	  completedata=mice::complete(imp,action=j)
	  colnames(completedata)=c(ynames,xnames,
	                          paste("lag",ynames,sep=''),
	                          paste("lag",xnames,sep=''),
	                          colnames(au))
	  
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
	  
	  colnames(newdata)=c("ID","Time","wp","hp","ca","cn")
	  
	  data <- dynr.data(newdata, id="ID", time="Time", 
	                    observed=c("wp","hp"),covariates=c("ca","cn"))
	  
	 modelnew=model
	 modelnew@data=data
	 #   model <- dynr.model(dynamics=model@dynamics, measurement=model@measurement,
	 #                      noise=model@noise, initial=model@initial, data=data, transform=model@transform,
	 #                     outfile=paste("trial4",i,".c",sep=""))
	  
	   
	   # model <- dynr.model(dynamics=dynm, measurement=meas,
	   #                     noise=mdcov, initial=initial, data=data,#transform=trans,
	   #                     outfile=paste("trial2",i,".c",sep=""))
	   
	   
	  trial <- dynr.cook(modelnew)  #names(trial) get names of the params
	  #summary(trial)
	  
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
