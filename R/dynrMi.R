##' Multiple Imputation of dynrModel objects
##' 
##' @param model dynrModel object
##' @param m number of multiple imputations
##' @param aux.variable names of auxiliary variables used in imputation
##' @param imp.obs logical. whether to impute the observed variables
##' @param imp.exo logical. whether to impute the exogenous variables
##' @param lag numeric. the number of lags to use
##' @param lag.variable names of variables to create lags on
##' @param leads logical. whether to use lags or leads
##' @param cook.save logical. whether to save dynr.cook object
##' 
##' @details
##' This function is in alpha-testing form.  Please do not use or rely on it for now. A full implementation is in progress.
dynr.mi <- function(model, m=5, aux.variable, imp.obs=FALSE, imp.exo=FALSE, lag, lag.variable, 
                    leads = FALSE, cook.save = FALSE){    #multiple lag; #factor  #get variable names
	
	data <- model$data$original.data
	k <- length(model$param.names)    # number of parameters estimated
	
	
	ynames <- model$data$observed.names
	xnames <- model$data$covariate.names
	y <- data[, colnames(data)==ynames]
	x <- data[, colnames(data)==xnames]
	au <- data[, colnames(data)==aux.variable] 
	ID <- model$data$id
	id <- unique(ID)   # number of subjects
	time <- model$data$time
	
	# data from dynr.model
	datanolag <- cbind(ID,y,x,au) 
	
	# variables to create lags on
	dataforlag <- subset(datanolag,select = c("ID", lag.variable))  
	
	# create a null data frame to store lag variables 
	datalag <- data.frame(ID, matrix(NA, nrow = length(ID), ncol = length(lag.variable)*lag))
	
	if(lag < 1){
	  warning("No lags/leads introduced.")
	  datalag <- dataforlag
	}
	
	for(i in id){
	  tmp <- dataforlag[dataforlag$ID == i,]
	  tmp <- as.matrix(subset(tmp, select = -ID))
	  P <- ncol(tmp)
	  nt <- nrow(tmp)
	  if (lag > nt-1){
	    lag <- nt-1
	    warning("The number of lags/leads should be smaller than the number of measurements.")
	  }
	  tmp1 <- matrix(NA,nrow = nt, ncol = P)
	  for(t in 1:lag){
	    if (leads == TRUE)  { tmp1[1:(nt-t),] <- tmp[(t+1):nt,] }
	    else  { tmp1[(t+1):nt,] <- tmp[1:(nt-t),] }
	    a <- P*(t-1)+2
	    b <- P*t+1
	    datalag[datalag$ID == i, ][, a:b] <- tmp1
	  }
	}

	# combine original variables and lag variables 
	dataformice <- data.frame(datanolag, subset(datalag, select = -ID))

	
	imp <- mice::mice(dataformice, m=m)
	
	
	pmcarqhat <- matrix(NA, nrow=m, ncol=k) #parameter estimates from each imputation
	pmcaru <- array(NA, dim=c(k,k,m)) #vcov of par estimates from each imputation
	
	for(j in 1:m){
		
		completedata <- mice::complete(imp, action=j) #obtain the jth imputation
		
		if(imp.obs==TRUE){
			imp.data.obs <- completedata[, colnames(completedata)==ynames]
		} else{
			imp.data.obs <- y
		}
		
		if(imp.exo==TRUE){
			imp.data.exo <- completedata[, colnames(completedata)==xnames]
		} else{
			imp.data.exo <- x
		}
		
		newdata <- cbind(ID, time, imp.data.obs, imp.data.exo)
		
		
		save(newdata,file="test.rdata")

		
		colnames(newdata) <- c("ID", "Time", ynames,xnames)
		
		data <- dynr.data(newdata, id="ID", time="Time",
		                observed=ynames, covariates=xnames)
		
		modelnew <- model
		modelnew@data <- data

		
		trial <- dynr.cook(modelnew)  #names(trial) get names of the params
		#summary(trial)
		
		if (cook.save == TRUE)
		  save(trial, file=paste0("cookresult",j,".rdata"))
		
		#getting parameter estimates
		pmcarqhat[j,] <- coef(trial)[1:k]
		pmcaru[, ,j] <- vcov(trial)[c(1:k),c(1:k)]
	}
	
	pqbarmcarimpute <- apply(pmcarqhat, 2, mean) 
	pubarmcarimpute <- apply(pmcaru, 1:2, mean)
	#ubar <- apply(u, c(2, 3), mean)
	pe.mcarimpute <- pmcarqhat - matrix(pqbarmcarimpute, nrow = m, ncol = k, byrow = TRUE)
	pb.mcarimpute <- (t(pe.mcarimpute) %*% pe.mcarimpute)/(m - 1)
	pvcovmcarimpute <- pubarmcarimpute + (1 + 1/m) * pb.mcarimpute #vcov for estimates
	psemcarimpute <- sqrt(diag(pvcovmcarimpute))
	
	t <- pqbarmcarimpute/psemcarimpute
	# TODO Don't just use 2 for the CIs !!!
	ci.upper <- pqbarmcarimpute + 2*psemcarimpute
	ci.lower <- pqbarmcarimpute - 2*psemcarimpute
	p <- pt( abs(t), df=nobs(model) - k, lower.tail=FALSE)
	
	result <- cbind(pqbarmcarimpute, psemcarimpute, t, ci.lower, ci.upper,p)
	
	colnames(result) <- c("Estimate", "Std. Error", "t value", "ci.lower", "ci.upper", #"",
	                     "Pr(>|t|)")
	row.names(result) <- names(coef(model))
	
	#TODO the current 'result' should be what is printed by summary() on the thing returned from dynr.mi()
	# This summary should use methods similar to summary.dynrCook and summary.lm
	# TODO return a 'mi' object similar to rest of mice
	return(result)
}
