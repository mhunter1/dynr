##' Multiple Imputation of dynrModel objects
##' 
##' @param model dynrModel object
##' @param aux.variable names of auxiliary variables used in imputation 
##' @param m number of multiple imputations
##' @param iter number of iterations in one imputation
##' @param imp.obs logical. whether to impute the observed variables
##' @param imp.exo logical. whether to impute the exogenous variables
##' @param lag numeric. the number of lags to use
##' @param lag.variable names of variables to create lags on
##' @param leads logical. whether to use lags or leads
##' @param diag logical. whether to use convergence diagnostics
##' @param cook.save logical. whether to save dynr.cook object
##' @param seed integer. a single value used to set seed in imputation
##' 
##' @details
##' This function is in alpha-testing form.  Please do not use or rely on it for now. A full implementation is in progress.
dynr.mi <- function(model, aux.variable, m=5, iter, imp.obs=FALSE, imp.exo=FALSE, lag, lag.variable, 
                    leads = FALSE, diag = TRUE, cook.save = FALSE, seed = NA){    #multiple lag; #factor  #get variable names
	
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
	
	# raw data 
	datanolag <- cbind(ID,y,x,au) 
	
	# select variables to create lags on
	dataforlag <- subset(datanolag,select = c("ID", lag.variable))  
	
	# create a null data frame to store lagged variables 
	datalag <- data.frame(ID, matrix(NA, nrow = length(ID), ncol = length(lag.variable)*lag))
	
	if(lag < 1){
	  warning("No lags/leads introduced.")
	  datalag <- dataforlag
	}
	
	for(i in id){
	  # tmp used to store original variables for each subject 
	  tmp <- dataforlag[dataforlag$ID == i,]
	  tmp <- as.matrix(subset(tmp, select = -ID))
	  P <- ncol(tmp)
	  nt <- nrow(tmp)
	  if (lag > nt-1){
	    lag <- nt-1
	    warning("The number of lags/leads should be smaller than the number of measurements.")
	  }
	  
	  for(t in 1:lag){
	    # tmp1 used to store lagged variables for each subject
	    tmp1 <- matrix(NA,nrow = nt, ncol = P)
	    if (leads == TRUE)  { tmp1[1:(nt-t),] <- tmp[(t+1):nt,] }
	    else  { tmp1[(t+1):nt,] <- tmp[1:(nt-t),] }
	    a <- P*(t-1)+2 
	    b <- P*t+1 
	    datalag[datalag$ID == i, ][, a:b] <- tmp1
	  }
	}

	# combine original variables and lagged variables 
	dataformice <- data.frame(datanolag, subset(datalag, select = -ID))

	
	imp <- mice::mice(dataformice, m=m, maxit = iter, seed = seed, printFlag = FALSE)
	
	
	# convergence diagnostics
	diag.mi = function(imp, nvariables, m,itermin,iter,burn){ #number of iterations should be more than itermin
	  
	  chains = m
	  
	  # reformat imp$chainMean to a matrix called coda
	  chainmean = imp$chainMean
	  chainmean2 = chainmean[!is.na(chainmean)]
	  coda = matrix(chainmean2, nrow = m*iter, byrow = T)
	  
	  # create a null matrix to store Rhats for all variables
	  Rhat = matrix(NA, nrow = iter, ncol = nvariables) 
	  
	  #cal Rhats for each iteration and each variable
	  for(i in itermin:iter){
	    
	    # select all values generated from i iterations
	    codai = coda[1:m*i,1:nvariables]
	    value = list()  #each list contains a chain 
	    for(k in 1:m){
	      value[[k]] = coda[(i*(k-1)+1+burn):(i*k),]  #drop the values of each chain in the beginning period
	    }
	    
	    # cal Rhats
	    iterations = i
	    for(j in 1:nvariables) {
	      
	      # create vectors to store means and variations of values in chains for each variable
	      chainmeans = rep(NA, chains)
	      chainvars = rep(NA, chains)
	      
	      for(k in 1:chains) {
	        sum = sum(value[[k]][,j])
	        var = var(value[[k]][,j])
	        mean = sum / iterations
	        
	        chainmeans[k] = mean
	        chainvars[k] = var
	      }
	      globalmean = sum(chainmeans) / chains
	      globalvar = sum(chainvars) / chains
	      
	      # Compute between- and within-variances and MPV
	      b = sum((chainmeans - globalmean)^2) * iterations / (chains - 1)
	      
	      varplus = (iterations - 1) * globalvar / iterations + b / iterations
	      
	      # Gelman-Rubin statistic (Rhat)
	      rhat = sqrt(varplus / globalvar)
	      Rhat[i,j] = rhat
	    }
	  }
	  return(Rhat)
	}
	
	
	if (diag == TRUE){
	  # trace plots from mice()
	  plot(imp, c(ynames,xnames)) 
	  
    # Rhat plots from diag.mi()
	  nvariables = length(c(ynames,xnames))
	  result = diag.mi(imp, nvariables, m, 2,iter,0)
	  
	  names =c(ynames,xnames)
	  for(j in 1:nvariables){
	    plot(2:iter, result[2:iter,j],type="l",ylim=c(min(na.omit(result)),max(na.omit(result))),ylab="Rhat",xlab="last iteration in chain",main=names[j])
	    abline(h=1.1,lty=2,col=2)  # criterion set as Rhat < 1.1
	  }
	}
	
	pmcarqhat <- matrix(NA, nrow=m, ncol=k) #parameter estimates from each imputation
	pmcaru <- array(NA, dim=c(k,k,m)) #vcov of par estimates from each imputation
	
	print("Cooking the dynr model to estimate free parameters")
	
	for(j in 1:m){
		
		completedata <- mice::complete(imp, action=j) #obtain the jth imputation
		
		if(imp.obs==TRUE){
			imp.data.obs <- completedata[, ynames]
		} else{
			imp.data.obs <- y
		}
		
		if(imp.exo==TRUE){
			imp.data.exo <- completedata[, xnames]
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

		
		trial <- dynr.cook(modelnew, verbose = FALSE)  #names(trial) get names of the params
		#summary(trial)
		
		# save dynr.cook results
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
