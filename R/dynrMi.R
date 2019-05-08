# convergence diagnostics
diag.mi = function(imp, nvariables, variables, m,itermin=2,iter,burn=0,Rhat){ 
  
  chains = m
  
  # reformat imp$chainMean to a matrix called coda
  chainmean = imp$chainMean
  chainmean2 = chainmean[!is.na(chainmean)]
  coda = matrix(chainmean2, nrow = chains*iter, byrow = T)
  
  # create a null matrix to store Rhats for all variables
  Rhatmatrix = matrix(NA, nrow = iter, ncol = nvariables) 
  
  #cal Rhats for each iteration and each variable
  for(i in itermin:iter){
    
    value = list()  #each list contains a chain 
    for(k in 1:chains){
      codak = coda[(iter*(k-1)+1):(iter*k),1:nvariables]  # select all values from each chain
      value[[k]] = codak[(1+burn):i,]  #drop the values of each chain in the beginning period
    }
    
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
      # globalvar: within-chain variances
      globalvar = sum(chainvars) / chains
      
      # b: between-chains variances
      b = sum((chainmeans - globalmean)^2) * iterations / (chains - 1)
      
      # varplus: pooled variance
      varplus = (iterations - 1) * globalvar / iterations + b / iterations
      
      # Gelman-Rubin statistic (Rhat)
      rhat = sqrt(varplus / globalvar)
      Rhatmatrix[i,j] = rhat
    } # close loop over variables
  } # close loop over iterations

  # trace plots from mice()
  p1 = plot(imp, variables)  
  
  # Rhat plot
  # prepare a long format data set for ggplot
  df = data.frame(matrix(ncol = 3, nrow = iter*nvariables))
  colnames(df) = c("variable","iteration","Rhatvalue")
  for(j in 1:nvariables){
    df$variable[(iter*(j-1)+1):(iter*j)] = rep(variables[j],iter)
    df$iteration[(iter*(j-1)+1):(iter*j)] = 1:iter
    df$Rhatvalue[(iter*(j-1)+1):(iter*j)] = Rhatmatrix[,j]
  }
  p2 = ggplot(subset(df, df$iteration!=1), aes(iteration,Rhatvalue)) +
    geom_line(aes(col = variable))+
    geom_hline(yintercept=Rhat,size=1) +
    labs(x="iteration", y="Rhat")+
    theme_classic()
  
  return(list(p1,p2))
}



##' Multiple Imputation of dynrModel objects
##' 
##' @param dynrModel dynrModel object. data and model setup
##' @param which.aux character. names of the auxiliary variables used in the imputation model
##' @param which.lag character. names of the variables to create lagged responses for imputation purposes
##' @param lag integer. number of lags of variables in the imputation model
##' @param which.lead character. names of the variables to create leading responses for imputation purposes
##' @param lead integer. number of leads of variables in the imputation model
##' @param m integer. number of multiple imputations
##' @param iter integer. number of MCMC iterations in each imputation
##' @param imp.obs logical. flag to impute the observed dependent variables
##' @param imp.exo logical. flag to impute the exogenous variables
##' @param diag logical. flag to use convergence diagnostics
##' @param Rhat numeric. value of the Rhat statistic used as the criterion in convergence diagnostics
##' @param conf.level numeric. confidence level used to generate confidence intervals
##' @param verbose logical. flag to print the intermediate output during the estimation process
##' @param seed integer. random number seed to be used in the MI procedure
##' 

dynr.mi <- function(dynrModel, which.aux=NULL, 
                    which.lag=NULL, lag=0,
                    which.lead=NULL, lead=0,
                    m=5, iter=5, 
                    imp.obs=FALSE, imp.exo=TRUE,
                    diag = TRUE, Rhat=1.1,
                    conf.level=0.95,
                    verbose=TRUE, seed=NA){    
	
	data <- dynrModel$data$original.data
	k <- length(dynrModel$param.names)    # number of parameters estimated
	
	
	ynames <- dynrModel$data$observed.names
	xnames <- dynrModel$data$covariate.names
	y <- data[, ynames]   # observed variables
	x <- data[, xnames]   # covariates
	au <- data[, which.aux]  # auxiliary variables
	ID <- dynrModel$data$id
	id <- unique(ID)   # a vector of IDs
	time <- dynrModel$data$time
	
	# raw data 
	datanolag <- cbind(ID,y,x,au) 
	
	# select variables to create lags on
	dataforlag <- subset(datanolag,select = c("ID", which.lag))  
	
	# create a null data frame to store lagged variables 
	datalag <- data.frame(ID, matrix(NA, nrow = length(ID), ncol = length(which.lag)*lag))
	
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
	dataformice <- data.frame(subset(datanolag,select=-ID), subset(datalag, select = -ID))

	# implement imputations
	imp <- mice::mice(dataformice, m=m, maxit = iter, seed = seed)
	
	# pool estimation results
	print("Cooking the dynr model to estimate free parameters")
	
	pmcarqhat <- matrix(NA, nrow=m, ncol=k) #parameter estimates from each imputation
	pmcaru <- array(NA, dim=c(k,k,m)) #vcov of par estimates from each imputation
	
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

		colnames(newdata) <- c("ID", "Time", ynames,xnames)
		
		data <- dynr.data(newdata, id="ID", time="Time",
		                observed=ynames, covariates=xnames)
		
		modelnew <- dynrModel
		modelnew@data <- data
		modelnew@compileLib = FALSE #avoid conflicts among concurrent output C scripts of model functions

		
		trial <- dynr.cook(modelnew, verbose = verbose)  
		
		
		#getting parameter estimates
		pmcarqhat[j,] <- coef(trial)[1:k]
		pmcaru[, ,j] <- vcov(trial)[c(1:k),c(1:k)]
	}
	
	pqbarmcarimpute <- apply(pmcarqhat, 2, mean) 
	pubarmcarimpute <- apply(pmcaru, 1:2, mean)

	pe.mcarimpute <- pmcarqhat - matrix(pqbarmcarimpute, nrow = m, ncol = k, byrow = TRUE)
	pb.mcarimpute <- (t(pe.mcarimpute) %*% pe.mcarimpute)/(m - 1)
	pvcovmcarimpute <- pubarmcarimpute + (1 + 1/m) * pb.mcarimpute #vcov for estimates
	psemcarimpute <- sqrt(diag(pvcovmcarimpute))
	
	t <- pqbarmcarimpute/psemcarimpute # t statistic
	df <- nobs(model)-k    # degree of freedom
	alpha <- 1-conf.level  # significance level

	ci.upper <- pqbarmcarimpute + qt(1-alpha/2,df)*psemcarimpute # upper CI
	ci.lower <- pqbarmcarimpute - qt(1-alpha/2,df)*psemcarimpute # lower CI
	p <- pt( abs(t), df, lower.tail=FALSE) # p value
	
	result <- cbind(pqbarmcarimpute, psemcarimpute, t, ci.lower, ci.upper,p)
	
	colnames(result) <- c("Estimate", "Std. Error", "t value", "ci.lower", "ci.upper",
	                     "Pr(>|t|)")
	row.names(result) <- names(coef(model))
	
	if(diag == TRUE){
	  # obtain diagnostic information from diag.mi()
	  nvariables = length(c(ynames,xnames))
	  variables =c(ynames,xnames)
	  diagplot = diag.mi(imp, nvariables, variables, m, iter, Rhat)
	  x11(); dev.off()  # avoid plot rendering errors
 	  return(list(diagplot,result))
	}
	else
	  return(result)
	
}
