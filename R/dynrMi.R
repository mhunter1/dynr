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
##' @return an object of `dynrMi' class
##' that is a list containing:
##' 1. the imputation information, including a data set 
##' containing structured lagged and leading variables and 
##' a `mids' object from mice() function;
##' 2. the diagnostic information, including trace plots,
##' an Rhat plot and a matrix containing Rhat values;
##' 3. the estimation results, including parameter estimates,
##' standard error estimates and confidence intervals.
##' 
##' @references
##' Ji, L., Chow, S-M., Schermerhorn, A.C., Jacobson, N.C., & Cummings, E.M. (2018). Handling 
##' Missing Data in the Modeling of Intensive Longitudinal Data. Structural Equation Modeling: 
##' A Multidisciplinary Journal, 1-22.
##' 
##' Yanling Li, Linying Ji, Zita Oravecz, Timothy R. Brick,
##' Michael D. Hunter, and Sy-Miin Chow. (2019).
##' dynr.mi: An R Program for Multiple Imputation in Dynamic Modeling.
##' International Journal of Computer, Electrical, Automation, Control
##' and Information Engineering, 13, 302-311.
##' 
##' @examples
##' # See the demo, MILinearDiscrete.R, for an illustrative example 
##' # of using dynr.mi to implement multiple imputation with 
##' # a vector autoregressive model
##' # dynrMi <- dynr.mi(dynrModel, which.aux=c("x1","x2"), 
##' # which.lag=c("wp","hp"), lag=1, which.lead=NULL, lead=0,
##' # m=5, iter=5, imp.obs=FALSE, imp.exo=TRUE,
##' # diag = TRUE, Rhat=1.1,
##' # conf.level=0.95, verbose=FALSE, seed=12345)

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
	y <- subset(data, select = ynames)   # observed variables
	x <- subset(data, select = xnames)   # covariates
	au <- subset(data, select = which.aux)  # auxiliary variables
	ID <- dynrModel$data$id
	id <- unique(ID)   # a vector of IDs
	time <- dynrModel$data$time
	
	# original data 
	datainit <- data.frame(ID,y,x,au) 
	
	# select variables to create lags on
	dataforlag <- subset(datainit,select = c("ID", which.lag))   
	
	# create a null data frame to store lagged variables 
	datalag <- data.frame(ID, matrix(NA, nrow = length(ID), ncol = length(which.lag)*lag))
	
	if(lag < 1){
	  warning("No lagged variables in the imputation model.")
	} else{
	  for(i in id){
	    # tmp used to store original variables for each subject 
	    tmp_lag <- dataforlag[dataforlag$ID == i, which.lag]
	    P_lag <- length(which.lag)  #number of variables to create lags on
	    nt <- nrow(tmp_lag) #number of time points
	    if (lag > nt-1){
	      lag <- nt-1
	      warning("Number of lags should be smaller than the number of time points.")
	    }
	    for(l in 1:lag){
	      # tmp1 used to store lagged variables for each subject
	      tmp1_lag <- data.frame(matrix(NA,nrow = nt, ncol = P_lag))
	      tmp1_lag[(l+1):nt,] <- tmp_lag[1:(nt-l),]
	      datalag[datalag$ID == i, ][, (P_lag*(l-1)+2):(P_lag*l+1)] <- tmp1_lag
	    }  #close loop over number of lags
	  } #close loop over subjects
	}

	# select variables to create leads on
	dataforlead <- subset(datainit,select = c("ID", which.lead))  
	
	# create a null data frame to store leading variables 
	datalead <- data.frame(ID, matrix(NA, nrow = length(ID), ncol = length(which.lead)*lead))
	
	if(lead > 0){
	  for(i in id){
	    # tmp used to store original variables for each subject 
	    tmp_lead <- dataforlead[dataforlead$ID == i, which.lead]
	    P_lead <- length(which.lead)  #number of variables to create lags on
	    nt <- nrow(tmp_lead) #number of time points
	    if (lead > nt-1){
	      lead <- nt-1
	      warning("Number of leads should be smaller than the number of time points.")
	    }
	    for(l in 1:lead){
	      # tmp1 used to store leading variables for each subject
	      tmp1_lead <- data.frame(matrix(NA,nrow = nt, ncol = P_lead))
	      tmp1_lead[1:(nt-l),] <- tmp_lead[(l+1):nt,]
	      datalead[datalead$ID == i, ][, (P_lead*(l-1)+2):(P_lead*l+1)] <- tmp1_lead
	    }  #close loop over number of leads
	  } #close loop over subjects
	}
	
	# combine original, lagged and leading variables 
	dataformice <- data.frame(subset(datainit,select=-ID), 
	                          subset(datalag, select = -ID),
	                          subset(datalead, select = -ID))

	print("Implementing imputations ...")
	
	imp <- mice::mice(dataformice, m=m, maxit = iter, seed = seed)
	
	
	# convergence diagnostics
	diag.mi = function(imp, nvariables,variables,m,itermin=2,iter,burn=0){ 
	  
	  chains = m
	  
	  # reformat imp$chainMean to a matrix called coda
	  chainmean = imp$chainMean[1:nvariables,,]
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
	  
	  Rhatmatrix = data.frame(Rhatmatrix)
	  colnames(Rhatmatrix) = variables
	  return(Rhatmatrix)
	  
	}
	
	pmcarqhat <- matrix(NA, nrow=m, ncol=k) #parameter estimates from each imputation
	pmcaru <- array(NA, dim=c(k,k,m)) #vcov of par estimates from each imputation
	
	print("Cooking and pooling of estimation results")
	
	for(j in 1:m){
		
		completedata <- mice::complete(imp, action=j) #obtain the jth imputation
		
		if(imp.obs==TRUE){
			imp.data.obs <- subset(completedata, select = ynames)
		} else{
			imp.data.obs <- y
		}
		
		if(imp.exo==TRUE){
			imp.data.exo <- subset(completedata, select = xnames)
		} else{
			imp.data.exo <- x
		}
		
		newdata <- cbind(ID, time, imp.data.obs, imp.data.exo)

		colnames(newdata) <- c("ID", "Time", ynames,xnames)
		
		if(length(xnames)==0){
		  data <- dynr.data(newdata, id="ID", time="Time",
		                    observed=ynames)
		}else{
		  data <- dynr.data(newdata, id="ID", time="Time",
		                    observed=ynames, covariates=xnames)
		}
		
		modelnew <- dynrModel
		modelnew@data <- data
		modelnew@compileLib = FALSE #avoid conflicts among concurrent output C scripts of model functions

		
		trial <- dynr.cook(modelnew, verbose = verbose)  
		
		
		# getting parameter estimates
		pmcarqhat[j,] <- coef(trial)[1:k]
		pmcaru[, ,j] <- vcov(trial)[c(1:k),c(1:k)]
	}
	
	pqbarmcarimpute <- apply(pmcarqhat, 2, mean) 
	pubarmcarimpute <- apply(pmcaru, 1:2, mean)

	pe.mcarimpute <- pmcarqhat - matrix(pqbarmcarimpute, nrow = m, ncol = k, byrow = TRUE)
	pb.mcarimpute <- (t(pe.mcarimpute) %*% pe.mcarimpute)/(m - 1)
	pvcovmcarimpute <- pubarmcarimpute + (1 + 1/m) * pb.mcarimpute #vcov for estimates
	psemcarimpute <- sqrt(diag(pvcovmcarimpute))
	
	t <- pqbarmcarimpute/psemcarimpute
	df <- nobs(dynrModel)-k
	alpha <- 1-conf.level

	ci.upper <- pqbarmcarimpute + qt(1-alpha/2,df)*psemcarimpute
	ci.lower <- pqbarmcarimpute - qt(1-alpha/2,df)*psemcarimpute
	p <- pt( abs(t), df, lower.tail=FALSE)
	
	result <- cbind(pqbarmcarimpute, psemcarimpute, t, ci.lower, ci.upper,p)
	
	colnames(result) <- c("Estimate", "Std. Error", "t value", "ci.lower", "ci.upper",
	                     "Pr(>|t|)")
	row.names(result) <- names(coef(dynrModel))
	
	
	if(diag == TRUE){
	  # obtain diagnostic information from diag.mi()
	  nvariables = length(c(ynames,xnames))
	  variables =c(ynames,xnames)
	  Rhatmatrix=diag.mi(imp, nvariables, variables,m=m, iter=iter)
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
	  p2 = ggplot(subset(df, df$iteration!=1), aes_string(x="iteration",y="Rhatvalue")) +
	    geom_line(aes_string(col = "variable"))+
	    geom_hline(yintercept=Rhat,size=1) +
	    labs(x="iteration", y="Rhat")+
	    theme_classic()
	  
	  #x11(); dev.off()  # avoid plot rendering errors
	  res = list(dataformice = dataformice,
	             imp = imp,
	             Rhat.matrix = Rhatmatrix,
	             trace.plot = p1, 
	             Rhat.plot = p2, 
	             parameters = pqbarmcarimpute,
	             standard.errors = psemcarimpute,
	             conf.intervals = cbind(ci.lower, ci.upper),
	             estimation.result = result)
	  class(res) <- "dynrMi"
	  invisible(res)
	}else{
	  res = list(dataformice = dataformice,
	             imp = imp,
	             parameters = pqbarmcarimpute,
	             standard.errors = psemcarimpute,
	             conf.intervals = cbind(ci.lower, ci.upper),
	             estimation.result = result)
	  class(res) <- "dynrMi"
	  invisible(res)
	}
	
}
