# Class definition for the dynrRun object that
#  stores all the output for a model that has 
#  been run.

# L = number of latent variables
# R = number of regimes
# T = number of time points
setClass(Class =  "dynrRun",
         representation = representation(
           fitted.parameters =  "numeric", #Can return
           transformed.parameters =  "numeric", #
           standard.errors =  "numeric",
           hessian =  "matrix",
           transformed.inv.hessian =  "matrix",
           conf.intervals = "matrix",
           exitflag = "numeric", #
           neg.log.likelihood = "numeric", #
           eta_regime_t = "array", # LxRxT #
           error_cov_regime_t  = "array", # LxLxRxT #
           innov_vec  = "array", # LxRxRxT #
           inverse_residual_cov = "array", # LxLxRxRxT #
           #Everything else from this point on
           pr_t_given_T  = "matrix", # RxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array", # LxLxT
           run.times = "numeric"
         )
)

# Initialize method for the new() function
setMethod("initialize", "dynrRun",
          function(.Object, x){
            .Object@fitted.parameters <- x$fitted.parameters
            .Object@transformed.parameters <- x$transformed.parameters
            .Object@standard.errors <- x$standard.errors
            .Object@hessian <- x$hessian.matrix
            .Object@transformed.inv.hessian <- x$transformed.inv.hessian
            .Object@conf.intervals <- x$conf.intervals
            .Object@exitflag <- x$exitflag
            .Object@neg.log.likelihood <- x$neg.log.likelihood
            .Object@eta_regime_t <- x$eta_regime_t
            .Object@error_cov_regime_t <- x$error_cov_regime_t
            .Object@innov_vec <- x$innov_vec
            .Object@inverse_residual_cov <- x$inverse_residual_cov
            .Object@pr_t_given_T <- x$pr_t_given_T
            .Object@eta_smooth_final <- x$eta_smooth_final
            .Object@error_cov_smooth_final <- x$error_cov_smooth_final
            return(.Object)
          }
)

setClass(Class =  "dynrDebug",
         representation = representation(
           fitted.parameters =  "numeric", #Can return
           transformed.parameters =  "numeric", #
           standard.errors =  "numeric",
           hessian =  "matrix",
           transformed.inv.hessian =  "matrix",
           conf.intervals = "matrix",
           exitflag = "numeric", #
           neg.log.likelihood = "numeric", #
           inverse.hessian = "matrix", 
           eta_regime_t = "array", # LxRxT #
           error_cov_regime_t  = "array", # LxLxRxT #
           eta_regime_regime_t_pred= "array", # LxRxRxT #
           error_cov_regime_regime_t_pred = "array", # LxLxRxRxT #
           eta_regime_regime_t_plus_1 = "array", # LxRxRxT #
           error_cov_regime_regime_t_plus_1 = "array", # LxLxRxRxT #
           innov_vec  = "array", # LxRxRxT #
           inverse_residual_cov = "array", # LxLxRxRxT #
           #Everything else from this point on
           pr_t_given_t   = "matrix", # RxT
           pr_t_given_t_less_1 = "matrix", # RxT
           pr_t_given_T  = "matrix", # RxT
           transprob_given_T  = "array", # RxRxT
           eta_regime_smooth = "array", # LxRxT
           error_cov_regime_smooth  = "array", # LxLxRxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array", # LxLxT
           run.times = "numeric"
         )
)

# Initialize method for the new() function
setMethod("initialize", "dynrDebug",
          function(.Object, x){
            .Object@fitted.parameters <- x$fitted.parameters
            .Object@transformed.parameters <- x$transformed.parameters
            .Object@standard.errors <- x$standard.errors
            .Object@hessian <- x$hessian.matrix
            .Object@transformed.inv.hessian <- x$transformed.inv.hessian
            .Object@conf.intervals <- x$conf.intervals
            .Object@exitflag <- x$exitflag
            .Object@neg.log.likelihood <- x$neg.log.likelihood
            .Object@inverse.hessian <- x$inverse.hessian.matrix
            .Object@eta_regime_t <- x$eta_regime_t
            .Object@error_cov_regime_t <- x$error_cov_regime_t
            .Object@eta_regime_regime_t_pred <- x$eta_regime_regime_t_pred
            .Object@error_cov_regime_regime_t_pred <- x$error_cov_regime_regime_t_pred
            .Object@eta_regime_regime_t_plus_1 <- x$eta_regime_regime_t_plus_1
            .Object@error_cov_regime_regime_t_plus_1 <- x$error_cov_regime_regime_t_plus_1
            .Object@innov_vec <- x$innov_vec
            .Object@inverse_residual_cov <- x$inverse_residual_cov
            .Object@pr_t_given_t <- x$pr_t_given_t
            .Object@pr_t_given_t_less_1 <- x$pr_t_given_t_less_1
            .Object@pr_t_given_T <- x$pr_t_given_T
            .Object@transprob_given_T <- x$transprob_given_T
            .Object@eta_regime_smooth <- x$eta_regime_smooth
            .Object@error_cov_regime_smooth <- x$error_cov_regime_smooth
            .Object@eta_smooth_final <- x$eta_smooth_final
            .Object@error_cov_smooth_final <- x$error_cov_smooth_final
            return(.Object)
          }
)

# Set the summary method of an object of class dynrRun
#  All this amounts to is writing a function that takes a
#  dynrRun object (and possibly other arguments)
#  and returns/does whatever we want.
# See Also the print method of summary.lm
#  getAnywhere(print.summary.lm)
summaryResults<-function(object){
             d <- data.frame(transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors)
             d$t_value<-ifelse(d$standard.errors==0, NA, d$transformed.parameters/d$standard.errors)
             d <-cbind(d,object@conf.intervals)
             neg2LL = -2*logLik(object)
             AIC = AIC(object)
             BIC = BIC(object)
             colnames(d) = c("Parameters","SE","t-value","CI.lower","CI.upper")
             cat("****************************SUMMARY******************************\n")
             print(d)
             cat(paste0("\n-2 log-likelihood value at convergence = ", round(neg2LL,2)))
             cat(paste0("\nAIC = ", round(AIC,2)))
             cat(paste0("\nBIC = ", round(BIC,2)))
             cat("\n*************************END OF SUMMARY**************************\n")
           }

setMethod( f = "summary",  signature = "dynrRun" ,
           definition = summaryResults)
setMethod( f = "summary",  signature = "dynrDebug" ,
		   definition = summaryResults)


displayDynrRun <- function(x){
  str(x)
  invisible(x)
}

setMethod("print", "dynrRun", function(x, ...) { 
  displayDynrRun(object) 
})

setMethod("show", "dynrRun", function(object) { 
  displayDynrRun(object) 
})


#------------------------------------------------------------------------------
# Some S3 methods
# coef
# logLik
# AIC
# BIC

coef.dynrRun <- function(object, ...){
  object@transformed.parameters
}


logLik.dynrRun <- function(object, ...){
  ans <- -object@neg.log.likelihood
  attr(ans, "df") <- length(object@fitted.parameters)
  attr(ans, "nobs") <- dim(object@eta_smooth_final)[2]
  class(ans) <- "logLik"
  return(ans)
}

# N.B. AIC() and BIC() are implicitly defined in terms
#  of logLik().

#------------------------------------------------------------------------------


dynr.run <- function(model, data,transformation, conf.level=.95, infile, verbose=TRUE,debug_flag=FALSE) {
	frontendStart <- Sys.time()
	if(missing(transformation)){
		transformation <- function(x){x}
	}
	model <- combineModelDataInformation(model, data)
	model <- preProcessModel(model)
	if(any(sapply(model$func_address, is.null.pointer))){
	    warning("Found null pointer(s) in 'func_address' list. (Re-)compiling your functions...")
	    if(missing(infile)){
	      stop("Cannot compile your functions because 'infile' argument is missing.")
	    }
	    model$func_address=dynr.funcaddress(isContinuousTime=model$isContinuousTime,infile=infile,verbose=verbose)
	}
	gc()
	backendStart <- Sys.time()
	output <- .Call("main_R", model, data, debug_flag,PACKAGE = "dynr")
	backendStop <- Sys.time()
	#gc()#garbage collection
	cat('Original exit flag: ', output$exitflag, '\n')
	output$exitflag <- output$exitflag+ifelse(output$exitflag<0,6,0)+ifelse(output$exitflag>0,5,0)
	cat('Modified exit flag: ', output$exitflag, '\n')
	cat(dynrExitFlags[output$exitflag], '\n')
	
	diagH = diag(output$hessian.matrix)
	diagH[diagH==0] = 10e-14
	diag(output$hessian.matrix) = diagH
	cat('Original fitted parameters: ', output$fitted.parameters, '\n', fill=TRUE)
	cat('Transformed fitted parameters: ', transformation(output$fitted.parameters), '\n', fill=TRUE)
	status = ifelse(any(!is.finite(output$hessian.matrix)) || !is.positive.definite(output$hessian.matrix), 0, 1)
	if (output$exitflag > 5 && status==1){
		output2 <- endProcessing(output, transformation, conf.level)
	}else{		
	  	output2 <- endProcessing2(output, transformation)
	  	cat('Hessian Matrix:',  '\n')
		print(output$hessian.matrix)
		cat('\n')
		#Print a message. Hessian matrix at convergence contains non-finite values or is
		#non-positive definite. ?
	}
	if (debug_flag){
		obj <- new("dynrDebug", output2)
	}else{
		obj <- new("dynrRun", output2)
	}
	
	frontendStop <- Sys.time()
	totalTime <- frontendStop-frontendStart
	backendTime <- backendStop-backendStart
	frontendTime <- totalTime-backendTime
	obj@run.times <- c(totalTime=totalTime, backendTime=backendTime, frontendTime=frontendTime)
	cat('Total Time:', totalTime, '\n')
	cat('Backend Time:', backendTime, '\n')
	rm(output2)
	rm(output)
	gc()
	return(obj)
}


endProcessing2 <- function(x, transformation){
  cat('Doing end processing for a failed trial\n')
  tParam <- transformation(x$fitted.parameters)
  x$transformed.parameters <- tParam
  
  nParam <- length(x$fitted.parameters)
  x$standard.errors <- rep(999,nParam)
  x$transformed.inv.hessian <-matrix(999, nrow=nParam,ncol=nParam)
  x$conf.intervals <-matrix(999, nrow=y,ncol=2)
  return(x)
}

endProcessing <- function(x, transformation, conf.level){
  #Sukruth to create another version of:
  #Do line 205, 207
	cat('Doing end processing\n')
	confx <- qnorm(1-(1-conf.level)/2)
	#Analytic Jacobian
	V1 = solve(x$hessian.matrix)
	J <- numDeriv::jacobian(func=transformation, x=x$fitted.parameters)
	iHess <- J %*% V1%*%t(J)
	tSE <- sqrt(diag(iHess))
	tParam <- transformation(x$fitted.parameters) #Can do
	CI <- c(tParam - tSE*confx, tParam + tSE*confx)
	x$transformed.parameters <- tParam #Can do
	x$standard.errors <- tSE
	x$transformed.inv.hessian <- iHess
	x$conf.intervals <- matrix(CI, ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
	return(x)
}

preProcessModel <- function(x){
	# Make sure that starting values are stored as double in R
	# so that C doesn't think they are integers (e.g. 1 vs 1.0)
	x$xstart <- as.double(x$xstart)
	x$ub <- as.double(x$ub)
	x$lb <- as.double(x$lb)
	return(x)
}

combineModelDataInformation <- function(model, data){
	# TODO add argument to dynrModel a la "usevars"
	# process usevars together with dynrData to drop
	# things from it that are not in usevars.
	model$num_sbj <- as.integer(length(unique(data[['id']])))
	model$dim_obs_var <- as.integer(ncol(data$observed))
	if ("covariates" %in% names(data)){
	  model$dim_co_variate <- as.integer(ncol(data$covariates))
	}else{
	  model$dim_co_variate <- as.integer(0)
	}
	return(model)
}


is.null.pointer <- function(pointer){
	a <- attributes(pointer)
	attributes(pointer) <- NULL
	out <- identical(pointer, new("externalptr"))
	attributes(pointer) <- a
	return(out)
}


# N.B. matrixcalc has a similar function based on the eigenvalues.
# On a  3x3 matrix ...
# If the matrix is PD, this check is 4x faster.
# If the matrix is not PD, this check is 1.25x slower.
# If the matrix has a zero eigenvalue, this check may not be useful.
is.positive.definite <- function(x){
	ret <- try(chol(x), silent=TRUE)
	ifelse(any(is(ret) %in% "try-error"), FALSE, TRUE)
}

dynrExitFlags <- c(
	'1'='Optimization halted because of a forced termination.',
	'2'='Optimization halted because of roundoff errors.',
	'3'='Optimization failed. Ran out of memory.',
	'4'='Optimization halted. Lower bounds are bigger than upper bounds.',
	'5'='Optimization failed. Check starting values.',
	'6'='Optimization terminated successfully',
	'7'='Optimization stopped because objective function reaches stopval',
	'8'='Optimization terminated successfully: ftol_rel or ftol_abs was reached.',
	'9'='Optimization terminated successfully: xtol_rel or xtol_abs was reached.',
	'10'='Maximum number of function evaluations reached.',
	'11'='Increase maxeval or change starting values.',
	'12'='Maximum optimization time reached.',
	'13'='Increase maxtime or change starting values.')


vec2mat<-function(vectr,dimension){
    n.vec=length(vectr)
    if (n.vec==dimension){
  	  #diagonal
  	  mat=diag(vectr)
    }else if (dimension==sqrt(2*n.vec+.25) - .5){
    	  #symmetric
  	  mat=matrix(0,dimension,dimension)
  	  diag(mat)<-vectr[1:dimension]
  	  for (i in 1:(dimension-1)){
  		  for (j in (i+1):dimension){
  		  	mat[i,j]<-vector[i+j+dimension-2]
  			mat[j,i]<-vector[i+j+dimension-2]
  		  }
  	  }
    }else{
  	  cat('Length of the vector does not match the matrix dimension!\n')
    }	
	return(mat)	
}

#transldl function for caluaclating the LDL values
transldl <- function(mat){
  L <- mat
  diag(L)<-1
  L[upper.tri(L)]<-0
  D<- diag(exp(diag(mat)))
  # final caluclation
  outldl <- L %*% D %*% t(L)
  return(outldl)
}
