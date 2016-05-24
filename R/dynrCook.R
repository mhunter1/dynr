# Class definition for the dynrCook object that
#  stores all the output for a model that has 
#  been run.

# L = number of latent variables
# R = number of regimes
# T = number of time points

setClass(Class =  "dynrCook",
         representation = representation(
           fitted.parameters =  "numeric", #Can return
           transformed.parameters =  "numeric", #
           standard.errors =  "numeric",
           hessian =  "matrix",
           transformed.inv.hessian =  "matrix",
           conf.intervals = "matrix",
           exitflag = "numeric", #
           neg.log.likelihood = "numeric", #
           #Everything else from this point on
           pr_t_given_T  = "matrix", # RxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array", # LxLxT
           run.times = "numeric",
           param.names = "character"
         )
)

# Initialize method for the new() function
setMethod("initialize", "dynrCook",
          function(.Object, x){
            .Object@fitted.parameters <- x$fitted.parameters
            .Object@transformed.parameters <- x$transformed.parameters
            .Object@standard.errors <- x$standard.errors
            .Object@hessian <- x$hessian.matrix
            .Object@transformed.inv.hessian <- x$transformed.inv.hessian
            .Object@conf.intervals <- x$conf.intervals
            .Object@exitflag <- x$exitflag
            .Object@neg.log.likelihood <- x$neg.log.likelihood
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
           eta_regime_t = "array", # LxRxT #
           error_cov_regime_t  = "array", # LxLxRxT #
           innov_vec  = "array", # LxRxRxT #
           inverse_residual_cov = "array", # LxLxRxRxT #
           #Everything else from this point on
           pr_t_given_T  = "matrix", # RxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array", # LxLxT
           run.times = "numeric",
           param.names = "character"
         ),
         contains = "dynrCook"
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


setClass(Class =  "dynrOutall",
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
           run.times = "numeric",
           param.names = "character"
         ),
         contains = "dynrCook"
)

# Initialize method for the new() function
setMethod("initialize", "dynrOutall",
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

# Set the summary method of an object of class dynrCook
#  All this amounts to is writing a function that takes a
#  dynrCook object (and possibly other arguments)
#  and returns/does whatever we want.
# See Also the print method of summary.lm
#  getAnywhere(print.summary.lm)
summaryResults<-function(object){
             d <- data.frame(names=object@param.names, transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors)
             d$t_value<-ifelse(d$standard.errors==0, NA, d$transformed.parameters/d$standard.errors)
             d <-cbind(d,object@conf.intervals)
             neg2LL = -2*logLik(object)
             AIC = AIC(object)
             BIC = BIC(object)
             colnames(d) = c("names","parameters","s.e.","t-value","ci.lower","ci.upper")
             print(d)
             cat(paste0("\n-2 log-likelihood value at convergence = ", signif(neg2LL,2)))
             cat(paste0("\nAIC = ", signif(AIC,2)))
             cat(paste0("\nBIC = ", signif(BIC,2)))
           }

setMethod( f = "summary",  signature = "dynrCook" ,
           definition = summaryResults)


display.dynrCook <- function(x, ...){
  str(x)
  invisible(x)
}

setMethod("print", "dynrCook", function(x, ...) { 
  display.dynrCook(x) 
})

setMethod("show", "dynrCook", function(object) { 
  display.dynrCook(object) 
})


#------------------------------------------------------------------------------
# Some S3 methods
# coef
# logLik
# AIC
# BIC

##' Extract fitted parameters from a dynrCook Object
##' 
##' @param object The dynrCook object for which the coefficients are desired
##' @param ... further named arguments, ignored for this method
##' 
##' @return A numeric vector of the fitted parameters.
##' 
##' @seealso Other S3 methods \code{\link{logLik.dynrCook}}
##' 
##' @examples
##' # Let cookedModel be the output from dynr.cook
##' #coef(cookedModel)
coef.dynrCook <- function(object, ...){
  object@transformed.parameters
}

#`coef<-.dynrCook` <- function(object, value){
#	object <- PopBackModel(object, value)
#	return(object)
#}

##' Extract the log likelihood from a dynrCook Object
##' 
##' @param object The dynrCook object for which the log likelihood is desired
##' @param ... further named arguments, ignored for this method
##' 
##' @details
##' The 'df' attribute for this object is the number of freely estimated parameters. The 'nobs' attribute is the total number of rows of data, adding up the number of time points for each person.
##' 
##' @return An object of class \code{logLik}.
##' 
##' @seealso Other S3 methods \code{\link{coef.dynrCook}}
##' 
##' @examples
##' # Let cookedModel be the output from dynr.cook
##' logLik(cookedModel)
logLik.dynrCook <- function(object, ...){
  ans <- -object@neg.log.likelihood
  attr(ans, "df") <- length(object@fitted.parameters)
  attr(ans, "nobs") <- dim(object@eta_smooth_final)[2]
  class(ans) <- "logLik"
  return(ans)
}

# N.B. AIC() and BIC() are implicitly defined in terms
#  of logLik().

##' Extract the Variance-Covariance Matrix of a dynrCook object
##' 
##' @param object The dynrCook object for which the variance-covariance matrix is desired
##' @param ... further named arguments, ignored by this method
##' 
##' @details
##' This is the inverse Hessian of the transformed parameters.
vcov.dynrCook <- function(object, ...){
	nm <- names(coef(object))
	rt <- object@transformed.inv.hessian
	dimnames(rt) <- list(nm, nm)
	return(rt)
}

#------------------------------------------------------------------------------

##' Cook a dynr model to estimate its free parameters
##' 
##' @param dynrModel a dynr model compiled using dynr.model, consisting of recipes for submodels, starting values, parameter names, and C code for each submodel
##' @param conf.level a cumulative proportion indicating the level of desired confidence intervals for the final parameter estimates (default is .95)
##' @param infile (not required for models specified through the recipe functions) the name of a file that has the C codes for all dynr submodels for those interested in specifying a model directly in C
##' @param verbose a flag (TRUE/FALSE) indicating whether more detailed intermediate output during the estimation process should be printed
##' @param debug_flag a flag (TRUE/FALSE) indicating whether users want additional dynr output that can be used for diagnostic purposes 
##' 
##' @details
##' TO BE COMPLETED: 
##' for "cooking dinner" -- namely, to start the estimation process
##' a description of things output when debug_flag = FALSE
##' a description of things output when debug_flag = TRUE
##' 
##' @seealso \code{\link{dynr.cook}}
##' 
##' @examples
##' #fitted.model <- dynr.cook(model, data)
dynr.cook <- function(dynrModel, conf.level=.95, infile, verbose=TRUE, debug_flag=FALSE) {
	outall_flag=FALSE#always set to FALSE except when a developer wants all the intermediate products from the C estimation algorithms.
	frontendStart <- Sys.time()
	transformation=dynrModel@transform@tfun
	data <- dynrModel$data

	#internalModelPrep convert dynrModel to a model list
	model <- internalModelPrep(
	  num_regime=dynrModel@num_regime,
	  dim_latent_var=dynrModel@dim_latent_var,
	  xstart=dynrModel@xstart,
	  ub=dynrModel@ub,
	  lb=dynrModel@lb,
	  options=dynrModel@options,
	  isContinuousTime=dynrModel@dynamics@isContinuousTime,
	  infile=dynrModel@outfile,
	  outfile=dynrModel@outfile,
	  compileLib=dynrModel@compileLib,
	  verbose=dynrModel@verbose
	)
	
	model <- combineModelDataInformation(model, data)
	model <- preProcessModel(model)
	if(any(sapply(model$func_address, is.null.pointer))){
	    warning("Found null pointer(s) in 'func_address' list. (Re-)compiling your functions...")
	    if(missing(infile)){
	      stop("Cannot compile your functions because 'infile' argument is missing.")
	    }
	    model$func_address=.C2funcaddress(isContinuousTime=model$isContinuousTime,infile=infile,verbose=verbose)
	}
	gc()
	backendStart <- Sys.time()
	output <- .Call(.Backend, model, data, debug_flag, outall_flag, PACKAGE = "dynr")
	backendStop <- Sys.time()
	#gc()#garbage collection
	cat('Original exit flag: ', output$exitflag, '\n')
	output$exitflag <- output$exitflag + ifelse(output$exitflag<0, 6, 0) + ifelse(output$exitflag>0, 5, 0)
	cat('Modified exit flag: ', output$exitflag, '\n')
	cat(.ExitFlags[output$exitflag], '\n')
	
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
	names(output2$transformed.parameters) <- dynrModel$param.names
	if (outall_flag){
		obj <- new("dynrOutall", output2)
	}else if(debug_flag){
		obj <- new("dynrDebug", output2)
	}else{
		obj <- new("dynrCook", output2)
	}
	
	obj@param.names <- dynrModel$param.names
	#populate transformed estimates to dynrModel
	#model<<-PopBackModel(model, obj@transformed.parameters)
  
	frontendStop <- Sys.time()
	totalTime <- frontendStop-frontendStart
	backendTime <- backendStop-backendStart
	frontendTime <- totalTime-backendTime
	obj@run.times <- as.numeric(c(totalTime=totalTime, backendTime=backendTime, frontendTime=frontendTime))
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
  x$conf.intervals <-matrix(999, nrow=nParam,ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
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

.ExitFlags <- c(
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

PopBackMatrix<-function(values.matrix, param.matrix, trans.parameters){
  if (class(values.matrix)=="list"){
    num_regime=length(values.matrix)
    if (num_regime>0){
      for (i in 1:num_regime){
         values.matrix[[i]][which(param.matrix[[i]]!=0,arr.ind = TRUE)]<-
        trans.parameters[param.matrix[[i]][which(param.matrix[[i]]!=0,arr.ind = TRUE)]]
      }
    }
  }else{
    values.matrix[which(param.matrix!=0,arr.ind = TRUE)]<-trans.parameters[param.matrix[which(param.matrix!=0,arr.ind = TRUE)]]
  }
  return(values.matrix)
}

PopBackFormula<- function(formula, paramnames, param.names, trans.parameters){
  string<-paste0(deparse(formula,width.cutoff = 500L),collapse="")
  for (i in 1:length(paramnames)){
    string<-gsub(paste0("param\\[", match(paramnames[i], param.names, nomatch=0)-1,"\\]"), trans.parameters[match(paramnames[i], param.names, nomatch=0)], string, perl = TRUE)
  }
  eval(parse(text=string))
}

PopBackModel<-function(dynrModel, trans.parameters){
  
  if (class(dynrModel$dynamics) == 'dynrDynamicsFormula'){
    dynrModel@dynamics@formula<-PopBackFormula(dynrModel@dynamics@formula,dynrModel@dynamics@paramnames,dynrModel@param.names,trans.parameters)
  }else{
    dynrModel@dynamics@values.dyn <- PopBackMatrix(dynrModel@dynamics@values.dyn, dynrModel@dynamics@params.dyn, trans.parameters)
    dynrModel@dynamics@values.exo <- PopBackMatrix(dynrModel@dynamics@values.exo, dynrModel@dynamics@params.exo, trans.parameters)
    dynrModel@dynamics@values.int <- PopBackMatrix(dynrModel@dynamics@values.int, dynrModel@dynamics@params.int, trans.parameters)
  }
  
    dynrModel@measurement@values.load<-PopBackMatrix(dynrModel@measurement@values.load, dynrModel@measurement@params.load, trans.parameters)
    dynrModel@measurement@values.exo<-PopBackMatrix(dynrModel@measurement@values.exo, dynrModel@measurement@params.exo, trans.parameters)
    dynrModel@measurement@values.int<-PopBackMatrix(dynrModel@measurement@values.int, dynrModel@measurement@params.int, trans.parameters)
  
  dynrModel@noise@values.latent<-PopBackMatrix(dynrModel@noise@values.latent, dynrModel@noise@params.latent, trans.parameters)
  dynrModel@noise@values.observed<-PopBackMatrix(dynrModel@noise@values.observed, dynrModel@noise@params.observed, trans.parameters)
  dynrModel@initial@values.inistate<-PopBackMatrix(dynrModel@initial@values.inistate, dynrModel@initial@params.inistate , trans.parameters) 
  dynrModel@initial@values.inicov<-PopBackMatrix(dynrModel@initial@values.inicov, dynrModel@initial@params.inicov, trans.parameters) 
  dynrModel@initial@values.regimep<-PopBackMatrix(dynrModel@initial@values.regimep, dynrModel@initial@params.regimep, trans.parameters)
  dynrModel@regimes@values<-PopBackMatrix(dynrModel@regimes@values, dynrModel@regimes@params, trans.parameters)
  
  return(dynrModel)
}
