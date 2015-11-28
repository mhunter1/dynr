
# Class definition for the dynrRun object that
#  stores all the output for a model that has 
#  been run.

# L = number of latent variables
# R = number of regimes
# T = number of time points
setClass(Class =  "dynrRun",
         representation = representation(
           fitted.parameters =  "numeric",
           transformed.parameters =  "numeric",
           standard.errors =  "numeric",
           hessian =  "matrix",
           transformed.inv.hessian =  "matrix",
           conf.intervals = "matrix",
           exitflag = "numeric",
           neg.log.likelihood = "numeric",
           inverse.hessian = "matrix",
           eta_regime_t = "array", # LxRxT
           error_cov_regime_t  = "array", # LxLxRxT
           eta_regime_regime_t_pred= "array", # LxRxRxT
           error_cov_regime_regime_t_pred = "array", # LxLxRxRxT
           eta_regime_regime_t_plus_1 = "array", # LxRxRxT
           error_cov_regime_regime_t_plus_1 = "array", # LxLxRxRxT
           innov_vec  = "array", # LxRxRxT
           inverse_residual_cov = "array", # LxLxRxRxT
           pr_t_given_t   = "matrix", # RxT
           pr_t_given_t_less_1 = "matrix", # RxT
           pr_t_given_T  = "matrix", # RxT
           transprob_given_T  = "array", # RxRxT
           eta_regime_smooth = "array", # LxRxT
           error_cov_regime_smooth  = "array", # LxLxRxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array" # LxLxT
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

setMethod( f = "summary",  signature = "dynrRun" ,
           definition = function(object){
             d <- data.frame(transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors)
             d$T_value<-ifelse(d$standard.errors==0, NA, d$transformed.parameters/d$standard.errors)
             d <-cbind(d,object@conf.intervals)
             Likfin = -1*object@neg.log.likelihood
             npar = length(object@fitted.parameters)
             AIC = 2*Likfin +2*npar
             InfDs.allT_1 = dim(object@eta_smooth_final)[2]
             BIC = 2*Likfin +npar*log(InfDs.allT_1)
             cat("********************************SUMMARY*****************************************\n")
             print(d)
             cat("\nAIC=");print(AIC)
             cat("\nBIC=");print(BIC)
             cat("**********END OF SUMMARY************\n")
           }
)


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
# logLike
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


dynr.run <- function(model, data, func_address, transformation, conf.level=.95) {
	frontendStart <- Sys.time()
	if(missing(transformation)){
		transformation <- function(x){x}
	}
	model <- combineModelDataInformation(model, data)
	model <- preProcessModel(model)
	backendStart <- Sys.time()
	output <- .Call("main_R", model, data, func_address, PACKAGE = "dynr")
	backendStop <- Sys.time()
	gc()#garbage collection
	print(output$exitflag)
	output$exitflag <- output$exitflag+ifelse(output$exitflag<0,6,0)+ifelse(output$exitflag>0,5,0)
	print(output$exitflag)
	switch(output$exitflag,
	       {	cat('Optimization halted because of a forced termination.','\n')
	       },
	       {	cat('Optimization halted because of roundoff errors.','\n')
	       },
	       {	cat('Optimization failed. Ran out of memory.','\n')
	       },
	       {	cat('Optimization halted. Lower bounds are bigger than upper bounds.','\n')
	       },
	       {	cat('Optimization failed. Check starting values.','\n')
	       },
	       {	cat('Optimization terminated successfully','\n')
	       },
	       {
	         cat('Optimization stopped because objective function reaches stopval','\n') 
	       },
	       {
	         cat('Optimization terminated successfully: ftol_rel or ftol_abs was reached.','\n')
	       },
	       {
	         cat('Optimization terminated successfully: xtol_rel or xtol_abs was reached.','\n')
	       },  
	       {
	         cat('Maximum number of function evaluations reached.','\n',
	             'Increase maxeval or change starting values.','\n')
	       },
	       {
	         cat('Maximum optimization time reached.','\n',
	             'Increase maxtime or change starting values.','\n')
	       }
	)
	if (output$exitflag > 5){
	output2 <- endProcessing(output, transformation, conf.level)
	obj <- new("dynrRun", output2)
	}else{
	obj <- NULL
	}
	frontendStop <- Sys.time()
	totalTime <- frontendStop-frontendStart
	backendTime <- backendStop-backendStart
	fronendTime <- totalTime-backendTime
	cat('Total Time:', totalTime, '\n')
	cat('Backend Time:', backendTime, '\n')
	# TODO add timing information to dynrRun object
	return(obj)
}

endProcessing <- function(x, transformation, conf.level){
  confx <- qnorm(1-(1-conf.level)/2)
  x$hessian.matrix[x$hessian.matrix==0] = 10e-14
  status = ifelse(abs(det(x$hessian.matrix)) <1e-6, 0, 1) #0 = singular hessian matrix
    
  if (x$exitflag < 0 || status == 0){
    tSE <- rep(NA,length(x$fitted.parameters))
    tHess <- matrix(NA,dim(x$inverse.hessian.matrix)[1],dim(x$inverse.hessian.matrix)[2])
    tParam <- rep(NA,length(x$fitted.parameters))
    CI <- matrix(NA,length(x$fitted.parameters),2)
  }else{
    #Analytic Jacobian
  V1 = solve(x$hessian.matrix);
  J <- numDeriv::jacobian(func=transformation, x=x$fitted.parameters)
	tHess <- J %*% V1%*%t(J)
	tSE <- sqrt(abs(diag(tHess)))
	tParam <- transformation(x$fitted.parameters)
	CI <- c(tParam - tSE*confx, tParam + tSE*confx)
}
	x$transformed.parameters <- tParam
	x$standard.errors <- tSE
	x$transformed.inv.hessian <- tHess
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
	# TODO remove these as.double calls
	# by changing how the arguments are processed in the backend.
	model$num_sbj <- as.double(length(unique(data[['id']])))
	model$dim_obs_var <- as.double(ncol(data$observed))
	model$dim_co_variate <- as.double(ncol(data$covariates))
	return(model)
}



