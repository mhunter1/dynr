

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
			transformed.hessian =  "matrix",
			conf.intervals = "matrix",
			exit.flag = "numeric",
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
		.Object@transformed.hessian <- x$transformed.hessian
		.Object@conf.intervals <- x$CI
		.Object@exit.flag <- x$exitflag
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
#setMethod("summary", "dynrRun",
#	function(object){
#		tval <- object@transformed.parameters/object@standard.errors
#		ret <- data.frame(transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors, object@conf.intervals, t.value=tval, p.value=2*pt(abs(tval), df=1, lower.tail=FALSE))
#		cat('\nSummary of dynrRun\n\n')
#		cat(paste('Log Likelihood: ', logLik(object), '\n', sep=''))
#		cat(paste('AIC: ', AIC(object), '\n', sep=''))
#		cat(paste('BIC: ', BIC(object), '\n', sep=''))
#		cat('Parameters:\n')
#		#print(ret)
#		return(ret)
#	}
#)
# See Also the print method of summary.lm
#  getAnywhere(print.summary.lm)

setMethod( f = "summary" ,  signature = "dynrRun" ,
           definition = function(object){
             d <- data.frame(transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors)
             d$T_value<-ifelse(d$standard.errors==0, NA, d$transformed.parameters/d$standard.errors)
             d <-cbind(d,object@conf.intervals)
             Likfin = -1*object@neg.log.likelihood
             npar = length(object@fitted.parameters)
             AIC = 2*Likfin +2*npar
             InfDs.allT_1 = (dim(object@pr_t_given_T)[2])
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
	attr(ans, "nobs") <- dim(object@eta_regime_t)[3]
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
	output <- endProcessing(output, transformation, conf.level)
	obj <- new("dynrRun", output)
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
	J <- numDeriv::jacobian(func=transformation, x=x$fitted.parameters)
	tHess <- J %*% x$inverse.hessian.matrix %*% t(J)
	tSE <- sqrt(abs(diag(tHess)))
	tParam <- transformation(x$fitted.parameters)
	x$transformed.parameters <- tParam
	x$standard.errors <- tSE
	x$transformed.hessian <- tHess
	confx <- qnorm(1-(1-conf.level)/2)
	x$CI <- matrix(c(tParam - tSE*confx, tParam + tSE*confx), ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
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



