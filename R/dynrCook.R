# Class definition for the dynrCook object that
#  stores all the output for a model that has 
#  been run.

# L = number of latent variables
# R = number of regimes
# T = number of time points

##' The dynrCook Class
##' 
##' @aliases
##' dynrDebug-class
##' $,dynrCook-method
##' print,dynrCook-method
##' show,dynrCook-method
##' 
##' @details
##' This is an internal class structure.  You should not use it directly.
##' Use \code{\link{dynr.cook}} instead.
setClass(Class =  "dynrCook", 
         representation = representation(
           fitted.parameters =  "numeric", #Can return
           transformed.parameters =  "numeric", #
           standard.errors =  "numeric",
           standard.errors.untransformed = "numeric",
           bad.standard.errors = "logical",
           hessian =  "matrix",
           transformed.inv.hessian =  "matrix",
           inv.hessian ='matrix',
           conf.intervals = "matrix",
           conf.intervals.endpoint.trans = "matrix",
           exitflag = "numeric", #
           neg.log.likelihood = "numeric", #
           pr_t_given_T  = "matrix", # RxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array", # LxLxT
           pr_t_given_t  = "matrix", # RxT
           eta_filtered = "matrix", # LxT
           error_cov_filtered = "array", # LxLxT
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
            .Object@standard.errors.untransformed <- x$standard.errors.untransformed
            .Object@bad.standard.errors <- x$bad.standard.errors
            .Object@hessian <- x$hessian.matrix
            .Object@transformed.inv.hessian <- x$transformed.inv.hessian
            .Object@inv.hessian <- x$inv.hessian
            .Object@conf.intervals <- x$conf.intervals
            .Object@conf.intervals.endpoint.trans <- x$conf.intervals.endpoint.trans
            .Object@exitflag <- x$exitflag
            .Object@neg.log.likelihood <- x$neg.log.likelihood
            .Object@pr_t_given_T <- x$pr_t_given_T
            .Object@eta_smooth_final <- x$eta_smooth_final
            .Object@error_cov_smooth_final <- x$error_cov_smooth_final
            .Object@pr_t_given_t <- x$pr_t_given_t
            .Object@eta_filtered <- x$eta_filtered
            .Object@error_cov_filtered <- x$error_cov_filtered
            return(.Object)
          }
)

setClass(Class =  "dynrDebug",
         representation = representation(
           fitted.parameters =  "numeric", #Can return
           transformed.parameters =  "numeric", #
           standard.errors =  "numeric",
           standard.errors.untransformed='numeric',
           bad.standard.errors = "logical",
           hessian =  "matrix",
           transformed.inv.hessian =  "matrix",
           inv.hessian =  "matrix",
           conf.intervals = "matrix",
           conf.intervals.endpoint.trans="matrix",
           exitflag = "numeric", #
           neg.log.likelihood = "numeric", #
           #Everything else from this point on
           pr_t_given_T  = "matrix", # RxT
           eta_smooth_final = "matrix", # LxT
           error_cov_smooth_final  = "array", # LxLxT
           pr_t_given_t  = "matrix", # RxT
           eta_filtered = "matrix", # LxT
           error_cov_filtered = "array", # LxLxT
           eta_predicted = "matrix", # LxT
           error_cov_predicted = "array", # LxLxT
           innov_vec = "matrix", # LxT
           residual_cov = "array", # LxLxT
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
            .Object@standard.errors.untransformed <- x$standard.errors.untransformed
            .Object@bad.standard.errors <- x$bad.standard.errors
            .Object@hessian <- x$hessian.matrix
            .Object@transformed.inv.hessian <- x$transformed.inv.hessian
            .Object@inv.hessian <- x$inv.hessian
            .Object@conf.intervals <- x$conf.intervals
            .Object@conf.intervals.endpoint.trans <- x$conf.intervals.endpoint.trans
            .Object@exitflag <- x$exitflag
            .Object@neg.log.likelihood <- x$neg.log.likelihood
            .Object@pr_t_given_T <- x$pr_t_given_T
            .Object@eta_smooth_final <- x$eta_smooth_final
            .Object@error_cov_smooth_final <- x$error_cov_smooth_final
            .Object@pr_t_given_t <- x$pr_t_given_t
            .Object@eta_filtered <- x$eta_filtered
            .Object@error_cov_filtered <- x$error_cov_filtered
            .Object@eta_predicted <- x$eta_predicted
            .Object@error_cov_predicted <- x$error_cov_predicted
            .Object@innov_vec <- x$innov_vec
            .Object@residual_cov <- x$residual_cov
            return(.Object)
          }
)

# Set the summary method of an object of class dynrCook
#  All this amounts to is writing a function that takes a
#  dynrCook object (and possibly other arguments)
#  and returns/does whatever we want.
# See Also the print method of summary.lm
#  getAnywhere(print.summary.lm)
summaryResults <- function(object, ...){
             ret <- list()
             d <- data.frame(transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors)
             row.names(d) <- names(coef(object))
             d$t_value <- ifelse(d$standard.errors==0, NA, d$transformed.parameters/d$standard.errors)
             d <- cbind(d, object@conf.intervals)
             #d$bad <- factor(c("", "!")[object@bad.standard.errors + 1])
             d$p <- pt( abs(d$t_value), df=nobs(object) - nrow(d), lower.tail=FALSE)
             neg2LL = -2*logLik(object)
             AIC = AIC(object)
             BIC = BIC(object)
             colnames(d) = c("Estimate", "Std. Error", "t value", "ci.lower", "ci.upper", #"",
               "Pr(>|t|)")
             ret$Coefficients <- d
             ret$neg2LL <- neg2LL
             ret$AIC <- AIC
             ret$BIC <- BIC
             class(ret) <- "summary.dynrCook"
             return(ret)
           }

print.summary.dynrCook <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...){
	cat("Coefficients:\n")
	printCoefmat(x$Coefficients, cs.ind=c(1L, 2L, 4L, 5L), tst.ind=3L, zap.ind=6L, digits = digits, signif.stars = signif.stars, ...)
	cat(paste0("\n-2 log-likelihood value at convergence = ", sprintf("%.02f", round(x$neg2LL,2))))
	cat(paste0("\nAIC = ", sprintf("%.02f", round(x$AIC,2))))
	cat(paste0("\nBIC = ", sprintf("%.02f", round(x$BIC,2))))
	cat("\n")
	invisible(x)
}

coef.summary.dynrCook <- function(object, ...){
	object$Coefficients
}

##' Get the summary of a dynrCook object
##' 
##' @param object The dynrCook object for which the summary is desired.
##' @param ... Further named arguments, passed to the print method (e.g., \code{digits} and \code{signif.stars}).
##' 
##' @details
##' The summary gives information on the free parameters estimated: names, parameter values, numerical Hessian-based standard errors, t-values (values divided by standard errors), and standard-error based confidence intervals.  Additionally, the likelihood, AIC, and BIC are provided.
##' 
##' Note that an exclamation point (!) in the final column of the summary table indicates that
##' the standard error and confidence interval for this parameter may not be trustworthy. The corresponding
##' element of the (transformed, inverse) Hessian was negative and an absolute value was taken to make it positive.
##' 
##' @method summary dynrCook
##' @return Object of class summary.dynrCook.  Primarily used for showing the results of a fitted model.
summary.dynrCook <- summaryResults


displayDynrCook <- function(x, ...){
	maxChar <- max(nchar(sn <- slotNames(x)))
	tfmt <- paste0("..$ %-", maxChar, "s:")
	for(aslot in sn){
		if( !(aslot %in% c("conf.intervals", "hessian", "transformed.inv.hessian", "neg.log.likelihood", "param.names")) ){
			cat(sprintf(tfmt, aslot))
			str(slot(x, aslot))
		}
	}
	return(invisible(x))
}

setMethod("print", "dynrCook", function(x, ...) { 
	displayDynrCook(x) 
})

setMethod("show", "dynrCook", function(object) { 
	displayDynrCook(object) 
})


#------------------------------------------------------------------------------
# Some S3 methods
# coef
# logLik
# AIC
# BIC

##' Extract fitted parameters from a dynrCook Object
##' 
##' aliases coef.dynrModel coef<- coef<-.dynrModel
##' 
##' @param object The dynrCook object for which the coefficients are desired
##' @param ... further named arguments, ignored for this method
##' 
##' @return A numeric vector of the fitted parameters.
##' 
##' @seealso Other S3 methods \code{\link{logLik.dynrCook}}
##' 
##' @examples
##' # Create a minimal cooked model called 'cook'
##' require(dynr)
##' 
##' meas <- prep.measurement(
##' 	values.load=matrix(c(1, 0), 1, 2),
##' 	params.load=matrix(c('fixed', 'fixed'), 1, 2),
##' 	state.names=c("Position","Velocity"),
##' 	obs.names=c("y1"))
##' 
##' ecov <- prep.noise(
##' 	values.latent=diag(c(0, 1), 2),
##' 	params.latent=diag(c('fixed', 'dnoise'), 2),
##' 	values.observed=diag(1.5, 1),
##' 	params.observed=diag('mnoise', 1))
##' 
##' initial <- prep.initial(
##' 	values.inistate=c(0, 1),
##' 	params.inistate=c('inipos', 'fixed'),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2))
##' 
##' dynamics <- prep.matrixDynamics(
##' 	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
##' 	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##' 	isContinuousTime=TRUE)
##' 
##' data(Oscillator)
##' data <- dynr.data(Oscillator, id="id", time="times", observed="y1")
##' 
##' model <- dynr.model(dynamics=dynamics, measurement=meas,
##' 	noise=ecov, initial=initial, data=data)
##' 
##' cook <- dynr.cook(model,
##' 	verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE)
##' 
##' # Now grab the coef!
##' coef(cook)
coef.dynrCook <- function(object, ...){
	object@transformed.parameters
}

#`coef<-.dynrCook` <- function(object, value){
#	object <- PopBackModel(object, value)
#	return(object)
#}

##' Extract the log likelihood from a dynrCook Object
##' 
##' @aliases deviance.dynrCook
##' 
##' @param object The dynrCook object for which the log likelihood is desired
##' @param ... further named arguments, ignored for this method
##' 
##' @details
##' The 'df' attribute for this object is the number of freely estimated parameters. The 'nobs' attribute is the total number of rows of data, adding up the number of time points for each person.
##' 
##' The \code{deviance} method returns minus two times the log likelihood.
##' 
##' @return In the case of \code{logLik}, an object of class \code{logLik}.
##' 
##' @seealso Other S3 methods \code{\link{coef.dynrCook}}
##' 
##' @examples
##' # Minimal model
##' require(dynr)
##' 
##' meas <- prep.measurement(
##' 	values.load=matrix(c(1, 0), 1, 2),
##' 	params.load=matrix(c('fixed', 'fixed'), 1, 2),
##' 	state.names=c("Position","Velocity"),
##' 	obs.names=c("y1"))
##' 
##' ecov <- prep.noise(
##' 	values.latent=diag(c(0, 1), 2),
##' 	params.latent=diag(c('fixed', 'dnoise'), 2),
##' 	values.observed=diag(1.5, 1),
##' 	params.observed=diag('mnoise', 1))
##' 
##' initial <- prep.initial(
##' 	values.inistate=c(0, 1),
##' 	params.inistate=c('inipos', 'fixed'),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2))
##' 
##' dynamics <- prep.matrixDynamics(
##' 	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
##' 	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##' 	isContinuousTime=TRUE)
##' 
##' data(Oscillator)
##' data <- dynr.data(Oscillator, id="id", time="times", observed="y1")
##' 
##' model <- dynr.model(dynamics=dynamics, measurement=meas,
##' 	noise=ecov, initial=initial, data=data)
##' 
##' cook <- dynr.cook(model,
##' 	verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE)
##' 
##' # Now get the log likelihood!
##' logLik(cook)
logLik.dynrCook <- function(object, ...){
	ans <- -object@neg.log.likelihood
	attr(ans, "df") <- length(object@fitted.parameters)
	attr(ans, "nobs") <- nobs(object) #dim(object@eta_smooth_final)[2]
	class(ans) <- "logLik"
	return(ans)
}

# N.B. AIC() and BIC() are implicitly defined in terms
#  of logLik().

##' @rdname logLik.dynrCook
deviance.dynrCook <- function(object, ...){
	as.numeric(-2*logLik(object))
}

##' Extract the number of observations for a dynrCook object
##' 
##' @param object A fitted model object
##' @param ... Further named arguments. Ignored.
##' 
##' @details
##' We return the total number of rows of data, adding up the number of time points for each person. For some purposes, you may want the mean number of observations per person or the number of people instead.  These are not currently supported via \code{nobs}.
##' 
##' @return
##' A single number. The total number of observations across all IDs.
##' 
##' @examples
##' # Minimal model
##' require(dynr)
##' 
##' meas <- prep.measurement(
##' 	values.load=matrix(c(1, 0), 1, 2),
##' 	params.load=matrix(c('fixed', 'fixed'), 1, 2),
##' 	state.names=c("Position","Velocity"),
##' 	obs.names=c("y1"))
##' 
##' ecov <- prep.noise(
##' 	values.latent=diag(c(0, 1), 2),
##' 	params.latent=diag(c('fixed', 'dnoise'), 2),
##' 	values.observed=diag(1.5, 1),
##' 	params.observed=diag('mnoise', 1))
##' 
##' initial <- prep.initial(
##' 	values.inistate=c(0, 1),
##' 	params.inistate=c('inipos', 'fixed'),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2))
##' 
##' dynamics <- prep.matrixDynamics(
##' 	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
##' 	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##' 	isContinuousTime=TRUE)
##' 
##' data(Oscillator)
##' data <- dynr.data(Oscillator, id="id", time="times", observed="y1")
##' 
##' model <- dynr.model(dynamics=dynamics, measurement=meas,
##' 	noise=ecov, initial=initial, data=data)
##' 
##' cook <- dynr.cook(model,
##' 	verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE)
##' 
##' # Now get the total number of observations
##' nobs(cook)
nobs.dynrCook <- function(object, ...){
	dim(object@eta_smooth_final)[2]
}
# TODO could give sample size for each individual through the ...
# arguments


##' Extract the Variance-Covariance Matrix of a dynrCook object
##' 
##' @param object The dynrCook object for which the variance-covariance matrix is desired
##' @param ... further named arguments, ignored by this method
##' 
##' @details
##' This is the inverse Hessian of the transformed parameters.
##'
##' @return matrix.  Asymptotic variance-covariance matrix of the transformed parameters.
vcov.dynrCook <- function(object, ...){
	nm <- names(coef(object))
	rt <- object@transformed.inv.hessian
	dimnames(rt) <- list(nm, nm)
	return(rt)
}
# TODO redefine this method as MDH proposed
#vcov.dynrCook <- function(object, transformed=TRUE, ...){
#	nm <- names(coef(object))
#	 if(transformed){
#  	    rt <- object@transformed.inv.hessian
#	 } else {
#	   rt <- object@inv.hessian}
#	dimnames(rt) <- list(nm, nm)
#	return(rt)
#}

##' Extract the free parameter names of a dynrCook object
##' 
##' @param x The dynrCook object from which the free parameter names are desired
setMethod("names", "dynrCook",
	function(x) {
		pnames <- names(coef(x))
		output <- c(pnames)
		output <- gsub("(\\w+\\W+.*)", "'\\1'", output)
		return(output)
	}
)

setMethod("$", "dynrCook",
          function(x, name){slot(x, name)}
)

.DollarNames.dynrCook <- function(x, pattern){
	if(missing(pattern)){
		pattern <- ''
	}
	output <- slotNames(x)
	output <- gsub("(\\w+\\W+.*)", "'\\1'", output)
	return(grep(pattern, output, value=TRUE))
}

##' Confidence Intervals for Model Parameters
##' 
##' @param object a fitted model object
##' @param parm which parameters are to be given confidence intervals
##' @param level the confidence level
##' @param type The type of confidence interval to compute. See details. Partial name matching is used.
##' @param transformation For \code{type='endpoint.transformation'} the transformation function used.
##' @param ... further named arguments. Ignored.
##' 
##' @details
##' The \code{parm} argument can be a numeric vector or a vector of names. If it is missing then it defaults to using all the parameters.
##' 
##' These are Wald-type confidence intervals based on the standard errors of the (transformed) parameters.  Wald-type confidence intervals are known to be inaccurate for variance parameters, particularly when the variance is near zero (See references for issues with Wald-type confidence intervals).
##' 
##' @return
##' A matrix with columns giving lower and upper confidence limits for each parameter. These will be labelled as (1-level)/2 and 1 - (1-level)/2 as a percentage (e.g. by default 2.5% and 97.5%).
##' 
##' @references
##' Pritikin, J.N., Rappaport, L.M. & Neale, M.C.  (In Press). Likelihood-Based Confidence Intervals for a Parameter With an Upper or Lower Bound.  Structural Equation Modeling.  DOI: 10.1080/10705511.2016.1275969
##' 
##' Neale, M. C. & Miller M. B. (1997). The use of likelihood based confidence intervals in genetic models. Behavior Genetics, 27(2), 113-120.
##' 
##' Pek, J. & Wu, H. (2015). Profile likelihood-based confidence intervals and regions for structural equation models. Psychometrica, 80(4), 1123-1145.
##' 
##' Wu, H. & Neale, M. C. (2012). Adjusted confidence intervals for a bounded parameter. Behavior genetics, 42(6), 886-898.
##' 
##' @examples
##' # Minimal model
##' require(dynr)
##' 
##' meas <- prep.measurement(
##' 	values.load=matrix(c(1, 0), 1, 2),
##' 	params.load=matrix(c('fixed', 'fixed'), 1, 2),
##' 	state.names=c("Position","Velocity"),
##' 	obs.names=c("y1"))
##' 
##' ecov <- prep.noise(
##' 	values.latent=diag(c(0, 1), 2),
##' 	params.latent=diag(c('fixed', 'dnoise'), 2),
##' 	values.observed=diag(1.5, 1),
##' 	params.observed=diag('mnoise', 1))
##' 
##' initial <- prep.initial(
##' 	values.inistate=c(0, 1),
##' 	params.inistate=c('inipos', 'fixed'),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2))
##' 
##' dynamics <- prep.matrixDynamics(
##' 	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
##' 	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##' 	isContinuousTime=TRUE)
##' 
##' data(Oscillator)
##' data <- dynr.data(Oscillator, id="id", time="times", observed="y1")
##' 
##' model <- dynr.model(dynamics=dynamics, measurement=meas,
##' 	noise=ecov, initial=initial, data=data)
##' 
##' cook <- dynr.cook(model,
##' 	verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE)
##' 
##' # Now get the confidence intervals
##' # But note that they are nonsense because we set hessian_flag=FALSE !!!!
##' confint(cook)
confint.dynrCook <- function(object, parm, level = 0.95, type = c("delta.method", "endpoint.transformation"), transformation =  NULL, ...){
	type <- match.arg(type)
	tlev <- (1-level)/2
	confx <- qnorm(1-tlev)
	
	vals <- coef(object)
	if(missing(parm)){
		parm <- names(vals)
	}
	vals <- vals[parm]
	
	if(type == "delta.method"){
		iHess <- vcov(object)[parm, parm, drop=FALSE]
		SE <- sqrt(diag(iHess))
		CI <- matrix(c(vals - SE*confx, vals + SE*confx), ncol=2)
	}
	if(type == "endpoint.transformation"){
		tSEalt <- sqrt(diag(object$inv.hessian))
		CI <- matrix(c(transformation(object$fitted.parameters - tSEalt*confx), transformation(object$fitted.parameters + tSEalt*confx)),ncol=2)
		# TODO Fix the above to use the vcov method to extract the inverse Hessian instead
		#  and especially give the User the ability to only get SOME of the free parameters via the 'param' argument
	}
	dimnames(CI) <- list(names(vals), c(paste(tlev*100, "%"), paste((1 - tlev)*100, "%")) )
	return(CI)
}

#------------------------------------------------------------------------------

##' Cook a dynr model to estimate its free parameters
##' 
##' @param dynrModel a dynr model compiled using dynr.model, consisting of recipes for submodels, 
##' starting values, parameter names, and C code for each submodel
##' @param conf.level a cumulative proportion indicating the level of desired confidence intervals for
##' the final parameter estimates (default is .95)
##' @param infile (not required for models specified through the recipe functions) the name of a file 
##' that has the C codes for all dynr submodels for those interested in specifying a model directly in C
##' @param optimization_flag a flag (TRUE/FALSE) indicating whether optimization is to be done.
##' @param hessian_flag a flag (TRUE/FALSE) indicating whether the Hessian matrix is to be calculated.
##' @param verbose a flag (TRUE/FALSE) indicating whether more detailed intermediate output during the 
##' estimation process should be printed
##' @param weight_flag a flag (TRUE/FALSE) indicating whether the negative log likelihood function should 
##' be weighted by the length of the time series for each individual
##' @param debug_flag a flag (TRUE/FALSE) indicating whether users want additional dynr output that can 
##' be used for diagnostic purposes
##' @param perturb_flag a flag (TRUE/FLASE) indicating whether to perturb the latent states during estimation. Only useful for ensemble forecasting.
##' 
##' @details
##' Free parameter estimation uses the SLSQP routine from NLOPT.
##' 
##' The typical items returned in the cooked model are the filtered and smoothed latent variable estimates. 
##' \code{eta_smooth_final}, \code{error_cov_smooth_final} and \code{pr_t_given_T} are respectively 
##' time-varying smoothed latent variable mean estimates, smoothed error covariance estimates, 
##' and smoothed regime probability. 
##' \code{eta_filtered}, \code{error_cov_filtered} and \code{pr_t_given_t} are respectively 
##' time-varying filtered latent variable mean estimates, filtered error covariance matrix estimates, 
##' and filtered regime probability.
##' Note that if \code{theta.formula} is provided in \code{dynrModel@dynamics}, this assumes  that random effects are present in the dynamic equation. This would call an internal function to insert the random effect components as additional state variables. In this case, the last set of elements (rows) in \code{eta_smooth_final} would contain the estimated random effect components.
##' 
##' When \code{debug_flag} is TRUE, then additional information is passed into the cooked model. 
##' \code{eta_predicted}, \code{error_cov_predicted}, \code{innov_vec}, and \code{residual_cov} are respectively 
##' time-varying predicted latent variable mean estimates, predicted error covariance matrix estimates, the error/residual estimates (innovation vector),
##' and the error/residual covariance matrix estimates.
##' 
##' The exit flag given after optimization has finished is from the SLSQP optimizer.  Generally, error codes have negative values and successful codes have positive values.  However, codes 5 and 6 do not indicate the model converged, but rather simply ran out of iterations or time, respectively.  A more full description of each code is available at \url{https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values} and is also listed in the table below.
##' 
##' \tabular{lcl}{
##' NLOPT Term \tab Numeric Code \tab Description \cr
##' SUCCESS \tab 1 \tab Generic success return value. \cr
##' STOPVAL_REACHED \tab 2 \tab Optimization stopped because stopval (above) was reached. \cr
##' FTOL_REACHED \tab 3 \tab Optimization stopped because ftol_rel or ftol_abs (above) was reached. \cr
##' XTOL_REACHED \tab 4 \tab Optimization stopped because xtol_rel or xtol_abs (above) was reached. \cr
##' MAXEVAL_REACHED \tab 5 \tab Optimization stopped because maxeval (above) was reached. \cr
##' MAXTIME_REACHED \tab 6 \tab Optimization stopped because maxtime (above) was reached. \cr
##' FAILURE \tab -1 \tab Generic failure code. \cr
##' INVALID_ARGS \tab -2 \tab Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera). \cr
##' OUT_OF_MEMORY \tab -3 \tab Ran out of memory. \cr
##' ROUNDOFF_LIMITED \tab -4 \tab Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.) \cr
##' FORCED_STOP \tab -5 \tab Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization's nlopt_opt object opt from the user's objective function or constraints. \cr
##' NONFINITE_FIT \tab -6 \tab Fit function is not finite (i.e., is NA, NaN, Inf or -Inf). \cr
##' }
##' The last row of this table corresponding to an exit code of -6, is not from NLOPT, but rather is specific to the \code{dynr} package.
##' 
##' @return Object of class dynrCook.
##' 
##' @seealso 
##' \code{\link{autoplot}}, \code{\link{coef}}, \code{\link{confint}},
##' \code{\link{deviance}}, \code{\link{initialize}}, \code{\link{logLik}},
##' \code{\link{names}}, \code{\link{nobs}}, \code{\link{plot}}, \code{\link{print}},
##' \code{\link{show}}, \code{\link{summary}}, \code{\link{vcov}}.          
##' 
##' @examples
##' # Minimal model
##' require(dynr)
##' 
##' meas <- prep.measurement(
##' 	values.load=matrix(c(1, 0), 1, 2),
##' 	params.load=matrix(c('fixed', 'fixed'), 1, 2),
##' 	state.names=c("Position","Velocity"),
##' 	obs.names=c("y1"))
##' 
##' ecov <- prep.noise(
##' 	values.latent=diag(c(0, 1), 2),
##' 	params.latent=diag(c('fixed', 'dnoise'), 2),
##' 	values.observed=diag(1.5, 1),
##' 	params.observed=diag('mnoise', 1))
##' 
##' initial <- prep.initial(
##' 	values.inistate=c(0, 1),
##' 	params.inistate=c('inipos', 'fixed'),
##' 	values.inicov=diag(1, 2),
##' 	params.inicov=diag('fixed', 2))
##' 
##' dynamics <- prep.matrixDynamics(
##' 	values.dyn=matrix(c(0, -0.1, 1, -0.2), 2, 2),
##' 	params.dyn=matrix(c('fixed', 'spring', 'fixed', 'friction'), 2, 2),
##' 	isContinuousTime=TRUE)
##' 
##' data(Oscillator)
##' data <- dynr.data(Oscillator, id="id", time="times", observed="y1")
##' 
##' model <- dynr.model(dynamics=dynamics, measurement=meas,
##' 	noise=ecov, initial=initial, data=data)
##' 
##' # Now cook the model!
##' cook <- dynr.cook(model,
##' 	verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE)
dynr.cook <- function(dynrModel, conf.level=.95, infile, optimization_flag=TRUE, hessian_flag = TRUE, verbose=TRUE, weight_flag=FALSE, debug_flag=FALSE, perturb_flag=FALSE, ...) {
    frontendStart <- Sys.time()
    transformation=dynrModel@transform@tfun
    data <- dynrModel$data
    if(xor(dynrModel@verbose, verbose)){ # If model@verbose does not agree with dynr.cook@verbose
	  if(verbose){
        message("'verbose' argument to dynr.cook() function did not agree with 'verbose' model slot.\nUsing function argument: verbose = TRUE\n")
			}
      dynrModel@verbose <- verbose
      # Always use 'verbose' function argument but only say so when they disagree and verbose=TRUE.
    }
    
    #browser()
	
	# if (.hasSlot(dynrModel$dynamics, 'theta.formula') && length(dynrModel$dynamics@theta.formula) > 0 && dynrModel$dynamics$saem==FALSE){
		# #get the initial values of b and startvars
		# model <- EstimateRandomAsLVModel(dynrModel, optimization_flag, hessian_flag, verbose, weight_flag, debug_flag)
		# fitted_model <- dynr.cook(model, optimization_flag=optimization_flag, hessian_flag = hessian_flag, verbose=verbose, weight_flag=weight_flag, debug_flag=debug_flag)
		# coefEst <- coef(fitted_model)
		# estimated.names <- intersect(names(dynrModel@xstart), names(coefEst))
		# dynrModel@xstart[estimated.names] <- coefEst[estimated.names]
		
		# return(fitted_model)
	# }	#internalModelPrep convert dynrModel to a model list
	
	#  ------- The following lines obtain the necessary components of SAEM processs ----------------------------
	dots <- list(...)
	
	if (class(dynrModel$dynamics) == "dynrDynamicsFormula"){
	  saem <- dynrModel$dynamics$saem 
	} else {
	  saem <- FALSE #SAEM only supports dynrDynamicsFormula for now
	}
	
	if(length(dots) > 0){
	  #obtaining saem parameters
	  if('saemp' %in% names(dots)){
	    saemp <-dots$saemp
	  }
	}
	
    if(saem==TRUE){
		sigmab.names <- unique(as.vector(dynrModel$random.params.inicov))
		sigmab.names <- sigmab.names[!sigmab.names %in% c('fixed', '0')]
		
		
		lambda.names <- unique(as.vector(dynrModel@measurement$params.load[[1]]))
		lambda.names <- lambda.names[!lambda.names %in% c('fixed', '0')]
		lambda.names <- dynrModel@param.names[lambda.names]
		
		noise.names <- unique(as.vector(dynrModel@noise@params.observed[[1]]))
		noise.names <- noise.names[!noise.names %in% c('fixed', '0')]
		noise.names <- dynrModel@param.names[noise.names]
		
		if(length(dynrModel@measurement@params.int) > 0 && length(dynrModel@measurement@params.int[[1]]) > 0){
		  mu.names <- dynrModel@measurement@params.int[[1]]
		  mu.names <- dynrModel@param.names[mu.names]
		  mu.values <- dynrModel@measurement@values.int[[1]]
		}
		else{
		  mu.names <- logical(0)
		  mu.values <- logical(0)
		}
		
		
		#browser()
		theta.variables <- extractVariablesfromFormula(dynrModel$dynamics@theta.formula)
	    startval.names <- names(dynrModel$dynamics@startval)
        r =formula2design( 
        dynrModel$dynamics@theta.formula,
        covariates=c(data$covariate.names, "1"),
        random.names=c(dynrModel$dynamics$random.names, 'dummy_random_variable_4895746'),
		beta.names=startval.names[startval.names %in% theta.variables == FALSE])    
        r$fixed= as.matrix(r$fixed[ ,colnames(r$fixed) != '0'])	
		r$random=r$random[, seq_len(length(dynrModel$dynamics$random.names)), drop=FALSE]
		
		
		
		# organize the lowerbound and upperbound vectors for saem
		param.names <- logical(0)
		if(length(dynrModel@dynamics@beta.names) > 0)
		  param.names <- c(param.names, dynrModel@dynamics@beta.names)
		if(length(mu.names) > 0)
		  param.names <- c(param.names, mu.names)
		if(length(lambda.names) > 0)
		  param.names <- c(param.names, lambda.names)
		if(length(noise.names) > 0)
		  param.names <- c(param.names, noise.names)
		#if(length(sigmab.names) > 0)
		#  param.names <- c(param.names, sigmab.names)		
		print('param.names')
		print(c(param.names, sigmab.names))
			
		
		#substitute the values in xstart into the expression 
		#the values in xstart is already reverse transformed
		dSigmaede<-matrix(sapply(dynrModel@dSigmaede, function(x){eval(x, as.list(dynrModel@xstart[dynrModel@noise@params.observed[[1]]]))}), nrow=nrow(dynrModel@dSigmaede), ncol=ncol(dynrModel@dSigmaede))
	    dSigmaede2<-matrix(sapply(dynrModel@dSigmaede2, function(x){eval(x, as.list(dynrModel@xstart[dynrModel@noise@params.observed[[1]]]))}), nrow=nrow(dynrModel@dSigmaede2), ncol=ncol(dynrModel@dSigmaede2))
		#dSigmaede2 <- t(dSigmaede2)
		
		#comment out [the values should be substituted in dynr.cook]
		dSigmabdb<-matrix(sapply(dynrModel@dSigmabdb, function(x){eval(x, dynrModel@known.vars)}), nrow=nrow(dynrModel@dSigmabdb), ncol=ncol(dynrModel@dSigmabdb))
		dSigmabdb2<-matrix(sapply(dynrModel@dSigmabdb2, function(x){eval(x, dynrModel@known.vars)}), nrow=nrow(dynrModel@dSigmabdb2), ncol=ncol(dynrModel@dSigmabdb2))
		#dSigmabdb2 <- t(dSigmabdb2)
		
		
		
		
		num.x <- length(dynrModel@initial$params.inistate[[1]])
		#num.subj <- length(unique(data$original.data[['id']]))
		num.subj <- length(unique(data$id))
		# ******examined (not extended)
		if(dynrModel@freeIC == FALSE){
		  random.names <- dynrModel@dynamics@random.names}
		else{
		  x.names <- dynrModel@measurement@state.names
		  random.names <- c(dynrModel@dynamics@random.names, paste0('b_',x.names))
		  #print (random.names)
		}


		#setting lower and upper bounds
		#lower_bound <- c(dynrModel@lb[param.names], rep(dynrModel@dynamics@random.lb, length(dynrModel@dynamics@random.values.inicov)))
		#upper_bound <- c(dynrModel@ub[param.names], rep(dynrModel@dynamics@random.ub, length(dynrModel@dynamics@random.values.inicov)))
		#browser()
		lower_bound <- c(dynrModel@lb[param.names], dynrModel@lb[sigmab.names])
		upper_bound <- c(dynrModel@ub[param.names], dynrModel@ub[sigmab.names])
		print('Lower and Uppder bounds:')
		print(lower_bound)
		print(upper_bound)
		
		
		# #temporarily commented out 
		#browser()
		# #get the initial values of startvars
		model <- ExpandRandomAsLVModel(dynrModel)
		fitted_model <- dynr.cook(model, optimization_flag=optimization_flag, hessian_flag = hessian_flag, verbose=verbose, weight_flag=weight_flag, debug_flag=debug_flag)
		coefEst <- coef(fitted_model)
		estimated.names <- intersect(names(dynrModel@xstart), names(coefEst))
		dynrModel@xstart[estimated.names] <- coefEst[estimated.names]
		#log transformation
		dynrModel@xstart[noise.names] <- log(coefEst[noise.names])
		dynrModel@xstart[sigmab.names] <- log(coefEst[sigmab.names])
	    #par_value <- c(dynrModel@xstart[param.names], log(dynrModel@dynamics@random.values.inicov))
		par_value <- c(dynrModel@xstart[param.names], dynrModel@xstart[sigmab.names])
		print('Starting values:')
		print(par_value)
		

		
		# #get the initial values of b
		# #temporarily commented out 
		b <- fitted_model@eta_smooth_final[(num.x + 1):nrow(fitted_model@eta_smooth_final),data$tstart[1:num.subj+1]]
		b[ b < dynrModel@dynamics@random.lb | b > dynrModel@dynamics@random.ub ] = 0
		b <- matrix(b, nrow=num.subj, ncol=length(random.names))

		# #browser()
		# #trueb <- data$trueb[data$tstart[1:num.subj+1], ]
		
		 
		# # obtain y0 form eta_smooth_final		
		y0 <- matrix(0, nrow=num.subj, ncol=num.x)
		# #temporarily set y0 to be the values from SAEM
		for(i in 1:num.subj){
		  if(length(dynrModel@initial$values.inistate[[1]]) > 0){
		    if(dynrModel@freeIC){
              #y0[i, ] <- dynrModel@initial$values.inistate[[1]] + b[i, (1:num.x)]
			  y0[i, ] <- fitted_model@eta_smooth_final[1:num.x, data$tstart[i]+1] + b[i, (1:num.x)]
			}
			else{
			  #y0[i, ] <- dynrModel@initial$values.inistate[[1]] 
			  y0[i, ] <- fitted_model@eta_smooth_final[1:num.x, data$tstart[i]+1]
			}
		  }
		}
		dynrModel@initial@y0 <- list(y0)

		
		
		# temprarily put the values from matlab in SAEMFixedIC070708 here
		#y0<- matrix(c(2.8943, 2.7370, 1.0120, 1.8056, 1.6517, 2.4414, 1.3905 ,-0.0999, 1.7179, 2.3390, 1.4471, 1.8978, 2.9108, 0.7028, 2.5545, 1.5768, 3.2847, 0.6709, 1.2588, 0.7579, 1.2183, 2.5117, 1.1434, 2.4538, 2.3197, 1.8457, 1.3066, 3.0685, 2.0675, 2.4461, 2.0649, 1.8809, 1.4026, 3.1053,-0.0348,-1.1215, 1.6681, 1.9565, 2.6411, 1.8611, 0.4022, 1.8196, 2.3910, 2.4407, 2.7521, 0.4220, 2.2558, 2.4675, 1.9894, 1.4271, 1.4571, 2.2556, 1.4189, 3.5124, 3.7355, 0.2686, 1.8448, 2.5865, 2.1420, 1.3040, 2.3083, 0.7862, 2.7918, 0.5926, 2.9407, 1.6763, 1.5685, 0.5470, 2.1088 , 1.8713, 0.5298, 3.1700, 1.3991, 1.5801, 1.6109, 2.3327, 2.7215, 2.8108, 3.0780,-0.0760, 2.0749, 2.6799, -0.2947, 2.2974, 2.0846, 3.7754, 2.6616, 1.9206, 2.2201, 1.0827, 2.3168, 1.7449, 2.0020, 1.6479, 2.6502, 2.4951, 1.9619, 1.2890, 0.9410, 1.5389, 3.0821, 1.4035, 1.1950, 0.2468, 0.9075, 3.5126, 2.4237, 0.7293, 2.5594, 1.5843, 2.3243, 1.8114, 2.0866, 2.1080, 2.6914, 2.4199, 2.3911, 2.4123, 0.6661, 1.6337, 1.6988, 1.6348, 0.6006, 2.5466, 2.5683, 2.7618, 1.6805, 3.3926, 1.3844, 1.6340, 1.4094, 2.6897, 1.2832, 1.0769, 2.5619, 2.9300, 1.7991, 1.9967, 2.2336 , -1.0434, 3.3407, 3.4599, 2.0019, 2.6400, 2.2386, 1.2927, 3.8334, 2.9341, 2.5831, 2.1816, 3.2356, 1.1322, 1.7038, 3.0659, 1.8767, 0.5955, 2.6554, 1.3435, 2.4426, 2.7236, 1.7719, 1.9509, 3.1603, 2.5185, 2.1783, 2.2093, 1.4967, 1.6336, 0.0459, 3.1868, 0.0687, 2.2556, 3.3849, 2.9790, 2.3975, 0.0534, 1.3006, 3.0543, 3.5471, 2.0166, 1.6686, 2.0387, 3.4958, 2.5447, 0.1896, 1.9407, 2.5985, 2.0634, 2.0356, 2.0175, 2.4723, 1.0182, 1.2379, 1.5167, 1.8446, 1.7766, 0.7159, 2.6900, 2.0392, 2.8739, 0.4624, 0.0706, 0.6249, 0.9172, 0.0421, 0.0793, 0.1179, 0.7781, 0.1407, 0.0826, 0.2892, 0.4860, 0.3644, 0.6703, 0.9860, 0.7820, 1.3539, 1.0307, 1.0159, 0.5982, 0.6500, 0.2560, 0.4808, 0.7606 , -0.3351 ,  -0.1885, 0.4104, 0.2887, 0.0800, 1.1323, 0.3694, 0.8045, 0.5759, 0.9073, 0.4047, 0.1989, 0.0260, 0.2530, 0.5172, 0.3613, 0.9084, 0.4301  , -0.0026, 0.3471, 1.0066, 0.6358, -0.1921, 0.1958, 0.4531, 0.7984, 0.4938, 0.2504, 0.4097, 0.7020, 0.8752, 0.6686, 0.2209, 0.6495, 0.3543, 1.0411, 0.0954, 0.3475, 0.5533, 0.9810, 0.7804, 0.7344, 0.7817, 0.0859, 0.3810, 0.3284 ,  -0.0828, 0.6099, 0.3803, 1.0668, 0.4722, 0.7875, 0.3675, 0.1021, 0.4357, -0.0533, 0.3793, 0.0354, 0.4712 ,  -0.2189, 0.3164, 0.5974, 0.5500, 0.4033, 0.5692, 0.4681, 1.1319, 0.6142, 0.0010, 0.4135, 0.4866, 0.6566, 0.5413, 0.6134, 1.0547, 0.3003, 0.8044, 0.9241, 0.5420, 0.8155, 1.1716, 0.5340, 0.4557, 0.2982, 0.4690, 0.9057, 0.4100, 1.0455, 0.5135, 0.5274, 0.2402, 0.4637, 0.1206, 0.5328, 0.6013, 0.0912, 0.6672, 0.1116, 0.5839, 0.2063, 0.8864, 0.1042, 0.4523, 0.9145, 0.4296, 0.3689, 0.7752, 0.8477, 0.5880, 0.3177, 0.4780, 0.7139, 0.7702, 0.8109, 0.2058, 0.8134, 0.7243, 0.1311, 0.8030, 0.7413, 0.4636, 0.8412, 0.5593, 0.3297, 0.3947, 0.4607, 0.6301, 0.6924, 1.3025, 0.8481, 0.8351, 0.8834, 0.4811, 0.6937 ,  -0.2270, 0.5713, 0.4597, 0.3622, 0.3677, 0.4871, 0.3591, 0.7494, 0.1330, 0.4007, 0.0113, 0.5606  , -0.1537, 0.7687, 0.2914, 0.5620, 0.3941, 0.6841, 0.2825, 0.1608, 0.4387, 0.8990, 0.5343, 0.4610, 0.3212, 0.2341, 0.3733, 0.3787, 1.0213, 0.3553, 0.4321, 0.3016, 0.1585, 0.6842, 0.4024, 0.6357  , -0.1361, 0.6505, 0.8996, 0.4596, -0.2569, 0.7014), nrow=num.subj, ncol=2)
		#par_value <- c(2.8419, 0.4093, 0.4394, 0.0036, -0.0014, 0.0059, 0.6886, 1.1834, -0.5696, -0.6346, -0.5124, -0.4153)
		#b <- matrix(c(1.1850, 0.6200, 0.8920, -0.0760, 0.5380, 0.9960, -0.1630, 0.5200, 0.6320, 0.0090, -0.8080, 1.2790, 1.5300, -0.2990, 1.1670, 0.4940, -1.2110, 0.9190, 0.7740, 0.7450, 2.3640, 0.8410, -0.5520, 0.5650, 0.0420, -0.0300, 0.9950, 0.7480, 0.6650, 0.1870, 0.6150, 1.3900, -0.1210, -0.6650, 1.0790, 1.1170, 0.5330, 0.4460, 0.0330, 0.2170, 0.8430, -0.1900, 0.7420, 0.1350, 0.3140, -0.0970, 0.9610, 1.4400, -0.6710, 1.2880, 1.2530, -0.0430, 0.0760, 0.4850, 1.4280, 0.5080, -0.1390, 1.1410, 1.4140, 0.3750, 0.4100, 1.0910, 0.7150, 1.2100, 0.6210, 0.4240, 0.6050, 1.2950, 1.3430, 0.6810, 0.6870, 0.4300, 1.2550, 0.0770, 0.0930, -0.7150, 1.0900, -0.5960, 0.8770, 0.5820, -0.1530, 1.0810, 0.7110, 0.5100, -0.4280, 1.0670, 0.5040, 1.3410, 0.7250, 0.8800, 1.0480, 0.4260, 0.7650, 0.7900, 0.5500, 0.8900, 0.7720, 1.5110, 0.6160, 0.4310, 0.1830, 1.6040, 1.1030, -0.5870, 1.0000, 0.8730, -0.2200, -0.0400, 0.2140, -0.9810, 0.7750, -0.4560, 0.5980, 0.0680, -0.4450, 1.4110, 1.1790, 0.1020, -0.2590, 1.0690, 0.5990, -0.2500, 1.0690, 0.3150, -0.1610, 0.5600, 0.6510, 1.2860, -0.1670, 0.7960, -0.1330, 1.3080, 0.4050, -0.3330, 0.8720, -0.1340, 0.9650, 0.0300, 1.0180, 0.0250, 0.4090, -0.4880, -0.8230, 0.0530, -0.2430, 0.0780, -0.1240, 0.2130, 1.6210, 0.3360, -0.1960, 1.1630, 1.6210, 0.7640, -0.4350, 1.5660, 0.1600, 0.9960, 0.4440, -1.0530, 0.0420, 0.5490, 0.4050, 0.6250, -0.0350, 0.0440, -0.2690, 0.1460, 1.0430, 0.1700, 0.3240, 1.3610, 0.2630, 0.6130, 0.4270, -0.3090, 0.9600, 0.1610, 0.3410, -0.5640, -0.0060, 1.0900, 1.0840, 0.8850, 0.2480, 0.0180, -0.6580, -0.0910, 1.4970, -0.0100, 0.9250, 0.0120, -0.8450, 1.7840, 0.2570, 1.0730, 2.2950, 0.6700, 1.2220, 0.6740),nrow=num.subj, ncol=length(random.names))
		#b <- matrix(0,nrow=num.subj, ncol=length(random.names))
		
		

        #browser() 
        model <- internalModelPrepSAEM(
            num_regime=dynrModel@num_regime,
            dim_latent_var=dynrModel@dim_latent_var,
            xstart=par_value,
            ub=upper_bound,
            lb=lower_bound,
            options=dynrModel@options,
            isContinuousTime=dynrModel@dynamics@isContinuousTime,
            infile=dynrModel@outfile,
            outfile=gsub(".c\\>","",dynrModel@outfile),
            compileLib=dynrModel@compileLib,
            verbose=dynrModel@verbose,
            num_theta=length(dynrModel@dynamics@theta.names),
            num_beta=length(dynrModel@dynamics@beta.names),
			total_t = nrow(dynrModel@data$original.data[dynrModel@data$original.data[['id']] == 1, ]),
            num_lambda=length(lambda.names),
            num_mu=length(mu.names),
            num_random=length(random.names),
            theta.formula=dynrModel@dynamics@theta.formula,
            random.names=random.names,
            p0=as.vector(dynrModel@initial@values.inicov[[1]]),
            lambda=dynrModel@measurement$values.load[[1]], #column-major
            b= b,
			#trueb=trueb,
			random.lb=dynrModel@dynamics@random.lb,
			random.ub=dynrModel@dynamics@random.ub,
			bAdaptParams=as.matrix(saemp@bAdaptParams),
            KKO=saemp@KKO,
			scaleb=saemp@scaleb,
			MAXGIB = as.integer(saemp@MAXGIB), 
			MAXITER = as.integer(saemp@MAXITER), 
			maxIterStage1 = as.integer(saemp@maxIterStage1), 
			setScaleb = as.integer(saemp@setScaleb),
			setAccept = saemp@setAccept,
			gainpara = saemp@gainpara, 
			gainparb = saemp@gainparb, 
			gainpara1 = saemp@gainpara1, 
			gainparb1 = saemp@gainparb1,
			num_bpar = length(sigmab.names),
			sigmab = dynrModel$random.values.inicov,
			sigmae = dynrModel@noise$values.observed[[1]],
			mu = mu.values,
			lower_bound=as.matrix(lower_bound),
			upper_bound=as.matrix(upper_bound),
			dmudparMu=dynrModel@dmudparMu,
			dmudparMu2=dynrModel@dmudparMu2, # to fit the input of SAEM
			dSigmaede=dSigmaede,
			dSigmaede2=dSigmaede2,
			dLambdparLamb=dynrModel@dLambdparLamb,
			dLambdparLamb2=dynrModel@dLambdparLamb2,
			dSigmabdb = dSigmabdb,
			dSigmabdb2 = dSigmabdb2,
			time_=dynrModel@data$time,
			freeIC=dynrModel@freeIC,
			par_value=as.matrix(par_value),
			seed=saemp@seed,
			y0=y0,
			r=r
        )

		
        
        #libname <- model@libname
        #model@libname <- NULL
        
        
        
        model <- combineModelDataInformationSAEM(model, data)
        model <- preProcessModel(model)
        
        #browser()
		addr <- .C2funcaddressSAEMRcpp(isContinuousTime=model$isContinuousTime, infile=model$infile, outfile=model$outfile, verbose=model$verbose, compileLib = model$compileLib)
		model$func_address <- addr$address
		libname <- addr$libname
		
        
		#call Rcpp_saem_interface
		output <- rcpp_saem_interface(model, data, weight_flag, debug_flag, optimization_flag, hessian_flag, verbose)
		

        return(output)
    }

    
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
		outfile=gsub(".c\\>","",dynrModel@outfile),
		compileLib=dynrModel@compileLib,
		verbose=dynrModel@verbose
	)
	libname <- model$libname
	model$libname <- NULL
	
	model <- combineModelDataInformation(model, data)
	model <- preProcessModel(model)
	if(any(sapply(model$func_address, is.null.pointer))){
		warning("Found null pointer(s) in 'func_address' list. (Re-)compiling your functions...")
		if(missing(infile)){
			stop("Cannot compile your functions because 'infile' argument is missing.")
		}
		addr <- .C2funcaddress(isContinuousTime=model$isContinuousTime, infile=infile, verbose=verbose)
		model$func_address <- addr$address
		libname <- addr$libname
	}
	seed <- sample(1073741824L, size=1)
	gc()
	backendStart <- Sys.time()
	output <- .Call(.Backend, model, data, weight_flag, debug_flag, optimization_flag, hessian_flag, verbose, perturb_flag, seed, PACKAGE = "dynr")
	backendStop <- Sys.time()
	dyn.unload(libname) # unload the compiled library
	# unlink(libname) # deletes the DLL
	#gc() # garbage collection
	cat('Original exit flag: ', output$exitflag, '\n')
	# Check to make sure likelihood is not NaN.
	output$exitflag <- ifelse(!is.finite(output$neg.log.likelihood), -6, output$exitflag)
	# Use lookup table for exit flags
	cat('Modified exit flag: ', output$exitflag, '\n')
	cat(.ExitFlags[as.character(output$exitflag)], '\n')
	cat('Original fitted parameters: ', output$fitted.parameters, '\n', fill=TRUE)
	cat('Transformed fitted parameters: ', transformation(output$fitted.parameters), '\n', fill=TRUE)
	output2 <- endProcessing(output, transformation, conf.level)
	if (output$exitflag > 0 && output2$hessian.status == 0 && sum(output2$bad.standard.errors) == 0){
		cat('Successful trial\n')
	} else {
		if(hessian_flag){
			msg <- paste0("These parameters may have untrustworthy standard errors: ", paste(dynrModel$param.names[output2$bad.standard.errors], collapse=", "), ".")
			warning(msg, call.=FALSE)
		}
	}
	names(output2$transformed.parameters) <- dynrModel$param.names
	if(debug_flag){
		obj <- new("dynrDebug", output2)
	}else{
		obj <- new("dynrCook", output2)
	}
	
	obj@param.names <- dynrModel$param.names
	#populate transformed estimates to dynrModel
	#model<<-PopBackModel(model, obj@transformed.parameters)
	
	finalEqualStart <- model$xstart == obj@fitted.parameters
	if(any(finalEqualStart) && optimization_flag){
		warning(paste0("Some parameters were left at their starting values.\nModel might not be identified, need bounds, or need different starting values.\nParameters that were unmoved: ", paste(obj@param.names[finalEqualStart], collapse=", ", sep="")), call.=FALSE)
	}
	
	frontendStop <- Sys.time()
	totalTime <- frontendStop-frontendStart
	backendTime <- backendStop-backendStart
	frontendTime <- totalTime-backendTime
	obj@run.times <- as.double(c(totalTime=totalTime, backendTime=backendTime, frontendTime=frontendTime), units="secs")
	cat('Total Time:', totalTime, '\n')
	cat('Backend Time:', backendTime, '\n')
	rm(output2)
	rm(output)
	gc()
	return(obj)
}


failedProcessing <- function(x, transformation){
	cat('Failed trial\n')
	tParam <- transformation(x$fitted.parameters)
	x$transformed.parameters <- tParam
	
	nParam <- length(x$fitted.parameters)
	x$standard.errors <- rep(as.numeric(NA), nParam)
	x$standard.errors.untransformed <- rep(as.numeric(NA), nParam)
	x$transformed.inv.hessian <- matrix(as.numeric(NA), nrow=nParam, ncol=nParam)
	x$inv.hessian <- matrix(as.numeric(NA), nrow=nParam, ncol=nParam)
	x$conf.intervals <- matrix(as.numeric(NA), nrow=nParam, ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
	x$conf.intervals.endpoint.trans <- matrix(as.numeric(NA), nrow=nParam,ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
	x$bad.standard.errors <- rep(TRUE, nParam)
	x$hessian.status <- 1
	return(x)
}


  ##Tried but not working broadly enough for non-positive definite hessians
  ##Calculating a pseudovariance matrix = V'V,
  ##where V = GCHOL(H-), GCHOL(.) is the generalized Cholesky, and H- is
  ##the generalized inverse of the Hessian
  ##See the pseudo-inverse Hessian described in the following chapter, without importance sampling:
  ##Gill, J., et al. (2003). Numerical Issues Involved in Inverting Hessian Matrices. Numerical Issues in Statistical Computing for the Social Scientist, John Wiley & Sons, Inc.: 143-176.
  #Hinv <- MASS::ginv(x) #Generalized inverse of x
  #V <- sechol(Hinv)   #V <- ifelse("try-error" %in% class(sechol(Hinv)), TRUE, FALSE)
  #V1 <- t(V) %*% V   

#J%*%(ginv(x$hessian))t(J) and flag the negative diagonal elements

endProcessing <- function(x, transformation, conf.level){
	cat('Doing end processing\n')
	x <- checkHessian(x, transformation)
	confx <- qnorm(1-(1-conf.level)/2)
	
	if(!all(is.na(x$inv.hessian))){
		#Numerical Jacobian
		J <- numDeriv::jacobian(func=transformation, x=x$fitted.parameters) # N.B. fitted.parameters has the untransformed/uncontrained free parameters (i.e. log variances that can be negative).
		iHess0 <- J %*% (MASS::ginv(x$hessian.matrix)) %*% t(J)
		bad.SE <- is.na(diag(iHess0)) | diag(iHess0) < 0
		
		iHess <- J %*% x$inv.hessian %*% t(J)
		tSE <- sqrt(diag(iHess))
		tParam <- transformation(x$fitted.parameters)
		CI <- c(tParam - tSE*confx, tParam + tSE*confx)
		
		# EndPoint Transformation
		tSEalt <- sqrt(diag(x$inv.hessian))
		x$standard.errors.untransformed <- tSEalt
		CIalt <- c(transformation(x$fitted.parameters - tSEalt*confx), transformation(x$fitted.parameters + tSEalt*confx))
		x$conf.intervals.endpoint.trans <- matrix(CIalt, ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
		
		
		x$transformed.parameters <- tParam
		x$standard.errors <- tSE
		x$transformed.inv.hessian <- iHess
		x$conf.intervals <- matrix(CI, ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
		x$bad.standard.errors <- bad.SE
	}
	return(x)
}

checkHessian <- function(x, transformation){
	# if any of these warnings get tripped
	#  increase value of hessian.status
	#  Any value greater than 0 is bad
	failHess <- 0
	# If the log lik is not finite replace the non-computed Hessian with all NAs
	if(!is.finite(x$neg.log.likelihood)){
		x$hessian.matrix <- matrix(NA, nrow(x$hessian.matrix), ncol(x$hessian.matrix))
		x <- failedProcessing(x, transformation)
		return(x)
	}
	# Warn if any of the Hessian elements are NA
	if(any(!is.finite(x$hessian.matrix))){
		warning("Found infinite, NaN, or missing values in Hessian.  Add that to the list of things that are problematic:\nproblematic <- c('mankind instead of humankind', 'claims of Cherokee heritage', 'gender income gap', ..., 'your Hessian')", call.=FALSE)
		x <- failedProcessing(x, transformation)
		return(x)
	}
	diagH <- diag(x$hessian.matrix)
	zdiag <- diagH == 0
	if(any(zdiag)){
		#warning("The following diagonal elements of the Hessian were zero")
		# No warning?
		diagH[zdiag] <- 10e-14
		diag(x$hessian.matrix) <- diagH
	}
	if (is.positive.definite(x$hessian.matrix) || is.positive.definite2(x$hessian.matrix)){ #N.B. We use both is.positive.definite() and is.positive.definite2().  Not sure why?
		useHess <- x$hessian.matrix
	} else{
		failHess <- failHess + 1
		msg <- "Hessian is not positive definite. The standard errors were computed using the nearest positive definite approximation to the Hessian matrix."
		warning(msg, call.=FALSE)
		useHess <- (Matrix::nearPD(x$hessian.matrix, conv.norm.type = "F"))$mat
	}
	V1 <- try(solve(useHess), silent=TRUE)
	if("try-error" %in% class(V1)){
		failHess <- failHess + 1
		warning("Hessian is not invertible; used pseudo-inverse.\nModel might not be identified or is not at an optimal solution.\nRegard standard errors suspiciously.", call.=FALSE)
		V1 <- MASS::ginv(useHess)
	}
	x$inv.hessian <- V1
	x$hessian.status <- failHess
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

combineModelDataInformationSAEM <- function(model, data){
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

    
    #model$max_t <- max(matrix(data[['time']], nrow = model$total_t, ncol= model$num_sbj)[model$total_t, ])
    model$max_t <- max(data[['time']])
	
	#all possible time points
	#InfDS.tspan
	model$tspan = sort(unique(data[['time']]))
	model$num_time = length(sort(unique(data[['time']])))
	#print ('num_time')
	#print (model$num_time)
	#print(model$tspan) 
	
    
    model$allT <- rep(0, model$num_sbj)
	#dif <- matrix(0, model$total_t, ncol= model$num_sbj)
	model$delt <- 2147493647
	model$tobs <- c()
	model$Y <- matrix(0,nrow=length(data$observed.names), ncol=length(data[['time']]))
	
	for(j in 1:length(data$observed.names)){
      model$Y[j,] <- data$observed[[j]]
    }
	
	#print(model$Y[3,11:20])
    for (i in sort(unique(data[['id']]))){
        table_i <- data$original.data[data$original.data[['id']] == i, ]
        model$allT[i] <- nrow(table_i)
		#InfDS.tobs{i}
		pos = sort(unlist(lapply(table_i[['time']], 
		  function(x, all_time){
            pos = match(x, all_time)
            return(pos)
          }, model$tspan)))
		model$tobs <- c(model$tobs, pos)
			
		for(j in 2:nrow(table_i)){
            #dif[j,i] <- table_i[j, 'time'] - table_i[j - 1, 'time']
			if(table_i[j, 'time'] - table_i[j - 1, 'time'] < model$delt)
				model$delt <- table_i[j, 'time'] - table_i[j - 1, 'time'] 
        }
    }
    #print(model$tobs)
	#print(model$allT) #correct here
	#print(model$delt) #correct here
	

	
	#print('H & Z')
	r = model$r
    Z= apply(r$random, 1, as.numeric)
    H = matrix(nrow=0, ncol=0)
    for (line in c(1: model$num_sbj)){
        U = c(1:length(data$covariate.names))
        for (u in c(1:length(data$covariate.names))){
            U[u] = data$original.data[((line-1)*model$total_t+1),data$covariate.names[u]]
        }
        
        
        temp <- r$fixed
        
        for (i in c(1:nrow(r$fixed))){
            for (j in c(1:ncol(r$fixed))){
                for (u in c(1:length(data$covariate.names))){
                    if(r$fixed[i,j] %in% data$covariate.names[u]){
                        r$fixed[i,j] <- U[u]
                    }
                }
            }
        }

        if(nrow(H) == 0 && ncol(H) == 0){
            H = t(apply(r$fixed, 1, as.numeric))
        } else{
            H=rbind(H,t(apply(r$fixed, 1, as.numeric)))
        }
        r$fixed <- temp
    }
    model$H <- as.matrix(H)
    model$Z <- as.matrix(Z)
	print(H)
    

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

is.positive.definite2 <- function(x) {
	"matrix" %in% class(try(MASS::ginv(x), silent=TRUE))
}

# From https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#Return_values
.ExitFlags <- c(
	'-5'='Optimization halted because of a forced termination.',
	'-4'='Optimization halted because of roundoff errors.',
	'-3'='Optimization failed. Ran out of memory.',
	'-2'='Optimization halted. Lower bounds are bigger than upper bounds.',
	'-1'='Optimization failed. Check starting values.',
	'0'='Optimization has been turned off.',
	'1'='Optimization terminated successfully',
	'2'='Optimization stopped because objective function reached stopval',
	'3'='Optimization terminated successfully: ftol_rel or ftol_abs was reached.',
	'4'='Optimization terminated successfully: xtol_rel or xtol_abs was reached.',
	'5'='Maximum number of function evaluations reached. Increase maxeval or change starting values.',
	'6'='Maximum optimization time reached. Increase maxtime or change starting values.',
	'-6'='Likelihood function is NaN and could not find a way out. Optimizer gave up but is not at a converged optimum.')

PopBackMatrix <- function(values.matrix, param.matrix, trans.parameters){
	if ("list" %in% class(values.matrix)){
		num_regime <- length(values.matrix)
		if (num_regime > 0){
			for (i in 1:num_regime){
				values.matrix[[i]][which(param.matrix[[i]] !=0 , arr.ind = TRUE)] <- trans.parameters[param.matrix[[i]][which(param.matrix[[i]] !=0 , arr.ind = TRUE)]]
			}
		}
	}else{
		values.matrix[which(param.matrix !=0 , arr.ind = TRUE)] <- trans.parameters[param.matrix[which(param.matrix != 0, arr.ind = TRUE)]]
	}
	return(values.matrix)
}

PopBackFormula <- function(formula, paramnames, param.names, trans.parameters){
	string <- paste0(deparse(formula, width.cutoff = 500L), collapse="")
	for (i in 1:length(paramnames)){
		string <- gsub(paste0("param\\[", match(paramnames[i], param.names, nomatch=0)-1,"\\]"), trans.parameters[match(paramnames[i], param.names, nomatch=0)], string, perl = TRUE)
	}
	eval(parse(text=string))
}

PopBackModel <- function(dynrModel, trans.parameters){
	if (class(dynrModel$dynamics) == 'dynrDynamicsFormula'){
		if (length(dynrModel@dynamics@paramnames) > 0){
			dynrModel@dynamics@formula <- PopBackFormula(dynrModel@dynamics@formula, dynrModel@dynamics@paramnames, dynrModel@param.names, trans.parameters)
		}
	}else{
		dynrModel@dynamics@values.dyn <- PopBackMatrix(dynrModel@dynamics@values.dyn, dynrModel@dynamics@params.dyn, trans.parameters)
		dynrModel@dynamics@values.exo <- PopBackMatrix(dynrModel@dynamics@values.exo, dynrModel@dynamics@params.exo, trans.parameters)
		dynrModel@dynamics@values.int <- PopBackMatrix(dynrModel@dynamics@values.int, dynrModel@dynamics@params.int, trans.parameters)
	}
	
	dynrModel@measurement@values.load <- PopBackMatrix(dynrModel@measurement@values.load, dynrModel@measurement@params.load, trans.parameters)
	dynrModel@measurement@values.exo <- PopBackMatrix(dynrModel@measurement@values.exo, dynrModel@measurement@params.exo, trans.parameters)
	dynrModel@measurement@values.int <- PopBackMatrix(dynrModel@measurement@values.int, dynrModel@measurement@params.int, trans.parameters)
	
	dynrModel@noise@values.latent <- PopBackMatrix(dynrModel@noise@values.latent, dynrModel@noise@params.latent, trans.parameters)
	dynrModel@noise@values.observed <- PopBackMatrix(dynrModel@noise@values.observed, dynrModel@noise@params.observed, trans.parameters)
	
	dynrModel@initial@values.inistate <- PopBackMatrix(dynrModel@initial@values.inistate, dynrModel@initial@params.inistate , trans.parameters) 
	dynrModel@initial@values.inicov <- PopBackMatrix(dynrModel@initial@values.inicov, dynrModel@initial@params.inicov, trans.parameters) 
	dynrModel@initial@values.regimep <- PopBackMatrix(dynrModel@initial@values.regimep, dynrModel@initial@params.regimep, trans.parameters)
	
	dynrModel@regimes@values<-PopBackMatrix(dynrModel@regimes@values, dynrModel@regimes@params, trans.parameters)
	
	# process model matrices to re-extract start values
	if(length(dynrModel$transform$inv.tfun.full) > 0 && is.numeric(trans.parameters)){
		trans.parameters <- dynrModel$transform$inv.tfun.full(trans.parameters)
		dynrModel@xstart <- trans.parameters
	}
	
	if(all(is.numeric(trans.parameters))){
		# Reset startval slots
		dynrModel@dynamics@startval <- trans.parameters[names(trans.parameters) %in% dynrModel@dynamics@paramnames]
		dynrModel@measurement@startval <- trans.parameters[names(trans.parameters) %in% dynrModel@measurement@paramnames]
		dynrModel@noise@startval <- trans.parameters[names(trans.parameters) %in% dynrModel@noise@paramnames]
		dynrModel@initial@startval <- trans.parameters[names(trans.parameters) %in% dynrModel@initial@paramnames]
		dynrModel@regimes@startval <- trans.parameters[names(trans.parameters) %in% dynrModel@regimes@paramnames]
		
		# Reset parameters numbers back to parameter names
		dynrModel@dynamics <- paramName2Number(dynrModel@dynamics, names(trans.parameters), invert=TRUE)
		dynrModel@measurement <- paramName2Number(dynrModel@measurement, names(trans.parameters), invert=TRUE)
		dynrModel@noise <- paramName2Number(dynrModel@noise, names(trans.parameters), invert=TRUE)
		dynrModel@initial <- paramName2Number(dynrModel@initial, names(trans.parameters), invert=TRUE)
		dynrModel@regimes <- paramName2Number(dynrModel@regimes, names(trans.parameters), invert=TRUE)
	}
	
	return(dynrModel)
}

#Computes the Sechol-Schnabel cholesky factorization
#See Sechol-Schnabel, R. B. and Eskow, E. 1990. "A New Modified Cholesky Factorization." SIAM Journal of Scientific Statistical Computing 11, 1136-58.
sechol <- function(A, tol = .Machine$double.eps, silent= TRUE ) {
	if (is.complex(A)) {
		warning("complex matrices not permitted at present")
		return(NULL)
	} else if (!is.numeric(A)) {
		warning("non-numeric argument to 'sechol'")
	return(NULL)
	}
	if (is.matrix(A)) {
		if (nrow(A) != ncol(A)) {
			warning("non-square matrix in 'sechol'")
			return(NULL)
		}
	} else {
		if (length(A) != 1) {
			warning("non-matrix argument to 'sechol'")
			return(NULL)
		}
		if (A>0) {
			return(as.matrix(sqrt(A)))
		} 
		warning("the leading minor of order 1 is not positive definite")
		return(NULL)
	}
	n <- nrow(A)
	L <- matrix(rep(0,n*n),ncol=ncol(A))
	tau <- tol ^(1/3) # made to match gauss
	gamm <- max(A)
	deltaprev <- 0
	Pprod <- diag(n)
	if (n > 2) {
		for (k in 1:(n-2)) {
			if( (min(diag(A[(k+1):n,(k+1):n]) - A[k,(k+1):n]^2/A[k,k]) < tau*gamm) 
					&& (min(svd(A[(k+1):n,(k+1):n])$d)) < 0) {
				dmax <- order(diag(A[k:n,k:n]))[(n-(k-1))]
				if (A[(k+dmax-1),(k+dmax-1)] > A[k,k]) {
					if (!silent) {
						print(paste("iteration:",k,"pivot on:",dmax,"with absolute:",(k+dmax-1)))
					}
					P <- diag(n)
					Ptemp <- P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
					A <- P%*%A%*%P
					L <- P%*%L%*%P
					Pprod <- P%*%Pprod
				}
				g <- rep(0,length=(n-(k-1)))
				for (i in k:n)	{
					if (i == 1) sum1 <- 0
					else sum1 <- sum(abs(A[i,k:(i-1)]))
					if (i == n) sum2 <- 0
					else sum2 <- sum(abs(A[(i+1):n,i]))
					g[i-(k-1)] <- A[i,i] - sum1 - sum2
				}
				gmax <- order(g)[length(g)]
				if (gmax != k)	{
					if (!silent) {
						print(paste("iteration:",k,
									"gerschgorin pivot on:", gmax, "with absolute:", (k+gmax-1)))
					}
					P <- diag(ncol(A))
					Ptemp <-	P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
					A <- P%*%A%*%P
					L <- P%*%L%*%P
					Pprod <- P%*%Pprod
				}
				normj <- sum(abs(A[(k+1):n,k]))
				delta <- max(0,deltaprev,-A[k,k]+max(normj,tau*gamm))
				if (delta > 0)	{
					A[k,k] <- A[k,k] + delta
					deltaprev <- delta
				}
			}
			
			L[k,k] <- A[k,k] <- sqrt(A[k,k])
			for (i in (k+1):n)	{
				L[i,k] <- A[i,k] <- A[i,k]/L[k,k]
				A[i,(k+1):i] <- A[i,(k+1):i] - L[i,k]*L[(k+1):i,k]
				if(A[i,i] < 0) A[i,i] <- 0
			}
		}
	}
	A[(n-1),n] <- A[n,(n-1)]
	eigvals <- eigen(A[(n-1):n,(n-1):n])$values
	delta <- max(0,deltaprev,
							 -min(eigvals)+tau*max((1/(1-tau))*(max(eigvals)-min(eigvals)),gamm))
	if (delta > 0)	{
		if (!silent) {
			print(paste("delta:",delta))
		}
		A[(n-1),(n-1)] <- A[(n-1),(n-1)] + delta
		A[n,n] <- A[n,n] + delta
		deltaprev <- delta
	}
	L[(n-1),(n-1)] <- A[(n-1),(n-1)] <- sqrt(A[(n-1),(n-1)])
	L[n,(n-1)] <- A[n,(n-1)] <- A[n,(n-1)]/L[(n-1),(n-1)]
	L[n,n] <- A[n,n] <- sqrt(A[n,n] - L[n,(n-1)]^2)
	
	r <- t(Pprod)%*%t(L)%*%t(Pprod)
	attr(r,"delta")=delta
	return(r)
}

