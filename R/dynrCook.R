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
##' # Let cookedModel be the output from dynr.cook
##' #coef(cookedModel)
coef.dynrCook <- function(object, ...){
    object@transformed.parameters
}

#`coef<-.dynrCook` <- function(object, value){
#   object <- PopBackModel(object, value)
#   return(object)
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
##' # Let cookedModel be the output from dynr.cook
##' #logLik(cookedModel)
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
##' # Let cookedModel be the output from dynr.cook
##' #nobs(cookedModel)
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
vcov.dynrCook <- function(object, ...){
    nm <- names(coef(object))
    rt <- object@transformed.inv.hessian
    dimnames(rt) <- list(nm, nm)
    return(rt)
}

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
##' @param ... further names arguments. Ignored.
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
##' # Let cookedModel be the output from dynr.cook
##' #confint(cookedModel)
confint.dynrCook <- function(object, parm, level = 0.95, type = c("delta.method","endpoint.transformation"), transformation= NULL, ...){
  type<-match.arg(type)
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
      tSEalt<-sqrt(diag(object$inv.hessian))
      CI <- matrix(c(transformation(object$fitted.parameters - tSEalt*confx), transformation(object$fitted.parameters + tSEalt*confx)),ncol=2)
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
##' 
##' When \code{debug_flag} is TRUE, then additional information is passed into the cooked model. 
##' \code{eta_predicted}, \code{error_cov_predicted}, \code{innov_vec}, and \code{residual_cov} are respectively 
##' time-varying predicted latent variable mean estimates, predicted error covariance matrix estimates, the error/residual estimates (innovation vector),
##' and the error/residual covariance matrix estimates.
##' 
##' The exit flag given after optimization has finished is from the SLSQP optimizer.  Generally, error codes have negative values and successful codes have positive values.  However, codes 5 and 6 do not indicate the model converged, but rather simply ran out of iterations or time, respectively.  A more full description of each code is available at \url{http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference#Return_values} and is also listed in the table below.
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
##' }
##' 
##' @seealso 
##' \code{\link{autoplot}}, \code{\link{coef}}, \code{\link{confint}},
##' \code{\link{deviance}}, \code{\link{initialize}}, \code{\link{logLik}},
##' \code{\link{names}}, \code{\link{nobs}}, \code{\link{plot}}, \code{\link{print}},
##' \code{\link{show}}, \code{\link{summary}}, \code{\link{vcov}}.          
##' 
##' @examples
##' #fitted.model <- dynr.cook(model)
dynr.cook <- function(dynrModel, conf.level=.95, infile, optimization_flag=TRUE, hessian_flag = TRUE, verbose=TRUE, weight_flag=FALSE, debug_flag=FALSE, saem=FALSE, ...) {

	
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
    
	dots <- list(...)
	if(length(dots) > 0){
	  saemp <-dots$saemp
	  #print(saemp)
	}
    if(saem==TRUE){
		#InfDS.Sigmab					
		#print(dynrModel$random.values.inicov)

		#InfDS.bpar
		sigmab.names <- unique(as.vector(dynrModel$random.params.inicov))
		sigmab.names <- sigmab.names[!sigmab.names %in% c('fixed', '0')]
		#print(sigmab.names)
		#num_bpar <- length(sigmab.names)
		
		lambda.names <- unique(as.vector(dynrModel@measurement$params.load[[1]]))
		lambda.names <- lambda.names[!lambda.names %in% c('fixed', '0')]
		
		noise.names <- unique(as.vector(dynrModel@noise@params.observed[[1]]))
		noise.names <- noise.names[!noise.names %in% c('fixed', '0')]
		
					
		# [TODO] to be generalized: bzeta			
        lb <- c(dynrModel@lb, 'b_zeta'=NA)
		ub <- c(dynrModel@ub, 'b_zeta'=NA)
		
		lower_bound <- as.vector(c(lb[dynrModel@dynamics@beta.names], lb[dynrModel@measurement@params.int[[1]]], lb[lambda.names], lb[noise.names], lb[sigmab.names]))
		upper_bound <- as.vector(c(ub[dynrModel@dynamics@beta.names], ub[dynrModel@measurement@params.int[[1]]], ub[lambda.names], ub[noise.names], ub[sigmab.names]))
		
		#print(lower_bound) 

		
		#[todo] put the b's initial code estimate
        #get the initial values of b
        b <- t(diag(3)%*%matrix(rnorm(600), nrow = 3, ncol=200))
        bEst <- getInitialVauleOfRandomEstimate(dynrModel) 
		
		
		#b, y0
		num.x <- length(model$initial$params.inistate[[1]])
		num.subj <- length(unique(data$original.data[['id']]))
		random.names <- model$dynamics@random.names
		y0 <- matrix(0, nrow=num.subj, ncol=num.x)
		b <- matrix(rnorm((num.subj*length(random.names))), nrow=num.subj, ncol=length(random.names))
		b_x <-  length(random.names) - num.x
		#print(b)
		for(i in 1:num.subj){
		  for (j in 1:length(random.names)){
			if(b[i, j] < model$dynamics@random.lb | b[i, j] > model$dynamics@random.ub){
			  b[i, j] <- 0
			}
		  }
		  for(j in 1:num.x){
			if(length(model$initial$values.inistate[[1]]) > 0){
			  y0[i, j] <- model$initial$values.inistate[[1]][j, 1] + b[i,(b_x+j)]
			}
		  }
		}
		model$initial@y0 <- list(y0)
		#print (inputs$initial@y0)

         
        model <- internalModelPrepSAEM(
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
            verbose=dynrModel@verbose,
            num_theta=length(dynrModel@dynamics@theta.names),
            num_beta=length(dynrModel@dynamics@beta.names),
            #total_t=length(unique(dynrModel@data$time)),
			total_t = nrow(dynrModel@data$original.data[dynrModel@data$original.data[['id']] == 1, ]),
            num_lambda=length(lambda.names),
            num_mu=length(dynrModel@measurement@params.int[[1]]),
            num_random=length(dynrModel@dynamics@random.names),
            theta.formula=dynrModel@dynamics@theta.formula,
            random.names=dynrModel@dynamics@random.names,
            p0=as.vector(dynrModel@initial@values.inicov[[1]]),
            lambda=as.vector(dynrModel@measurement$values.load[[1]]), #column-major
            b= b,
			random.lb=dynrModel@dynamics@random.lb,
			random.ub=dynrModel@dynamics@random.ub,
			bAdaptParams=saemp@bAdaptParams,
            KKO=saemp@KKO,
			MAXGIB = as.integer(saemp@MAXGIB), 
			MAXITER = as.integer(saemp@MAXITER), 
			maxIterStage1 = as.integer(saemp@maxIterStage1), 
			gainpara = saemp@gainpara, 
			gainparb = saemp@gainparb, 
			gainpara1 = saemp@gainpara1, 
			gainparb1 = saemp@gainparb1,
			num_bpar = length(sigmab.names),
			sigmab = dynrModel$random.values.inicov,
			mu = model@measurement@values.int[[1]],
			lower_bound=lower_bound,
			upper_bound=upper_bound,
			dmudparMu=model@measurement@dmudparMu,
			dmudparMu2=model@measurement@dmudparMu2
        )
        #print(model$MAXGIB)
		#print(model$gainpara)
        
        libname <- model$libname
        model$libname <- NULL
        
        
        
        model <- combineModelDataInformationSAEM(model, data)
        model <- preProcessModel(model)
        
        #print(model$b)
        
        # the following code compiles the generated file of dynr.model, comment out now for testing
        # some related comments also be commented in dynrModelInternal.R
        # if(any(sapply(model$func_address, is.null.pointer))){
            # warning("Found null pointer(s) in 'func_address' list. (Re-)compiling your functions...")
            # if(missing(infile)){
                # stop("Cannot compile your functions because 'infile' argument is missing.")
            # }
            # addr <- .C2funcaddressSAEM(isContinuousTime=model$isContinuousTime, infile=infile, verbose=verbose)
            # print('C2funcaddressSAEM done')
            # model$func_address <- addr$address
            # libname <- addr$libname
        # }
        # gc()
        # backendStart <- Sys.time()
        
        
        output <- .Call(.BackendS, model, data, weight_flag, debug_flag, optimization_flag, hessian_flag, verbose, PACKAGE = "dynr")

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
    gc()
    backendStart <- Sys.time()
    
    
    output <- .Call(.Backend, model, data, weight_flag, debug_flag, optimization_flag, hessian_flag, verbose, PACKAGE = "dynr")
    
    
    backendStop <- Sys.time()
    dyn.unload(libname) # unload the compiled library
    # unlink(libname) # deletes the DLL
    #gc() # garbage collection
    cat('Original exit flag: ', output$exitflag, '\n')
    # Check to make sure likelihood is not NaN.
    output$exitflag <- ifelse(is.na(output$neg.log.likelihood), -6, output$exitflag)
    # Use lookup table for exit flags
    cat('Modified exit flag: ', output$exitflag, '\n')
    cat(.ExitFlags[as.character(output$exitflag)], '\n')
    
    diagH = diag(output$hessian.matrix)
    diagH[diagH==0] = 10e-14
    diag(output$hessian.matrix) = diagH
    cat('Original fitted parameters: ', output$fitted.parameters, '\n', fill=TRUE)
    cat('Transformed fitted parameters: ', transformation(output$fitted.parameters), '\n', fill=TRUE)
    # if any of the Hessian elements are non-finite or the overall Hessian is non-positive definite
    #  set status = 0, otherwise it is 1
    nonfiniteH <- any(!is.finite(output$hessian.matrix))
    nonpdH <- !is.positive.definite2(output$hessian.matrix)
    status = ifelse(nonfiniteH || nonpdH, 0, 1)
    output2 <- endProcessing(output, transformation, conf.level)
    if (output$exitflag > 0 && status==1 && length(dynrModel$param.names[output2$bad.standard.errors])==0){
        cat('Successful trial\n')
    } else {
        #cat('Check the hessian matrix from your dynr output. \n')
        #cat('Hessian Matrix:',  '\n')
        #print(output$hessian.matrix)
        cat('\n')
        if(nonfiniteH){
            msg <- paste("Non-finite values in the Hessian matrix.")
            warning(msg)
        } else if(nonpdH || length(dynrModel$param.names[output2$bad.standard.errors]) > 0){
            msg <- "Hessian is not positive definite. The standard errors were computed using the nearest positive definite approximation to the Hessian matrix."
            if (length(dynrModel$param.names[output2$bad.standard.errors]) > 0){
                msg <- paste(c(msg,"These parameters may have untrustworthy standard errors: ", paste(dynrModel$param.names[output2$bad.standard.errors],collapse=", "),"."),collapse="")  
            }
            warning(msg)
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
        warning(paste0("Some parameters were left at their starting values.\nModel might not be identified, need bounds, or need different starting values.\nParameters that were unmoved: ", paste(obj@param.names[finalEqualStart], collapse=", ", sep="")))
    }
    
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
    x$transformed.inv.hessian <- matrix(999, nrow=nParam,ncol=nParam)
    x$conf.intervals <- matrix(999, nrow=nParam,ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
    return(x)
}


  ##Tried but not working broadly enough for non-positive definite hessians
  ##Calculating a pseudovariance matrix = V'V,
  ##where V = GCHOL(H-), GCHOL(.) is the generalized Cholesky, and H- is
  ##the generalized inverse of the Hessian
  ##See the pseudo-inverse Hessian described in the following chapter, without importance sampling:
  ##Gill, J., et al. (2003). Numerical Issues Involved in Inverting Hessian Matrices. Numerical Issues in Statistical Computing for the Social Scientist, John Wiley & Sons, Inc.: 143-176.
  #Hinv <- MASS::ginv(x) #Generalized inverse of x
  #V <- sechol(Hinv)   #V <- ifelse(class(sechol(Hinv)) == "try-error", TRUE, FALSE)
  #V1 <- t(V) %*% V   

#J%*%(ginv(x$hessian))t(J) and flag the negative diagonal elements

# TODO adjust endProcessing logic/workflow.  We're computing some things twice, and sometimes in different ways.

endProcessing <- function(x, transformation, conf.level){
    cat('Doing end processing\n')
    confx <- qnorm(1-(1-conf.level)/2)
    if (is.positive.definite(x$hessian.matrix)){ #N.B. We use is.positive.definite() here, but is.positive.definite2() above.  Why?
        useHess <- x$hessian.matrix
    }
    else{
        useHess <- (Matrix::nearPD(x$hessian.matrix, conv.norm.type = "F"))$mat
    }
    V1 <- try(solve(useHess))
    if(class(V1) == "try-error"){
        warning("Hessian is not invertible; used pseudo-inverse.\nModel might not be identified or is not at an optimal solution.\nRegard standard errors suspiciously.")
        V1 <- MASS::ginv(useHess)
    }
    
    #Identifies too many problematic parameters
    #evector <- eigen(x$hessian.matrix)$vectors
    #bad.values <- eigen(x$hessian.matrix)$values < 0
    #bad.evec <-  evector[,bad.values]
    #bad.evec <- apply(bad.evec,2,function(x){abs(x/sum(x))}) #normalize within column
    #bad.evecj <- apply(bad.evec,2,function(x){x > .5})#For each column find the substantial row (param) entries  
    #bad.SE <- apply(bad.evecj,1,function(x){ifelse(length(x[x=="TRUE"]) > 0, TRUE,FALSE)}) #Flag parameters that have been identified as problematic at least once
    
    #Numerical Jacobian
    J <- numDeriv::jacobian(func=transformation, x=x$fitted.parameters) # N.B. fitted.parameters has the untransformed/uncontrained free parameters (i.e. log variances that can be negative).
    iHess0 <- J %*% (MASS::ginv(x$hessian)) %*% t(J)
    bad.SE <- diag(iHess0) < 0
    
    iHess <- J %*% V1 %*% t(J)
    tSE <- sqrt(diag(iHess))
    tParam <- transformation(x$fitted.parameters) #Can do
    CI <- c(tParam - tSE*confx, tParam + tSE*confx)
    
    # EndPoint Transformation
    tSEalt<-sqrt(diag(V1))
    x$standard.errors.untransformed <- tSEalt
    CIalt <- c(transformation(x$fitted.parameters - tSEalt*confx), transformation(x$fitted.parameters + tSEalt*confx))
    x$conf.intervals.endpoint.trans <- matrix(CIalt, ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
    
    
    x$transformed.parameters <- tParam #Can do
    x$standard.errors <- tSE
    x$inv.hessian <- V1
    x$transformed.inv.hessian <- iHess
    x$conf.intervals <- matrix(CI, ncol=2, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
    x$bad.standard.errors <- bad.SE
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
	
    
    model$allT <- rep(0, model$num_sbj)
	#dif <- matrix(0, model$total_t, ncol= model$num_sbj)
	model$delt <- 2147493647
    for (i in unique(data[['id']])){
        table_i <- data$original.data[data$original.data[['id']] == i, ]
        model$allT[i] <- nrow(table_i)
		#print(nrow(table_i))
		for(j in 2:nrow(table_i)){
            #dif[j,i] <- table_i[j, 'time'] - table_i[j - 1, 'time']
			if(table_i[j, 'time'] - table_i[j - 1, 'time'] < model$delt)
				model$delt <- table_i[j, 'time'] - table_i[j - 1, 'time'] 
        }
    }
	#print(dif)
	#model$delt <- min(dif[2:model$total_t, ], na.omit = TRUE)
    #print(model$allT) #correct here
	#print(model$delt) #correct here

    
    #H & Z 
    r =formula2design( 
        model$theta.formula,
        #model$theta.formula[[2]],
        #model$theta.formula[[3]],
        covariates=c(data$covariate.names, "1"),
        random.names=model$random.names)    
    r$fixed= r$fixed[ ,colnames(r$fixed) != '0']
    
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
    model$H <- H
    model$Z <- Z
    model$tspan <- unique(data[['time']])
	#print(H)
    

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
  class(try(MASS::ginv(x),silent=TRUE))=="matrix"
}

# From http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference#Return_values
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
    if (length(dynrModel@dynamics@paramnames) > 0){
        dynrModel@dynamics@formula<-PopBackFormula(dynrModel@dynamics@formula,dynrModel@dynamics@paramnames,dynrModel@param.names,trans.parameters) 
    }
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
  
  # process model matrices to re-extract start values
  if(length(dynrModel$transform$inv.tfun.full) > 0 && is.numeric(trans.parameters)){
    trans.parameters <- dynrModel$transform$inv.tfun.full(trans.parameters)
    dynrModel@xstart <- trans.parameters
  }
  return(dynrModel)
}

#Computes the Sechol-Schnabel cholesky factorization
#See Sechol-Schnabel, R. B. and Eskow, E. 1990. "A New Modified Cholesky Factorization." SIAM Journal of Scientific Statistical Computing 11, 1136-58.
sechol <- function(A, tol = .Machine$double.eps, silent= TRUE )  {
  if (is.complex(A))  {
    warning("complex matrices not permitted at present")
    return(NULL)
  } else if (!is.numeric(A))  {
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
  tau <- tol ^(1/3)  # made to match gauss
  gamm <- max(A)
  deltaprev <- 0
  Pprod <- diag(n)
  if (n > 2)  {
    for (k in 1:(n-2))  {
      if( (min(diag(A[(k+1):n,(k+1):n]) - A[k,(k+1):n]^2/A[k,k]) < tau*gamm) 
          && (min(svd(A[(k+1):n,(k+1):n])$d)) < 0) {
        dmax <- order(diag(A[k:n,k:n]))[(n-(k-1))]
        if (A[(k+dmax-1),(k+dmax-1)] > A[k,k])  {
          if (!silent) {
            print(paste("iteration:",k,"pivot on:",dmax,"with absolute:",(k+dmax-1)))
          }
          P <- diag(n)
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        g <- rep(0,length=(n-(k-1)))
        for (i in k:n)  {
          if (i == 1) sum1 <- 0
          else sum1 <- sum(abs(A[i,k:(i-1)]))
          if (i == n) sum2 <- 0
          else sum2 <- sum(abs(A[(i+1):n,i]))
          g[i-(k-1)] <- A[i,i] - sum1 - sum2
        }
        gmax <- order(g)[length(g)]
        if (gmax != k)  {
          if (!silent) {
            print(paste("iteration:",k,
                        "gerschgorin pivot on:",gmax,"with absolute:",(k+gmax-1)))
          }
          P <- diag(ncol(A))
          Ptemp <-  P[k,]; P[k,] <- P[(k+dmax-1),]; P[(k+dmax-1),] = Ptemp
          A <- P%*%A%*%P
          L <- P%*%L%*%P
          Pprod <- P%*%Pprod
        }
        normj <- sum(abs(A[(k+1):n,k]))
        delta <- max(0,deltaprev,-A[k,k]+max(normj,tau*gamm))
        if (delta > 0)  {
          A[k,k] <- A[k,k] + delta
          deltaprev <- delta
        }
      }
      
      L[k,k] <- A[k,k] <- sqrt(A[k,k])
      for (i in (k+1):n)  {
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
  if (delta > 0)  {
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
  
  r = t(Pprod)%*%t(L)%*%t(Pprod)
  attr(r,"delta")=delta
  return(r)
}

getInitialVauleOfRandomEstimate<- function(dynrModel){
  mdcov <- prep.noise(
    values.latent=dynrModel@noise@values.latent,
    params.latent=dynrModel@noise@params.latent,
	values.observed=dynrModel@noise@values.observed,
	params.observed=dynrModel@noise@params.observed
  )

  meas <- prep.measurement(
    values.load = dynrModel@measurement@values.load,
    params.load = dynrModel@measurement@params.load,
    obs.names = dynrModel@measurement@obs.names,
    state.names = dynrModel@measurement@state.names
  )

  initial <- prep.initial(
    values.inistate=dynrModel@initial@values.inistate,
    params.inistate=dynrModel@initial@params.inistate,
    values.inicov=dynrModel@initial@values.inicov, 
    params.inicov=dynrModel@initial@params.inicov
  )

  formula <- unlist(dynrModel@dynamics@formula2)[1:length(dynrModel@measurement@state.names)]
  print(formula)
  #print(prep.thetaFormula(formula[1], dynrModel@measurement@params.int, dynrModel@dynamics@random.names))
  #print(prep.thetaFormula(formula[2], dynrModel@measurement@params.int, dynrModel@dynamics@random.names))
  
  for(i in 1:length(formula))
    formula[i] <- prep.thetaFormula(formula[i], dynrModel@measurement@params.int, dynrModel@dynamics@random.names)
  #formula2 <- lapply(formula, function(x){prep.thetaFormula(x, dynrModel@measurement@params.int, dynrModel@dynamics@random.names)})
  print(formula)
  print("...")
  dynm<-prep.formulaDynamics(
    formula=formula,
	#formula = formula,
    startval=dynrModel@dynamics@startval,
    isContinuousTime=dynrModel@dynamics@isContinuousTime,
    beta.names=names(dynrModel@dynamics@startval)
  )

  model <- dynr.model(dynamics=dynm, measurement=meas,
                    noise=mdcov, initial=initial, data=dynrModel@data, saem=FALSE,
                    outfile="VanDerPol_.c")

  fitted_model <- dynr.cook(model, optimization_flag = TRUE, hessian_flag = FALSE, saem=FALSE)

}
