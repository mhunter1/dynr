# Author: Michael D. Hunter
# Date: 2020-06-20 10:26:21
# Filename: dynrPredict.R
# Purpose: Write a predict method for dynr
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# TODO Use of filtered versus predicted?
# TODO Use of time in row 1 or row 2.  Off by 1 issue?
# TODO Wrap everything in a function.  Function takes dynrModel at dynrCook point estimates?
# TODO Generalize to multiple IDs.  Via 'newdata' argument.

##' \code{predict} method for \code{dynrModel} objects
##' 
##' @param object a dynrModel object from which predictions are desired
##' @param newdata an optional \code{data.frame} or \code{ts} object. See details.
##' @param interval character indicating what kind of intervals are desired.  'none' gives no intervals, 'confidence', gives confidence intervals, 'prediction' gives prediction intervals.
##' @param method character the method used to create the forecasts.  See details.
##' @param level the confidence or predictions level, ignored if not using intervals
##' @param type character the type of thing you want predicted: latent variables or manifest variables.
##' @param ... further named arguments, e.g., \code{size} for the ensemble size when using the ensemble prediction
##' 
##' @return A list of the prediction estimates, intervals, and ensemble members.
##' 
##' @details
##' The \code{newdata} argument is either a \code{data.frame} or \code{ts} object.  It passed as the \code{dataframe} argument of \code{dynr.data} and must accept the same further arguments as the data in the model passed in the \code{object} argument (e.g., same \code{id}, \code{time}, \code{observed}, and \code{covariates} arguments).
##' 
##' The available methods for prediction are 'kalman' and 'ensemble'.  The 'kalman' method uses the Kalman filter to create predictions.  The 'ensemble' method simulates a set of initial conditions and lets those run forward in time.  The distribution of this ensemble provides the predictions.  The mean is the value predicted.  The quantiles of the distribution provide the intervals.
predict.dynrModel <- function(object,
                              newdata=NULL,
                              interval=c('none', 'confidence', 'prediction'),
                              method=c('kalman', 'ensemble'),
                              level=.95,
                              type=c('latent', 'observed'),
                              ...
                             ){
	method <- match.arg(arg=tolower(method), choices=c('kalman', 'ensemble'))
	interval <- match.arg(arg=tolower(interval), choices=c('none', 'confidence', 'prediction'))
	type <- match.arg(arg=tolower(type), choices=c('latent', 'observed'))
	if(!is.null(newdata)){
		ddat <- dynr.data(newdata, id=object$data$idVar, time=object$data$timeVar, observed=object$data$observed.names)
	} else {
		ddat <- object$data
	}
	# Process ... arg
	dots <- list(...)
	size <- dots[['size']]
	# Create model
	model0 <- dynr.model(dynamics=object$dynamics, measurement=object$measurement,
		noise=object$noise, initial=object$initial, data=ddat, outfile='forecast.c')
	# Run model to get filtered estimates
	cook0 <- dynr.cook(model0, debug_flag=TRUE, verbose=FALSE,
		optimization_flag=FALSE, hessian_flag=FALSE)
	# Set initial conditions to final filtered estimates ???
	# But only if we're forecasting ahead?
	# TODO determine details of the needed structure of newdata
	# TODO check/generalize to multiple latent variables and multiple people
	# Need to have a reasonable return structure
	if(method == 'kalman'){
		predKal <- cook0$eta_predicted
		kalSE <- apply(cook0$error_cov_predicted, 3, function(x){qnorm(1-(1-level)/2)*sqrt(diag(x))})
		kalCI <- matrix(c(predKal - kalSE, predKal + kalSE), nrow=2, byrow=TRUE)
		return(list(estimate=predKal, CI=kalCI))
	} else if(method == 'ensemble'){
		numEns <- size
		predEnsK <- array(NA, dim=c(nrow(ddat$observed), model0$dim_latent_var, numEns))
		for(n in 1:numEns){
			if(n == 1){model0@compileLib <- TRUE} else{model0@compileLib=FALSE}
			cook1 <- dynr.cook(model0, debug_flag=TRUE,
				verbose=FALSE, optimization_flag=FALSE, hessian_flag=FALSE,
				perturb_flag=TRUE)
			predEnsK[,,n] <- cook1$eta_predicted
		}
		# Make quantile-based 95% CIs for ensemble forecast
		ensCI <- apply(predEnsK, c(1, 2), quantile, probs=c((1-level)/2, 1-(1-level)/2))
		# Return ensemble mean, quantile CI, and ensemble members instead
		return(list(estimate=apply(predEnsK, c(1, 2), mean), CI=ensCI, members=predEnsK))
	} else {
		stop(paste('Unknown method', method, 'when trying to predict'))
	}
}



# Proposed syntax example
#myNewData <- data.frame(time=seq(1+dt, 1+dt+T, length.out=numMissing), ID=1)
#predict(resOU, newdata=myNewData, method='ensemble', intervals='confidence', level=.80)
#predict(resOU, newdata=myNewData) #method=kalman, no intervals

#d <- modelOU$data
#bob <- predict(modelOU2)
#kal <- predict(modelOU2, method='ensemble', size=100)
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



#------------------------------------------------------------------------------



