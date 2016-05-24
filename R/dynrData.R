
#setClass(Class = "dynrData",
#	representation = representation(
#	data
#	)
#)

#setMethod("initialize", "dynrData",
#	function(.Object){
#		return(.Object)
#	}
#)

#setMethod("summary", "dynrData",
#	function(.Object){
#		summary(.Object@data)
#	}
#)

#setMethod("print", "dynrData", function(x,...) { 
#	displayDynrData(x) 
#})

#setMethod("show", "dynrData", function(x,...) { 
#	displayDynrData(x) 
#})

##' Create a list of data for parameter estimation (cooking dynr) using \code{\link{dynr.cook}}
##'
##' @param dataframe the data in the format of a data frame in R that contain a column of subject ID numbers 
##' (i.e., an ID variable), a column indicating subject-specific measurement occasions
##' (i.e., a TIME variable), at least one column of observed values, and any number of covariates. 
##' If the data are fit to a discrete-time model, the TIME variable should contain subject-specific sequences 
##' of consecutively equally spaced numbers (e.g, 1, 2, 3, ...).
##' If the data are fit to a continuous-time model, the TIME varibles can contain subject-specific increasing sequences 
##' of irregularly spaced real numbers.  
##' Missing values in the observed variables shoud be indicated by NA. Missing values in covariates are not allowed. 
##' That is, missing values in the covariates, if there are any, should be imputed first. 
##' @param id a character string of the name of the ID variable in the data.
##' @param time a character string of the name of the TIME variable in the data.
##' @param observed a vector of character strings of the names of the observed variables in the data.
##' @param covariates a vector of character strings of the names of the covariates in the data, which can be missing.
##' 
##' @examples
##' data(EMG)
##' ds <- EMG
##' ds$ID <- rep(1, nrow(EMG))
##' ds$t <- 1:nrow(EMG)
##' dd <- dynr.data(ds, id='ID', time='t', observed='EMG', covariates='self')
dynr.data <- function(dataframe, id, time, observed, covariates){
	ids <- unique(dataframe[ , id])
	tstart <- c(
		sapply(1:length(ids),
			function(i){
				min(which(dataframe[,id] %in% ids[i]))
			}
		)-1,
		dim(dataframe)[1])
	data.object <- list(
		id=dataframe[,id],
		tstart=as.integer(tstart),
		time=as.double(dataframe[,time]),
		observed=data.frame(apply(dataframe[ , observed, drop=FALSE], 2, as.double)),
		observed.names=observed
		)
	if(!missing(covariates)){
	  data.object$covariate.names <- covariates
		data.object$covariates <- data.frame(apply(dataframe[,covariates, drop=FALSE], 2, as.double))
		names(data.object$covariates) <- paste0("covar", 1:length(covariates))
	}
	names(data.object$observed) <- paste0("obs", 1:length(observed))
	return(data.object)
}



