
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
##' @param dataframe either a ``ts'' class object of time series data for a single subject or 
##' a data frame object of data for potentially multiple subjects that contain a column of subject ID numbers 
##' (i.e., an ID variable), a column indicating subject-specific measurement occasions
##' (i.e., a TIME variable), at least one column of observed values, and any number of covariates. 
##' If the data are fit to a discrete-time model, the TIME variable should contain subject-specific sequences 
##' of (subsets of) consecutively equally spaced numbers (e.g, 1, 2, 3, ...). That is, the
##' program assumes that the input data.frame is equally spaced with potential missingness. If
##' the measurement occasions for a subject are a subset of an arithmetic sequence but are not
##' consecutive, NAs will be inserted automatically to create an equally spaced data set before
##' estimation.
##' If the data are fit to a continuous-time model, the TIME varibles can contain subject-specific increasing sequences 
##' of irregularly spaced real numbers.  
##' Missing values in the observed variables shoud be indicated by NA. Missing values in covariates are not allowed. 
##' That is, missing values in the covariates, if there are any, should be imputed first. 
##' @param id a character string of the name of the ID variable in the data. Optional for a ``ts'' class object. 
##' @param time a character string of the name of the TIME variable in the data. Optional for a ``ts'' class object.
##' @param observed a vector of character strings of the names of the observed variables in the data. 
##' Optional for a ``ts'' class object.
##' @param covariates (optional) a vector of character strings of the names of the covariates in the data,
##'  which can be missing.
##' 
##' @examples
##' data(EMGsim)
##' dd <- dynr.data(EMGsim, id = 'id', time = 'time', observed = 'EMG', covariates = 'self')
##' 
##' z <- ts(matrix(rnorm(300), 100, 3), start = c(1961, 1), frequency = 12)
##' dz <- dynr.data(z)
dynr.data <- function(dataframe, id = 'id', time = 'time', observed, covariates){
	if (is.ts(dataframe)){
		# ts class 
		# single or multivariate time series
		# one subject
		tsp=attributes(dataframe)$tsp
		dataframe = as.data.frame(dataframe)
		if (missing(observed)){
			observed = colnames(dataframe)
		}
		dataframe[ , id] <- 1
		dataframe[ , time] <- seq(from = tsp[1], to = tsp[2], length.out = nrow(dataframe))
	}
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
		observed.names=observed,
		original.data=dataframe,
		idVar=id,
		timeVar=time
		)
	if(!missing(covariates)){
	  data.object$covariate.names <- covariates
		data.object$covariates <- data.frame(apply(dataframe[,covariates, drop=FALSE], 2, as.double))
		names(data.object$covariates) <- paste0("covar", 1:length(covariates))
	}
	names(data.object$observed) <- paste0("obs", 1:length(observed))
	fid <- as.factor(data.object$id)
	time.split <- split(data.object$time, fid)
	diff.npos <- sapply(time.split, function(x){d <- diff(x); any(d <= 0)})
	if(any(diff.npos)){
		msg <- paste0("Time steps are not all increasing.\nFound zero or negative time differences for IDs: ", paste(levels(fid)[diff.npos], sep="", collapse=", "))
		stop(msg)
	}
	return(data.object)
}



