
# TODO Add elements to the class representation for other things passed
#  back from backend.
# Currently everything not named here is stored in 'misc'
setClass(Class = "dynrRun",
	representation = representation(
	fitted.parameters = "numeric",
	transformed.parameters = "numeric",
	standard.errors = "numeric",
	hessian = "matrix",
	transformed.hessian = "matrix",
	conf.intervals = "matrix",
	misc = "list"
	)
)

# TODO Add population of class slots for other things added to representation
setMethod("initialize", "dynrRun",
	function(.Object, x){
		.Object@fitted.parameters <- x$fitted.parameters
		.Object@transformed.parameters <- x$transformed.parameters
		.Object@standard.errors <- x$standard.errors
		.Object@hessian <- x$hessian.matrix
		.Object@transformed.hessian <- x$transformed.hessian
		.Object@conf.intervals <- x$CI
		x[c('fitted.parameters', 'transformed.parameters', 'standard.errors', 'hessian.matrix', 'transformed.hessian', 'CI')] <- NULL
		.Object@misc <- x
		return(.Object)
	}
)

#setMethod("summary", "dynrRun",
#	function(.Object){
#		#something nice
#	}
#)

#setMethod("print", "dynrRun", function(x,...) { 
#	displayDynrRun(x) 
#})

#setMethod("show", "dynrRun", function(x,...) { 
#	displayDynrRun(x) 
#})

#setMethod("plot", "dynrRun",
#	function(.Object){
#		#something nice
#	}
#)

dynr.run <- function(model, data, transformation, conf.level=.95) {
	if(missing(transformation)){
		transformation <- function(x){x}
	}
	frontendStart <- Sys.time()
	model <- preProcessModel(model)
	backendStart <- Sys.time()
	output <- .Call("main_R", model, data, PACKAGE = "dynr")
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
	tHess <- J %*% solve(x$hessian.matrix) %*% t(J)
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


