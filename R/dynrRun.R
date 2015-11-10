
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
	#output2 <- endProcessing(output, transformation, conf.level)
	#obj <- new("dynrRun", output2)
	frontendStop <- Sys.time()
	totalTime <- frontendStop-frontendStart
	backendTime <- backendStop-backendStart
	fronendTime <- totalTime-backendTime
	cat('Total Time:', totalTime, '\n')
	cat('Backend Time:', backendTime, '\n')
	# TODO add timing information to dynrRun object
	return(output)
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
	CI <- cbind(tParam - tSE*confx, tParam + tSE*confx)
}
	x$transformed.parameters <- tParam
	x$standard.errors <- tSE
	x$transformed.inv.hessian <- tHess
	x$conf.intervals <- matrix(CI, ncol=2, byrow=T, dimnames=list(NULL, c('ci.lower', 'ci.upper')))
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



