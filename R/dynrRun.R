

dynr.run <- function(model, data, transformation, conf.level=.95) {
	if(missing(transformation)){
		transformation <- function(x){x}
	}
	tmp <- .Call("main_R", model, data, PACKAGE = "dynr")
	tmp <- endProcessing(tmp, transformation, conf.level)
	return(tmp)
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
	x$CI <- matrix(c(tParam - tSE*confx, tParam + tSE*confx), ncol=2, dimnames=list(NULL, c('lower', 'upper')))
	return(x)
}


