# TODO Add elements to the class representation for other things passed
#  back from backend.
# Currently everything not named here is stored in 'misc'

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

#[1] "exitflag"                         "neg.log.likelihood"              
#[3] "fitted.parameters"                "hessian.matrix"                  
#[5] "inverse.hessian.matrix"           "eta_regime_t"                    
#[7] "error_cov_regime_t"               "eta_regime_regime_t_pred"        
#[9] "error_cov_regime_regime_t_pred"   "eta_regime_regime_t_plus_1"      
#[11] "error_cov_regime_regime_t_plus_1" "innov_vec"                       
#[13] "inverse_residual_cov"             "pr_t_given_t"                    
#[15] "pr_t_given_t_less_1"              "pr_t_given_T"                    
#[17] "transprob_given_T"                "eta_regime_smooth"               
#[19] "error_cov_regime_smooth"          "eta_smooth_final"                
#[21] "error_cov_smooth_final" 

#tval = d$transformed.parameters/d$standard.errors to be added
#cbind(d$transformed.parameters,d$standard.errors,) to be added


# TODO Add population of class slots for other things added to representation
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
		.Object@error_cov_regime_regime_t_plus_1  <- x$error_cov_regime_regime_t_plus_1
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
		#x[c('fitted.parameters', 'transformed.parameters', 'standard.errors', 'hessian.matrix', 'transformed.hessian', 'CI')] <- NULL
		#.Object@misc <- x
		return(.Object)
	}
)

# Set the summary method of an object of class dynrRun
#  All this amounts to is writing a function that takes a
#  dynrRun object (and possibly other arguments)
#  and returns/does whatever we want.
setMethod("summary", "dynrRun",
	function(object){
		ret <- data.frame(transformed.parameters=object@transformed.parameters, standard.errors=object@standard.errors, CI = object@conf.intervals)
		return(ret)
	}
)

displayDynrRun <- function(x){
	str(x)
}

setMethod("print", "dynrRun", function(x, ...) { 
	displayDynrRun(object) 
})

setMethod("show", "dynrRun", function(object) { 
	displayDynrRun(object) 
})

findR = function(y){
  Rindex2 = 1:model$num_regime
  Rtoret = Rindex2[y==max(y)]
  return(Rtoret)
}

Mode <- function(y) {
  uy <- unique(y)
  uy[which.max(tabulate(match(y, uy)))]
}

setMethod("plot", "dynrRun",
	function(x, y=NULL, data, graphingPar=par(no.readonly = TRUE), ylab = "Smoothed state values", xlab = "Time", numSubjDemo=2, legend.cex=1.2){
		opar = par(no.readonly = TRUE)
		par(graphingPar)
		thesmooth = data.frame(t(x@eta_smooth_final))  
		colnames(thesmooth) = paste0("state",1:model$dim_latent_var)
		ID = data[["id"]]
		rowIndex = 1:length(ID)
		thes = sort(sample(1:model$"num_sbj",numSubjDemo))
		par(mfrow=c(ifelse(numSubjDemo%%2 > 0,
			numSubjDemo,numSubjDemo/2),
			ifelse(numSubjDemo%%2 > 0,
			1,2)))
		if (model$num_regime > 1){
			thePr = t(x@pr_t_given_T)
			mostLikelyRegime = apply(thePr,1,findR)  
			theR = Mode(mostLikelyRegime)
		}
		for (s in 1:numSubjDemo){
			T = length(ID[ID==thes[s]])
			therow = rowIndex[ID==thes[s]] 
	 		plot(c(1,T),c(min(thesmooth)-1,max(thesmooth)+quantile(unlist(thesmooth),.1)),
				ylab=ifelse(exists("ylab"),ylab,"State values"), 
				xlab=ifelse(exists("xlab"),xlab,"Time"),
				main=ifelse(exists("main"),main,""),
				type='n')
			if (model$num_regime > 1){
				times=1:T
				thepri = mostLikelyRegime[therow]
				rect(times[thepri==theR]-.5,quantile(unlist(thesmooth),.01),times[thepri==theR]+.5,quantile(unlist(thesmooth),.99),col="yellow",density=30)
				#legend('topleft',paste0("Regime",theR), 	
				#bty="n",cex=1.4,col=c(NA),fill=c("yellow"),
				#density=c(100))
			}
			time2 = if(T>500){
				seq(1,T,9)
			}else {1:T}
			
			for (j in 1:model$dim_latent_var){
				lines(1:T,thesmooth[therow,j],lty=1,lwd=1,col=j)
				points(time2,thesmooth[therow,j][time2],pch=as.character(j),col=j)
			}
			stateNames = NULL
			for(j in 1:model$dim_latent_var){
				stateNames = c(stateNames, paste0('State ',j))
			}
			legend('topright',
			paste0("State",as.character(1:model$dim_latent_var)),
				lty=1:model$dim_latent_var,
				lwd=2,col=1:model$dim_latent_var,bty="n",
				pch=as.character(1:model$dim_latent_var),cex=legend.cex)
		}
		on.exit(par(opar))
	}
)


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


