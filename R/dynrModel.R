# TODO
# create dynrModel objects
# They should output a list like this
#model <- list(num_sbj=20,
#              dim_latent_var=2,
#              dim_obs_var=2,
#              dim_co_variate=2, 
#              num_regime=2,
#              isDiscretTime=0,
#              xstart=c(rep(log(.1), 4), log(10.0), log(10.0), -3.0, 9.0, -1.5, -0.5, 95.0,-.3,-.3),
#              num_func_param=13,
#              ub=c(rep(10, 6), rep(20, 4), 1000, 20, 20),
#              lb=c(rep(-10, 6), rep(-20, 4), 0, -20, -20)
#)

# Really all the model needs right now is
# num_regime, xstart, ub, lb
# Everything else can be gathered from the data.
# Importantly, the thing handed to the backend must
# remain a list exactly like the above.


default.model.options <- list(xtol_rel=1e-7, stopval=-9999, ftol_rel=1e-10, 
                              ftol_abs=-1, maxeval=500, maxtime=-1)
#' Create the model for dynr
#' 
#' @param num_regime An integer number of the regimes.
#' @param dim_latent_var An integer number of the latent variables.
#' @param xstart The starting values for parameter estimation.
#' @param ub The upper bounds of the estimated parameters.
#' @param lb The lower bounds of the estimated parameters.
#' @param isDiscreteTime A binary flag indicating whether the model is a discrete-time model (0 = no; 1 = yes)
#' @param options A list of NLopt estimation options. By default, xtol_rel=1e-7, stopval=-9999, ftol_rel=-1, ftol_abs=-1, maxeval=as.integer(-1), and maxtime=-1.
#' @return A list of model statements to be passed to dynr.run().
dynr.model <- function(num_regime=1, dim_latent_var, xstart, ub, lb, isDiscreteTime=0, options=default.model.options){
	if(!is.list(options)){
		stop("'options' argument to dynr.model function must be a list.")
	}
	options <- processModelOptionsArgument(options)
	xlen <- length(xstart)
	ulen <- length(ub)
	llen <- length(lb)
	if( (xlen != ulen) || (xlen != llen) || (ulen != llen)){
		stop("Length of 'xstart', 'ub', and 'lb' must match.")
	}
	if( (length(num_regime) != 1) || (round(num_regime) != num_regime) ){
		stop("Number of regimes (num_regime) must be a single integer.")
	}
	if( (length(dim_latent_var) != 1) || (round(dim_latent_var) != dim_latent_var) ){
		stop("Number of latent variables (dim_latent_var) must be a single integer.")
	}
	# TODO remove these as.double calls
	# by changing how the arguments are processed in the backend.
	return(list(num_regime=num_regime, dim_latent_var=dim_latent_var, xstart=xstart, ub=ub, lb=lb, isDiscreteTime=isDiscreteTime, num_func_param=as.double(length(xstart)), options=options))
}

processModelOptionsArgument <- function(opt){
	if(!identical(opt, default.model.options)){
		nameMatch <- names(opt) %in% names(default.model.options)
		if( any(!nameMatch) ){
			msg <- paste("Tried to set invalid option(s): ", paste(names(opt)[!nameMatch], collapse=', '))
			stop(msg)
		}
		newopt <- default.model.options
		for(i in 1:length(opt)){
			newopt[[names(opt)[i]]] <- opt[[i]]
		}
		newopt$maxeval <- as.integer(newopt$maxeval)
		return(newopt)
	}else{
		return(opt)
	}
}

