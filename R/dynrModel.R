# TODO
# create dynrModel objects
# They should output a list like this
#model <- list(num_sbj=20,
#              dim_latent_var=2,
#              dim_obs_var=2,
#              dim_co_variate=2, 
#              num_regime=2,
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

default.model.options <- list(xtol_rel=1e-7, stopval=-9999, ftol_rel=-1, ftol_abs=-1, maxeval=as.integer(-1), maxtime=-1)

dynr.model <- function(num_regime=1, dim_latent_var, xstart, ub, lb, options=default.model.options){
	if(!is.list(options)){
		stop("'options' argument to dynr.model function must be a list.")
	}
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
	return(list(num_regime=num_regime, dim_latent_var=dim_latent_var, xstart=xstart, ub=ub, lb=lb, num_func_param=as.double(length(xstart)), options=options))
}

