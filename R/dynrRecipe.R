# A Recipe is a helper function that takes user input
#  and produces a C function definition that can be
#  compiled by dynrFuncaddress.

# TODO add check that there is no cross loading on one of the ID variables

dynr.loadings <- function(map, params, idvar){
	if(missing(idvar)){
		idvar <- sapply(map, '[', 1)
	}
	
	allVars <- unique(unlist(map))
	
	nx <- length(allVars)
	ne <- length(map)
	
	if(!all(idvar %in% c(names(map), unlist(map)))){
		stop("The 'idvar' must all be either in the names of the 'map' argument or parts of the part 'map' argument.")
	}
	
	paramsNeeded <- length(unlist(map)) - sum(idvar %in% allVars)
	if(length(params) != paramsNeeded){
		stop(paste0("Number of free parameters provided (", length(params), ") does not match number needed (", paramsNeeded, ")."))
	}
	
	SOME_CONDITION_MET <- TRUE
	SOMETHING_WITH_PARAMNUMS <- 1
	
	k <- 1
	ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
	for(j in 1:ne){
		for(i in 1:nx){
			if(allVars[i] %in% map[[j]] & !(allVars[i] %in% idvar)){
				ret <- paste(ret,
					'    gsl_matrix_set(Ht, ', i-1, ', ', j-1,
					', param[', params[k]-1, ']);\n', sep='')
				k <- k+1
			}
		}
	}
	ret <- paste(ret, "\n}\n\n")
	cat(ret)
}

# Examples
# Single factor model with one latent variable
#dynr.loadings( list(eta1=paste0('y', 1:4)), 4:6)

# Two factor model with simple structure
#dynr.loadings( list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:7)), c(4:6, 1:2))

# Two factor model with repeated use of a free parameter
#dynr.loadings( list(eta1=paste0('y', 1:4), eta2=paste0('y', 5:8)), c(4:6, 1:2, 4))

# Two factor model with a cross loading
#dynr.loadings( list(eta1=paste0('y', 1:4), eta2=c('y5', 'y2', 'y6')), c(4:6, 1:2))

