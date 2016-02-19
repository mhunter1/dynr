# A Recipe is a helper function that takes user input
#  and produces a C function definition that can be
#  compiled by dynrFuncaddress.

# TODO add check that there is no cross loading on one of the ID variables

# TODO add assignment of first loading to 1 in C
# TODO rename function to more like fastLoadings
# TODO create new mid-level function (for which this is a wrapper) that takes three
#  matrices (parameter staring values, free-fixed, and parameter numbers)
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
	
	
	k <- 1
	valuesMat <- matrix(0, nx, ne)
	paramsMat <- matrix(0, nx, ne)
	ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
	for(j in 1:ne){
		for(i in 1:nx){
			if(allVars[i] %in% map[[j]] & !(allVars[i] %in% idvar)){
				ret <- paste(ret,
					'    gsl_matrix_set(Ht, ', i-1, ', ', j-1,
					', param[', params[k]-1, ']);\n', sep='')
				paramsMat[i, j] <- params[k]
				valuesMat[i, j] <- .8
				k <- k+1
			} else if(allVars[i] %in% map[[j]] & allVars[i] %in% idvar){
				valuesMat[i, j] <- 1
			}
		}
	}
	ret <- paste(ret, "\n}\n\n")
	cat(ret)
	return(list(values=valuesMat, params=paramsMat, C=ret))
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


# values, free, and params are all MxN matrices
dynr.matrixLoadings <- function(values, params){
	ne <- ncol(values)
	nx <- nrow(values)

	ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
	for(j in 1:ne){
		for(i in 1:nx){
			if(params[i, j] > 0){
				ret <- paste(ret,
					'    gsl_matrix_set(Ht, ', i-1, ', ', j-1,
					', param[', params[i, j] - 1, ']);\n', sep='')
			} else if(values[i, j] != 0){
				ret <- paste(ret,
					'    gsl_matrix_set(Ht, ', i-1, ', ', j-1,
					', ', values[i, j], ');\n', sep='')
			}
		}
	}
	ret <- paste(ret, "\n}\n\n")
	cat(ret)

}


# Examples
# dynr.matrixLoadings(diag(1, 5), diag(1:5))
# dynr.matrixLoadings(matrix(1, 5, 5), diag(1:5))


# Error covariance matrix---
dynr.error_cov <- function(map, params, idvar)
{
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
  
  nx <- length(allVars)
  ne <- length(map)
  
  k <- 1
  ret <- "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){\n\n"
  for(j in 1:ne){
    for(i in 1:nx){
      if(allVars[i] %in% map[[j]] & !(allVars[i] %in% idvar)){
        ret <- paste(ret,
                     '    gsl_matrix_set(y_noise_cov, ', i-1, ', ', j-1,
                     ', param[', (params[k]-1), ']);\n', sep='')
        k <- k+1
      }
    }
  }
  ret <- paste(ret, "\n}\n\n")
  cat(ret)
}


