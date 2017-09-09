#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2017-09-09
# Filename: dynrTogo.R
# Purpose: Non-matrix helper functions for all the prep.*() recipe functions.
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

.covTypes <- c('identity', 'symmetric', 'diagonal', 'scalar', 'zero', 'ar1')

#------------------------------------------------------------------------------

##' Order some noise to-go
##' 
##' @param latent.type character. type of the latent covariance matrices (dynamic noise)
##' @param latent.dim numeric. dimension of the latent process
##' @param latent.group grouping structure of free parameters for latent covariance across regimes
##' @param observed.type character. type of the observed covariance matrices (measurement noise)
##' @param observed.dim numeric. dimension of the observed process
##' @param observed.group grouping structure of free parameters for observed covariance across regimes
##' 
##' @details
##' Just like food that you order to-go, this function quickly and easily gets you something tasty, but you have fewer options and less control than when preparing recipes from scratch.
##'  
##' Currently supported matrix types are 'identity', 'symmetric', 'diagonal', 'scalar', and 'zero'.
##' 
##' @seealso
##' \code{\link{dynr.noise}}
##' 
##' @examples
##' 
##' library(dynr)
##' 
##' togo.noise('identity', 1, '', 'symmetric', 4, '')
togo.noise <- function(latent.type, latent.dim, latent.group,
	observed.type, observed.dim, observed.group){
	# TODO match partial strings for 'type' and be case insensitive
	#  probably via grep or match.arg
	mltype <- match.arg(tolower(latent.type), .covTypes, several.ok=TRUE)
	# Throw error if length of mltype and latent.type differ.  This means some were not matched.
	motype <- match.arg(tolower(observed.type), .covTypes, several.ok=TRUE)
	
	# Check that lengths of type arguments match or are 1
	# Check that lengths of dim arguments match or are 1
	# Check that lengths of group arguments match or are 1
	# Check that lengths of type, dim, group arguments match or are 1
	# Write function that returns TRUE when arguments have same length or one of them is one.
	# Then check all 6*5/2 = 15 combos for this
	
	
	# Fill in default groupings
	# every regimes has completely separate free params by default
	# But check for any regimes that share params must be of the same type
	
	# Build matrices
	vlat <- buildValuesMatrix(latent.type, latent.dim, 'dnoise')
	plat <- buildParamsMatrix(latent.type, latent.dim, latent.group, 'dnoise')
	vobs <- buildValuesMatrix(observed.type, observed.dim, 'mnoise')
	pobs <- buildParamsMatrix(observed.type, observed.dim, observed.group, 'mnoise')
	
	obj <- prep.noise(values.latent=vlat, params.latent=plat, values.observed=vobs, params.observed=pobs)
	return(obj)
}

#------------------------------------------------------------------------------

buildValuesMatrix <- function(type, dim, recipe){
	if(length(dim) == 1) dim <- c(dim, dim)
	if(type=='diagonal' || type == 'scalar'){
		if(recipe=='dnoise'){
			value <- 1.1
		} else if(recipe=='mnoise'){
			value <- 0.2
		}
		mat <- diag(value, dim[1], dim[2])
	} else if(type=='symmetric'){
		if(recipe=='dnoise'){
			mat <- matrix(0.4, dim[1], dim[2])
			diag(mat) <- 1.1
		} else if(recipe=='mnoise'){
			mat <- matrix(0.2, dim[1], dim[2])
			diag(mat) <- 0.5
		}
	} else if(type=='identity' || type=='zero'){
		value <- ifelse(type == 'identity', 1, 0)
		mat <- diag(value, dim[1], dim[2])
	} else if(type=='ar1'){
		stop("Sorry, but Mike is a bonehead and hasn't yet implemented this.")
	}
	return(mat)
}

#buildValuesMatrix('diagonal', 4, 'dnoise')
#buildValuesMatrix('diagonal', 4, 'mnoise')
#buildValuesMatrix('scalar', 4, 'dnoise')
#buildValuesMatrix('scalar', 4, 'mnoise')
#buildValuesMatrix('identity', 4, 'dnoise')
#buildValuesMatrix('identity', 4, 'mnoise')
#buildValuesMatrix('symmetric', 4, 'dnoise')
#buildValuesMatrix('symmetric', 4, 'mnoise')
#buildValuesMatrix('zero', 4, 'dnoise')
#buildValuesMatrix('zero', 4, 'mnoise')
#buildValuesMatrix('ar1', 4, 'dnoise')
#buildValuesMatrix('ar1', 4, 'mnoise')


#------------------------------------------------------------------------------

buildParamsMatrix <- function(type, dim, group, recipe){
	if(length(dim) == 1) dim <- c(dim, dim)
	if(type=='diagonal'){
		if(recipe=='dnoise'){
			value <- paste0('zeta_', 1:dim[1], group)
		} else if(recipe=='mnoise'){
			value <- paste0('epsilon_', 1:dim[1], group)
		}
		mat <- diag(value, dim[1], dim[2])
	} else if(type == 'scalar'){
		if(recipe=='dnoise'){
			value <- paste0('zeta', ifelse(group=='', '', '_'), group)
		} else if(recipe=='mnoise'){
			value <- paste0('epsilon', ifelse(group=='', '', '_'), group)
		}
		mat <- diag(value, dim[1], dim[2])
	} else if(type=='symmetric'){
		if(recipe=='dnoise'){
			mat <- outer(paste0('zeta_', 1:dim[1]), paste0(1:dim[2], group), paste0)
		} else if(recipe=='mnoise'){
			mat <- outer(paste0('epsilon_', 1:dim[1]), paste0(1:dim[2], group), paste0)
		}
		mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
	} else if(type=='identity' || type=='zero'){
		value <- 'fixed'
		mat <- diag(value, dim[1], dim[2])
	} else if(type=='ar1'){
		stop("Sorry, but Mike is a bonehead and hasn't yet implemented this.")
	}
	return(mat)
}



#buildParamsMatrix('diagonal', 4, 'dnoise', group='')
#buildParamsMatrix('diagonal', 4, 'mnoise', group='')
#buildParamsMatrix('scalar', 4, 'dnoise', group='')
#buildParamsMatrix('scalar', 4, 'mnoise', group='')
#buildParamsMatrix('identity', 4, 'dnoise', group='')
#buildParamsMatrix('identity', 4, 'mnoise', group='')
#buildParamsMatrix('symmetric', 4, 'dnoise', group='')
#buildParamsMatrix('symmetric', 4, 'mnoise', group='')
#buildParamsMatrix('zero', 4, 'dnoise', group='')
#buildParamsMatrix('zero', 4, 'mnoise', group='')
#buildParamsMatrix('ar1', 4, 'dnoise', group='')
#buildParamsMatrix('ar1', 4, 'mnoise', group='')


#------------------------------------------------------------------------------



#------------------------------------------------------------------------------




