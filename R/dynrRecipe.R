# A Recipe is a helper function that takes user input
#  and produces a C function definition that can be
#  compiled by dynrFuncaddress.

# TODO add check that there is no cross loading on one of the ID variables

# TODO add assignment of first loading to 1 in C
# TODO rename function to more like fastLoadings
# TODO create new mid-level function (for which this is a wrapper) that takes three
#  matrices (parameter staring values, free-fixed, and parameter numbers)


#------------------------------------------------------------------------------
# Create recipe for function measurement

# TODO add ability to use covariates in these functions

#--------------------------------------
# brief input version

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
	for(j in 1:ne){
		for(i in 1:nx){
			if(allVars[i] %in% map[[j]] & !(allVars[i] %in% idvar)){
				paramsMat[i, j] <- params[k]
				valuesMat[i, j] <- .8
				k <- k+1
			} else if(allVars[i] %in% map[[j]] & allVars[i] %in% idvar){
				valuesMat[i, j] <- 1
			}
		}
	}
	return(list(values=valuesMat, params=paramsMat))
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


#--------------------------------------
# matrix input version

# values, and params are all MxN matrices
# a zero param is taken to be fixed.

dynr.matrixLoadings <- function(values, params){
	ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
	ret <- paste(ret, setGslMatrixElements(values, params, "Ht"), sep="\n")
	ret <- paste(ret, "\n    gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);\n")
	ret <- paste(ret, "\n}\n\n")

	return(ret)
}


# Examples
# a <- dynr.loadings( list(eta1=paste0('y', 1:4), eta2=c('y5', 'y2', 'y6')), c(4:6, 1:2))
# dynr.matrixLoadings(a$values, a$params)
#
# dynr.matrixLoadings(diag(1, 5), diag(1:5))
# dynr.matrixLoadings(matrix(1, 5, 5), diag(1:5))
# dynr.matrixLoadings(diag(1, 5), diag(0, 5)) #identity measurement model


#------------------------------------------------------------------------------
# Error covariance matrix
# N.B. This function produces BOTH the latent and observed error covariance matrices.

# TODO break below into two functions that populate into the same C function
#  at cook time.

#--------------------------------------
# matrix input version
dynr.matrixErrorCov <- function(values.latent, params.latent, values.observed, params.observed){
	ret <- "void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){\n\n"
	ret <- paste(ret, setGslMatrixElements(values.latent, params.latent, "eta_noise_cov"), sep="\n")
	ret <- paste(ret, setGslMatrixElements(values.observed, params.observed, "y_noise_cov"), sep="\n")
	ret <- paste(ret, "\n}\n\n")

	return(ret)
}


#------------------------------------------------------------------------------
# Regime switching matrix/function

dynr.regimes <- function(){
	
}


#------------------------------------------------------------------------------
# "Dynamics" functions

# produces "drift" or "state-transition" matrices in continuous or discrete time, respectively.
# second input is for the jacobian of the drift/state-transition


dynr.dynamics <- function(){
	
}


dynr.linearDynamics <- function(){
	
}


# nonlinear functions for lookup
# exp, log, +, -, *, ^, **, sqrt,
# sin, cos, tan, asin, acos, atan,
# sinh, cosh, tanh, asinh, acosh, atanh


#------------------------------------------------------------------------------
dynr.dP_dt <- "/**\n * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end\n * but user does not need to modify it or care about it.\n */\nvoid mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert matrix to vector*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));\n\t\t\t/*printf(\"%lu\",i+j+nx-1);}*/\n\t\t}\n\t}\n}\nvoid mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){\n\tsize_t i,j;\n\tsize_t nx=mat->size1;\n\t/*convert vector to matrix*/\n\tfor(i=0; i<nx; i++){\n\t\tgsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));\n\t\tfor (j=i+1;j<nx;j++){\n\t\t\tgsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));\n\t\t\tgsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));\n\t\t}\n\t}\n}\nvoid function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){\n\t\n\tsize_t nx;\n\tnx = (size_t) floor(sqrt(2*(double) p->size));\n\tgsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(p,P_mat);\n\tgsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);\n\tfunction_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);\n\tgsl_matrix *dFP=gsl_matrix_calloc(nx,nx);\n\tgsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);\n\tgsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);\n\tgsl_matrix_transpose_memcpy(dP_dt, dFP);\n\tgsl_matrix_add(dP_dt, dFP);\n\tsize_t n_Q_vec=(1+nx)*nx/2;\n\tgsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);\n\tsize_t i;\n\tfor(i=1;i<=n_Q_vec;i++){\n\t\t\tgsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);\n\t}\n\tgsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);\n\tmathfunction_vec_to_mat(Q_vec,Q_mat);\n\tgsl_matrix_add(dP_dt, Q_mat);\n\tmathfunction_mat_to_vec(dP_dt, F_dP_dt);\n\tgsl_matrix_free(P_mat);\n\tgsl_matrix_free(F_dx_dt_dx);\n\tgsl_matrix_free(dFP);\n\tgsl_matrix_free(dP_dt);\n\tgsl_vector_free(Q_vec);\n\tgsl_matrix_free(Q_mat);\n}\n"


#------------------------------------------------------------------------------
#  Utility function written by Lu as a wrapper to take 
function_gsl_matrix_set(Pattern, StartVal, Fit, MatrixName){
	code=""
	if (Pattern=="Symm"){
		for (index_col in 1:ncol){
			if (fit){
				code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col, param[",index_col*(index_col+1)+index_col,"]);"),""), collapse="\n")
			}else{
				code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col,", StartVal[index_col*(index_col+1)+index_col],");"),""), collapse="\n")
			}
			for (index_row in (index_col+1):nrow){
				if (fit){
					code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col, param[",index_row*(index_col+1)+index_col,"]);"),""), collapse="\n")
				}else{
					code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col,", StartVal[index_row*(index_col+1)+index_col],");"),""), collapse="\n")
					code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_col, index_row,", StartVal[index_row*(index_col+1)+index_col],");"),""), collapse="\n")
				}
			}
		}
	}else if (Pattern=="Diag"){
		for (index_col in 1:ncol){
			if (fit){
				code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col, param[",index_col,"]);"),""), collapse="\n")
			}else{
				code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col,", StartVal[index_col],");"),""), collapse="\n")
			}
		}
	}else{#Full
		for (index_col in 1:ncol){
			for (index_row in 1:nrow){
				if (fit){
					code <- paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col, param[",index_row*ncol+index_col,"]);"),""), collapse="\n")
				}else{
					code <-paste(c(code,paste0("\t\t\tgsl_matrix_set(", MatrixName, ", index_row, index_col,", StartVal[index_row*ncol+index_col],");"),""), collapse="\n")
				}
			}
		}
	}
	return(code)
}



setGslMatrixElements <- function(values, params, name){
	ret <- ""
	numRow <- nrow(values)
	numCol <- ncol(values)
	for(j in 1:numCol){
		for(i in 1:numRow){
			if(params[i, j] > 0){
				ret <- paste(ret,
					'    gsl_matrix_set(', name, ', ', i-1, ', ', j-1,
					', param[', params[i, j] - 1, ']);\n', sep='')
			} else if(values[i, j] != 0){
				ret <- paste(ret,
					'    gsl_matrix_set(', name, ', ', i-1, ', ', j-1,
					', ', values[i, j], ');\n', sep='')
			}
		}
	}
	return(ret)
}



