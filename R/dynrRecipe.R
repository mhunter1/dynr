# A Recipe is a helper function that takes user input
#  and produces a C function definition that can be
#  compiled by dynrFuncaddress.

dynr.loadings <- function(map, paramNums){
	nx <- length(unlist(map))
	ne <- length(map)
	
	SOME_CONDITION_MET <- TRUE
	SOMETHING_WITH_PARAMNUMS <- 1
	
	ret <- "void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){\n\n"
	for(i in 1:nx){
		for(j in 1:ne){
			if(SOME_CONDITION_MET){
				ret <- paste(ret,
					'    gsl_matrix_set(Ht, ', i-1, ', ', j-1,
					', param[', SOMETHING_WITH_PARAMNUMS, ']);\n', sep='')
			}
		}
	}
	ret <- paste(ret, "\n}\n\n")
	cat(ret)
}

#dynr.loadings( list(eta1=paste0('y', 1:4)) )

