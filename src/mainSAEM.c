/*
Authors: Hui-Ju Hung
Date: 2019-05-20
Filename: mainSAEM.c
Purpose: Hello World for integrating SAEM and dynr 
*/

#include <math.h>/*sqrt(double)*/
/*#include <stdio.h>*/
#include <string.h>
#include "nlopt.h"
#include "math_function.h"
#include "ekf.h"
#include "data_structure.h"
#include "brekfis.h"
#include "adaodesolver.h"
#include "model.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "wrappernegloglike.h"
#include "numeric_derivatives.h"
#include "estimation.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "print_function.h"

/* get the list element named str, or return NULL */
/*
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    size_t i;
    for (i = 0; i < length(list); i++){
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}
	return elmt;
}*/
/**
 * The gateway function for the R interface
 * @param model_list is a list in R of all model specifications.
 * @param data_list is a list in R of the outputs prepared by dynr.data()
 * @param weight_flag_in a flag for weighting the neg loglike function by individual data length
 * @param debug_flag_in a flag for returning a longer list of outputs for debugging purposes
 * @param optimization_flag_in a flag for running optimization
 * @param hessian_flag_in a flag for calculating hessian matrix
 * @param verbose_flag_in a flag of whether or not to print debugging statements before and during estimation.
 */
SEXP main_SAEM(SEXP model_list, SEXP data_list, SEXP weight_flag_in, SEXP debug_flag_in, SEXP optimization_flag_in, SEXP hessian_flag_in, SEXP verbose_flag_in)
{
	bool debug_flag=*LOGICAL(PROTECT(debug_flag_in));
	bool optimization_flag=*LOGICAL(PROTECT(optimization_flag_in));
	bool hessian_flag=*LOGICAL(PROTECT(hessian_flag_in));
	bool verbose_flag=*LOGICAL(PROTECT(verbose_flag_in));
	bool weight_flag=*LOGICAL(PROTECT(weight_flag_in));
	UNPROTECT(5);
    /** =======================Interface : Start to Set up the data and the model========================= **/
	

    /* From the SEXP called model_list, get the list element named "num_sbj" */
	static Data_and_Model data_model;
	data_model.pc.verbose_flag = (bool) verbose_flag;

    /* From the SEXP called model_list, get the list element named "num_sbj" */
	/*number of subjects*/
	SEXP num_sbj_sexp = PROTECT(getListElement(model_list, "num_sbj"));
	data_model.pc.num_sbj=(size_t) *INTEGER(num_sbj_sexp);
	DYNRPRINT(verbose_flag, "num_sbj: %lu\n", (long unsigned int) data_model.pc.num_sbj);
	
	/*number of function parameters*/
	SEXP num_func_param_sexp = PROTECT(getListElement(model_list, "num_func_param"));
	data_model.pc.num_func_param=(size_t) *INTEGER(num_func_param_sexp);
	DYNRPRINT(verbose_flag, "num_func_param: %lu\n", (long unsigned int) data_model.pc.num_func_param);
	
	/*number of latent variables*/
	SEXP dim_latent_var_sexp = PROTECT(getListElement(model_list, "dim_latent_var"));
	data_model.pc.dim_latent_var=(size_t) *INTEGER(dim_latent_var_sexp);
	DYNRPRINT(verbose_flag, "dim_latent_var: %lu\n", (long unsigned int) data_model.pc.dim_latent_var);
	
	/*number of observed variables*/
	SEXP dim_obs_var_sexp = PROTECT(getListElement(model_list, "dim_obs_var"));
	data_model.pc.dim_obs_var=(size_t) *INTEGER(dim_obs_var_sexp);
	DYNRPRINT(verbose_flag, "dim_obs_var: %lu\n", (long unsigned int) data_model.pc.dim_obs_var);
	
	/*number of covariates*/
	SEXP dim_co_variate_sexp = PROTECT(getListElement(model_list, "dim_co_variate"));
	data_model.pc.dim_co_variate=(size_t) *INTEGER(dim_co_variate_sexp);
	DYNRPRINT(verbose_flag, "dim_co_variate: %lu\n", (long unsigned int) data_model.pc.dim_co_variate);
	
	/*number of regimes*/
	SEXP num_regime_sexp = PROTECT(getListElement(model_list, "num_regime"));
	data_model.pc.num_regime=(size_t) *INTEGER(num_regime_sexp);
	DYNRPRINT(verbose_flag, "num_regime: %lu\n", (long unsigned int) data_model.pc.num_regime);
	UNPROTECT(6);
	
	
	
	
	SEXP out = PROTECT(allocVector(REALSXP, 3));
	for (int i = 0; i < 3; i++) {
		REAL(out)[i] = 1 + i;
	}
	UNPROTECT(1);
	
	
	return out;

}




