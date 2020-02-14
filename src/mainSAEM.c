/*
Authors: Hui-Ju Hung
Date: 2019-05-20
Filename: mainSAEM.c
Purpose: Hello World for integrating SAEM and dynr 
*/

#include <math.h> /*sqrt(double)*/
/*#include <stdio.h>*/
#include <stdbool.h> /*for bool type*/
#include <string.h>
/* #include "nlopt.h"
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
#include "estimation.h" */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "print_function.h" 
/*#include "example1.h"*/
/*#include "example2.h"*/

/*
[ASK] How to compile the cpp file
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
*/
/*#include "structure_prototype.h"*/

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
}
*/
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
/*extern "C" SEXP main_SAEM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);*/

SEXP main_SAEM(SEXP model_list, SEXP data_list, SEXP weight_flag_in, SEXP debug_flag_in, SEXP optimization_flag_in, SEXP hessian_flag_in, SEXP verbose_flag_in)
{

	/*has error ask*/
	int total_time_all_subj_int = 60000;
	int i;
	
	double *timeDiscrete;
	printf("total_time_all_subj %d\n", total_time_all_subj_int);
	if (total_time_all_subj_int  > 0){
		
		timeDiscrete = (double *)malloc((total_time_all_subj_int + 1)* sizeof(double));
		printf("timeDiscrete(1):\n");
		//timeDiscrete = REAL(PROTECT(getListElement(model_list, "time_")));
		timeDiscrete = REAL(PROTECT(getListElement(data_list, "time")));
		UNPROTECT(1);
		printf("timeDiscrete(2):\n");
		
		for(i = 0; i < 10;i++){
			printf(" %lf", timeDiscrete[i]);
			
		}
		printf("\n");
		
		

    }else{
        timeDiscrete = NULL;
    }
	
	
	/*Inconsistent variables*/
	//printf("Nbeta %d NLambda %d\n", Nbeta, NLambda);

	
	
	
	printf("SAEM process starts\n");
	/*saem_interface(100, Nsubj, NxState, Ny, Nu, Ntheta, Nbeta, totalT, NLambda, Nmu, Nb, delt, U1, b, H, Z, maxT, allT, y0, lb, ub, MAXGIB, MAXITER, maxIterStage1, gainpara, gainparb, gainpara1, gainparb1, bAdaptParams, Nbpar, mu, tspan, lower_bound, upper_bound, Lambda, dmudparMu, dmudparMu2, num_time, Y, tobs, timeDiscrete, Sigmab);
	*/
	
	SEXP out = PROTECT(allocVector(REALSXP, 3));
	for (int i = 0; i < 3; i++) {
		REAL(out)[i] = 1.5;
	}
	UNPROTECT(1);
	
	
	

	return out;

}




