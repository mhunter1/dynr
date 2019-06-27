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

	/*static Data_and_Model data_model; */
	/*data_model.pc.verbose_flag = (bool) verbose_flag;*/
	/* C_INFDS InfDS;*/

    /* From the SEXP called model_list, get the list element named "num_sbj" */
	/*number of subjects: Nsubj*/
	SEXP num_sbj_sexp = PROTECT(getListElement(model_list, "num_sbj"));
	int Nsubj=(size_t) *INTEGER(num_sbj_sexp);
	DYNRPRINT(verbose_flag, "Nsubj: %lu\n", (long unsigned int) Nsubj);
		
	/*number of latent variables: NxState*/
	SEXP dim_latent_var_sexp = PROTECT(getListElement(model_list, "dim_latent_var"));
	int NxState=(size_t) *INTEGER(dim_latent_var_sexp);
	DYNRPRINT(verbose_flag, "NxState: %lu\n", (long unsigned int) NxState);
	
	
	/*number of observed variables: Ny*/
	SEXP dim_obs_var_sexp = PROTECT(getListElement(model_list, "dim_obs_var"));
	int Ny=(size_t) *INTEGER(dim_obs_var_sexp);
	DYNRPRINT(verbose_flag, "Ny: %lu\n", (long unsigned int) Ny);
	
	
	/*number of covariates: Nu*/
	SEXP dim_co_variate_sexp = PROTECT(getListElement(model_list, "dim_co_variate"));
	int Nu =(size_t) *INTEGER(dim_co_variate_sexp);
	DYNRPRINT(verbose_flag, "Nu: %lu\n", (long unsigned int) Nu);
	
	
	/*number of regimes: always 1 in SAEM*/
	SEXP num_regime_sexp = PROTECT(getListElement(model_list, "num_regime"));
	int num_regime =(size_t) *INTEGER(num_regime_sexp);
	DYNRPRINT(verbose_flag, "num_regime: %lu\n", (long unsigned int) num_regime);
	
	
	/*Ntheta*/
	SEXP num_theta_sexp = PROTECT(getListElement(model_list, "num_theta"));
	int Ntheta=(size_t) *INTEGER(num_theta_sexp);
	DYNRPRINT(verbose_flag, "Ntheta: %lu\n", (long unsigned int) Ntheta);
	
	/*Nbeta*/
	SEXP num_beta_sexp = PROTECT(getListElement(model_list, "num_beta"));
	int Nbeta =(size_t) *INTEGER(num_beta_sexp);
	DYNRPRINT(verbose_flag, "Nbeta: %lu\n", (long unsigned int) Nbeta);
	
	/*totalT*/
	SEXP total_t_sexp = PROTECT(getListElement(model_list, "total_t"));
	int totalT =(size_t) *INTEGER(total_t_sexp);
	DYNRPRINT(verbose_flag, "totalT: %lu\n", (long unsigned int) totalT);
	
	
	/*maxT*/
	SEXP max_t_sexp = PROTECT(getListElement(model_list, "max_t"));
	double maxT =*REAL(max_t_sexp);
	DYNRPRINT(verbose_flag, "maxT: %lf\n", maxT);

	
	/*NLambda*/
	SEXP num_lambda_sexp = PROTECT(getListElement(model_list, "num_lambda"));
	int NLambda =(size_t) *INTEGER(num_lambda_sexp);
	DYNRPRINT(verbose_flag, "NLambda: %lu\n", (long unsigned int) NLambda);
	
	
	/*maxT*/
	SEXP delt_sexp = PROTECT(getListElement(model_list, "delt"));
	double delt=*REAL(delt_sexp);
	DYNRPRINT(verbose_flag, "delt: %lf\n", delt);
	
	/*Nmu*/
	SEXP num_mu_sexp = PROTECT(getListElement(model_list, "num_mu"));
	int Nmu = (size_t) *INTEGER(num_mu_sexp);
	DYNRPRINT(verbose_flag, "Nmu: %lu\n", (long unsigned int) Nmu);
	
	/*Nb*/
	SEXP num_random_sexp = PROTECT(getListElement(model_list, "num_random"));
	int Nb = (size_t) *INTEGER(num_mu_sexp);
	DYNRPRINT(verbose_flag, "Nb: %lu\n", (long unsigned int) Nb);
	
	/*maxT*/
	SEXP kko_sexp = PROTECT(getListElement(model_list, "KKO"));
	double KKO=*REAL(kko_sexp);
	DYNRPRINT(verbose_flag, "KKO: %lf\n", KKO);
	
	UNPROTECT(14);
	/*----------*/
	
	/*U1: covariate matrix*/ 
	double **U1;
	/*
	if (Nu > 0){
        for(t=0; t<data_model.pc.total_obs; t++){
            data_model.co_variate[t]=gsl_vector_calloc(data_model.pc.dim_co_variate);
        }
        for(index=0;index<data_model.pc.dim_co_variate;index++){
            sprintf(str_number, "%lu", (long unsigned int) index+1);
            sprintf(str_name, "%s", "covar");
	        ptr_index=REAL(PROTECT(getListElement(covariates_sexp, strncat(str_name, str_number, strlen(str_number)))));
            for(t=0; t<data_model.pc.total_obs; t++){
                gsl_vector_set(data_model.co_variate[t],index, ptr_index[t]);
            }
			UNPROTECT(1);
        }
		
    }else{
        data_model.co_variate=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
        
        for(t=0; t<data_model.pc.total_obs; t++){
            data_model.co_variate[t]=NULL;
        }
    }
	*/
	
	int row, col;
	double *temp, **P0;
	char str_number[64], str_name[64];
	

	if (NxState > 0){
		temp = (double *)malloc(NxState * NxState* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"p0")));
		UNPROTECT(1);
		
		P0 = (double **)malloc((NxState + 1)* sizeof(double *));
		for(row = 0;row < NxState; row++){
			for(col = 0;col < NxState; col++){
				if(col == 0){
					P0[row] = (double *)malloc((NxState + 1)* sizeof(double));
				}
				P0[row][col] = temp[col * NxState + row];
			}
		}
		
		printf("P0:\n");
		for(row = 0; row < NxState; row++){
			for(col = 0;col < NxState; col++){
				printf(" %lf", P0[row][col]);
			}
			printf("\n");
		}	
    }
	else{
		P0 = NULL;
    }
	
	double *bAdaptParams;
	bAdaptParams = (double *)malloc(4* sizeof(double));
	bAdaptParams = REAL(PROTECT(getListElement(model_list,"bAdaptParams")));
	UNPROTECT(1);
	printf("bAdaptParams: %lf %lf %lf\n",bAdaptParams[0], bAdaptParams[1], bAdaptParams[2]);
	
	
	SEXP out = PROTECT(allocVector(REALSXP, 3));
	for (int i = 0; i < 3; i++) {
		REAL(out)[i] = 1.5;
	}
	UNPROTECT(1);
	
	//printf("here");
	return out;

}



