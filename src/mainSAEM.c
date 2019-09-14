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
	/*DYNRPRINT(verbose_flag, "Nsubj: %lu\n", (long unsigned int) Nsubj);*/
		
	/*number of latent variables: NxState*/
	SEXP dim_latent_var_sexp = PROTECT(getListElement(model_list, "dim_latent_var"));
	int NxState=(size_t) *INTEGER(dim_latent_var_sexp);
	/*DYNRPRINT(verbose_flag, "NxState: %lu\n", (long unsigned int) NxState);*/
	
	
	/*number of observed variables: Ny*/
	SEXP dim_obs_var_sexp = PROTECT(getListElement(model_list, "dim_obs_var"));
	int Ny=(size_t) *INTEGER(dim_obs_var_sexp);
	/*DYNRPRINT(verbose_flag, "Ny: %lu\n", (long unsigned int) Ny);*/
	
	
	/*number of covariates: Nu*/
	SEXP dim_co_variate_sexp = PROTECT(getListElement(model_list, "dim_co_variate"));
	int Nu =(size_t) *INTEGER(dim_co_variate_sexp);
	/*DYNRPRINT(verbose_flag, "Nu: %lu\n", (long unsigned int) Nu);*/
	
	
	/*number of regimes: always 1 in SAEM*/
	SEXP num_regime_sexp = PROTECT(getListElement(model_list, "num_regime"));
	int num_regime =(size_t) *INTEGER(num_regime_sexp);
	/*DYNRPRINT(verbose_flag, "num_regime: %lu\n", (long unsigned int) num_regime);*/
	
	
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
	/*DYNRPRINT(verbose_flag, "maxT: %lf\n", maxT);*/

	
	/*NLambda*/
	SEXP num_lambda_sexp = PROTECT(getListElement(model_list, "num_lambda"));
	int NLambda =(size_t) *INTEGER(num_lambda_sexp);
	DYNRPRINT(verbose_flag, "NLambda: %lu\n", (long unsigned int) NLambda);
	
	
	/*delt*/
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
	
	/*KKO*/ /*not feed*/
	SEXP kko_sexp = PROTECT(getListElement(model_list, "KKO"));
	double KKO=*REAL(kko_sexp);
	/*DYNRPRINT(verbose_flag, "KKO: %lf\n", KKO);*/
	
	UNPROTECT(14);
	/*----------*/
	
	/*U1: covariate matrix*/ 
	double **U1;
	int row, col;
	double *temp, **P0, **Lamdba, **Y, **b, **H, **Z, **allT, *tspan;
	char str_number[64], str_name[64];
	
	
	if (Nu > 0){
		U1 = (double **)malloc((Nsubj + 1)* sizeof(double *));
		for(int row = 0; row < Nsubj; row++){
			U1[row] = (double *)malloc((Nu+1)* sizeof(double));
		}
		
		temp = (double *)malloc(Nsubj * totalT * 2* sizeof(double));
		SEXP covariates_sexp = PROTECT(getListElement(data_list, "covariates"));
		UNPROTECT(1);
		for(int u = 0;u < Nu; u++){
			sprintf(str_name, "covar%u", (long unsigned int) u+1);
			temp = REAL(PROTECT(getListElement(covariates_sexp, str_name)));
			UNPROTECT(1);
			for(int i = 0; i < Nsubj;i++){
				U1[i][u] = temp[i *totalT];
				/*printf("i =%d, u = %d %lf\n", i, u, U1[i][u]);
				printf("temp[%d] = %lf %lf\n", i, temp[i *totalT], U1[i][u]);*/
			}
		}
		/*
		printf("U1:\n");
		for(int i = 0; i < Nsubj;i++){
			for(int u = 0;u < Nu; u++){
				printf(" %6lf", U1[i][u]);
			}
			printf("\n");
		}
		*/
    }else{
        U1 = NULL;
    }
	
	
	
	
	/*not feed*/
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
		/*
		printf("P0:\n");
		for(row = 0; row < NxState; row++){
			for(col = 0;col < NxState; col++){
				printf(" %lf", P0[row][col]);
			}
			printf("\n");
		}	
		*/
    }
	else{
		P0 = NULL;
    }
	
	if (NxState > 0 && Ny > 0){
		temp = (double *)malloc(NxState * Ny* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"lambda")));
		UNPROTECT(1);
		
		Lamdba = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny; row++){
			for(col = 0;col < NxState; col++){
				if(col == 0){
					Lamdba[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				Lamdba[row][col] = temp[col * Ny + row];
			}
		}
		/*
		printf("Lamdba:\n");
		for(row = 0; row < Ny; row++){
			for(col = 0;col < NxState; col++){
				printf(" %lf", Lamdba[row][col]);
			}
			printf("\n");
		}	
		*/
    }
	else{
		Lamdba = NULL;
    }
	
	double *bAdaptParams;
	bAdaptParams = (double *)malloc(4* sizeof(double));
	bAdaptParams = REAL(PROTECT(getListElement(model_list,"bAdaptParams")));
	UNPROTECT(1);
	printf("bAdaptParams: %lf %lf %lf\n",bAdaptParams[0], bAdaptParams[1], bAdaptParams[2]);
	
	if (Ny > 0){
		Y = (double **)malloc((Ny + 1)* sizeof(double *));
		for(int row = 0; row < Ny; row++){
			Y[row] = (double *)malloc((Nsubj*totalT+1)* sizeof(double));
		}
		
		temp = (double *)malloc(Nsubj * totalT * 2* sizeof(double));
		SEXP observed_sexp = PROTECT(getListElement(data_list, "observed"));
		UNPROTECT(1);
		for(int u = 0;u < Ny; u++){
			sprintf(str_name, "obs%u", (long unsigned int) u+1);
			temp = REAL(PROTECT(getListElement(observed_sexp, str_name)));
			UNPROTECT(1);
			for(int i = 0; i < Nsubj * totalT;i++){
				Y[u][i] = temp[i];
				/*printf("i =%d, u = %d %lf\n", i, u, U1[i][u]);
				printf("temp[%d] = %lf %lf\n", i, temp[i *totalT], U1[i][u]);*/
			}
		}
		/*
		printf("Y:\n");
		for(int u = 0;u < Ny; u++){
			for(int i = 0; i < 10;i++){
				printf(" %6lf", Y[u][i]);
			}
			printf("\n");
		}*/
		
    }else{
        Y = NULL;
    }
	
	
	if (Nb > 0 && Nsubj > 0){
		temp = (double *)malloc(Nsubj * Nb* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"b")));
		UNPROTECT(1);
		
		b = (double **)malloc((Nsubj + 1)* sizeof(double *));
		for(row = 0;row < Nsubj; row++){
			for(col = 0;col < Nb; col++){
				if(col == 0){
					b[row] = (double *)malloc((Nb + 1)* sizeof(double));
				}
				b[row][col] = temp[col * Nsubj + row];
			}
		}
		/*
		printf("b:\n");
		for(int i = 0; i < Nsubj;i++){
			for(int u = 0;u < Nb; u++){
				printf(" %6lf", b[i][u]);
			}
			printf("\n");
		}
		*/
		
    }else{
        b = NULL;
    }
	
	//H
	/*Nbetax*/
	int Nbetax = 5;
	
	if (Ntheta > 0 && Nsubj > 0 && Nbetax > 0){
		temp = (double *)malloc(Ntheta * Nsubj * Nbetax * sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"H")));
		UNPROTECT(1);
		
		H = (double **)malloc((Ntheta * Nsubj  + 1)* sizeof(double *));
		for(row = 0;row < Ntheta * Nsubj; row++){
			for(col = 0;col < Nbetax; col++){
				if(col == 0){
					H[row] = (double *)malloc((Nbetax + 1)* sizeof(double));
				}
				H[row][col] = temp[col * Ntheta * Nsubj + row];
			}
		}
		
		/*
		printf("H:\n");
		for(int i = 0; i < Ntheta * Nsubj;i++){
			for(int u = 0;u < Nbetax; u++){
				printf(" %6lf", H[i][u]);
			}
			printf("\n");
		}
		*/
		
		
    }else{
        H = NULL;
    }
	
	//Z
	if (Nb > 0 && Ntheta > 0){
		temp = (double *)malloc(Ntheta * Nb* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"Z")));
		UNPROTECT(1);
		
		Z = (double **)malloc((Ntheta + 1)* sizeof(double *));
		for(row = 0;row < Ntheta; row++){
			for(col = 0;col < Nb; col++){
				if(col == 0){
					Z[row] = (double *)malloc((Nb + 1)* sizeof(double));
				}
				Z[row][col] = temp[col * Ntheta + row];
			}
		}
		
		/*
		printf("Z:\n");
		for(int i = 0; i < Ntheta;i++){
			for(int u = 0;u < Nb; u++){
				printf(" %6lf", Z[i][u]);
			}
			printf("\n");
		}
		*/
    }else{
        Z = NULL;
    }
	
	if (Nsubj > 0){
		allT = (double *)malloc((Nsubj+1)* sizeof(double));
		allT = REAL(PROTECT(getListElement(model_list,"allT")));
		UNPROTECT(1);
		
		/*
		printf("allT:\n");
		for(row = 0; row < Nsubj; row++)
			printf("%lf\n", allT[row]);
		*/
		
    }

	/*totalT*/
	if (totalT > 0){
		tspan = (double *)malloc((totalT+1)* sizeof(double));
		tspan = REAL(PROTECT(getListElement(model_list,"tspan")));
		UNPROTECT(1);
		
		/*
		printf("tspan:\n");
		for(row = 0; row < totalT; row++)
			printf("%lf ", tspan[row]);
		printf("\n");
		*/
		
		
    }
	
	/*Inconsistent variables*/
	Nbeta = 0;	
	NLambda = 2;
	Nbetax = 5;
	printf("Nbeta %d NLambda %d\n", Nbeta, NLambda);
	
	
	/*example1();*/
	printf("start to call MainUseThis %d %d %d\n", Nsubj, NxState, Ny);
	//interface(100, Nsubj, NxState, Ny, Nu, Ntheta, Nbeta, totalT, NLambda, Nmu, Nb, delt, U1, b, H, Z, maxT, allT, y0);
	
	
	SEXP out = PROTECT(allocVector(REALSXP, 3));
	for (int i = 0; i < 3; i++) {
		REAL(out)[i] = 1.5;
	}
	UNPROTECT(1);
	
	
	

	return out;

}




