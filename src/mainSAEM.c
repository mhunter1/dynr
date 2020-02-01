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
	/*Nbetax -- seems to be the same as Nbeta here*/
	int Nbetax = Nbeta;
	
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
	int Nb = (size_t) *INTEGER(num_random_sexp);
	DYNRPRINT(verbose_flag, "Nb: %lu\n", (long unsigned int) Nb);
	
	/*KKO*/ /*not feed*/
	SEXP kko_sexp = PROTECT(getListElement(model_list, "KKO"));
	double KKO=*REAL(kko_sexp);
	DYNRPRINT(verbose_flag, "KKO: %lf\n", KKO);
	
	/*Nbpar*/
	SEXP num_bpar_sexp = PROTECT(getListElement(model_list, "num_bpar"));
	int Nbpar = (size_t) *INTEGER(num_bpar_sexp);
	DYNRPRINT(verbose_flag, "Nbpar: %lu\n", (long unsigned int) Nbpar);
	
	/*num_time: length(tspan)*/
	SEXP num_time_sexp = PROTECT(getListElement(model_list, "num_time"));
	int num_time=*INTEGER(num_time_sexp);
	DYNRPRINT(verbose_flag, "num_time: %lu\n", num_time);
	
	UNPROTECT(16);
	
	/*----------*/
	SEXP lb_sexp = PROTECT(getListElement(model_list, "random.lb"));
	double lb = *REAL(lb_sexp);
	SEXP ub_sexp = PROTECT(getListElement(model_list, "random.ub"));
	double ub = *REAL(ub_sexp);
	
	DYNRPRINT(verbose_flag, "lb: %lf\n", lb);
	DYNRPRINT(verbose_flag, "ub: %lf\n", ub);
	
	UNPROTECT(2);
	
	/*----------*/
	SEXP maxgib_sexp = PROTECT(getListElement(model_list, "MAXGIB"));
	int MAXGIB =(size_t) *INTEGER(maxgib_sexp);
	
	SEXP maxiter_sexp = PROTECT(getListElement(model_list, "MAXITER"));
	int MAXITER =(size_t) *INTEGER(maxiter_sexp);
	
	SEXP maxiterstage1_sexp = PROTECT(getListElement(model_list, "maxIterStage1"));
	int maxIterStage1 =(size_t) *INTEGER(maxiterstage1_sexp);
	
	SEXP gainpara_sexp = PROTECT(getListElement(model_list, "gainpara"));
	double gainpara = *REAL(gainpara_sexp);
	
	SEXP gainparb_sexp = PROTECT(getListElement(model_list, "gainparb"));
	double gainparb = *REAL(gainparb_sexp);
	
	SEXP gainpara1_sexp = PROTECT(getListElement(model_list, "gainpara1"));
	double gainpara1 = *REAL(gainpara1_sexp);
	
	SEXP gainparb1_sexp = PROTECT(getListElement(model_list, "gainparb1"));
	double gainparb1 = *REAL(gainparb1_sexp);
	
	DYNRPRINT(verbose_flag, "[SAEM Parameters] %d %d %d %lf %lf %lf %lf\n", MAXGIB, MAXITER, maxIterStage1, gainpara, gainparb, gainpara1, gainparb1);
	
	UNPROTECT(7);
	
	/*----------*/
	
	/*U1: covariate matrix*/ 
	double **U1;
	int row, col;
	int i, j, y, u;
	double *temp, **P0, **Lambda, **Y, **b, **H, **Z, **dmudparMu, **dmudparMu2;
	double *lower_bound, *upper_bound, *tspan, *mu;
	double *allT;
	char str_number[64], str_name[64];
	double total_time_all_subj;
	
	/*U1: covariate matrix*/ 
	if (Nu > 0){
		U1 = (double **)malloc((Nsubj + 1)* sizeof(double *));
		for(row = 0; row < Nsubj; row++){
			U1[row] = (double *)malloc((Nu+1)* sizeof(double));
		}
		
		temp = (double *)malloc(Nsubj * totalT * 2* sizeof(double));
		SEXP covariates_sexp = PROTECT(getListElement(data_list, "covariates"));
		UNPROTECT(1);
		for(u = 0;u < Nu; u++){
			sprintf(str_name, "covar%u", (long unsigned int) u+1);
			temp = REAL(PROTECT(getListElement(covariates_sexp, str_name)));
			UNPROTECT(1);
			for(i = 0; i < Nsubj;i++){
				U1[i][u] = temp[i *totalT];
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
		temp = (double *)malloc((NxState * Ny + 1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"lambda")));
		UNPROTECT(1);
		
		Lambda = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny; row++){
			for(col = 0;col < NxState; col++){
				if(col == 0){
					Lambda[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				Lambda[row][col] = temp[col * Ny + row];
			}
		}
		
		/*
		printf("lambda:\n");
		for(row = 0; row < Ny; row++){
			for(col = 0;col < NxState; col++){
				printf(" %lf", Lambda[row][col]);
			}
			printf("\n");
		}
		*/
		
    }
	else{
		Lambda = NULL;
    }
	
	double *bAdaptParams;
	bAdaptParams = (double *)malloc(4* sizeof(double));
	bAdaptParams = REAL(PROTECT(getListElement(model_list,"bAdaptParams")));
	UNPROTECT(1);
	//printf("bAdaptParams: %lf %lf %lf\n",bAdaptParams[0], bAdaptParams[1], bAdaptParams[2]);
	
	
	
	
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
		
		
		printf("Z:\n");
		for(int i = 0; i < Ntheta;i++){
			for(int u = 0;u < Nb; u++){
				printf(" %6lf", Z[i][u]);
			}
			printf("\n");
		}
		
    }else{
        Z = NULL;
    }
	
	//printf("allT\n");
	int total_time_all_subj_int;
	if (Nsubj > 0){
		allT = (double *)malloc((Nsubj+1)* sizeof(double));
		allT = REAL(PROTECT(getListElement(model_list,"allT")));
		UNPROTECT(1);
		
		
		//printf("allT:\n");
		total_time_all_subj = 0;
		for(row = 0; row < Nsubj; row++){
			//printf("%lf  ", allT[row]);
			total_time_all_subj += allT[row];
			//printf("total_time_all_subj %lf\n", total_time_all_subj);
		}
		//printf("\n");
		
		total_time_all_subj_int = (int)(total_time_all_subj + 0.001);
		printf("total_time_all_subj %d\n", total_time_all_subj_int);
    }

	printf("here %d\n", totalT);
	/*tspan*/
	if (totalT > 0){
		tspan = (double *)malloc((num_time+1)* sizeof(double));
		tspan = REAL(PROTECT(getListElement(model_list,"tspan")));
		UNPROTECT(1);
		
		/*
		printf("tspan:\n");
		for(row = 0; row < num_time; row++)
			printf("%lf ", tspan[row]);
		printf("\n");
		*/
		
    }
	else{
		tspan = NULL;
    }
	
	//printf("mu\n");
	/*mu*/
	if (Nmu > 0){
		mu = (double *)malloc((Nmu+1)* sizeof(double));
		mu = REAL(PROTECT(getListElement(model_list,"mu")));
		UNPROTECT(1);
		
		/*
		printf("mu:\n");
		for(row = 0; row < Nmu; row++)
			printf("%lf ", mu[row]);
		printf("\n");
		*/
    }
	else{
		mu = NULL;
    }
	
	
	if(Nmu > 0 || Ny > 0 || NLambda > 0 || Nbeta > 0 || Nbpar > 0){
		int Npar = Nmu + Ny + NLambda + Nbeta + Nbpar;
		
		lower_bound = (double *)malloc((Npar+1)* sizeof(double));
		lower_bound = REAL(PROTECT(getListElement(model_list,"lower_bound")));
		UNPROTECT(1);
		
		upper_bound = (double *)malloc((Npar+1)* sizeof(double));
		upper_bound = REAL(PROTECT(getListElement(model_list,"upper_bound")));
		UNPROTECT(1);
		
		
		printf("lower and upper bound:\n");
		for(row = 0; row < Npar; row++)
			printf("%6lf %6lf\n", lower_bound[row], upper_bound[row]);
		printf("\n");
	}
	
	
	
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny  + 1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"dmudparMu")));
		UNPROTECT(1);
		
		dmudparMu = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny ; row++){
			for(col = 0;col < Ny; col++){
				if(col == 0){
					dmudparMu[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				dmudparMu[row][col] = temp[col * (Ny) + row];
			}
		}
		
		/*
		printf("dmudparMu:\n");
		for(row = 0; row < Ny; row++){
			for(col = 0;col < Ny; col++){
				printf(" %lf", dmudparMu[row][col]);
			}
			printf("\n");
		}
		*/		
	
    }
	else{
		dmudparMu = NULL;
    }
	
	
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny * Ny + 1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"dmudparMu2")));
		UNPROTECT(1);
		
		dmudparMu2 = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny * Ny; row++){
			for(col = 0;col < Ny; col++){
				if(col == 0){
					dmudparMu2[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				dmudparMu2[row][col] = temp[col * (Ny * Ny) + row];
			}
		}
		
		/*
		printf("dmudparMu2:\n");
		for(row = 0;row < Ny * Ny; row++){
			for(col = 0;col < Ny; col++){
				printf(" %lf", dmudparMu2[row][col]);
			}
			printf("\n");
		}
		*/		
	
    }
	else{
		dmudparMu2 = NULL;
    }
	
	
	//Z
	double **y0;
	if (Nsubj > 0 && NxState > 0){
		temp = (double *)malloc(Nsubj * NxState* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"y0")));
		UNPROTECT(1);
		
		y0 = (double **)malloc((Nsubj + 1)* sizeof(double *));
		for(row = 0;row < Nsubj; row++){
			for(col = 0;col < NxState; col++){
				if(col == 0){
					y0[row] = (double *)malloc((NxState + 1)* sizeof(double));
				}
				y0[row][col] = temp[col * Nsubj + row];
			}
		}
		
		
		printf("y0:\n");
		for(i = 0; i < Nsubj;i++){
			for(u = 0;u < NxState; u++){
				printf(" %6lf", y0[i][u]);
			}
			printf("\n");
		}
		
		
    }else{
        y0 = NULL;
    }
	

	printf("time %d\n", total_time_all_subj_int);
	
	int *tobs;
	printf("total_time_all_subj %d\n", total_time_all_subj_int);
	if (total_time_all_subj_int  > 0){		
		tobs = (int *)malloc((total_time_all_subj_int + 1)* sizeof(int));
		tobs = INTEGER(PROTECT(getListElement(model_list,"tobs")));
		UNPROTECT(1);
	
		
		printf("tobs:\n");
		for(i = 59990; i < total_time_all_subj_int;i++){
			printf(" %d", tobs[i]);
		}
		printf("\n");
		
    }else{
        tobs = NULL;
    }
	
	
	if(Ny > 0 && total_time_all_subj_int > 0){
		Y = (double **)malloc((Ny + 1)* sizeof(double));
		SEXP observed_sexp = PROTECT(getListElement(data_list, "observed"));
		UNPROTECT(1);
		for(i = 0; i < Ny; i++){
			Y[i] = (double *)malloc((total_time_all_subj_int + 1)* sizeof(double));
			sprintf(str_name, "obs%u", (long unsigned int) i+1);
			Y[i] = REAL(PROTECT(getListElement(observed_sexp, str_name)));
			UNPROTECT(1);
		}
		
		
		printf("Y:\n");
		for(u = 0;u < Ny; u++){
			for(i = 0; i < 10;i++){
				printf(" %6lf", Y[u][i]);
			}
			printf("\n");
		}
		
		
	}
	else{
		Y = NULL;
	}
	
	
	
	double **Sigmab;
	if (Nb > 0){
		temp = (double *)malloc((Nb * Nb +  1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"sigmab")));
		UNPROTECT(1);
		
		Sigmab = (double **)malloc((Nb + 1)* sizeof(double *));
		for(row = 0;row < Nb; row++){
			for(col = 0;col < Nb; col++){
				if(col == 0){
					Sigmab[row] = (double *)malloc((Nb + 1)* sizeof(double));
				}
				Sigmab[row][col] = temp[col * Nb + row];
			}
		}
 		printf("Sigmab:\n");
 		for(row = 0; row < Ny; row++){
 			for(col = 0;col < Ny; col++){
 				printf(" %lf", Sigmab[row][col]);
 			}
 		printf("\n");
		}
	}
	else{
		Sigmab = NULL;
	}
	
	
	double **Sigmae;
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny +  1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"sigmae")));
		UNPROTECT(1);
		
		Sigmae = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny; row++){
			for(col = 0;col < Ny; col++){
				if(col == 0){
					Sigmab[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				Sigmab[row][col] = temp[col * Nb + row];
			}
		}
		
		
 		printf("Sigmae:\n");
 		for(row = 0; row < Ny; row++){
 			for(col = 0;col < Ny; col++){
 				printf(" %lf", Sigmab[row][col]);
 			}
 		printf("\n");
		}	
 																			
	}
	else{
		Sigmae = NULL;
	}

	/*----*/
	
	double **dSigmaede;
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny * Ny +  1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"dSigmaede")));
		UNPROTECT(1);
		
		dSigmaede = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny; row++){
			for(col = 0;col < Ny *  Ny; col++){
				if(col == 0){
					dSigmaede[row] = (double *)malloc((Ny *  Ny + 1)* sizeof(double));
				}
				dSigmaede[row][col] = temp[col * Ny + row];
			}
		}
		
		
 		printf("dSigmaede:\n");
 		for(row = 0; row < Ny; row++){
 			for(col = 0;col < Ny * Ny; col++){
 				printf(" %lf", dSigmaede[row][col]);
 			}
 		printf("\n");
		}	
 																			
	}
	else{
		dSigmaede = NULL;
	}
	
	double **dSigmaede2;
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny * Ny * Ny +  1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"dSigmaede2")));
		UNPROTECT(1);
		
		dSigmaede2 = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny * Ny * Ny; row++){
			for(col = 0;col < Ny; col++){
				if(col == 0){
					dSigmaede2[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				dSigmaede2[row][col] = temp[col * Ny + row];
			}
		}
		
		
 		printf("dSigmaede2:\n");
 		for(row = 0; row < Ny * Ny * Ny; row++){
 			for(col = 0;col < Ny; col++){
 				printf(" %lf", dSigmaede2[row][col]);
 			}
 		printf("\n");
		}	
 																			
	}
	else{
		dSigmaede2 = NULL;
	}
	
	/*
	double **dSigmabdb;
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny +  1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"dSigmabdb")));
		UNPROTECT(1);
		
		dSigmabdb = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny; row++){
			for(col = 0;col < Ny; col++){
				if(col == 0){
					dSigmabdb[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				dSigmabdb[row][col] = temp[col * Nb + row];
			}
		}
		
		
 		printf("dSigmabdb:\n");
 		for(row = 0; row < Ny; row++){
 			for(col = 0;col < Ny; col++){
 				printf(" %lf", dSigmabdb[row][col]);
 			}
 		printf("\n");
		}	
 																			
	}
	else{
		dSigmabdb = NULL;
	}
	
	double **dSigmabdb2;
	if (Ny > 0){
		temp = (double *)malloc((Ny * Ny +  1)* sizeof(double));
		temp = REAL(PROTECT(getListElement(model_list,"dSigmabdb2")));
		UNPROTECT(1);
		
		dSigmabdb2 = (double **)malloc((Ny + 1)* sizeof(double *));
		for(row = 0;row < Ny; row++){
			for(col = 0;col < Ny; col++){
				if(col == 0){
					dSigmabdb2[row] = (double *)malloc((Ny + 1)* sizeof(double));
				}
				dSigmabdb2[row][col] = temp[col * Nb + row];
			}
		}
		
		
 		printf("dSigmabdb2:\n");
 		for(row = 0; row < Ny; row++){
 			for(col = 0;col < Ny; col++){
 				printf(" %lf", dSigmabdb2[row][col]);
 			}
 		printf("\n");
		}	
 																			
	}
	else{
		dSigmabdb2 = NULL;
	}
	*/

	/*has error ask*/
	
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




