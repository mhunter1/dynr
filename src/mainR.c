/*
Authors: Lu Ou, Mike Hunter, Sy-Miin Chow
Date: 2015-07-30
Filename: estimation_nlopt.c
Purpose: Obtain parameters estimates by minimizing the negative log likelihood function
Run with
   ./runnlopt_onair.sh
   type in the command: gsl-config --cflags to find the compiler flag
   type in the command: gsl-config --libs to find the flag
Note
  You must have NLOpt installed on your machine to run this.
  Installing NLOpt basically entailed
    1.  Download the tar.gz
    2.  Extract the tar bar
    3.  cd into trunk of folder
    4.  ./configure
    5.  make
    6.  sudo make install
  You may also need to change -I/dir and -L/dir to the locations
    where NLOpt was installed.
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
/**
 * The gateway function for the R interface
 * @param model_list is a list in R of all model specifications.
 * @param data_list is a list in R of the outputs prepared by dynr.data()
 * @param weight_flag_in a flag for weighting the neg loglike function by individual data length
 * @param debug_flag_in a flag for returning a longer list of outputs for debugging purposes
 * @param outall_flag_in a flag for returning all possible outputs
 * @param verbose_flag_in a flag of whether or not to print debugging statements before and during estimation.
 */
SEXP main_R(SEXP model_list, SEXP data_list, SEXP weight_flag_in, SEXP debug_flag_in, SEXP outall_flag_in, SEXP verbose_flag_in)
{
    size_t index,index_col,index_row;
    bool debug_flag=*LOGICAL(PROTECT(debug_flag_in));
	bool outall_flag=*LOGICAL(PROTECT(outall_flag_in));
	bool verbose_flag=*LOGICAL(PROTECT(verbose_flag_in));
	bool weight_flag=*LOGICAL(PROTECT(weight_flag_in));
    /** =======================Interface : Start to Set up the data and the model========================= **/

    static Data_and_Model data_model;

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

    /*function specifications*/
    SEXP func_address_list = PROTECT(getListElement(model_list, "func_address"));
	
	SEXP f_measure_sexp = PROTECT(getListElement(func_address_list, "f_measure"));
	SEXP f_regime_switch_sexp = PROTECT(getListElement(func_address_list, "f_regime_switch"));
	SEXP f_noise_cov_sexp = PROTECT(getListElement(func_address_list, "f_noise_cov"));
	SEXP f_initial_condition_sexp = PROTECT(getListElement(func_address_list, "f_initial_condition"));
	SEXP f_transform_sexp = PROTECT(getListElement(func_address_list, "f_transform"));
    *(void **) (&data_model.pc.func_measure) = R_ExternalPtrAddr(f_measure_sexp);
    *(void **) (&data_model.pc.func_regime_switch) = R_ExternalPtrAddr(f_regime_switch_sexp);
    *(void **) (&data_model.pc.func_noise_cov) = R_ExternalPtrAddr(f_noise_cov_sexp);
    *(void **) (&data_model.pc.func_initial_condition) = R_ExternalPtrAddr(f_initial_condition_sexp);
    *(void **) (&data_model.pc.func_transform) = R_ExternalPtrAddr(f_transform_sexp);

/*
 *   data_model.pc.func_dx_dt=function_dx_dt;
 *   data_model.pc.func_dP_dt=function_dP_dt;
 *   data_model.pc.func_initial_condition=function_initial_condition;
 *   data_model.pc.func_regime_switch=function_regime_switch;
 *   data_model.pc.func_noise_cov=function_noise_cov;
 */
	/*whether a continuous-time model is used*/
	SEXP isContinuousTime_sexp = PROTECT(getListElement(model_list, "isContinuousTime"));
	data_model.pc.isContinuousTime=*LOGICAL(isContinuousTime_sexp);
	DYNRPRINT(verbose_flag, "isContinuousTime: %s\n", data_model.pc.isContinuousTime? "true" : "false");
	
    if (data_model.pc.isContinuousTime){
		SEXP f_dx_dt_sexp = PROTECT(getListElement(func_address_list, "f_dx_dt"));
		SEXP f_dF_dx_sexp = PROTECT(getListElement(func_address_list, "f_dF_dx"));
		SEXP f_dP_dt_sexp = PROTECT(getListElement(func_address_list, "f_dP_dt"));
	    *(void **) (&data_model.pc.func_dx_dt) = R_ExternalPtrAddr(f_dx_dt_sexp);
	    *(void **) (&data_model.pc.func_dF_dx) = R_ExternalPtrAddr(f_dF_dx_sexp);
	    *(void **) (&data_model.pc.func_dP_dt) = R_ExternalPtrAddr(f_dP_dt_sexp);
	    data_model.pc.adaodesolver=false;/*true: use adapative ode solver; false: RK4*/
	    if (data_model.pc.adaodesolver){
	        data_model.pc.func_dynam=function_dynam_ada;
	    }else{
	        data_model.pc.func_dynam=rk4_odesolver;
	    }
		data_model.pc.func_jacob_dynam=function_jacob_dynam_rk4;
    }else{
	    data_model.pc.func_dx_dt=NULL;
	    data_model.pc.func_dF_dx=NULL;
	    data_model.pc.func_dP_dt=NULL;
		SEXP f_dynamic_sexp = PROTECT(getListElement(func_address_list, "f_dynamic"));
		SEXP f_jacob_dynamic_sexp = PROTECT(getListElement(func_address_list, "f_jacob_dynamic"));
    	*(void **) (&data_model.pc.func_dynam) = R_ExternalPtrAddr(f_dynamic_sexp);
    	*(void **) (&data_model.pc.func_jacob_dynam) = R_ExternalPtrAddr(f_jacob_dynamic_sexp);
    }
	
    data_model.pc.isnegloglikeweightedbyT=weight_flag;
    data_model.pc.second_order=false;

    /*specify the start position for each subject: User always need to provide a txt file called tStart.txt*/
    /*for example, 500 time points for each sbj, specify 0 500 1000 ... 10000 also the end point*/
    /*n subjects -> n+1 indices*/
    data_model.pc.index_sbj=(size_t *)malloc((data_model.pc.num_sbj+1)*sizeof(size_t *));

    double *ptr_index;/*used for multiple times*/
	int *ptr_index_int;
    ptr_index_int=INTEGER(PROTECT(getListElement(data_list, "tstart")));
    for(index=0;index<=data_model.pc.num_sbj;index++){
        data_model.pc.index_sbj[index]= ptr_index_int[index];
    }
    /*DYNRPRINT(verbose_flag, "index_sbj 2: %lu\n", (long unsigned int) data_model.pc.index_sbj[1]);*/

    data_model.pc.total_obs=*(data_model.pc.index_sbj+data_model.pc.num_sbj);/*total observations for all subjects*/
    DYNRPRINT(verbose_flag, "total_obs: %lu\n", (long unsigned int) data_model.pc.total_obs);

    /** read in the data**/
	/*observed data*/
	SEXP observed_sexp = PROTECT(getListElement(data_list,"observed")); 
	/*covariates*/
	SEXP covariates_sexp = PROTECT(getListElement(data_list,"covariates"));
		
    data_model.y=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
    size_t t;
    for(t=0; t<data_model.pc.total_obs; t++){
        data_model.y[t]=gsl_vector_calloc(data_model.pc.dim_obs_var);/*y[t] corresponds to y(),which is a gsl_vector; loop through total_obj*/

    }

    size_t enough_length=(ceil(log10((double)data_model.pc.dim_obs_var>data_model.pc.dim_co_variate?data_model.pc.dim_obs_var:data_model.pc.dim_co_variate))+1)*sizeof(char);
    char *str_number=(char *)malloc(enough_length);
    char str_name[enough_length+7];

    for(index=0;index<data_model.pc.dim_obs_var;index++){
        sprintf(str_number, "%lu", (long unsigned int) index+1);
        sprintf(str_name, "%s", "obs");
        /*DYNRPRINT(verbose_flag, "The str_number is %s\n",str_number);
        DYNRPRINT(verbose_flag, "The str_name length is %lu\n",strlen(str_name));*/
        ptr_index=REAL(PROTECT(getListElement(observed_sexp, strncat(str_name, str_number, strlen(str_number)))));
        for(t=0; t<data_model.pc.total_obs; t++){
            gsl_vector_set(data_model.y[t],index, ptr_index[t]);
        }
		UNPROTECT(1);
    }


    if (data_model.pc.dim_co_variate > 0){
    data_model.co_variate=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));

        for(t=0; t<data_model.pc.total_obs; t++){
            data_model.co_variate[t]=gsl_vector_calloc(data_model.pc.dim_co_variate);
        }


        for(index=0;index<data_model.pc.dim_co_variate;index++){
            sprintf(str_number, "%lu", (long unsigned int) index+1);
            sprintf(str_name, "%s", "covar");
            /*DYNRPRINT(verbose_flag, "The str_number is %s\n",str_number);
            DYNRPRINT(verbose_flag, "The str_name length is %lu\n",strlen(str_name));*/
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

    data_model.y_time=(double *)malloc(data_model.pc.total_obs*sizeof(double));
        memcpy(data_model.y_time,REAL(PROTECT(getListElement(data_list, "time"))),data_model.pc.total_obs*sizeof(double));

    /*DYNRPRINT(verbose_flag, "In main_R:\n");
    print_vector(data_model.y[0]);
    DYNRPRINT(verbose_flag, "\n");
    print_vector(data_model.co_variate[0]);
    DYNRPRINT(verbose_flag, "\n");
    DYNRPRINT(verbose_flag, "y_time_1 is %lf\n",data_model.y_time[0]);
    */

    /** Optimization options **/
    SEXP option_list = PROTECT(getListElement(model_list, "options"));
        SEXP xtol_rel_sexp = PROTECT(getListElement(option_list, "xtol_rel"));
        SEXP stopval_sexp = PROTECT(getListElement(option_list, "stopval"));
        SEXP ftol_rel_sexp = PROTECT(getListElement(option_list, "ftol_rel"));
        SEXP ftol_abs_sexp = PROTECT(getListElement(option_list, "ftol_abs"));
        SEXP maxeval_sexp = PROTECT(getListElement(option_list, "maxeval"));
        SEXP maxtime_sexp = PROTECT(getListElement(option_list, "maxtime"));
    double *xtol_rel = REAL(xtol_rel_sexp);
    double *stopval = REAL(stopval_sexp);
    double *ftol_rel = REAL(ftol_rel_sexp);
    double *ftol_abs = REAL(ftol_abs_sexp);
    int *maxeval = INTEGER(maxeval_sexp);
    double *maxtime = REAL(maxtime_sexp);

    /** Optimization bounds and starting values **/

    double params[data_model.pc.num_func_param];
    	memcpy(params,REAL(PROTECT(getListElement(model_list, "xstart"))),sizeof(params));
    /*DYNRPRINT(verbose_flag, "Array paramvec allocated.\n");*/
    /*print_array(params,data_model.pc.num_func_param);*/
    /*DYNRPRINT(verbose_flag, "\n");*/
    double fittedpar[data_model.pc.num_func_param];
    	memcpy(fittedpar, params,sizeof(params));
    /*DYNRPRINT(verbose_flag, "Array params copied.\n");*/
    /*print_array(fittedpar,data_model.pc.num_func_param);*/
    /*DYNRPRINT(verbose_flag, "\n");*/
    double ub[data_model.pc.num_func_param];
    double lb[data_model.pc.num_func_param];
    memcpy(ub,REAL(PROTECT(getListElement(model_list, "ub"))),sizeof(ub));
    memcpy(lb,REAL(PROTECT(getListElement(model_list, "lb"))),sizeof(lb));
    int h;
    /*DYNRPRINT(verbose_flag, "ub value h %f\n", ub[1]);*/
    for (h=0; h < data_model.pc.num_func_param; h++){
      if (ub[h]==9999 || !R_FINITE(ub[h])){ub[h]=HUGE_VAL;}
      if (lb[h]==9999 || !R_FINITE(lb[h])){lb[h]=-HUGE_VAL;}
    }
    /*DYNRPRINT(verbose_flag, "Arrays allocated.\n");*/
    /*double ub[6] = {4, 4, 4, 4, 4, 4};
    double lb[6] = {-4,-4,-4,-4,-12, -12}; */
    /*double params[]={log(1),log(2),0,0,-10,-10};*//* some initial guess*/
    /*log(1.2)=0.1823216,log(1.8)=0.5877867,-0.5,-0.5,log(0.0001)=-9.21034,log(0.0001)=-9.21034*//* lower bounds */


    /** =======================Interface: Model and data set up========================= **/

    /** =================Optimization: start======================**/

    double minf; /* the minimum objective value, upon return */

    gsl_matrix *Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    gsl_matrix *inv_Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    int status=opt_nlopt(&data_model,data_model.pc.num_func_param,ub,lb,&minf,fittedpar,Hessian_mat,inv_Hessian_mat,xtol_rel,stopval,ftol_rel,ftol_abs,maxeval, maxtime);


	/*DYNRPRINT(verbose_flag, "Optimization done.\n");*/
    /** =================Optimization: done======================**/
    
    /** =================Extended Kim Filter and Smoother: start======================**/

    /*declaritions for C program*/

    size_t index_sbj_t,regime_j,regime_k;

	    /**initialization**/
	    /*input and output of filter & input of smooth: eta^k_it|t*/
	    gsl_vector ***eta_regime_j_t=(gsl_vector ***)malloc(data_model.pc.total_obs*sizeof(gsl_vector **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		eta_regime_j_t[index_sbj_t]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    eta_regime_j_t[index_sbj_t][regime_j]=gsl_vector_calloc(data_model.pc.dim_latent_var);
		}
	    }
	    /*input and output of filter & input of smooth: error_cov^k_it|t*/
	    gsl_matrix ***error_cov_regime_j_t=(gsl_matrix ***)malloc(data_model.pc.total_obs*sizeof(gsl_matrix **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_regime_j_t[index_sbj_t]=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    error_cov_regime_j_t[index_sbj_t][regime_j]=gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
		}
	    }

	    /*output of filter and input of smooth: eta^regime_jk_it|t-1*/
	    gsl_vector ****eta_regime_jk_pred=(gsl_vector ****)malloc(data_model.pc.total_obs*sizeof(gsl_vector ***));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		eta_regime_jk_pred[index_sbj_t]=(gsl_vector ***)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    eta_regime_jk_pred[index_sbj_t][regime_j]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
		}
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
		    eta_regime_jk_pred[index_sbj_t][regime_j][regime_k]=gsl_vector_calloc(data_model.pc.dim_latent_var);
		}
	    }
	    /*output of filter and input of smooth: error_cov^regime_jk_it|t-1*/
	    gsl_matrix ****error_cov_regime_jk_pred=(gsl_matrix ****)malloc(data_model.pc.total_obs*sizeof(gsl_matrix ***));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_regime_jk_pred[index_sbj_t]=(gsl_matrix ***)malloc(data_model.pc.num_regime*sizeof(gsl_matrix **));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    error_cov_regime_jk_pred[index_sbj_t][regime_j]=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
		}
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
			error_cov_regime_jk_pred[index_sbj_t][regime_j][regime_k]=gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
		}
	    }
	    /*output of filter and input of smooth: Pr(S_it=k|Y_it)*/
	    gsl_vector **pr_t=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		pr_t[index_sbj_t]=gsl_vector_calloc(data_model.pc.num_regime);
	    }
	    /*output of filter and input of smooth: Pr(S_it=k|Y_i,t-1)*/
	    gsl_vector **pr_t_given_t_minus_1=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		pr_t_given_t_minus_1[index_sbj_t]=gsl_vector_calloc(data_model.pc.num_regime);
	    }

	    /*output of smooth: eta_it|T*/
	    gsl_vector **eta_smooth=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		eta_smooth[index_sbj_t]=gsl_vector_calloc(data_model.pc.dim_latent_var);
	    }
	    /*output of smooth: error_cov_it|T*/
	    gsl_matrix **error_cov_smooth=(gsl_matrix **)malloc(data_model.pc.total_obs*sizeof(gsl_matrix *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_smooth[index_sbj_t]=gsl_matrix_calloc(data_model.pc.dim_latent_var,data_model.pc.dim_latent_var);
	    }
	    /*output of smooth: Pr(S_it=k|Y_iT)*/
	    gsl_vector **pr_T=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		pr_T[index_sbj_t]=gsl_vector_calloc(data_model.pc.num_regime);
	    }
	    /*output of smooth: Pr(S_i,t+1=h, S_it=k|Y_iT)*/
	    gsl_vector ***transprob_T=(gsl_vector ***)malloc(data_model.pc.total_obs*sizeof(gsl_vector **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		transprob_T[index_sbj_t]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    transprob_T[index_sbj_t][regime_j]=gsl_vector_calloc(data_model.pc.num_regime);
		}
	    }
		
	    /*output of filter: eta_it|t*/
	    gsl_vector **eta_t=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		eta_t[index_sbj_t] = gsl_vector_calloc(data_model.pc.dim_latent_var);
	    }
		/*output of filter: error_cov_it|t*/
	    gsl_matrix **error_cov_t=(gsl_matrix **)malloc(data_model.pc.total_obs*sizeof(gsl_matrix *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_t[index_sbj_t] = gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
		}
	    

	    /*output of filter: innovation vector*/
	    gsl_vector ****innov_v=(gsl_vector ****)malloc(data_model.pc.total_obs*sizeof(gsl_vector ***));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		innov_v[index_sbj_t]=(gsl_vector ***)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    innov_v[index_sbj_t][regime_j]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
		}
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
		    innov_v[index_sbj_t][regime_j][regime_k]=gsl_vector_calloc(data_model.pc.dim_obs_var);
		}
	    }
	    /*output of filter: inverse_residual_cov*/
	    gsl_matrix ****inv_residual_cov=(gsl_matrix ****)malloc(data_model.pc.total_obs*sizeof(gsl_matrix ***));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		inv_residual_cov[index_sbj_t]=(gsl_matrix ***)malloc(data_model.pc.num_regime*sizeof(gsl_matrix **));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    inv_residual_cov[index_sbj_t][regime_j]=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
		}
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
			inv_residual_cov[index_sbj_t][regime_j][regime_k]=gsl_matrix_calloc(data_model.pc.dim_obs_var, data_model.pc.dim_obs_var);
		}
	    }



	    /**Actual Functions**/
	    /** initialize regime parameter**/
	    ParamInit pi;
	    pi.eta_0=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    for(index=0;index<data_model.pc.num_regime;index++){
		(pi.eta_0)[index]=gsl_vector_calloc(data_model.pc.num_sbj*data_model.pc.dim_latent_var);
	    }

	    pi.error_cov_0=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
	    for(index=0;index<data_model.pc.num_regime;index++){
	    	    (pi.error_cov_0)[index]=gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
	    }
	    pi.pr_0=gsl_vector_calloc(data_model.pc.num_regime);

	    /** set parameter **/
	    Param par;

	    /*function parameters*/
	    par.func_param=(double *)malloc(data_model.pc.num_func_param*sizeof(double));
            size_t i;
	    for(i=0;i<data_model.pc.num_func_param;i++){
	    	    par.func_param[i]=fittedpar[i];
	    }


	    data_model.pc.func_initial_condition(par.func_param, data_model.co_variate, pi.pr_0, pi.eta_0, pi.error_cov_0, data_model.pc.index_sbj);

	    par.eta_noise_cov=gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
	    par.y_noise_cov=gsl_matrix_calloc(data_model.pc.dim_obs_var, data_model.pc.dim_obs_var);
	    par.regime_switch_mat=gsl_matrix_calloc(data_model.pc.num_regime, data_model.pc.num_regime);

	    /* calculate the log_like */
	    data_model.pc.func_transform(par.func_param);

	    model_constraint_init(&(data_model.pc), &pi);

	    /*print_matrix(par.y_noise_cov);
	    DYNRPRINT(verbose_flag, "\n");*/

        double neg_log_like=EKimFilter(data_model.y, data_model.co_variate, data_model.y_time, &(data_model.pc), &pi, &par,
	    	eta_regime_j_t, error_cov_regime_j_t,
			eta_regime_jk_pred, error_cov_regime_jk_pred,
	    	pr_t, pr_t_given_t_minus_1,
			innov_v, inv_residual_cov, 
			eta_t, error_cov_t);

	    if ( status< 0) {
			MYPRINT("nlopt failed!\n");
	    }else{
			MYPRINT("Starting Hessian calculation ...\n");
		    data_model.pc.isnegloglikeweightedbyT=false;
			hessianRichardson(fittedpar, &data_model, function_neg_log_like, neg_log_like, Hessian_mat); /*information matrix*/
		    data_model.pc.isnegloglikeweightedbyT=weight_flag;
			MYPRINT("Finished Hessian calculation.\n");
			/* mathfunction_inv_matrix(Hessian_mat, inv_Hessian_mat); */ /*variance*/
		}


	    EKimSmoother(data_model.y_time, data_model.co_variate, &data_model.pc, &par, 
			pr_t_given_t_minus_1, pr_t, 
			eta_regime_jk_pred, error_cov_regime_jk_pred, 
			eta_regime_j_t, error_cov_regime_j_t,
			eta_smooth, error_cov_smooth, pr_T, transprob_T);

    /** =================Extended Kim Filter and Smoother: done======================**/

    /** =================Interface: SEXP Output====================== **/
	DYNRPRINT(verbose_flag, "Creating and allocating R output ... \n");
	SEXP res_list;
	SEXP res_names;
	/*TODO Lu Delete the outall_flag*/
	if (outall_flag){
	    res_list=PROTECT(allocVector(VECSXP,21));
	    res_names=PROTECT(allocVector(STRSXP, 21));
	}else if (debug_flag){
	    res_list=PROTECT(allocVector(VECSXP,9));
	    res_names=PROTECT(allocVector(STRSXP, 9));
	}else{
	    res_list=PROTECT(allocVector(VECSXP,10));
	    res_names=PROTECT(allocVector(STRSXP, 10));
	}


    SEXP exitflag=PROTECT(allocVector(INTSXP,1));
    *INTEGER(exitflag)=status;
    DYNRPRINT(verbose_flag, "exitflag created and copied.\n");
    SEXP negloglike=PROTECT(allocVector(REALSXP,1));
    *REAL(negloglike)=neg_log_like;/*should be the same as minf*/
    DYNRPRINT(verbose_flag, "negloglike created and copied.\n");
    SEXP fittedout=PROTECT(allocVector(REALSXP, data_model.pc.num_func_param));
	/*DYNRPRINT(verbose_flag, "fittedout created.\n");*/
	memcpy(REAL(fittedout),fittedpar,sizeof(fittedpar));
	/*DYNRPRINT(verbose_flag, "fittedout copied.\n");
	print_array(REAL(fittedout),data_model.pc.num_func_param);
	DYNRPRINT(verbose_flag, "\n");*/
	DYNRPRINT(verbose_flag, "fittedout created and copied.\n");

    SEXP hessian=PROTECT(allocMatrix(REALSXP, data_model.pc.num_func_param, data_model.pc.num_func_param));
         ptr_index=REAL(hessian);
         double tmp;
         for (index=0;index<data_model.pc.num_func_param;index++){
             ptr_index[index+data_model.pc.num_func_param*index]=gsl_matrix_get(Hessian_mat,index,index);
             for (index_col=index+1;index_col<data_model.pc.num_func_param;index_col++){
                 tmp=gsl_matrix_get(Hessian_mat,index,index_col);
                 ptr_index[index+data_model.pc.num_func_param*index_col]=tmp;
                 ptr_index[index_col+data_model.pc.num_func_param*index]=tmp;
             }
         }
    DYNRPRINT(verbose_flag, "hessian created and copied.\n");
    /*eta_regime_t: input and output of filter & input of smooth: eta^k_it|t*/
    /*error_cov_regime_t:input and output of filter & input of smooth: error_cov^k_it|t*/
    /*eta_regime_regime_t_pred: output of filter and input of smooth: eta^regime_jk_it|t-1*/
    /*error_cov_regime_regime_t_pred: output of filter and input of smooth: error_cov^regime_jk_it|t-1*/
    /*eta_regime_regime_t_plus_1: output of filter: eta^regime_jk_it|t*/
    /*error_cov_regime_regime_t_plus_1:output of filter: error_cov^regime_jk_it|t*/
    /*innov_vec:output of filter: innovation vector*/
    /*inverse_residual_cov:output of filter: inverse_residual_cov*/
    /*pr_t_given_t:output of filter and input of smooth: Pr(S_it=k|Y_it)*/
    /*pr_t_given_t_less_1:output of filter and input of smooth: Pr(S_it=k|Y_i,t-1)*/
    /*pr_t_given_T: output of smooth: Pr(S_it=k|Y_iT)*/
    /*transprob_given_T: output of smooth: Pr(S_i,t+1=h, S_it=k|Y_iT)*/
    /*eta_regime_smooth:output of smooth: eta^k_it|T*/
    /*error_cov_regime_smooth:output of smooth: error_cov^k_it|T*/
    /*eta_smooth_final: output of smooth: eta_it|T*/
    /*error_cov_smooth_final:output of smooth: error_cov_it|T*/
	     SEXP dims_eta_smooth_final=PROTECT(allocVector(INTSXP,2));
	     memcpy(INTEGER(dims_eta_smooth_final), ((int[]){data_model.pc.dim_latent_var, data_model.pc.total_obs}),2*sizeof(int));
	     SEXP eta_smooth_final = PROTECT(Rf_allocArray(REALSXP,dims_eta_smooth_final));
	     index=0;
	     ptr_index=REAL(eta_smooth_final);
	     for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
	         for(index_col=0; index_col<data_model.pc.dim_latent_var; index_col++){
	             ptr_index[index]=gsl_vector_get(eta_smooth[index_sbj_t],index_col);
	             index++;
	         }
	     }
	     DYNRPRINT(verbose_flag, "eta_smooth_final created and copied.\n");

	     SEXP dims_error_cov_smooth_final=PROTECT(allocVector(INTSXP,3));
	     memcpy(INTEGER(dims_error_cov_smooth_final), ((int[]){data_model.pc.dim_latent_var,  data_model.pc.dim_latent_var,  data_model.pc.total_obs}),3*sizeof(int));
	     SEXP error_cov_smooth_final = PROTECT(Rf_allocArray(REALSXP,dims_error_cov_smooth_final));
	     index=0;
	     ptr_index=REAL(error_cov_smooth_final);
	     for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
	         for(index_col=0; index_col<data_model.pc.dim_latent_var; index_col++){
	             for(index_row=0; index_row<data_model.pc.dim_latent_var; index_row++){
	                 ptr_index[index]=gsl_matrix_get(error_cov_smooth[index_sbj_t],index_row, index_col);
	                 index++;
	             }
	         }
	     }
	     DYNRPRINT(verbose_flag, "error_cov_smooth_final created and copied.\n");

		 SEXP dims_pr_t_given_T=PROTECT(allocVector(INTSXP,2));
		 memcpy(INTEGER(dims_pr_t_given_T), ((int[]){data_model.pc.num_regime,  data_model.pc.total_obs}),2*sizeof(int));
		 SEXP pr_t_given_T = PROTECT(Rf_allocArray(REALSXP,dims_pr_t_given_T));
		 index=0;
		 ptr_index=REAL(pr_t_given_T);
		 for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
			 for(index_col=0; index_col<data_model.pc.num_regime; index_col++){
				 ptr_index[index]=gsl_vector_get(pr_T[index_sbj_t],index_col);
				 index++;
			 }
		 }
		 DYNRPRINT(verbose_flag, "pr_t_given_T created and copied.\n");
		 
		 SEXP dims_pr_t_given_t=PROTECT(allocVector(INTSXP,2));
		 memcpy(INTEGER(dims_pr_t_given_t), ((int[]){data_model.pc.num_regime,  data_model.pc.total_obs}),2*sizeof(int));
		 SEXP pr_t_given_t = PROTECT(Rf_allocArray(REALSXP,dims_pr_t_given_t));
		 index=0;
		 ptr_index=REAL(pr_t_given_t);
		 for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
			 for(index_col=0; index_col<data_model.pc.num_regime; index_col++){
				 ptr_index[index]=gsl_vector_get(pr_t[index_sbj_t],index_col);
				 index++;
			 }
		 }
		 DYNRPRINT(verbose_flag, "pr_t_given_t created and copied.\n");

		 	 
 		/*filtered state estimate*/
 		SEXP dims_eta_filtered=PROTECT(allocVector(INTSXP, 2));
 		memcpy(INTEGER(dims_eta_filtered), ((int[]){data_model.pc.dim_latent_var, data_model.pc.total_obs}),2*sizeof(INTEGER(dims_eta_filtered)));
 		SEXP eta_filtered = PROTECT(Rf_allocArray(REALSXP, dims_eta_filtered));
 		index = 0;
 		ptr_index = REAL(eta_filtered);
 		for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
 				for(index_col=0; index_col<data_model.pc.dim_latent_var; index_col++){
 					ptr_index[index]=gsl_vector_get(eta_t[index_sbj_t],index_col);
 					index++;
 				}
 		}
  		DYNRPRINT(verbose_flag, "eta_filtered created and copied.\n");
		
 		/*filtered error covariance estimate*/
 		SEXP dims_error_cov_filtered=PROTECT(allocVector(INTSXP, 3));
 		memcpy(INTEGER(dims_error_cov_filtered), ((int[]){data_model.pc.dim_latent_var,  data_model.pc.dim_latent_var,  data_model.pc.total_obs}),3*sizeof(int));
 		SEXP error_cov_filtered = PROTECT(Rf_allocArray(REALSXP, dims_error_cov_filtered));
 		index = 0;
 		ptr_index=REAL(error_cov_filtered);
 		for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
 				for(index_col=0; index_col<data_model.pc.dim_latent_var; index_col++){
 					for(index_row=0; index_row<data_model.pc.dim_latent_var; index_row++){
 						ptr_index[index]=gsl_matrix_get(error_cov_t[index_sbj_t],index_row, index_col);
 						index++;
 					}
 				}
 		}
 	 	DYNRPRINT(verbose_flag, "error_cov_filtered created and copied.\n");

	     SET_STRING_ELT(res_names, 0, mkChar("exitflag"));
	     SET_VECTOR_ELT(res_list, 0, exitflag);
	     SET_STRING_ELT(res_names, 1, mkChar("neg.log.likelihood"));
	     SET_VECTOR_ELT(res_list, 1, negloglike);

	     SET_STRING_ELT(res_names, 2, mkChar("fitted.parameters"));
	     SET_VECTOR_ELT(res_list, 2, fittedout);
	     SET_STRING_ELT(res_names, 3, mkChar("hessian.matrix"));
	     SET_VECTOR_ELT(res_list, 3, hessian);
    
	 	 SET_STRING_ELT(res_names, 4, mkChar("eta_smooth_final"));
	     SET_VECTOR_ELT(res_list, 4, eta_smooth_final);
	     SET_STRING_ELT(res_names, 5, mkChar("error_cov_smooth_final"));
	     SET_VECTOR_ELT(res_list, 5, error_cov_smooth_final);
	     SET_STRING_ELT(res_names, 6, mkChar("pr_t_given_T"));
	     SET_VECTOR_ELT(res_list, 6, pr_t_given_T);

		 SET_STRING_ELT(res_names, 7, mkChar("eta_filtered"));
		 SET_VECTOR_ELT(res_list, 7, eta_filtered);
		 SET_STRING_ELT(res_names, 8, mkChar("error_cov_filtered"));
		 SET_VECTOR_ELT(res_list, 8, error_cov_filtered);
		 SET_STRING_ELT(res_names, 9, mkChar("pr_t_given_t"));
		 SET_VECTOR_ELT(res_list, 9, pr_t_given_t);
		 
		 
	if (outall_flag|debug_flag){
		
		/*innov vec*/
		SEXP dims_innov_vec=PROTECT(allocVector(INTSXP,4));
		memcpy(INTEGER(dims_innov_vec), ((int[]){data_model.pc.dim_latent_var,  data_model.pc.num_regime,  data_model.pc.num_regime,  data_model.pc.total_obs}),4*sizeof(int));
		SEXP innov_vec = PROTECT(Rf_allocArray(REALSXP,dims_innov_vec));
		index=0;
		ptr_index=REAL(innov_vec);
		for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
			for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++){
				for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
					for(index_col=0; index_col<data_model.pc.dim_obs_var; index_col++){
						ptr_index[index]=gsl_vector_get(innov_v[index_sbj_t][regime_j][regime_k],index_col);
						index++;
					}
				}
			}
		}
		SET_STRING_ELT(res_names, 10, mkChar("innov_vec"));
		SET_VECTOR_ELT(res_list, 10, innov_vec);
		UNPROTECT(2);
		DYNRPRINT(verbose_flag, "innov_vec created and copied.\n");
		
		/*TODO output the residual covariance*/
		/*inverse residual cov*/
		SEXP dims_inverse_residual_cov=PROTECT(allocVector(INTSXP,5));
		memcpy(INTEGER(dims_inverse_residual_cov), ((int[]){data_model.pc.dim_latent_var,  data_model.pc.dim_latent_var,  data_model.pc.num_regime,data_model.pc.num_regime,  data_model.pc.total_obs}),5*sizeof(int));
		SEXP inverse_residual_cov = PROTECT(Rf_allocArray(REALSXP,dims_inverse_residual_cov));
		index=0;
		ptr_index=REAL(inverse_residual_cov);
		for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
			for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++){
				for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
					for(index_col=0; index_col<data_model.pc.dim_obs_var; index_col++){
						for(index_row=0; index_row<data_model.pc.dim_obs_var; index_row++){
							ptr_index[index]=gsl_matrix_get(inv_residual_cov[index_sbj_t][regime_j][regime_k],index_row, index_col);
							index++;
						}
					}
				}
			}
		}
		SET_STRING_ELT(res_names, 11, mkChar("inverse_residual_cov"));
		SET_VECTOR_ELT(res_list, 11, inverse_residual_cov);
		UNPROTECT(2);
		DYNRPRINT(verbose_flag, "inverse_residual_cov created and copied.\n");
		
		/*TODO output predicted states and covariance estimates */
	}
	
	/*TODO clean this part*/
	if (outall_flag){
		
		SEXP invhessian=PROTECT(allocMatrix(REALSXP, data_model.pc.num_func_param, data_model.pc.num_func_param));
		ptr_index=REAL(invhessian);
		for (index=0;index<data_model.pc.num_func_param;index++){
			ptr_index[index+data_model.pc.num_func_param*index]=gsl_matrix_get(inv_Hessian_mat,index,index);
			for (index_col=index+1;index_col<data_model.pc.num_func_param;index_col++){
				tmp=gsl_matrix_get(inv_Hessian_mat,index,index_col);
				ptr_index[index+data_model.pc.num_func_param*index_col]=tmp;
				ptr_index[index_col+data_model.pc.num_func_param*index]=tmp;
			}
		}
		SET_STRING_ELT(res_names, 16, mkChar("inverse.hessian.matrix"));
		SET_VECTOR_ELT(res_list, 16, invhessian);
		UNPROTECT(1);
		DYNRPRINT(verbose_flag, "invhessian created and copied.\n");
				
		
		SEXP dims_pr_t_given_t_less_1=PROTECT(allocVector(INTSXP,2));
		memcpy(INTEGER(dims_pr_t_given_t_less_1), ((int[]){data_model.pc.num_regime,  data_model.pc.total_obs}),2*sizeof(int));
		SEXP pr_t_given_t_less_1 = PROTECT(Rf_allocArray(REALSXP,dims_pr_t_given_t_less_1));
		index=0;
		ptr_index=REAL(pr_t_given_t_less_1);
		for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
			for(index_col=0; index_col<data_model.pc.num_regime; index_col++){
				ptr_index[index]=gsl_vector_get(pr_t_given_t_minus_1[index_sbj_t],index_col);
				index++;
			}
		}
		SET_STRING_ELT(res_names, 17, mkChar("pr_t_given_t_less_1"));
		SET_VECTOR_ELT(res_list, 17, pr_t_given_t_less_1);
		UNPROTECT(2);
		DYNRPRINT(verbose_flag, "pr_t_given_t_less_1 created and copied.\n");
		
		
		SEXP dims_transprob_given_T=PROTECT(allocVector(INTSXP,3));
		memcpy(INTEGER(dims_transprob_given_T), ((int[]){data_model.pc.num_regime,  data_model.pc.num_regime,  data_model.pc.total_obs}),3*sizeof(int));
		SEXP transprob_given_T = PROTECT(Rf_allocArray(REALSXP,dims_transprob_given_T));
		index=0;
		ptr_index=REAL(transprob_given_T);
		for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
			for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
				for(index_col=0; index_col<data_model.pc.num_regime; index_col++){
					ptr_index[index]=gsl_vector_get(transprob_T[index_sbj_t][regime_j],index_col);
					index++;
				}
			}
		}
		SET_STRING_ELT(res_names, 18, mkChar("transprob_given_T"));
		SET_VECTOR_ELT(res_list, 18, transprob_given_T);
		UNPROTECT(2);
		DYNRPRINT(verbose_flag, "transprob_given_T created and copied.\n");	
	
	}
	
	
	setAttrib(res_list, R_NamesSymbol, res_names);
	
    DYNRPRINT(verbose_flag, "R return list completed.\n");
    /** =================Interface: Output done====================== **/

    /** =================Free Allocated space====================== **/
	DYNRPRINT(verbose_flag, "Freeing objects before return ... \n");
    if (data_model.pc.isContinuousTime){
		if (outall_flag){
			UNPROTECT(90-17-19-8 + 6);
		}else if (debug_flag){
			UNPROTECT(90-17-19-8 + 6);
		}else{
			UNPROTECT(90-17-19-8 + 6);
		}

	}else{
		if (outall_flag){
			UNPROTECT(90-17-19-8-1 + 6);
		}else if (debug_flag){
			UNPROTECT(90-17-19-8-1 + 6);
		}else{
			UNPROTECT(90-17-19-8-1 + 6);
		}
	}/*unprotect objects: find all PROTECT in the script, then -2*2Cancel-6ENDUNP-4outputflag-2/3CTflag -1 comment=-17*/

    free(str_number);
    free(data_model.pc.index_sbj);

    for(index=0; index<data_model.pc.total_obs; index++){
    gsl_vector_free(data_model.y[index]);}
    free(data_model.y);

    for(index=0; index<data_model.pc.total_obs; index++){
        gsl_vector_free(data_model.co_variate[index]);
    }
    free(data_model.co_variate);


    free(data_model.y_time);

    gsl_matrix_free(Hessian_mat);
    gsl_matrix_free(inv_Hessian_mat);

    gsl_vector_free(pi.pr_0);

    for(index=0;index<data_model.pc.num_regime;index++){
        gsl_vector_free((pi.eta_0)[index]);
    }
    free(pi.eta_0);

    for(index=0;index<data_model.pc.num_regime;index++){
        gsl_matrix_free((pi.error_cov_0)[index]);
    }
    free(pi.error_cov_0);

    gsl_matrix_free(par.regime_switch_mat);

    gsl_matrix_free(par.eta_noise_cov);

    gsl_matrix_free(par.y_noise_cov);

    free(par.func_param);


    /*input and output of filter & input of smooth: eta^k_it|t*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            gsl_vector_free(eta_regime_j_t[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(eta_regime_j_t[index_sbj_t]);
    }
    free(eta_regime_j_t);

    /*input and output of filter & input of smooth: error_cov^k_it|t*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            gsl_matrix_free(error_cov_regime_j_t[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(error_cov_regime_j_t[index_sbj_t]);
    }
    free(error_cov_regime_j_t);

    /*output of filter and input of smooth: eta^regime_jk_it|t-1*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
            gsl_vector_free(eta_regime_jk_pred[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            free(eta_regime_jk_pred[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(eta_regime_jk_pred[index_sbj_t]);
    }
    free(eta_regime_jk_pred);

    /*output of filter and input of smooth: error_cov^regime_jk_it|t-1*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
                gsl_matrix_free(error_cov_regime_jk_pred[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            free(error_cov_regime_jk_pred[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(error_cov_regime_jk_pred[index_sbj_t]);
    }
    free(error_cov_regime_jk_pred);

    /*output of filter and input of smooth: Pr(S_it=k|Y_it)*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        gsl_vector_free(pr_t[index_sbj_t]);
    }
    free(pr_t);

    /*output of filter and input of smooth: Pr(S_it=k|Y_i,t-1)*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        gsl_vector_free(pr_t_given_t_minus_1[index_sbj_t]);
    }
    free(pr_t_given_t_minus_1);

    /*output of smooth: eta_it|T*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        gsl_vector_free(eta_smooth[index_sbj_t]);
    }
    free(eta_smooth);

    /*output of smooth: error_cov_it|T*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        gsl_matrix_free(error_cov_smooth[index_sbj_t]);
    }
    free(error_cov_smooth);

    /*output of smooth: Pr(S_it=k|Y_iT)*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        gsl_vector_free(pr_T[index_sbj_t]);
    }
    free(pr_T);

    /*output of smooth: Pr(S_i,t+1=h, S_it=k|Y_iT)*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            gsl_vector_free(transprob_T[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(transprob_T[index_sbj_t]);
    }
    free(transprob_T);
	
    /*output of filter: eta_it|t*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		gsl_vector_free(eta_t[index_sbj_t]);
    }
	free(eta_t);
	/*output of filter: error_cov_it|t*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		gsl_matrix_free(error_cov_t[index_sbj_t]);
	}
	free(error_cov_t);
	
    /*output of filter: innovation vector*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
            gsl_vector_free(innov_v[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            free(innov_v[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(innov_v[index_sbj_t]);
    }
    free(innov_v);

    /*output of filter: inverse_residual_cov*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
                gsl_matrix_free(inv_residual_cov[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            free(inv_residual_cov[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(inv_residual_cov[index_sbj_t]);
    }
    free(inv_residual_cov);


    return res_list;
}




