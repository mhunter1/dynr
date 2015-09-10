/*
Author: Lu Ou, Mike Hunter
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
#include <stdio.h>
#include <string.h>
#include <nlopt.h>
#include "math_function.h"
#include "cdaekf.h"
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


/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    size_t i;
    for (i = 0; i < length(list); i++)
	if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	   elmt = VECTOR_ELT(list, i);
	   break;
	}
    return elmt;
}
/**
 * The getway function for the R interface
 * @param model_list is a list in R of all model specifications.
 * @param paramvec is a vector in R of the parameter starting values
 * @param ubvec is a vector in R of the upper bounds of search region
 * @param lbvec is a vecotr in R of the lower bounds of the search region
 */
SEXP main_R(SEXP model_list,SEXP data_list)
{
    size_t index,index_col;

    /** =======================Interface : Start to Set up the data and the model========================= **/    
     
    static Data_and_Model data_model;
     
    /* From the SEXP called model_list, get the list element named "num_sbj" */
    data_model.pc.num_sbj=(size_t) *REAL(getListElement(model_list, "num_sbj"));/*number of subjects*/
    printf("num_sbj: %lu\n",data_model.pc.num_sbj);
    
    data_model.pc.num_func_param=(size_t) *REAL(getListElement(model_list, "num_func_param")); /*number of function parameters*/
    printf("num_func_param: %lu\n",data_model.pc.num_func_param);
    
    data_model.pc.dim_latent_var=(size_t) *REAL(getListElement(model_list, "dim_latent_var"));/*number of latent variables*/
    printf("dim_latent_var: %lu\n",data_model.pc.dim_latent_var);
    
    data_model.pc.dim_obs_var=(size_t) *REAL(getListElement(model_list, "dim_obs_var")); /*number of observed variables*/
    printf("dim_obs_var: %lu\n",data_model.pc.dim_obs_var);

    data_model.pc.dim_co_variate=(size_t) *REAL(getListElement(model_list, "dim_co_variate"));/*number of covariates*/   
    printf("dim_co_variate: %lu\n",data_model.pc.dim_co_variate);
    
    data_model.pc.num_regime=(size_t) *REAL(getListElement(model_list, "num_regime")); /*number of regimes*/
    printf("num_regime: %lu\n",data_model.pc.num_regime);

    


    /*function specifications*/
    data_model.pc.func_measure=function_measurement;
    data_model.pc.func_dF_dx=function_dF_dx;
    data_model.pc.func_jacobdynamic=function_jacobdynamic;
    data_model.pc.func_dx_dt=function_dx_dt;
    data_model.pc.func_dP_dt=function_dP_dt;
    data_model.pc.func_initial_condition=function_initial_condition;
    data_model.pc.func_regime_switch=function_regime_switch;
    data_model.pc.func_noise_cov=function_noise_cov;
    data_model.pc.isnegloglikeweightedbyT=false;
    data_model.pc.second_order=false;/*true;*/
    data_model.pc.adaodesolver=false;/*true: use adapative ode solver; false: RK4*/
    if (data_model.pc.adaodesolver){
        data_model.pc.func_dynam=function_dynam_ada;
    }else{
        data_model.pc.func_dynam=rk4_odesolver;
    }
    

    /*specify the start position for each subject: User always need to provide a txt file called tStart.txt*/
    /*for example, 500 time points for each sbj, specify 0 500 1000 ... 10000 also the end point*/
    /*n subjects -> n+1 indices*/   
    data_model.pc.index_sbj=(size_t *)malloc((data_model.pc.num_sbj+1)*sizeof(size_t *));  
    
    /*printf("tstart 0: %lu\n",(size_t) *REAL(getListElement(data_list, "tstart")));
    printf("tstart 1: %lu\n",(size_t) REAL(getListElement(data_list, "tstart"))[1]);
    printf("tstart 2: %lu\n",(size_t) REAL(getListElement(data_list, "tstart"))[2]);*/
    
    double *ptr_index;/*used for multiple times*/
    
    ptr_index=REAL(getListElement(data_list, "tstart"));
    for(index=0;index<=data_model.pc.num_sbj;index++){    
        data_model.pc.index_sbj[index]=(size_t) ptr_index[index];
    }
    printf("index_sbj 2: %lu\n",data_model.pc.index_sbj[1]);
    
    data_model.pc.total_obs=*(data_model.pc.index_sbj+data_model.pc.num_sbj);/*total observations for all subjects*/
    printf("total_obs: %lu\n",data_model.pc.total_obs);
    
    /** read in the data**/
    
    data_model.y=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
    size_t t;
    for(t=0; t<data_model.pc.total_obs; t++){
        data_model.y[t]=gsl_vector_calloc(data_model.pc.dim_obs_var);/*y[t] corresponds to y(),which is a gsl_vector; loop through total_obj*/

    }
    
    size_t enough_length=(ceil(log10((double)data_model.pc.dim_obs_var>data_model.pc.dim_co_variate?data_model.pc.dim_obs_var:data_model.pc.dim_co_variate))+1)*sizeof(char);
    char *str_number=(char *)malloc(enough_length);
    char str_name[enough_length+7];
     
    for(index=0;index<data_model.pc.dim_obs_var;index++){
        sprintf(str_number, "%lu", index+1);
        sprintf(str_name, "%s", "obs");
        /*printf("The str_number is %s\n",str_number);
        printf("The str_name length is %lu\n",strlen(str_name));*/      
        ptr_index=REAL(getListElement(getListElement(data_list,"observed"), strncat(str_name, str_number, strlen(str_number))));
        
        for(t=0; t<data_model.pc.total_obs; t++){   
            gsl_vector_set(data_model.y[t],index, ptr_index[t]);
        }
    }

    data_model.co_variate=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
    for(t=0; t<data_model.pc.total_obs; t++){
        data_model.co_variate[t]=gsl_vector_calloc(data_model.pc.dim_co_variate);
    }
    
    for(index=0;index<data_model.pc.dim_co_variate;index++){
        sprintf(str_number, "%lu", index+1);
        sprintf(str_name, "%s", "covar");
        /*printf("The str_number is %s\n",str_number);
        printf("The str_name length is %lu\n",strlen(str_name));*/
        ptr_index=REAL(getListElement(getListElement(data_list,"covariates"), strncat(str_name, str_number, strlen(str_number))));
        
        for(t=0; t<data_model.pc.total_obs; t++){   
            gsl_vector_set(data_model.co_variate[t],index, ptr_index[t]);
        }
    }

    
    data_model.y_time=(double *)malloc(data_model.pc.total_obs*sizeof(double));
        memcpy(data_model.y_time,REAL(getListElement(data_list, "time")),data_model.pc.total_obs*sizeof(double));


    printf("In main_R:\n");
    print_vector(data_model.y[0]);
    printf("\n"); 
    print_vector(data_model.co_variate[0]);
    printf("\n");
    printf("y_time_1 is %lf\n",data_model.y_time[0]);
    
    
    /** Optimization options **/
    
    double params[data_model.pc.num_func_param];
    /*printf("xstart 0: %lu\n",(size_t) *REAL(getListElement(model_list, "xstart")));*/
    	memcpy(params,REAL(getListElement(model_list, "xstart")),sizeof(params));  
    /*printf("Array paramvec allocated.\n");*/
    /*print_array(params,data_model.pc.num_func_param);*/
    /*printf("\n");*/
    double fittedpar[data_model.pc.num_func_param];
    	memcpy(fittedpar, params,sizeof(params));
    printf("Array paramvec copied.\n");
    print_array(fittedpar,data_model.pc.num_func_param);
    printf("\n");
    double ub[data_model.pc.num_func_param];
    double lb[data_model.pc.num_func_param];
    	memcpy(ub,REAL(getListElement(model_list, "ub")),sizeof(ub));
    	memcpy(lb,REAL(getListElement(model_list, "lb")),sizeof(lb));
    /*printf("Arrays allocated.\n");*/
    /*double ub[6] = {4, 4, 4, 4, 4, 4};
    double lb[6] = {-4,-4,-4,-4,-12, -12}; */
    /*double params[]={log(1),log(2),0,0,-10,-10};*//* some initial guess*/
    /*log(1.2)=0.1823216,log(1.8)=0.5877867,-0.5,-0.5,log(0.0001)=-9.21034,log(0.0001)=-9.21034*//* lower bounds */

    
    /** =======================Interface: Model and data set up========================= **/
     
    /** =================Optimization: start======================**/	
    
    double minf; /* the minimum objective value, upon return */
	
    gsl_matrix *Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    gsl_matrix *inv_Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
	
	if (opt_nlopt(&data_model,data_model.pc.num_func_param,ub,lb,&minf,fittedpar,Hessian_mat,inv_Hessian_mat,1e-7) < 0) {
		/*printf("nlopt failed!\n");*/
	}
	else {
		/*printf("found minimum at \n");
		print_array(fittedpar,data_model.pc.num_func_param);
		printf("\n f = %0.10g\n", minf);*/
		printf("The hessian matrix is \n");
		print_matrix(Hessian_mat);
		printf("\n");
		printf("The inverse hessian matrix is \n");
		print_matrix(inv_Hessian_mat);
		printf("\n");
	}
	
	
	/*printf("Optimization done.\n");*/
    /** =================Optimization: done======================**/
	
    /** =================Interface: SEXP Output====================== **/
    SEXP res_list=PROTECT(allocVector(VECSXP,2));
    SEXP res_names=PROTECT(allocVector(STRSXP, 2));
    
    SEXP fittedout=PROTECT(allocVector(REALSXP,data_model.pc.num_func_param));
	/*printf("fittedout created.\n");*/
	memcpy(REAL(fittedout),fittedpar,sizeof(fittedpar));
	/*printf("fittedout copied.\n");
	print_array(REAL(fittedout),data_model.pc.num_func_param);
	printf("\n");*/
    
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
    SET_STRING_ELT(res_names, 0, mkChar("fitted.parameters")); 
    SET_VECTOR_ELT(res_list, 0, fittedout);
    SET_STRING_ELT(res_names, 1, mkChar("hessian.matrix")); 
    SET_VECTOR_ELT(res_list, 1, hessian);
    
    setAttrib(res_list, R_NamesSymbol, res_names);

	
	/*printf("Done.\n");*/
    /** =================Interface: Output done====================== **/
    
    /** =================Free Allocated space====================== **/
    UNPROTECT(4);/*unprotect 3 objects*/
    
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
    
    return res_list;
}




