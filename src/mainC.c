/*
Author: Lu Ou
Date: 2015-09-10
Filename: main.c
Purpose: the main function called in C to obtain parameters estimates by minimizing the negative log likelihood function;
        good for debugging purposes.
Run with
   ./runnlopt_onair.sh
   type in the command: gsl-config --cflags to find the compiler flag
   type in the command: gsl-config --libs to find the flag
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


int main()
{    
    size_t index;

    /** =======================Interface : Start to Set up the data and the model========================= **/    
     
    static Data_and_Model data_model;
    data_model.pc.num_sbj=217;/*number of subjects*/

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


    data_model.pc.dim_latent_var=4;/*number of latent variables*/
    data_model.pc.dim_obs_var=2;/*number of observed variables*/
    data_model.pc.dim_co_variate=1; /*number of covariates*/
    data_model.pc.num_func_param=6; /*number of function parameters*/
    data_model.pc.num_regime=1;/*number of regimes*/


    data_model.pc.index_sbj=(size_t *)malloc((data_model.pc.num_sbj+1)*sizeof(size_t *));
    size_t i;

    /*specify the start position for each subject: User always need to provide a txt file called tStart.txt*/
    /*for example, 500 time points for each sbj, specify 0 500 1000 ... 10000 also the end point*/
    /*n subjects -> n+1 indices*/   
    FILE *file_tstart=fopen("../data/tStartPANAsim.txt","r");
    if (file_tstart == NULL) {
        perror("fopen");
        printf("-1");
    }

    int errorcheck;
    for(i=0;i<=data_model.pc.num_sbj;i++){
           errorcheck=fscanf(file_tstart,"%lu",data_model.pc.index_sbj+i);
            if (errorcheck == EOF) {
                if (ferror(file_tstart)) {
                    perror("fscanf");
                }
                else {
                    fprintf(stderr, "Error: fscanf reached end of file, no matching characters, no matching failure\n");
                }
                printf("-1");
            }
            else if (errorcheck != 1) {
                fprintf(stderr, "Error: fscanf successfully matched and assigned %i input items\n", errorcheck);
                printf("-1");
            }
        }
    if (fclose(file_tstart) == EOF) {
        perror("fclose");
        printf("-1");
    }
    
    data_model.pc.total_obs=*(data_model.pc.index_sbj+data_model.pc.num_sbj);/*total observations for all subjects*/
 
    
    /** read in the data**/
    data_model.y=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
    size_t t;
    for(t=0; t<data_model.pc.total_obs; t++){
        data_model.y[t]=gsl_vector_calloc(data_model.pc.dim_obs_var);/*y[t] corresponds to y(),which is a gsl_vector; loop through total_obj*/

    }
    
    data_model.co_variate=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
    for(t=0; t<data_model.pc.total_obs; t++){
        data_model.co_variate[t]=gsl_vector_calloc(data_model.pc.dim_co_variate);
    }

    data_model.y_time=(double *)malloc(data_model.pc.total_obs*sizeof(double));

    FILE *file_data=fopen("../data/dataPANAsim.txt","r");
    if (file_data == NULL) {
        perror("fopen");
    }
    
    for(t=0;t<data_model.pc.total_obs;t++){
        errorcheck=0;
        errorcheck+=fscanf(file_data,"%*lu %lf",data_model.y_time+t);
        for(index=0;index<data_model.pc.dim_obs_var;index++){errorcheck+=fscanf(file_data," %lf",data_model.y[t]->data+index);}
        for(index=0;index<data_model.pc.dim_co_variate;index++){errorcheck+=fscanf(file_data," %lf",data_model.co_variate[t]->data+index);}
        
        if (errorcheck == EOF) {
            if (ferror(file_data)) {
                perror("fscanf");
            }
            else {
                fprintf(stderr, "Error: fscanf reached end of file, no matching characters, no matching failure\n");
            }
        }
        else if (errorcheck != (data_model.pc.dim_obs_var+data_model.pc.dim_co_variate+1)) {
            fprintf(stderr, "Error: fscanf successfully matched and assigned %i input items\n", errorcheck);
        }
        
    }
    if (fclose(file_data) == EOF) {
        perror("fclose");
    }


    printf("In main_R:\n");
    print_vector(data_model.y[0]);
    printf("\n"); 
    print_vector(data_model.co_variate[0]);
    printf("\n");
    printf("y_time_1 is %lf\n",data_model.y_time[0]);
    
    
    /** Optimization options **/
    double ub[] = {5, 5, 5, 5, 5, 5};
    double lb[] = {-5,-5,-5,-5,-15, -15}; /* lower bounds */
    double params[]={log(1),log(2),0,0,-10,-10};/* some initial guess*/
    /*log(1.2)=0.1823216,log(1.8)=0.5877867,-0.5,-0.5,log(0.0001)=-9.21034,log(0.0001)=-9.21034*/
    double fittedpar[data_model.pc.num_func_param];
    	memcpy(fittedpar, params,sizeof(params));
    printf("Array paramvec copied.\n");
    print_array(fittedpar,data_model.pc.num_func_param);
    printf("\n");
    /*printf("Arrays allocated.\n");*/
    
    /** =======================Interface: Model and data set up========================= **/
     
    /** =================Optimization: start======================**/	
    
    double minf; /* the minimum objective value, upon return */
	
    gsl_matrix *Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    gsl_matrix *inv_Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    /*int status=opt_nlopt(&data_model,data_model.pc.num_func_param,ub,lb,&minf,fittedpar,Hessian_mat,inv_Hessian_mat,1e-7);*/
    int status=-1;
    
	if ( status< 0) {
		printf("nlopt failed!\n");
	}
	else {
		printf("found minimum at \n");
		print_array(fittedpar,data_model.pc.num_func_param);
		printf("\n f = %0.10g\n", minf);
		printf("The hessian matrix is \n");
		print_matrix(Hessian_mat);
		printf("\n");
		printf("The inverse hessian matrix is \n");
		print_matrix(inv_Hessian_mat);
		printf("\n");
	}
	
	
	/*printf("Optimization done.\n");*/
    /** =================Optimization: done======================**/
	
    /** =================Free Allocated space====================== **/

    free(data_model.pc.index_sbj);
    data_model.pc.index_sbj=NULL;
    
    for(index=0; index<data_model.pc.total_obs; index++){
        gsl_vector_free(data_model.y[index]);
        data_model.y[index]=NULL;
    }
    free(data_model.y);
    data_model.y=NULL;

    for(index=0; index<data_model.pc.total_obs; index++){
        gsl_vector_free(data_model.co_variate[index]);
        data_model.co_variate[index]=NULL;
    }
    free(data_model.co_variate);
    data_model.co_variate=NULL;

    free(data_model.y_time);
    data_model.y_time=NULL;
    
    gsl_matrix_free(Hessian_mat);
    Hessian_mat=NULL;
    gsl_matrix_free(inv_Hessian_mat);
    inv_Hessian_mat=NULL;
    
	
	return 0;
}


