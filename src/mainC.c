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
#include "nlopt.h"
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
           /*errorcheck=fscanf(file_tstart,"%lu", (long unsigned int *) (data_model.pc.index_sbj+i));*/
            errorcheck=fscanf(file_tstart,"%lu", data_model.pc.index_sbj+i);
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
    
    double minf=0; /* the minimum objective value, upon return */
	
    gsl_matrix *Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    gsl_matrix *inv_Hessian_mat=gsl_matrix_calloc(data_model.pc.num_func_param,data_model.pc.num_func_param);
    /*double ub[] = {5, 5, 5, 5, 5, 5};
    double lb[] = {-5,-5,-5,-5,-15, -15};*/ /* lower bounds */
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
		    eta_regime_j_t[index_sbj_t][regime_j]=gsl_vector_alloc(data_model.pc.dim_latent_var);
		}
	    }
	    /*input and output of filter & input of smooth: error_cov^k_it|t*/
	    gsl_matrix ***error_cov_regime_j_t=(gsl_matrix ***)malloc(data_model.pc.total_obs*sizeof(gsl_matrix **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_regime_j_t[index_sbj_t]=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    error_cov_regime_j_t[index_sbj_t][regime_j]=gsl_matrix_alloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
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
		    eta_regime_jk_pred[index_sbj_t][regime_j][regime_k]=gsl_vector_alloc(data_model.pc.dim_latent_var);
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
			error_cov_regime_jk_pred[index_sbj_t][regime_j][regime_k]=gsl_matrix_alloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
		}
	    }
	    /*output of filter: eta^regime_jk_it|t*/
	    gsl_vector ****eta_regime_jk_t_plus_1=(gsl_vector ****)malloc(data_model.pc.total_obs*sizeof(gsl_vector ***));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		eta_regime_jk_t_plus_1[index_sbj_t]=(gsl_vector ***)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    eta_regime_jk_t_plus_1[index_sbj_t][regime_j]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
		}
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
		    eta_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]=gsl_vector_alloc(data_model.pc.dim_latent_var);
		}
	    }
	    /*output of filter: error_cov^regime_jk_it|t*/
	    gsl_matrix ****error_cov_regime_jk_t_plus_1=(gsl_matrix ****)malloc(data_model.pc.total_obs*sizeof(gsl_matrix ***));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_regime_jk_t_plus_1[index_sbj_t]=(gsl_matrix ***)malloc(data_model.pc.num_regime*sizeof(gsl_matrix **));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j]=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
		}
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
			error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]=gsl_matrix_alloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
		}
	    }
	    /*output of filter and input of smooth: Pr(S_it=k|Y_it)*/
	    gsl_vector **pr_t=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		pr_t[index_sbj_t]=gsl_vector_alloc(data_model.pc.num_regime);
	    }
	    /*output of filter and input of smooth: Pr(S_it=k|Y_i,t-1)*/
	    gsl_vector **pr_t_given_t_minus_1=(gsl_vector **)malloc(data_model.pc.total_obs*sizeof(gsl_vector *));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		pr_t_given_t_minus_1[index_sbj_t]=gsl_vector_calloc(data_model.pc.num_regime);
	    }
	
	    /*output of smooth: eta^k_it|T*/
	    gsl_vector ***eta_regime_j_smooth=(gsl_vector ***)malloc(data_model.pc.total_obs*sizeof(gsl_vector **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		eta_regime_j_smooth[index_sbj_t]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    eta_regime_j_smooth[index_sbj_t][regime_j]=gsl_vector_calloc(data_model.pc.dim_latent_var);
		}
	    }
	    /*output of smooth: error_cov^k_it|T*/
	    gsl_matrix ***error_cov_regime_j_smooth=(gsl_matrix ***)malloc(data_model.pc.total_obs*sizeof(gsl_matrix **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		error_cov_regime_j_smooth[index_sbj_t]=(gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    error_cov_regime_j_smooth[index_sbj_t][regime_j]=gsl_matrix_alloc(data_model.pc.dim_latent_var,data_model.pc.dim_latent_var);
		}
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
		pr_T[index_sbj_t]=gsl_vector_alloc(data_model.pc.num_regime);
	    }
	    /*output of smooth: Pr(S_i,t+1=h, S_it=k|Y_iT)*/
	    gsl_vector ***transprob_T=(gsl_vector ***)malloc(data_model.pc.total_obs*sizeof(gsl_vector **));
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		transprob_T[index_sbj_t]=(gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	    }
	    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
		for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
		    transprob_T[index_sbj_t][regime_j]=gsl_vector_alloc(data_model.pc.num_regime);
		}
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
		    innov_v[index_sbj_t][regime_j][regime_k]=gsl_vector_alloc(data_model.pc.dim_obs_var);
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
			inv_residual_cov[index_sbj_t][regime_j][regime_k]=gsl_matrix_alloc(data_model.pc.dim_obs_var, data_model.pc.dim_obs_var);
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
		/*Note that i was previously declared as size_t*/
	    for(i=0;i<data_model.pc.num_func_param;i++){
	    	    par.func_param[i]=fittedpar[i];
	    }


	    function_initial_condition(par.func_param, data_model.co_variate, pi.pr_0, pi.eta_0, pi.error_cov_0);

	    par.eta_noise_cov=gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
	    par.y_noise_cov=gsl_matrix_calloc(data_model.pc.dim_obs_var, data_model.pc.dim_obs_var);
	    par.regime_switch_mat=gsl_matrix_alloc(data_model.pc.num_regime, data_model.pc.num_regime);
 
	    /* calculate the log_like */
	    function_transform(&(data_model.pc), &pi, &par);

	    model_constraint_init(&(data_model.pc), &pi);

	    /*print_matrix(par.y_noise_cov);
	    printf("\n");*/

            double neg_log_like;

	    neg_log_like=EKimFilter(data_model.y, data_model.co_variate, data_model.y_time, &(data_model.pc), &pi, &par,
	    eta_regime_j_t,error_cov_regime_j_t,eta_regime_jk_pred,error_cov_regime_jk_pred,eta_regime_jk_t_plus_1,error_cov_regime_jk_t_plus_1,
	    pr_t, pr_t_given_t_minus_1,innov_v,inv_residual_cov);
	    
            
	    EKimSmoother(data_model.y_time, data_model.co_variate, &data_model.pc, &par, pr_t_given_t_minus_1, pr_t, eta_regime_jk_pred,error_cov_regime_jk_pred,eta_regime_j_t,error_cov_regime_j_t,
	    	    eta_regime_j_smooth,error_cov_regime_j_smooth,eta_smooth,error_cov_smooth,pr_T,transprob_T);

    /** =================Extended Kim Filter and Smoother: done======================**/
	
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


    /*output of filter: eta^regime_jk_it|t*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
            gsl_vector_free(eta_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            free(eta_regime_jk_t_plus_1[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(eta_regime_jk_t_plus_1[index_sbj_t]);
    }
    free(eta_regime_jk_t_plus_1);

    /*output of filter: error_cov^regime_jk_it|t*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            for(regime_k=0; regime_k<data_model.pc.num_regime; regime_k++)
                gsl_matrix_free(error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j][regime_k]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            free(error_cov_regime_jk_t_plus_1[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(error_cov_regime_jk_t_plus_1[index_sbj_t]);
    }
    free(error_cov_regime_jk_t_plus_1);

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

    /*output of smooth: eta^k_it|T*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            gsl_vector_free(eta_regime_j_smooth[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(eta_regime_j_smooth[index_sbj_t]);
    }
    free(eta_regime_j_smooth);


    /*output of smooth: error_cov^k_it|T*/
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        for(regime_j=0; regime_j<data_model.pc.num_regime; regime_j++){
            gsl_matrix_free(error_cov_regime_j_smooth[index_sbj_t][regime_j]);
        }
    }
    for(index_sbj_t=0;index_sbj_t<data_model.pc.total_obs;index_sbj_t++){
        free(error_cov_regime_j_smooth[index_sbj_t]);
    }
    free(error_cov_regime_j_smooth);


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
    

	
	return 0;
}


