/*=================================================================
 * Author: Lu Ou
 * Date: 7/30/2015
 * wrappernegloglike.C
 * This function is a wrapper function to get the negative log-likelihood, which is to be optimized by the optimizer.
 *
 *=================================================================*/


#include "wrappernegloglike.h"


double function_neg_log_like(const double *params, void *data){
	double neg_log_like;
	size_t index;
	
	/** model configuration **/
	Data_and_Model data_model=*((Data_and_Model *)data);/*dereference the void pointer*/
	
	/*MYPRINT("In function_neg_log_like:\n");
	print_vector(data_model.y[0]);
	MYPRINT("\n"); 
	print_vector(data_model.y[1]);
	MYPRINT("\n");
	print_vector(data_model.co_variate[0]);
	MYPRINT("\n");
	print_vector(data_model.co_variate[1]);
	MYPRINT("\n");
	MYPRINT("y_time_1 is %lf\n",data_model.y_time[0]);
	MYPRINT("y_time_2 is %lf\n",data_model.y_time[1]);
	exit(0);*/
	
	/** initialize regime parameter**/
	ParamInit pi;/*, fin_pi;*/
	
	// Allocate initial means
	pi.eta_0 = (gsl_vector **)malloc(data_model.pc.num_regime*sizeof(gsl_vector *));
	for(index=0; index < data_model.pc.num_regime; index++){
		(pi.eta_0)[index] = gsl_vector_calloc(data_model.pc.num_sbj*data_model.pc.dim_latent_var);
	}
	
	// Allocate initial covariance
	pi.error_cov_0 = (gsl_matrix **)malloc(data_model.pc.num_regime*sizeof(gsl_matrix *));
	for(index=0; index < data_model.pc.num_regime; index++){
		(pi.error_cov_0)[index] = gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
	}
	
	// Allocate initial probability
	pi.pr_0 = (gsl_vector **)malloc(data_model.pc.num_sbj*sizeof(gsl_vector *));
	for(index=0; index < data_model.pc.num_sbj; index++){
		(pi.pr_0)[index] = gsl_vector_calloc(data_model.pc.num_regime);
	}
	
	
	/** set parameter **/
	Param par;
	
	/*function parameters*/
	par.func_param = (double *)malloc(data_model.pc.num_func_param*sizeof(double));
	
	size_t i;
	for(i=0; i < data_model.pc.num_func_param; i++){
		par.func_param[i] = params[i];
	}
	print_array(par.func_param, data_model.pc.num_func_param);
	
	// Set initial conditions
	data_model.pc.func_initial_condition(par.func_param, data_model.co_variate, pi.pr_0, pi.eta_0, pi.error_cov_0, data_model.pc.index_sbj);
	
	/* Allocate noise covariances and regime switching matrix*/
	par.eta_noise_cov = gsl_matrix_calloc(data_model.pc.dim_latent_var, data_model.pc.dim_latent_var);
	par.y_noise_cov = gsl_matrix_calloc(data_model.pc.dim_obs_var, data_model.pc.dim_obs_var);
	par.regime_switch_mat = gsl_matrix_calloc(data_model.pc.num_regime, data_model.pc.num_regime);
	
	
	/** calculate the log_like **/
	data_model.pc.func_transform(par.func_param);
	model_constraint_init(&data_model.pc, &pi);
	
	neg_log_like = brekfis(data_model.y, data_model.co_variate, data_model.pc.total_obs, data_model.y_time, &data_model.pc, &pi, &par);
	MYPRINT("%lf\n", neg_log_like);
	
	
	/** Free allocated space **/
	for(index=0; index < data_model.pc.num_sbj; index++){
		gsl_vector_free((pi.pr_0)[index]);
	}
	free(pi.pr_0);
	
	for(index=0; index<data_model.pc.num_regime; index++){
		gsl_vector_free((pi.eta_0)[index]);
	}
	free(pi.eta_0);
	
	for(index=0; index < data_model.pc.num_regime; index++){
		gsl_matrix_free((pi.error_cov_0)[index]);
	}
	free(pi.error_cov_0);
	
	gsl_matrix_free(par.regime_switch_mat);
	gsl_matrix_free(par.eta_noise_cov);
	gsl_matrix_free(par.y_noise_cov);
	free(par.func_param);
	
	return neg_log_like;
}



