#### No-covariates
# Single regime, no covariates
initial <- prep.initial(
	values.inistate=c(0, 1),
	params.inistate=c('inipos', 'fixed'), #initial position is free parameter 5, initial slope is fixed at 1
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2)) #initial covariance is fixed to a diagonal matrix of 1s.

#/******************C Code generated********************/#
#void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0, size_t *index_sbj){
#	
#	gsl_vector *Pvector = gsl_vector_calloc(1);
#	gsl_vector_set(Pvector, 0, 1);
#	mathfunction_softmax(Pvector, pr_0);
#	
#	size_t num_regime=pr_0->size;
#	size_t dim_latent_var=error_cov_0[0]->size1;
#	size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
#	size_t i;
#	size_t regime;
#	for(regime=0; regime < num_regime; regime++){
#		for(i=0; i < num_sbj; i++){
#			gsl_vector_set((eta_0)[regime], i*dim_latent_var+0, param[4]);
#			gsl_vector_set((eta_0)[regime], i*dim_latent_var+1, 1);
#		}
#	}
#	gsl_vector_free(Pvector);
#}


#### One covariate
# Single regime, one covariate on the inital mean
initial <- prep.initial(
	values.inistate=matrix(
		c(0, .5,
		  1,  0), byrow=TRUE,
		nrow=2, ncol=2), #nrow = numLatentState, ncol = numCovariates + 1
	params.inistate=matrix(
		c('iniPosInt', 'iniPosSlopeU1',
		'fixed', 'fixed'), byrow=TRUE,
		nrow=2, ncol=2), #initial position has free intercept and free u1 effect, initial slope is fixed at 1
	values.inicov=diag(1, 2),
	params.inicov=diag('fixed', 2), #initial covariance is fixed to a diagonal matrix of 1s.
	covariates='u1')

#eta_0 = [ iniPosInt + iniPosSlopeU1*u1;
#		  1         + 0*u1]

# TODO: update C code to create eta_0 as above
#/******************C Code generated********************/#
#void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0, size_t *index_sbj){
#	
#	gsl_vector *Pvector = gsl_vector_calloc(1);
#	gsl_vector_set(Pvector, 0, 1);
#	mathfunction_softmax(Pvector, pr_0);
#	
#	size_t num_regime=pr_0->size;
#	size_t dim_latent_var=error_cov_0[0]->size1;
#	size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
#	size_t i;
#	size_t regime;
#	for(regime=0; regime < num_regime; regime++){
#		for(i=0; i < num_sbj; i++){
#			gsl_vector_set((eta_0)[regime], i*dim_latent_var+0, param[4]);
#			gsl_vector_set((eta_0)[regime], i*dim_latent_var+1, 1);
#		}
#	}
#	gsl_vector_free(Pvector);
#}
