#ifndef BREKFIS_H_INCLUDED
#define BREKFIS_H_INCLUDED

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include "data_structure.h"


double brekfis(gsl_vector ** y, gsl_vector **co_variate, size_t total_time, double *y_time, const ParamConfig *config,  ParamInit *init, Param *param);

/**
 * run brekfis. It will call gsl_multimin_fminimizer to minize the brekfis_obj() function with nmsimplex2 algorithm.
 * More details about the algorithm is available at: http://www.gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-without-Derivatives.html#Multimin-Algorithms-without-Derivatives
 * INPUT:
 * @param y the observation
 * @param co_variate the observed co-variates, similar to controlling signal of the ekf.
 * @param total_time total number of time points
 * @param pc model configuration
 * @param pi initial value for some parameter
 * @param par the model and user-defined function parameters
 * OUTPUT:
 * @param fin_pi the final parameter initial values
 * @param fin_par the learned model and user-defined function parameters.
 */


void model_constraint_par(const ParamConfig *pc, Param *par);
void model_constraint_init(const ParamConfig *pc, ParamInit *pi);


/****************************Extended Kim Filter************************/
/**
* This function implements the extended Kim Filter
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* y -- the data
* y_time -- continous time points in the data
* co_variate -- covariates
* config -- model configuration,
* param -- parameters

*
* eta_regime_j_t -- eta^k_it|t*
* error_cov_regime_j_t -- error_cov^k_it|t*
* *
* Output*
* *
* **>>Output via using pointers<<**
* eta_regime_j_t -- eta^k_it|t*
* error_cov_regime_j_t -- error_cov^k_it|t*
* pr_t_given_t_minus_1 -- Pr(S_it=k|Y_i,t-1) *
* pr_t -- Pr(S_it=k|Y_it) *
* eta_regime_jk_pred -- eta^regime_jk_it|t-1 *
* error_cov_regime_jk_pred -- error_cov^regime_jk_it|t-1 *
* eta_regime_jk_t_plus_1 -- filtered regime specific state estimate *
* error_cov_regime_jk_t_plus_1 -- filtered regime specific error covariance estimate *
* innov_v -- innovation vector
* inv_residual_cov -- inverse of the residual covariance
* eta_t -- filtered state estimate
* error_cov_t -- filtered state estimate
**/

double EKimFilter(gsl_vector ** y, gsl_vector **co_variate, double *y_time, const ParamConfig *config, ParamInit *init, Param *param,
    gsl_vector ***eta_regime_j_t,gsl_matrix ***error_cov_regime_j_t,gsl_vector ****eta_regime_jk_pred,gsl_matrix ****error_cov_regime_jk_pred,gsl_vector ****eta_regime_jk_t_plus_1,gsl_matrix ****error_cov_regime_jk_t_plus_1,
    gsl_vector **pr_t, gsl_vector **pr_t_given_t_minus_1,gsl_vector ****innov_v,gsl_matrix ****inv_residual_cov,
	gsl_vector **eta_t, gsl_matrix **error_cov_t);

/****************************Extended Kim Smoother************************/
/**
* This function implements the extended Kim Smoother
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* y_time -- continous time points in the data
* co_variate -- covariates
* config -- model configuration,
* param -- parameters
* pr_t_given_t_minus_1 -- Pr(S_it=k|Y_i,t-1) *
* pr_t -- Pr(S_it=k|Y_it) *
* eta_regime_jk_pred -- eta^regime_jk_it|t-1 *
* error_cov_regime_jk_pred -- error_cov^regime_jk_it|t-1 *
* eta_regime_j_t -- eta^k_it|t*
* error_cov_regime_j_t -- error_cov^k_it|t*
* *
* Output*
* *
* **>>Output via using pointers<<**
* transprob_T -- Pr(S_i,t+1=h, S_it=k|Y_iT)*
* pr_T -- Pr(S_it=k|Y_iT)*
* eta_smooth -- eta_it|T *
* error_cov_smooth -- error_cov_it|T *
* *
**/

void EKimSmoother(double *y_time, gsl_vector **co_variate, const ParamConfig *config, const Param *param,
    gsl_vector **pr_t_given_t_minus_1, gsl_vector **pr_t, 
	gsl_vector ****eta_regime_jk_pred,gsl_matrix ****error_cov_regime_jk_pred,gsl_vector ***eta_regime_j_t,gsl_matrix ***error_cov_regime_j_t,
    gsl_vector **eta_smooth,gsl_matrix **error_cov_smooth,gsl_vector **pr_T,gsl_vector ***transprob_T);

#endif
