/**
 * This file contains the definition of data structures used in the project.
 * @author Peifeng Yin
 * @created Feb. 11, 2014
 */

#ifndef DATA_STRUCTURE_H_INCLUDED
#define DATA_STRUCTURE_H_INCLUDED
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdbool.h>

/**
 * configuration of the model
 */
typedef struct ParamConfig{
    size_t num_regime; /** number of regimes */
    size_t dim_latent_var; /** dimension of latent variable */
    size_t dim_obs_var; /** dimension of observed variable */
    size_t num_func_param; /** the number parameters for user-defined functions **/
    bool second_order; /** whether second-order ekf is used **/
    bool adaodesolver; /** whether adaptive ode solver is used **/
    bool isnegloglikeweightedbyT;/** whether the negative loglikelihood is weighted by individual T**/
    bool isAnalytic; //Whether to use analytic gradient
    size_t dim_co_variate;
    size_t num_sbj; /** number of subjects **/
    size_t *index_sbj;
    size_t total_obs;
    bool isContinuousTime; /** Flag for continuous-time model: 1 = yes; 0 = no**/
    bool verbose_flag; /** Flag for printing verbose output, including every function evaluation; 1 = yes; 0 = no**/
	bool is_cov_formula;

    /** time, regime, parameter, eta_t, co_variate, Hk, y_t **/
    void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *);
    /** double t, size_t regime, const gsl_vector *x,double *param, co_variate, gsl_vector *F_dx_dt**/
    void (*func_dx_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *);
    /**double t, size_t regime, const gsl_vector *x,double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx**/
    void (*func_dF_dx)(double, size_t, double *, const gsl_vector *, gsl_matrix *);
    /**const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *param, size_t num_func_param, const gsl_vector *co_variate,
        void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
	gsl_matrix *Jx**/
    void (*func_jacob_dynam)(const double, const double, size_t, const gsl_vector *,
        double *, size_t,const gsl_vector *,
        void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
		gsl_matrix *);
    /** time, const gsl_vector *p,double *param, gsl_vector *F_dP_dt**/
    void (*func_dP_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *);
    /**double *param, gsl_vector **co_variate, gsl_vector **pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0, size_t index_sbj**/
    void (*func_initial_condition)(double *, gsl_vector **, gsl_vector **, gsl_vector **, gsl_matrix **, size_t *);
    /**size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat**/
    void (*func_regime_switch)(size_t, size_t, double *, const gsl_vector *, gsl_matrix *);
    /**size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov**/
    void (*func_noise_cov)(size_t, size_t, double *, gsl_matrix *, gsl_matrix *,const gsl_vector *);
    /**double *param**/
    void (*func_transform)(double *);
    /** tstart, tend, regime, xstart,gparameters, n_gparam,co_variate, (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),x_tend **/
    void (*func_dynam)(const double, const double, size_t, const gsl_vector *,
        double *, size_t, const gsl_vector *, void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *), gsl_vector *);
 } ParamConfig;

/**
 * The model and the data
 */
typedef struct Data_and_Model{
    ParamConfig pc; /** model configuraiton */
    gsl_vector **y; /** observed variables */
    gsl_vector **co_variate; /** covariates */
    double *y_time; /** observed real times**/

 } Data_and_Model;

typedef struct ParamInit{
    gsl_vector **eta_0; /** initial value for latent variable in different regimes **/
    gsl_matrix **error_cov_0; /** initial value for covariance matrix of latent variable in different regimes **/
    gsl_vector **pr_0; /** the probability of staying at each regime for the start **/
} ParamInit;

/**
 * A collection of parameters: after constraint functions, eta_noise_cov will hold the process noise covariance matrix, not elements in the L and D in LDL' decomposition.
 */
typedef struct Param{
    gsl_matrix *regime_switch_mat; /** regime switch probability **/
    gsl_matrix *eta_noise_cov; /** Q: noise covariance matrix for latent variable **/
    gsl_matrix *y_noise_cov; /** R: noise covariance matrix for observation **/
    double *func_param;
} Param;

#endif
