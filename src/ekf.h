#ifndef EKF_H_INCLUDED
#define EKF_H_INCLUDED
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "adaodesolver.h"
/******************************************************************************
* Discrete/continuous-discrete extended kalman filter (EKF)
* *
* calculate the filtered estimates and error covariance matrix using CDA-EKF
* store the innovation vector and innovation covariance matrix for calculating log-likelihood
* *
* Parameters/Input *
* *
* **>>Parameters/Input Constants specified by user<<**
* t -- the previous time t
* t_plus_1 -- the current time t+1
* regime -- which regime is used
* eta_t -- a vector of the filtered estimates at time t
* error_cov_t -- a vector of the elements in the filtered error covariance matrix at time t: P_t[1,1],P_t[2,2],P_t[3,3],P_t[1,2],P_t[1,3],P_t[2,3]
* y_t_plus_1 -- data y_t+1
* co_variate -- covariates
* eta_noise_cov -- process noise covariance
* y_noise_cov -- measurement error covariance i.e. Rk or R
* H_t_plus_1 -- Lambda matrix in the measurement function/model, can be time-varying
* params -- array of parameters
* func_measure -- a pointer to the measurement function without measurement error
* func_dx_dt -- a pointer to a function. the derivative function, dx/dt, *dx_dt is the function *
* func_dP_dt -- a pointer to a function. the derivative function, dP/dt, *dP_dt is the function *
* *
* **>>Parameters/Input Pointers<<**
* eta_t_plus_1 -- a vector of filtered new estimates at time t+1
* error_cov_t_plus_1 -- a vector of the elements in the filtered error covariance matrix at time t: P_t[1,1],P_t[2,2],P_t[3,3],P_t[1,2],P_t[1,3],P_t[2,3]
* innov_v-- innovation vector
* innov_cov-- innovation covariance matrix
* *
* Output*
* *
* **>>Output via using pointers<<**
* eta_t_plus_1 -- a vector of filtered new estimates at time t+1
* error_cov_t_plus_1 -- a vector of the elements in the filtered error covariance matrix at time t: P_t[1,1],P_t[2,2],P_t[3,3],P_t[1,2],P_t[1,3],P_t[2,3]
* innov_v-- innovation vector
* innov_cov-- innovation covariance matrix
* *
*********************************************/

double ext_kalmanfilter(size_t t, size_t regime,
        gsl_vector *eta_t,  gsl_matrix *error_cov_t,
	const gsl_vector *y_t_plus_1,const gsl_vector *co_variate, const double *y_time,
	const gsl_matrix *eta_noise_cov, const gsl_matrix *y_noise_cov,
        double *params,size_t num_func_param,
		bool isContinuousTime,
        void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *),
        void (*func_dx_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
        void (*func_dP_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
        void (*func_dynam)(const double, const double, size_t, const gsl_vector *,double *, size_t, const gsl_vector *, void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),gsl_vector *),
        gsl_vector *eta_t_plus_1, gsl_matrix *error_cov_t_plus_1, gsl_vector *innov_v, gsl_matrix *inv_innov_cov);

double ext_kalmanfilter_updateonly(size_t t, size_t regime,
     gsl_vector *eta_t,  gsl_matrix *error_cov_t,
	const gsl_vector *y_t_plus_1,const gsl_vector *co_variate, const double *y_time,
	const gsl_matrix *eta_noise_cov, const gsl_matrix *y_noise_cov,
        double *params,
        void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *),
        gsl_vector *eta_t_plus_1, gsl_matrix *error_cov_t_plus_1, gsl_vector *innov_v, gsl_matrix *inv_innov_cov);
size_t find_miss_data(const gsl_vector *y, gsl_vector *non_miss);
double ext_kalmanfilter_smoother(size_t t, size_t regime,
        gsl_vector *eta_t,  gsl_matrix *error_cov_t,
	const gsl_vector *y_t_plus_1,const gsl_vector *co_variate, const double *y_time,
	const gsl_matrix *eta_noise_cov, const gsl_matrix *y_noise_cov,
        double *params,size_t num_func_param,
		bool isContinuousTime,
        void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *),
        void (*func_dx_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
        void (*func_dP_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
        void (*func_dynam)(const double, const double, size_t, const gsl_vector *,double *, size_t, const gsl_vector *, void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),gsl_vector *),
        gsl_vector *eta_pred, gsl_matrix *error_cov_pred, gsl_vector *eta_t_plus_1, gsl_matrix *error_cov_t_plus_1, gsl_vector *innov_v, gsl_matrix *inv_innov_cov);
double ext_kalmanfilter_updateonly_smoother(size_t t, size_t regime,
     gsl_vector *eta_t,  gsl_matrix *error_cov_t,
	const gsl_vector *y_t_plus_1,const gsl_vector *co_variate, const double *y_time,
	const gsl_matrix *eta_noise_cov, const gsl_matrix *y_noise_cov,
    double *params,
    void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *),
    gsl_vector *eta_pred, gsl_matrix *error_cov_pred, gsl_vector *eta_t_plus_1, gsl_matrix *error_cov_t_plus_1, gsl_vector *innov_v, gsl_matrix *inv_innov_cov);
#endif
