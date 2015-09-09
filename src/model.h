#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED
#include "data_structure.h"
void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht,gsl_vector *y);
void function_dx_dt(double t, size_t regime, const gsl_vector *x,double *param, const gsl_vector *, gsl_vector *F_dx_dt);
void function_dynam_ada(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
                        double *gparameters,const gsl_vector *co_variate,
                        void (*g)(double, size_t, const gsl_vector *, double *, const gsl_vector *, gsl_vector *),
                        gsl_vector *x_tend);
void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx);
void function_jacobdynamic(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
                           double *param, size_t num_func_param, const gsl_vector *co_variate,
                           void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
                           gsl_matrix *Jx);
void function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, const gsl_vector *, gsl_vector *F_dP_dt);
void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0);
void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat);
void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov);
void function_transform(const ParamConfig *pc, ParamInit *pi, Param *par);
ParamConfig model_configure();
void run_model();

#endif
