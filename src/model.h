#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

void function_dynam_ada(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
                        double *gparameters,const gsl_vector *co_variate,
                        void (*g)(double, size_t, const gsl_vector *, double *, const gsl_vector *, gsl_vector *),
                        gsl_vector *x_tend);

void function_jacobdynamic(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
                           double *param, size_t num_func_param, const gsl_vector *co_variate,
                           void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
                           gsl_matrix *Jx);
#endif
