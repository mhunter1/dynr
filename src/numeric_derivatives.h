#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include "wrappernegloglike.h"
#include <string.h>
#include <math.h>/*sqrt(double),pow*/

void forward_diff_grad(double *grad_approx, double ref_fit, const double *x, void * data, double (*func_obj)(const double *, void *));

void hessian(const double *x,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian);

double myfunc_wrapper(unsigned n, const double *x, double *grad, void *my_func_data);

void hessianR(const double *x,void *data,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian);

