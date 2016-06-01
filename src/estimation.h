#ifndef ESTIMATION_H_INCLUDED
#define ESTIMATION_H_INCLUDED

#include <math.h>/*sqrt(double)*/

#include <string.h>

/*#include <stdio.h>*/
#include "nlopt.h"
#include "math_function.h"
#include "numeric_derivatives.h"
#include <gsl/gsl_matrix.h>



int opt_nlopt(void *my_func_data,size_t num_func_param,double *ub,double *lb,double *minf,double *fittedpar,gsl_matrix *Hessian_mat,gsl_matrix *inv_Hessian_mat,double *xtol_rel, double *stopval, double *ftol_rel, double *ftol_abs, int *maxeval, double *maxtime);
#endif

