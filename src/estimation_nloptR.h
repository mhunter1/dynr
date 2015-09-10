
#include <math.h>/*sqrt(double)*/
#include <stdio.h>
#include <string.h>
#include <nlopt.h>
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
#include "PANAmodel.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

int opt_nlopt(void *my_func_data,size_t num_func_param,double *ub,double *lb,double *minf,double *fittedpar,gsl_matrix *Hessian_mat,gsl_matrix *inv_Hessian_mat,double x_tol);

SEXP getListElement(SEXP list, const char *str);

SEXP main_R(SEXP model_list,SEXP data_list);



