
#include <math.h>/*sqrt(double)*/
/*#include <stdio.h>*/
#include <string.h>
#include "nlopt.h"
#include "math_function.h"
#include "ekf.h"
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
#include "estimation.h"
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

SEXP getListElement(SEXP list, const char *str);

SEXP main_R(SEXP model_list, SEXP data_list, SEXP weight_flag_in, SEXP debug_flag_in, SEXP optimization_flag_in, SEXP hessian_flag_in, SEXP verbose_flag_in);



