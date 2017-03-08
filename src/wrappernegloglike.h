#ifndef WRAPPERNEGLOGLIKE_H_INCLUDED
#define WRAPPERNEGLOGLIKE_H_INCLUDED

#include "brekfis.h"
#include "ekf.h"
#include "data_structure.h"
#include "math_function.h"
#include "model.h"
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include "print_function.h"
double function_neg_log_like(const double *params, void *data);
#endif
