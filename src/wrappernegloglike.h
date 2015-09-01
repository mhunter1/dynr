#ifndef WRAPPERNEGLOGLIKE_H_INCLUDED
#define WRAPPERNEGLOGLIKE_H_INCLUDED

#include "headers/brekfis.h"
#include "headers/cdaekf.h"
#include "headers/data_structure.h"
#include "headers/math_function.h"
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include "headers/model.h"
double function_neg_log_like(const double *params, void *data);
#endif
