/*
Author: Lu Ou
Date: 2015-07-30
Filename: estimation_nlopt.c
Purpose: Obtain parameters estimates by minimizing the negative log likelihood function
Run with
   ./runnlopt_onair.sh
   type in the command: gsl-config --cflags to find the compiler flag
   type in the command: gsl-config --libs to find the flag
Note
  You must have NLOpt installed on your machine to run this.
  Installing NLOpt basically entailed
    1.  Download the tar.gz
    2.  Extract the tar bar
    3.  cd into trunk of folder
    4.  ./configure
    5.  make
    6.  sudo make install
  You may also need to change -I/dir and -L/dir to the locations
    where NLOpt was installed.
*/

#include <math.h>/*sqrt(double)*/
#include <stdio.h>
#include <string.h>
#include "nlopt.h"
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

int opt_nlopt(void *my_func_data,size_t num_func_param,double *ub,double *lb,double *minf,double *fittedpar,gsl_matrix *Hessian_mat,gsl_matrix *inv_Hessian_mat,double *xtol_rel, double *stopval, double *ftol_rel, double *ftol_abs, int *maxeval, double *maxtime)
{
    printf("Optimize function called.\n");
    nlopt_opt opt;

    
	
    /*opt = nlopt_create(NLOPT_LD_MMA, 2); */
    opt = nlopt_create(NLOPT_LD_SLSQP, num_func_param); /* algorithm and dimensionality */
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, myfunc_wrapper,my_func_data);	
    nlopt_set_xtol_rel(opt, * xtol_rel);
    if(*stopval==-9999){
    	nlopt_set_stopval(opt, -HUGE_VAL);
    } else{
		nlopt_set_stopval(opt, * stopval);
	}
    nlopt_set_ftol_rel(opt, *ftol_rel);
    nlopt_set_ftol_abs(opt, *ftol_abs);
    nlopt_set_maxeval(opt, *maxeval);
    nlopt_set_maxtime(opt, *maxtime);

	
    int status=nlopt_optimize(opt, fittedpar, minf);
    if ( status< 0) {
		/*printf("nlopt failed!\n");*/
    }else{
		/*printf("found minimum at \n");
		print_array(fittedpar,num_func_param);
		printf("\n f = %0.10g\n", *minf);*/
		hessianR(fittedpar,my_func_data,function_neg_log_like, *minf, Hessian_mat);/*information matrix*/
		mathfunction_inv_matrix(Hessian_mat, inv_Hessian_mat);/*variance*/
		/*printf("The hessian matrix is \n");
		print_matrix(Hessian_mat);
		printf("\n");
		printf("The inverse hessian matrix is \n");
		print_matrix(inv_Hessian_mat);
		printf("\n");*/
	}

	
	nlopt_destroy(opt);

	

	/*printf("Done.\n");*/
	return status;
}



