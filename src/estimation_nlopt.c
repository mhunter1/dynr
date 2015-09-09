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


int main()
{
    size_t num_func_param=6;
    double ub[6] = {4, 4, 4, 4, 4, 4};
    double lb[6] = {-4,-4,-4,-4,-12, -12}; /* lower bounds */
    nlopt_opt opt;
    size_t index;
	
    	/*opt = nlopt_create(NLOPT_LD_MMA, 2); */
	opt = nlopt_create(NLOPT_LD_SLSQP, num_func_param); /* algorithm and dimensionality */
        nlopt_set_upper_bounds(opt, ub);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective(opt, myfunc_wrapper, NULL);
	
	nlopt_set_xtol_rel(opt, 1e-5);
	
	double params[]={log(1),log(2),0,0,-10,-10};/* some initial guess*/
    /*log(1.2)=0.1823216,log(1.8)=0.5877867,-0.5,-0.5,log(0.0001)=-9.21034,log(0.0001)=-9.21034*/
	double minf; /* the minimum objective value, upon return */
    
    gsl_matrix *Hessian_mat=gsl_matrix_calloc(num_func_param,num_func_param);
    gsl_matrix *inv_Hessian_mat=gsl_matrix_calloc(num_func_param,num_func_param);
	
	if (nlopt_optimize(opt, params, &minf) < 0) {
		printf("nlopt failed!\n");
	}
	else {
		printf("found minimum at \n");
        print_array(params,num_func_param);
        printf("\n f = %0.10g\n", minf);
        hessian(params,function_neg_log_like, minf, Hessian_mat);/*information matrix*/
        mathfunction_inv_matrix(Hessian_mat, inv_Hessian_mat);/*variance*/
        printf("The hessian matrix is \n");
        print_matrix(Hessian_mat);
        printf("\n");
        printf("The inverse hessian matrix is \n");
        print_matrix(inv_Hessian_mat);
        printf("\n");
	}
	
	nlopt_destroy(opt);
    gsl_matrix_free(Hessian_mat);
    gsl_matrix_free(inv_Hessian_mat);
	
	return 0;
}


