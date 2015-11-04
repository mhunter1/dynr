
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "math_function.h"
#include "brekfis.h"
#include "adaodesolver.h"
#include <gsl/gsl_matrix.h>



void function_dynam_ada(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *gparameters,const gsl_vector *co_variate,
        void (*g)(double, size_t, const gsl_vector *, double *, const gsl_vector *, gsl_vector *),
        gsl_vector *x_tend){

        double tau_max=(tend-tstart)/10;/*specify tau_max*/
        double global_error_limit=10;/*specify global error limit*/
        adaptive_ode_kf(tstart, tend, xstart,tau_max,global_error_limit, regime,  gparameters, co_variate, (*g),x_tend);

}


void function_jacobdynamic(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *param, size_t num_func_param, const gsl_vector *co_variate,
        void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
	gsl_matrix *Jx){

        size_t np=xstart->size;
        gsl_matrix *k1=gsl_matrix_alloc(np,np);
        gsl_vector *diag=gsl_vector_alloc(np);
        gsl_matrix *k2=gsl_matrix_alloc(np,np);
        gsl_matrix *k3=gsl_matrix_alloc(np,np);
        gsl_matrix *k4=gsl_matrix_alloc(np,np);
        gsl_vector *x1=gsl_vector_alloc(np);
        gsl_vector *x2=gsl_vector_alloc(np);
        gsl_vector *x3=gsl_vector_alloc(np);
        size_t i;
        double params_aug[num_func_param+np];
        for (i=0;i<num_func_param;i++)
            params_aug[i]=param[i];

        /*printf("xstart:");
        print_vector(xstart);*/

	double delta=tend-tstart;
	for (i=0;i<np;i++)
           params_aug[num_func_param+i]=gsl_vector_get(xstart,i);
	(*g)(tstart, regime, params_aug, co_variate, k1);
	mathfunction_diagout_scale(k1,delta/2,diag);
	gsl_vector_memcpy(x1, xstart);
	gsl_vector_add(x1, diag);/*x1<-xstart+delta/2*diag(k1)*/
	gsl_matrix_scale(k1, delta/6);/*k1<-delta/6*k1*/

	for (i=0;i<np;i++)
           params_aug[num_func_param+i]=gsl_vector_get(x1,i);
	(*g)(tstart, regime, params_aug, co_variate, k2);
        mathfunction_diagout_scale(k2,delta/2,diag);
	gsl_vector_memcpy(x2, xstart);
	gsl_vector_add(x2, diag);/*x1<-xstart+delta/2*diagk2*/
	gsl_matrix_scale(k2, delta/3);/*k2<-delta/3*k2*/
	gsl_matrix_add(k1,k2);/*k1<-delta/6*k1+delta/3*k2*/
        for (i=0;i<np;i++)
           params_aug[num_func_param+i]=gsl_vector_get(x2,i);
        (*g)(tstart, regime, params_aug, co_variate, k3);
        mathfunction_diagout_scale(k3,delta,diag);
	gsl_vector_memcpy(x3, xstart);
	gsl_vector_add(x3, diag);/*x3<-xstart+delta*k3*/
        gsl_matrix_scale(k3, delta/3);/*k3<-delta*k3*/
	gsl_matrix_add(k1,k3);/*k1<-delta/6*k1+delta/3*k2+delta/3*k3*/
        for (i=0;i<np;i++)
           params_aug[num_func_param+i]=gsl_vector_get(x3,i);
	(*g)(tstart, regime, params_aug, co_variate, k4);
        gsl_matrix_scale(k4,delta/6);/*k4<-delta/6*k4*/
	gsl_matrix_add(k1,k4);/*k1<-delta/6*k1+delta/3*k2+delta/3*k3+delta/6*k4*/

	gsl_matrix_set_identity(Jx);
	gsl_matrix_add(Jx,k1);

	/*printf("xend:");
	print_vector(x_tend);*/
	gsl_matrix_free(k1);
	gsl_matrix_free(k2);
	gsl_matrix_free(k3);
	gsl_matrix_free(k4);
	gsl_vector_free(x1);
	gsl_vector_free(x2);
	gsl_vector_free(x3);
	gsl_vector_free(diag);

}


