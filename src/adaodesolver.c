 /********************************************************
* The Adaptive ODE solver *
* @Author: Lu Ou.                       *
* @created Fall 2014
* Purpose:    *
* Usage:                                        *
* References: A                                   *
* File formats:                                 *
* Restrictions:*
* Revision history:*
* Error handling:*
* Notes:*
********************************************************/
#include <stdio.h>/* standard i/o, printf */
#include <string.h> /*srting functions*/
#include <stdlib.h>
#include "math_function.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>/*sqrt(double)*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "print_function.h"
/******************************************************************************
* adaptive_ode
* *
* Solves an ode by adaptively adjusting the step-size with global and local error control *
* *
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* x -- a p*T matrix of the state estimates*
* t -- a vector of time points*
* tau -- a vector of the step-sizes *
* global_error -- a p*T matrix of the step-sizes *
* global_error_norm -- a vector of measures of the global error *
* *
* **>>Parameters/Input Constants specified by user<<**
* xstart -- a vector of the starting values for vector x *
* tstart -- the staring time *
* tend -- the staring time *
* g -- a pointer to a function, the derivative function, dx/dt, of x and t, *g is the function *
* gparameters -- an array of the parameters in the function g*
* *
* **>>Parameters/Input Constants specified by default but can be changed<<**
* tau_max -- the limit of the step-size*
* global_error_limit -- the parameters in the function g *
* *
* Output*
* *
* **>>Output via using pointers<<**
* x -- a p*T matrix of the state estimates*
* t -- a vector of time points*
* tau -- a vector of the step-sizes *
* global_error -- a vector of the step-sizes *
* global_error_norm -- a vector of measures of the global error *
* *
* **>>Returns<<** Leave it here for now
* x(tend)*
*********************************************/
/*change the variables and values here and the iterative functions in adaptive_ode() according to the pairs you select*/
const double a1=0.8849002;/*(9+4*sqrt(3))/18;automatically double 3.0*/
const double a2=0.1150998;/*(9-4*sqrt(3))/18;*/
const double d1=0.1314459;/*(3+sqrt(3))/36;*/
const double d2=0.03522081;/*(3-sqrt(3))/36;*/
const double c1=0.2113249;/*(3-sqrt(3))/6;*/
const double c2=0.7886751;/*(3+sqrt(3))/6;*/
const double b=0.5;
const double highorder=4.0;
const size_t n_alloc=5000;
/*const double tau_max=0.1;*/
/*const double global_error_limit=0.05;*/

void adaptive_ode(const double tstart, const double tend,
	const gsl_vector *xstart,
	const double tau_max,const double global_error_limit,
	gsl_vector *t, gsl_vector *tau, gsl_vector *global_error_norm,
	gsl_matrix *x, gsl_matrix *global_error, size_t regime,
	double *gparameters, void (*g)(double, size_t, const gsl_vector *,double *, gsl_vector *)
	){

    int M=1;
    int i;
    int np=xstart->size;
    double local_error_limit=pow(global_error_limit,(highorder-1)/(highorder-2));
    double local_error_norm;
    double taul_new;

    gsl_vector *x_l12=gsl_vector_alloc(np);/*stage values*/
    gsl_vector *x_l22=gsl_vector_alloc(np);
    gsl_vector *local_error=gsl_vector_alloc(np);
    gsl_vector *x_column_lplus1=gsl_vector_alloc(np);
    gsl_vector *g_l=gsl_vector_calloc(np);
    gsl_vector *g_lplus1=gsl_vector_calloc(np);
    gsl_vector *g_integral=gsl_vector_alloc(np);
    gsl_vector *g_xl12=gsl_vector_calloc(np);
    gsl_vector *g_xl22=gsl_vector_calloc(np);
    (*g)(tstart,regime,xstart,gparameters,g_l);

    size_t index=0;
    for(index=0; index<xstart->size; index++){
    gsl_matrix_set(x, index, 0, gsl_vector_get(xstart,index));
    gsl_matrix_set(global_error, index, 0, 0);
    }

    gsl_vector_set(t, 0, tstart);
    gsl_vector_set(tau, 0, tau_max);
    gsl_vector_set(global_error_norm, 0, 0);



    while(M>0){
    	size_t l=0;
    	M=0;
    	double sum=0;

    	gsl_vector_set(global_error_norm,0,0);

        while((gsl_vector_get(t,l)<tend)&&(gsl_vector_get(global_error_norm,l)<=100*global_error_limit)){

            gsl_vector_set(t,l+1,gsl_vector_get(t,l)+gsl_vector_get(tau,l));
            for(index=0; index<np; index++)
                gsl_matrix_set(x,index,l+1,gsl_matrix_get(x,index,l));

            for(i=3;i!=0;i--){
                gsl_matrix_get_col(x_column_lplus1,x,l+1);
                (*g)(gsl_vector_get(t,l+1),regime,x_column_lplus1,gparameters,g_lplus1);

                for(index=0; index<np; index++){
                    gsl_vector_set(x_l12,index,a1*gsl_matrix_get(x,index,l)+a2*gsl_matrix_get(x,index,l+1)+gsl_vector_get(tau,l)*(d1*gsl_vector_get(g_l,index)-d2*gsl_vector_get(g_lplus1,index)));
                    gsl_vector_set(x_l22,index,a2*gsl_matrix_get(x,index,l)+a1*gsl_matrix_get(x,index,l+1)+gsl_vector_get(tau,l)*(d2*gsl_vector_get(g_l,index)-d1*gsl_vector_get(g_lplus1,index)));
                }
                (*g)(gsl_vector_get(t,l)+c1*gsl_vector_get(tau,l),regime,x_l12,gparameters,g_xl12);
                (*g)(gsl_vector_get(t,l)+c2*gsl_vector_get(tau,l),regime,x_l22,gparameters,g_xl22);
                for(index=0; index<np; index++)
                    gsl_vector_set(g_integral, index,gsl_vector_get(tau,l)*b*(gsl_vector_get(g_xl12,index)+gsl_vector_get(g_xl22,index)));
                for(index=0; index<np; index++)
                    gsl_matrix_set(x,index,l+1,gsl_matrix_get(x,index,l)+gsl_vector_get(g_integral,index));
            }



            for(index=0; index<np; index++)
            	gsl_vector_set(local_error,index,gsl_vector_get(tau,l)*b*(gsl_vector_get(g_l,index)+gsl_vector_get(g_lplus1,index))-gsl_vector_get(g_integral,index));
            local_error_norm=gsl_blas_dasum(local_error);

            taul_new=0.8*gsl_vector_get(tau,l)*pow(local_error_limit/local_error_norm,1/(highorder-1));

            if (local_error_norm>local_error_limit){
                gsl_vector_set(tau,l,taul_new);
            }else {
                for(index=0; index<np; index++)
                    gsl_matrix_set(global_error,index,l+1,gsl_matrix_get(global_error,index,l)+gsl_vector_get(local_error,index));
                sum=0;
                for(index=0; index<np; index++)
                    sum+=fabs(gsl_matrix_get(global_error,index,l+1));
                gsl_vector_set(global_error_norm,l+1,sum);
                if (gsl_vector_get(global_error_norm,l+1)>global_error_limit)
                    M=1;
                gsl_vector_set(tau,l+1,mathfunction_min(taul_new, tend-gsl_vector_get(t,l+1),tau_max));

                l=l+1;
                gsl_vector_memcpy(g_l,g_lplus1);
            }
        }
        if (M>0){
            local_error_limit=local_error_limit*pow(0.8*global_error_limit/gsl_vector_max(global_error_norm),(highorder-1)/(highorder-2));
        }
    }



    gsl_vector_free(x_l12);
    gsl_vector_free(x_l22);
    gsl_vector_free(local_error);
    gsl_vector_free(g_l);
    gsl_vector_free(g_lplus1);
    gsl_vector_free(g_integral);
    gsl_vector_free(x_column_lplus1);
    gsl_vector_free(g_xl12);
    gsl_vector_free(g_xl22);

}

void function_F(double t, size_t regime, const gsl_vector *x,double *gparameters, gsl_vector *F){
    gsl_vector_set(F,0,gsl_vector_get(x,1));
    gsl_vector_set(F,1,gparameters[0]*(1-pow(gsl_vector_get(x,0),2))*gsl_vector_get(x,1)-gsl_vector_get(x,0));
}

void debug_adaptive_ode(){
    size_t np=2;
    size_t regime=1;
    gsl_vector *xstart=gsl_vector_alloc(np);
    gsl_vector *t=gsl_vector_alloc(n_alloc);
    gsl_vector *tau=gsl_vector_alloc(n_alloc);
    gsl_vector *global_error_norm=gsl_vector_alloc(n_alloc);
    double gparameters[]={3.0};
    gsl_matrix *x=gsl_matrix_alloc(np, n_alloc);
    gsl_matrix *global_error=gsl_matrix_alloc(np, n_alloc);


    gsl_vector_set(xstart,0,0.5);
    gsl_vector_set(xstart,1,0.5);


    adaptive_ode(0.0,  10.0,xstart,0.1,0.05,t,tau,global_error_norm,
    	    x,global_error,regime,
    	    gparameters,function_F);




    {
     /*FILE * f = fopen("test.dat", "w");
     print_matrix(x);
     print_vector(t);
     print_vector(tau);
     fclose (f);*/
    }

    gsl_vector_free(xstart);
    gsl_vector_free(t);
    gsl_vector_free(tau);
    gsl_vector_free(global_error_norm);
    gsl_matrix_free(x);
    gsl_matrix_free(global_error);


}
/******************************************************************************
* adaptive_ode_kf (a wrapper of the function adaptive_ode)
* *
* Solves an ode by adaptively adjusting the step-size with global and local error control *
* Returns x(tend)
* *
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* x_tend -- a vector of x(tend)*
* *
* **>>Parameters/Input Constants specified by user<<**
* tstart -- the staring time *
* tend -- the staring time *
* regime -- the regime *
* xstart -- a vector of the starting values for vector x *
* gparameters -- an array of the parameters in the function g *
* co_variate -- a vector of covariates *
* g -- a pointer to a function, the derivative function, dx/dt, of x and t, *g is the function *
* tau_max -- the limit of the step-size*
* global_error_limit -- the parameters in the function g *
* *
* Output*
* *
* **>>Output via using pointers<<**
* x_tend -- a vector of x(tend)*
* *
*********************************************/
void adaptive_ode_kf(const double tstart, const double tend,
	const gsl_vector *xstart,
	const double tau_max,const double global_error_limit, size_t regime,
	double *gparameters, size_t n_gparam, const gsl_vector *co_variate, void (*g)(double, size_t,const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
	gsl_vector *x_tend
	){

    int M=1;
    int i;
    size_t np=xstart->size;
    double local_error_limit=pow(global_error_limit,(highorder-1)/(highorder-2));
    double local_error_norm;
    double taul_new;
    size_t l;

    gsl_vector *t=gsl_vector_alloc(n_alloc);
    gsl_vector *tau=gsl_vector_alloc(n_alloc);
    gsl_vector *global_error_norm=gsl_vector_alloc(n_alloc);
    gsl_matrix *x=gsl_matrix_alloc(np, n_alloc);
    gsl_matrix *global_error=gsl_matrix_alloc(np, n_alloc);

    gsl_vector *x_l12=gsl_vector_alloc(np);/*stage values*/
    gsl_vector *x_l22=gsl_vector_alloc(np);
    gsl_vector *local_error=gsl_vector_alloc(np);
    gsl_vector *x_column_lplus1=gsl_vector_alloc(np);
    gsl_vector *g_l=gsl_vector_calloc(np);
    gsl_vector *g_lplus1=gsl_vector_calloc(np);
    gsl_vector *g_integral=gsl_vector_alloc(np);
    gsl_vector *g_xl12=gsl_vector_calloc(np);
    gsl_vector *g_xl22=gsl_vector_calloc(np);
    (*g)(tstart,regime,xstart, gparameters, n_gparam,co_variate,g_l);

    size_t index=0;
    for(index=0; index<xstart->size; index++){
    gsl_matrix_set(x, index, 0, gsl_vector_get(xstart,index));
    gsl_matrix_set(global_error, index, 0, 0);
    }

    gsl_vector_set(t, 0, tstart);
    gsl_vector_set(tau, 0, tau_max);
    gsl_vector_set(global_error_norm, 0, 0);



    while(M>0){
    	l=0;
    	M=0;
    	double sum=0;

    	gsl_vector_set(global_error_norm,0,0);

        while((gsl_vector_get(t,l)<tend)&&(gsl_vector_get(global_error_norm,l)<=100*global_error_limit)){
            /*printf("%f ",gsl_vector_get(t,l));
            printf("%f ",gsl_vector_get(tau,l));
            printf("%lu",l);
            printf("\n");*/

            gsl_vector_set(t,l+1,gsl_vector_get(t,l)+gsl_vector_get(tau,l));
            for(index=0; index<np; index++)
                gsl_matrix_set(x,index,l+1,gsl_matrix_get(x,index,l));

            for(i=3;i!=0;i--){
                gsl_matrix_get_col(x_column_lplus1,x,l+1);
                (*g)(gsl_vector_get(t,l+1),regime, x_column_lplus1,gparameters, n_gparam, co_variate,g_lplus1);

                for(index=0; index<np; index++){
                    gsl_vector_set(x_l12,index,a1*gsl_matrix_get(x,index,l)+a2*gsl_matrix_get(x,index,l+1)+gsl_vector_get(tau,l)*(d1*gsl_vector_get(g_l,index)-d2*gsl_vector_get(g_lplus1,index)));
                    gsl_vector_set(x_l22,index,a2*gsl_matrix_get(x,index,l)+a1*gsl_matrix_get(x,index,l+1)+gsl_vector_get(tau,l)*(d2*gsl_vector_get(g_l,index)-d1*gsl_vector_get(g_lplus1,index)));
                }
                (*g)(gsl_vector_get(t,l)+c1*gsl_vector_get(tau,l),regime,x_l12,gparameters, n_gparam, co_variate,g_xl12);
                (*g)(gsl_vector_get(t,l)+c2*gsl_vector_get(tau,l),regime,x_l22,gparameters, n_gparam, co_variate,g_xl22);
                for(index=0; index<np; index++)
                    gsl_vector_set(g_integral, index,gsl_vector_get(tau,l)*b*(gsl_vector_get(g_xl12,index)+gsl_vector_get(g_xl22,index)));
                for(index=0; index<np; index++)
                    gsl_matrix_set(x,index,l+1,gsl_matrix_get(x,index,l)+gsl_vector_get(g_integral,index));
            }



            for(index=0; index<np; index++)
            	gsl_vector_set(local_error,index,gsl_vector_get(tau,l)*b*(gsl_vector_get(g_l,index)+gsl_vector_get(g_lplus1,index))-gsl_vector_get(g_integral,index));
            local_error_norm=gsl_blas_dasum(local_error);

            taul_new=0.8*gsl_vector_get(tau,l)*pow(local_error_limit/local_error_norm,1/(highorder-1));
            
            /*added by Lu on Jul., 13, 2015*/
            /*if (taul_new<(tend-tstart)/n_alloc){taul_new=(tend-tstart)/n_alloc;}*/
            /*added by Lu on Jul., 13, 2015*/
            if (local_error_norm>local_error_limit){
                gsl_vector_set(tau,l,taul_new);
            }else {
                for(index=0; index<np; index++)
                    gsl_matrix_set(global_error,index,l+1,gsl_matrix_get(global_error,index,l)+gsl_vector_get(local_error,index));
                sum=0;
                for(index=0; index<np; index++)
                    sum+=fabs(gsl_matrix_get(global_error,index,l+1));
                gsl_vector_set(global_error_norm,l+1,sum);
                if (gsl_vector_get(global_error_norm,l+1)>global_error_limit)
                    M=1;
                gsl_vector_set(tau,l+1,mathfunction_min(taul_new, tend-gsl_vector_get(t,l+1),tau_max));

                l=l+1;
                gsl_vector_memcpy(g_l,g_lplus1);
            }
        }
        if (M>0){
            local_error_limit=local_error_limit*pow(0.8*global_error_limit/gsl_vector_max(global_error_norm),(highorder-1)/(highorder-2));
        }
    }
    for(index=0; index<np; index++)
    	    gsl_vector_set(x_tend,index,gsl_matrix_get(x,index,l));

    /*print_vector(t);*/
    /*print_vector(global_error_norm);*/
    /*print_matrix(x);*/


    gsl_vector_free(x_l12);
    gsl_vector_free(x_l22);
    gsl_vector_free(local_error);
    gsl_vector_free(g_l);
    gsl_vector_free(g_lplus1);
    gsl_vector_free(g_integral);
    gsl_vector_free(x_column_lplus1);
    gsl_vector_free(g_xl12);
    gsl_vector_free(g_xl22);


    gsl_vector_free(t);
    gsl_vector_free(tau);
    gsl_vector_free(global_error_norm);

    gsl_matrix_free(x);
    gsl_matrix_free(global_error);

}

void function_F_debug(double t, size_t regime, const gsl_vector *x,double *gparameters, size_t n_gparam, const gsl_vector *co_variate, gsl_vector *F){
    gsl_vector_set(F,0,gsl_vector_get(x,1));
    gsl_vector_set(F,1,gparameters[0]*(1-pow(gsl_vector_get(x,0),2))*gsl_vector_get(x,1)-gsl_vector_get(x,0));
}
void debug_adaptive_ode_kf(){
    size_t np=2;
    size_t regime=1;
    gsl_vector *xstart=gsl_vector_alloc(np);
    gsl_vector *x_tend=gsl_vector_alloc(np);
    double gparameters[]={3.0};
	size_t n_gparam=1;

    gsl_vector_set(xstart,0,0.5);
    gsl_vector_set(xstart,1,0.5);


    adaptive_ode_kf(0, 10,xstart,0.1,0.05,regime, gparameters, n_gparam, NULL, function_F_debug,x_tend);

     print_vector(x_tend);


    gsl_vector_free(xstart);
    gsl_vector_free(x_tend);



}

/******************************************************************************
* rk4_odesolver
* *
* Solves the dynamic functions (an ode) by by means of fourth-order *
* Returns x(tend)
* *
* Parameters/Input *
* *
* **>>Parameters/Input Pointers<<**
* x_tend -- a vector of x(tend)*
* *
* **>>Parameters/Input Constants specified by user<<**

* tstart -- the staring time *
* tend -- the staring time *
* regime -- the regime *
* xstart -- a vector of the starting values for vector x *
* gparameters -- an array of the parameters in the function g *
* co_variate -- a vector of covariates *
* g -- a pointer to a function, the derivative function, dx/dt, of x and t, *g is the function *
* *
* Output*
* *
* **>>Output via using pointers<<**
* x_tend -- a vector of x(tend)*
* *
*********************************************/
void rk4_odesolver(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *gparameters, size_t n_gparam,const gsl_vector *co_variate,
        void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
	gsl_vector *x_tend){

        size_t np=xstart->size;
        gsl_vector *k1=gsl_vector_calloc(np);
        gsl_vector *k2=gsl_vector_calloc(np);
        gsl_vector *k3=gsl_vector_calloc(np);
        gsl_vector *k4=gsl_vector_calloc(np);
        gsl_vector *x1=gsl_vector_alloc(np);
        gsl_vector *x2=gsl_vector_alloc(np);
        gsl_vector *x3=gsl_vector_alloc(np);
        
        /*printf("xstart:");
        print_vector(xstart);*/

	double delta=tend-tstart;
	(*g)(tstart, regime, xstart, gparameters, n_gparam, co_variate, k1);
	gsl_vector_scale(k1, delta/2);/*k1<-delta/2*k1*/
	gsl_vector_memcpy(x1, xstart);
	gsl_vector_add(x1, k1);/*x1<-xstart+delta/2*k1*/
	gsl_vector_scale(k1,1.0/3);/*k1<-delta/6*k1*/
        (*g)(tstart, regime, x1, gparameters, n_gparam, co_variate, k2);
        gsl_vector_scale(k2, delta/2);/*k2<-delta/2*k2*/
	gsl_vector_memcpy(x2, xstart);
	gsl_vector_add(x2, k2);/*x1<-xstart+delta/2*k2*/
	gsl_vector_scale(k2,2.0/3);/*k2<-delta/3*k2*/
	gsl_vector_add(k1,k2);/*k1<-delta/6*k1+delta/3*k2*/
        (*g)(tstart, regime, x2, gparameters, n_gparam, co_variate, k3);
        gsl_vector_scale(k3, delta);/*k3<-delta*k3*/
	gsl_vector_memcpy(x3, xstart);
	gsl_vector_add(x3, k3);/*x3<-xstart+delta*k3*/
        gsl_vector_scale(k3,1.0/3);/*k3<-delta/3*k3*/
	gsl_vector_add(k1,k3);/*k1<-delta/6*k1+delta/3*k2+delta/3*k3*/
	(*g)(tstart, regime, x3, gparameters, n_gparam, co_variate, k4);
        gsl_vector_scale(k4,delta/6);/*k4<-delta/6*k4*/
	gsl_vector_add(k1,k4);/*k1<-delta/6*k1+delta/3*k2+delta/3*k3+delta/6*k4*/

	gsl_vector_memcpy(x_tend, xstart);
	gsl_vector_add(x_tend, k1);
	
	/*printf("xend:");
	print_vector(x_tend);*/
	gsl_vector_free(k1);
	gsl_vector_free(k2);
	gsl_vector_free(k3);
	gsl_vector_free(k4);
	gsl_vector_free(x1);
	gsl_vector_free(x2);
	gsl_vector_free(x3);

}
