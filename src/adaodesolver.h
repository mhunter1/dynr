#ifndef ADAODESOLVER_H_INCLUDED
#define ADAODESOLVER_H_INCLUDED

/*#include <stdio.h>*//* standard i/o, printf */
#include <string.h> /*srting functions*/
#include <stdlib.h>
#include "math_function.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>/*sqrt(double)*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

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
static const double a1=0.8849002;/*(9+4*sqrt(3))/18;automatically double 3.0*/
static const double a2=0.1150998;/*(9-4*sqrt(3))/18;*/
static const double d1=0.1314459;/*(3+sqrt(3))/36;*/
static const double d2=0.03522081;/*(3-sqrt(3))/36;*/
static const double c1=0.2113249;/*(3-sqrt(3))/6;*/
static const double c2=0.7886751;/*(3+sqrt(3))/6;*/
static const double b=0.5;
static const double highorder=4.0;
static const size_t n_alloc=5000;
/*const double tau_max=0.1;*/
/*const double global_error_limit=0.05;*/

void adaptive_ode(const double tstart, const double tend,
	const gsl_vector *xstart,
	const double tau_max,const double global_error_limit,
	gsl_vector *t, gsl_vector *tau, gsl_vector *global_error_norm,
	gsl_matrix *x, gsl_matrix *global_error, size_t regime,
	double *gparameters, size_t n_gparam, void (*g)(double,size_t,const gsl_vector *,double *, gsl_vector *)
	);
void function_F(size_t t, const gsl_vector *x,double *gparameters, gsl_vector *F);
void debug_adaptive_ode();
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
	gsl_vector *x_tend);
void debug_adaptive_ode_kf();
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
	gsl_vector *x_tend);
#endif
