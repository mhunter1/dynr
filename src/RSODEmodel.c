/**
 * This file implements the Van der pol Oscillator Model
 * ny=3, nf=3;
 There is no process noise in this ODE model
 search for eta_noise_cov
 they should be commmented out

 dx(t)=F(t,x(t))dt+G(t)dw(t)

 dx_t/dt=d[x_i1_t x_i2_t b_izeta]'/dt
 =F(t,x)=[x_i2_t | -rho*x_i1_t+zeta_i*(1-x_i1_t^2)*x_i2_t | 0 ], rho=(2*Pi/tau)^2
 =[x_i2_t | -rho*x_i1_t+zeta_i*(1-x_i1_t^2)*x_i2_t | 0 ], zeta_i=Z_0+gamma_1*Z_1+gamma_2*Z_2+b_izeta
 =[x_i2_t | -rho*x_i1_t+(Z_0+gamma_1*Z_1+gamma_2*Z_2)*x_i2_t+b_izeta*x_i2_t-(Z_0+gamma_1*Z_1+gamma_2*Z_2)*x_i2_t*x_i1_t^2-b_izeta*x_i2_t*x_i1_t^2|0]

 dF/dx=[
 0 1 0;
 -rho-2*(Z_0+gamma_1*Z_1+gamma_2*Z_2+b_izeta)*x_i2_t*x_i1_t |(Z_0+gamma_1*Z_1+gamma_2*Z_2+b_izeta)*(1-x_i1_t^2) |(1-x_i1_t^2)*x_i2_t;
 0 0 0]

 z_k=mu+H_kx_k+v_k,v_k \sim N(0,R_k)

 z_k=[y_i1_tij | y_i2_tij | y_i3_tij]'
 =[mu_1 mu_2 mu_3]'+ [1 lambda_21 lambda_31; 0 0 0; 0 0 0]'[x_i1_t x_i2_t b_izeta]'+[e_i1_tij | e_i2_tij | e_i3_tij]'
 =mu+H_kx_k+v_k

 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "data_structure.h"
#include "math_function.h"
#include "brekfis.h"
#include "adaodesolver.h"
#include <gsl/gsl_matrix.h>

/**
 * The measurement function
 * @param param mu1 mu2 mu3 lambda21 lambda31
 */
void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht,gsl_vector *y){
    /*printf("t=%d\nparam: ", t);
    print_array(param, 6);
    printf("\n");
    printf("eta: ");
    print_vector(eta);
    printf("\n");*/


    gsl_matrix_set(Ht,0,0,1.0);
    gsl_matrix_set(Ht,1,1,1.0);


    gsl_vector_set(y, 0, gsl_vector_get(eta, 0));
    gsl_vector_set(y, 1, gsl_vector_get(eta, 1));


    /*printf("y: ");
    print_vector(y);
    printf("\n\n");*/
}

/**
 * The dx/dt function
 */
/**
 * jacobian for dynamic function
 * @param x is state variable
 */
void function_dx_dt(double t, size_t regime, const gsl_vector *x,double *param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){

    double r1, r2, a12, a21;

    /* Specify the ODEs for the RS ODE model */
    /*gparameters[]={-.8,-1,-0.6931,2.3026,2.3026,-4.0000,7.0000,-0.5000,-0.5000,90.0000};*/

    switch (regime) {
        case 0:
            r1 = param[0];
            r2 = param[1];
            a12 = 0.0;
            a21 = 0.0;
            break;
        case 1:
            r1 = 0.0;
            r2 = 0.0;
            a12 = param[2];
            a21 = param[3];
            break;
    }

    /*dx[0] = -r1 * x[0] + a12 * (x[1] - x[0]);
     dx[1] = -r2 * (x[1] - par[9]) - a21 * (x[1] - x[0]);
     printf("\nTesting testing in functionF!");*/
    gsl_vector_set(F_dx_dt,0,-r1 * gsl_vector_get(x,0) + a12 * (gsl_vector_get(x,1) - gsl_vector_get(x,0)));
    gsl_vector_set(F_dx_dt,1,-r2 * (gsl_vector_get(x,1) - param[10]) - a21 * (gsl_vector_get(x,1) - gsl_vector_get(x,0)));
}
void function_dynam_ada(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *gparameters,const gsl_vector *co_variate,
        void (*g)(double, size_t, const gsl_vector *, double *, const gsl_vector *, gsl_vector *),
        gsl_vector *x_tend){

        double tau_max=(tend-tstart)/10;/*specify tau_max*/
        double global_error_limit=10;/*specify global error limit*/
        adaptive_ode_kf(tstart, tend, xstart,tau_max,global_error_limit, regime,  gparameters, co_variate, (*g),x_tend);

}
/**
 * The dF/dx function
 * The partial derivative of the jacobian of the dynamic function with respect to the variable x
 * @param param includes at the end the current state estimates in the same order as the states following the model parameters
 */

void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){


    /*Supply the Jacobian matrix for the ODEs
      ODE functions go down the rows; latent states go across columns*/
    double r1, r2, a12, a21;
    switch (regime) {
        case 0:
            r1 = param[0];
            r2 = param[1];
            a12 = 0.0;
            a21 = 0.0;
            break;
        case 1:
            r1 = 0.0;
            r2 = 0.0;
            a12 = param[2];
            a21 = param[3];
            break;
    }
    gsl_matrix_set(F_dx_dt_dx,0,0,-r1-a12);
    gsl_matrix_set(F_dx_dt_dx,0,1,a12);
    gsl_matrix_set(F_dx_dt_dx,1,0,a21);
    gsl_matrix_set(F_dx_dt_dx,1,1,-r2-a21);

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
/**
 * The dP/dt function
 */

void function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){

    size_t nx;
    nx = (size_t) floor(sqrt(2*(double) p->size));

    gsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);
    mathfunction_vec_to_mat(p,P_mat);

    gsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);

    function_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);


    gsl_matrix *dFP=gsl_matrix_calloc(nx,nx);
    gsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);

    /*printf("P_mat:\n");
    print_matrix(P_mat);
    print_matrix(F_dx_dt_dx);
    print_matrix(dFP);*/

    gsl_matrix_transpose_memcpy(dP_dt, dFP);
    gsl_matrix_add(dP_dt, dFP);

    /*if 0, add small amont so that Ppred invertible*/
    size_t i;
    for(i=0;i<nx;i++){
      gsl_matrix_set(dP_dt,i,i,gsl_matrix_get(dP_dt,i,i)+1e-4);
    }

    /*print_matrix(dP_dt);*/
    mathfunction_mat_to_vec(dP_dt, F_dP_dt);
    /*print_vector(F_dP_dt);
    exit(0);*/

    gsl_matrix_free(P_mat);
    gsl_matrix_free(F_dx_dt_dx);
    gsl_matrix_free(dFP);
    gsl_matrix_free(dP_dt);

}

/**
 * Set the initial condition
 */

void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0){

    gsl_vector_set(pr_0,0,1);

    size_t num_regime=pr_0->size;
    size_t dim_latent_var=error_cov_0[0]->size1;
    size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
    /*printf("%lu %lu %lu\n",num_regime,dim_latent_var,num_sbj);*/
    size_t i,j;
    for(j=0;j<num_regime;j++){
        for(i=0;i<num_sbj;i++){
            /*printf("%lu %lu\n",i,j);*/
            gsl_vector_set((eta_0)[j],i*dim_latent_var,70.0);
            gsl_vector_set((eta_0)[j],i*dim_latent_var+1,40.0);
        }/*statevar_1_p1 statevar_2_p1 statevar_1_p2 statevar_2_p2 ..., eta_0[] with a length of num_sbj*dim_latent_var*/

        gsl_matrix_set((error_cov_0)[j],0,0,log(225.0));
        gsl_matrix_set((error_cov_0)[j],1,1,log(100.0));
    }
}

/**
 * Set the regime-switch transition probability matrix
 */

void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){

    double p11, p12, p21, p22;
    switch (type) {
        case 0:
            p21 = 0;
            p11 = 1;
            break;
        case 1:
            /**p11 = (exp(param[6] + param[7] + param[9]*gsl_vector_get(co_variate,1)))/(exp(0)+(exp(param[6] + param[7] + param[9]*gsl_vector_get(co_variate,1))));
            p21 = (exp(param[6] + param[8]*gsl_vector_get(co_variate,1)))/(exp(0)+(exp(param[6] + param[8]*gsl_vector_get(co_variate,1))));**//*co_variate[t](2) is time-varying*/
            
            p11 = (exp(param[6] + param[7] + param[9]*gsl_vector_get(co_variate,1)+param[12]*gsl_vector_get(co_variate,0)))/(exp(0)+(exp(param[6] + param[7] + param[9]*gsl_vector_get(co_variate,1)+param[12]*gsl_vector_get(co_variate,0))));
            p21 = (exp(param[6] + param[8]*gsl_vector_get(co_variate,1)+param[11]*gsl_vector_get(co_variate,0)))/(exp(0)+(exp(param[6] + param[8]*gsl_vector_get(co_variate,1)+param[11]*gsl_vector_get(co_variate,0))));/*co_variate[t](2) is time-varying*/

            break;
    }
    p12 = 1-p11;
    p22 = 1-p21;

    gsl_matrix_set(regime_switch_mat,0,0,p11);
    gsl_matrix_set(regime_switch_mat,0,1,p12);
    gsl_matrix_set(regime_switch_mat,1,0,p21);
    gsl_matrix_set(regime_switch_mat,1,1,p22);

}

/**
 * Set the noise covariance matrix
 *They are to be LDL' transformed.
 *e.g., [a b
         b c]
 *-->LDL', L=[1 0;b 1], D=diag(a,c)
 */

void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){
    size_t i;
    for (i=0;i<eta_noise_cov->size1;i++){
        gsl_matrix_set(eta_noise_cov,i,i,-20);
    }
    /*gsl_matrix_set_zero(par.y_noise_cov);*/
    gsl_matrix_set(y_noise_cov,0,0, param[4]);
    gsl_matrix_set(y_noise_cov,1,1, param[5]);

}

/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 */
void function_transform(const ParamConfig *pc, ParamInit *pi, Param *par){
    size_t i;
    for (i=0;i<4;i++){
      par->func_param[i]=exp(par->func_param[i]);
    }
}



