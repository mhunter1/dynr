
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/**
 * The measurement function
 */
void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht,gsl_vector *y){

    gsl_matrix_set(Ht,0,0,1.0);
    gsl_matrix_set(Ht,1,1,1.0);


    gsl_vector_set(y, 0, gsl_vector_get(eta, 0));
    gsl_vector_set(y, 1, gsl_vector_get(eta, 1));

}

/**
 * The dx/dt function
 */

void function_dx_dt(double t, size_t regime, const gsl_vector *x,double *param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){

    double r1=0, r2=0, a12=0, a21=0;

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


/**
 * The dF/dx function
 * The partial derivative of the jacobian of the dynamic function with respect to the variable x
 * @param param includes at the end the current state estimates in the same order as the states following the model parameters
 */

void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){


    /*Supply the Jacobian matrix for the ODEs
      ODE functions go down the rows; latent states go across columns*/
    double r1=0, r2=0, a12=0, a21=0;
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
void function_transform(double *param){
    size_t i;
    for (i=0;i<4;i++){
      param[i]=exp(param[i]);
    }
}


/**
 * The dP/dt function: depend on function_dF_dx, No need to change
 */
void mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec){
    size_t i,j;
    size_t nx=mat->size1;
    /*convert matrix to vector*/
    for(i=0; i<nx; i++){
        gsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));
        for (j=i+1;j<nx;j++){
            gsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));
            /*printf("%lu",i+j+nx-1);}*/
        }
    }
}

void mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat){
    size_t i,j;
    size_t nx=mat->size1;
    /*convert vector to matrix*/
    for(i=0; i<nx; i++){
        gsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));
        for (j=i+1;j<nx;j++){
            gsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));
            gsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));
        }
    }
}

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



