
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>

/**
 * The measurement function
 */
void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht,gsl_vector *y){

    gsl_matrix_set(Ht,0,0,1.0);
    gsl_matrix_set(Ht,0,1,0.0);
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);

}

/**
 * The dx/dt function
 */

void function_dx_dt(double t, size_t regime, const gsl_vector *x, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dx_dt){

	/*if linear*/
	size_t nrow=F_dx_dt->size, ncol=x->size;
	/*size_t index_col,index_row;*/
	
	gsl_matrix *A=gsl_matrix_calloc(nrow,ncol);
	gsl_matrix_set(A,0,0,0.0);
	gsl_matrix_set(A,0,1,1.0);
	gsl_matrix_set(A,1,0,param[0]);
	gsl_matrix_set(A,1,1,param[1]);
	/*in R, use the code below to generate the above lines*/
    /*for(index_col=0; index_col<ncol; index_col++){
        for(index_row=0; index_row<nrow; index_row++){
			gsl_matrix_set(A,index_row, index_col,param[index_row*ncol+index_col]);
		}
	}*/
	gsl_blas_dgemv(CblasNoTrans, 1.0, A, x, 0.0, F_dx_dt);
	gsl_matrix_free(A);
	
	if (co_variate){
		ncol=co_variate->size;
		gsl_matrix *B=gsl_matrix_calloc(nrow,ncol);
		/*in R, use the code below to generate the above lines*/
	    /*for(index_col=0; index_col<ncol; index_col++){
	        for(index_row=0; index_row<nrow; index_row++){
				gsl_matrix_set(B,index_row, index_col,param[index_row*ncol+index_col]);
			}
		}*/
		gsl_blas_dgemv(CblasNoTrans, 1.0, B, co_variate, 1.0, F_dx_dt);
		gsl_matrix_free(B);
	}
	
    
}


/**
 * The dF/dx function
 * The partial derivative of the jacobian of the DE function with respect to the variable x
 * @param param includes at the end the current state estimates in the same order as the states following the model parameters
 */

void function_dF_dx(double t, size_t regime, double *param, const gsl_vector *co_variate, gsl_matrix *F_dx_dt_dx){


    /*Supply the Jacobian matrix for the ODEs
      ODE functions go down the rows; latent states go across columns*/
	
	gsl_matrix_set(F_dx_dt_dx,0,0,0.0);
	gsl_matrix_set(F_dx_dt_dx,0,1,1.0);
	gsl_matrix_set(F_dx_dt_dx,1,0,param[0]);
	gsl_matrix_set(F_dx_dt_dx,1,1,param[1]);
	
	/*When the model is linear*/
    /*in R, use the code below to generate the above lines*/
	/*size_t nrow=F_dx_dt_dx->size1, ncol=F_dx_dt_dx->size2,index_col,index_row;
	
    for(index_col=0; index_col<ncol; index_col++){
        for(index_row=0; index_row<nrow; index_row++){
			gsl_matrix_set(F_dx_dt_dx,index_row, index_col,param[index_row*ncol+index_col]);
		}
	}*/


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
			gsl_vector_set((eta_0)[j],i*dim_latent_var,param[4]);
			gsl_vector_set((eta_0)[j],i*dim_latent_var+1,1.0);
			
			/*in R, use code like below to produce the lines */
			/*for(index_col=0; index_col<dim_latent_var; index_col++){
					gsl_vector_set((eta_0)[j],i*dim_latent_var+index_col,param[index_col]);
			}*/
			
        }/*statevar_1_p1 statevar_2_p1 statevar_1_p2 statevar_2_p2 ..., eta_0[] with a length of num_sbj*dim_latent_var*/

		
		gsl_matrix_set((error_cov_0)[j],0, 0, log(1.0));
		gsl_matrix_set((error_cov_0)[j],1, 1, log(1.0));
		
		/*in R, use code like below to generate the above lines */
		/*if Diagonal*/
		/*
		for(index_col=0; index_col<dim_latent_var; index_col++){
				gsl_matrix_set((error_cov_0)[j],index_col, index_col,param[index_col]);
		}*/
		
		/*if Symmetric*/
		/*
		for(index_col=0; index_col<dim_latent_var; index_col++){
			gsl_matrix_set((error_cov_0)[j],index_col, index_col,param[index_col*(index_col+1)+index_col]);
	        for(index_row=0; index_row<index_col; index_row++){
				gsl_matrix_set((error_cov_0)[j],index_row, index_col,param[index_row*(index_col+1)+index_col]);
				gsl_matrix_set((error_cov_0)[j],index_col, index_row,param[index_row*(index_col+1)+index_col]);
			}
		}	
		*/

    }

}

/**
 * Set the regime-switch transition probability matrix
 */

void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){
     /*if there is only one regime*/
	gsl_matrix_set_identity(regime_switch_mat);
	
}

/**
 * Set the noise covariance matrix
 *They are to be LDL' transformed.
 *e.g., [a b
         b c]
 *-->LDL', L=[1 0;b 1], D=diag(a,c)
 */

void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){

    gsl_matrix_set(eta_noise_cov,0,0,-10.0);
	gsl_matrix_set(eta_noise_cov,1,1,param[2]);
	gsl_matrix_set(y_noise_cov,0,0, param[3]);
	/*in R use for loops to generate the lines above for symmetric and diagonal matrices*/

}

/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 * Do not include parameters in noise_cov matrices 
 */
void function_transform(double *param){
}


/**
 * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end
 * but user does not need to modify it or care about it.
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

void function_dP_dt(double t, size_t regime, const gsl_vector *p, double *param, size_t n_param, const gsl_vector *co_variate, gsl_vector *F_dP_dt){
    
    size_t nx; 
    nx = (size_t) floor(sqrt(2*(double) p->size));
    
    gsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);
    mathfunction_vec_to_mat(p,P_mat);
    
    gsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);
    
    function_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);
    
    
    gsl_matrix *dFP=gsl_matrix_calloc(nx,nx);
    gsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);
    

    gsl_matrix_transpose_memcpy(dP_dt, dFP);
    gsl_matrix_add(dP_dt, dFP);
    
	size_t n_Q_vec=(1+nx)*nx/2;
	gsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);
	size_t i;
    for(i=1;i<=n_Q_vec;i++){
        gsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);
    }

	gsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);
	mathfunction_vec_to_mat(Q_vec,Q_mat);

	gsl_matrix_add(dP_dt, Q_mat);

    mathfunction_mat_to_vec(dP_dt, F_dP_dt);

    
    gsl_matrix_free(P_mat);
    gsl_matrix_free(F_dx_dt_dx);
    gsl_matrix_free(dFP);
    gsl_matrix_free(dP_dt);
	gsl_vector_free(Q_vec);
	gsl_matrix_free(Q_mat);
    
}



