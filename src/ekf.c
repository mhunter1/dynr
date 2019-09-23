 /*********************************************************************
* The Discrete/Continous-Discrete Extended Kalman Filter              *
* @Author: Lu Ou.                                                     *
* @created Fall 2014                                                  *
* Purpose:                                                            *
* Usage:                                                              *
* References: A                                                       *
* File formats:                                                       *
* Restrictions:                                                       *
* Revision history:                                                   *
* Error handling:                                                     *
* Notes:                                                              *
**********************************************************************/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "math_function.h"
#include "ekf.h"
#include "adaodesolver.h"
#include <R.h>
#include <Rinternals.h>
#include "print_function.h"


/******************************************************************************
* extended kalman filter (EKF)
* *
* calculate the filtered estimates and error covariance matrix using EKF
* store the innovation vector and innovation covariance matrix for calculating log-likelihood
* *
* Parameters/Input *
* *
* **>>Parameters/Input Constants specified by user<<**
* t -- the previous time t
* t_plus_1 -- the current time t+1
* regime -- which regime is used
* eta_t -- a vector of the filtered estimates at time t
* error_cov_t -- a vector of the elements in the filtered error covariance matrix at time t: P_t[1,1],P_t[2,2],P_t[3,3],P_t[1,2],P_t[1,3],P_t[2,3]
* y_t_plus_1 -- data y_t+1
* co_variate -- covariates
* eta_noise_cov -- process noise covariance
* y_noise_cov -- measurement error covariance i.e. Rk or R
* H_t_plus_1 -- Lambda matrix in the measurement function/model, can be time-varying
* params -- array of parameters
* func_measure -- a pointer to the measurement function without measurement error
* func_dx_dt -- a pointer to a function. the derivative function, dx/dt, *dx_dt is the function *
* func_dP_dt -- a pointer to a function. the derivative function, dP/dt, *dP_dt is the function *
* *
* **>>Parameters/Input Pointers<<**
* eta_t_plus_1 -- a vector of filtered new estimates at time t+1
* error_cov_t_plus_1 -- a vector of the elements in the filtered error covariance matrix at time t: P_t[1,1],P_t[2,2],P_t[3,3],P_t[1,2],P_t[1,3],P_t[2,3]
* innov_v-- innovation vector
* innov_cov-- innovation covariance matrix
* *
* Output*
* *
* **>>Output via using pointers<<**
* eta_t_plus_1 -- a vector of filtered new estimates at time t+1
* error_cov_t_plus_1 -- a vector of the elements in the filtered error covariance matrix at time t: P_t[1,1],P_t[2,2],P_t[3,3],P_t[1,2],P_t[1,3],P_t[2,3]
* innov_v-- innovation vector
* innov_cov-- innovation covariance matrix
* *
*********************************************/

double ext_kalmanfilter(size_t t,
		size_t regime,
		gsl_vector *eta_t,
		gsl_matrix *error_cov_t,
		const gsl_vector *y_t_plus_1,
		const gsl_vector *co_variate,
		const double *y_time,
		const gsl_matrix *eta_noise_cov,
		const gsl_matrix *y_noise_cov,
		double *params,
		size_t num_func_param,
		bool isContinuousTime,
		void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *),
		void (*func_dx_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
		void (*func_dP_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
		void (*func_dF_dx)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
		void (*func_dynam)(const double, const double, size_t, const gsl_vector *,double *, size_t, const gsl_vector *,
			void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *), gsl_vector *),
		void (*func_jacob_dynam)(const double, const double, size_t, const gsl_vector *, double *, size_t,const gsl_vector *,
			void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *), gsl_matrix *),
		gsl_vector *eta_t_plus_1, gsl_matrix *error_cov_t_plus_1,
		gsl_vector *innov_v,
		gsl_matrix *inv_innov_cov,
		bool isFirstTime){
	
	int DEBUG_EKF = 0; /*0=false/no; 1=true/yes*/
	if(DEBUG_EKF){
		MYPRINT("Called ext_kalmanfilter\n");
	}
	
	double det = 0;
	size_t nx=eta_t->size;
	
	size_t i,j;
	
	/*NB ext_kalmanfilter_smoother() has innov_cov passed to it rather than allocating it*/
	
	
	/** handling missing data **/
	gsl_vector *y_non_miss_index = gsl_vector_alloc(y_t_plus_1->size);
	size_t num_non_miss = find_miss_data(y_t_plus_1, y_non_miss_index);
	gsl_vector *y_small = gsl_vector_alloc(num_non_miss);
	filter_vector(y_t_plus_1, y_non_miss_index, y_small);
	gsl_vector *innov_v_small = gsl_vector_alloc(num_non_miss);
	
	/* Allocate matrix pointers */
	gsl_matrix *H_t_plus_1=gsl_matrix_calloc(y_t_plus_1->size, nx);
	gsl_matrix *H_small=gsl_matrix_calloc(num_non_miss, nx);
	gsl_matrix *ph=gsl_matrix_calloc(eta_t->size, num_non_miss); /* P*H' - error_cov*jacob'*/
	gsl_matrix *innov_cov=gsl_matrix_calloc(num_non_miss, num_non_miss);
	gsl_matrix *y_noise_cov_small=gsl_matrix_calloc(num_non_miss, num_non_miss);
	gsl_matrix *kalman_gain=gsl_matrix_calloc(eta_t->size, num_non_miss);
	gsl_matrix *inv_innov_cov_small=gsl_matrix_calloc(num_non_miss, num_non_miss);
	
	filter_matrix_rows_cols(y_noise_cov, y_non_miss_index, y_noise_cov_small);
	
	
	/*------------------------------------------------------*\
	* update xk *
	\*------------------------------------------------------*/
	if(DEBUG_EKF){
		MYPRINT("y_time:\n");
		MYPRINT("%f ",y_time[t-1]);
		MYPRINT("%f",y_time[t]);
		MYPRINT("\n");
		MYPRINT("regime: %lu\n",regime);
		MYPRINT("eta_previous:\n");
		print_vector(eta_t);
		MYPRINT("\n");
		MYPRINT("parameters:\n");
		print_array(params,num_func_param);
		MYPRINT("\n");
	}
	
	/*NB ext_kalmanfilter_smoother() also memcpy's eta_t_plus_1
	  Add if(isFirstTime || isSmoother ){} */
	
	if(isFirstTime){
		gsl_vector_memcpy(eta_t_plus_1,eta_t);
	} else{
		func_dynam(y_time[t-1], y_time[t], regime, eta_t, params, num_func_param, co_variate, func_dx_dt, eta_t_plus_1); /** y_time - observed time**/
	}
	
	/*MYPRINT("eta_pred:\n");
	print_vector(eta_t_plus_1);
	MYPRINT("\n");*/
	
	
	/*if(t_plus_1==0 || t_plus_1==1){
		MYPRINT("eta(%d):",t_plus_1-1);
		print_vector(eta_t);
		MYPRINT("\n");
	}
	if(t_plus_1==0 || t_plus_1==1){
		MYPRINT("eta(%d):", t_plus_1);
		print_vector(eta_t_plus_1);
		MYPRINT("\n");
	}*/
	
	
	/*------------------------------------------------------*\
	* update P *
	\*------------------------------------------------------*/
	if (!isFirstTime & isContinuousTime){ /*NB ((!isFirstTime || isSmoother) & isContinuousTime)*/
		
		gsl_vector *Pnewvec = gsl_vector_calloc(nx*(nx+1)/2);
		gsl_vector *error_cov_t_vec = gsl_vector_calloc(nx*(nx+1)/2);
		
		/*------------------------------------------------------*\
		* update P *
		\*------------------------------------------------------*/
		/*print_matrix(error_cov_t);
		MYPRINT("\n");
		exit(0);*/
		/* error_cov_t_vec = vech(error_cov_t) */
		for(i=0; i<nx; i++){
			gsl_vector_set(error_cov_t_vec, i, gsl_matrix_get(error_cov_t, i, i));
			for (j=i+1; j<nx; j++){
				gsl_vector_set(error_cov_t_vec, i+j+nx-1, gsl_matrix_get(error_cov_t, i, j));
				/*MYPRINT("%lu",i+j+nx-1);}*/
			}
		}
		/*print_vector(error_cov_t_vec);
		MYPRINT("\n");*/
		
		/*dpparams include params, eta_t, and eta_noise_cov_vec*/
		size_t n_dpparams=num_func_param+nx+(nx+1)*nx/2;
		double dpparams[n_dpparams];
		for (i=0; i<num_func_param; i++)
			dpparams[i]=params[i];
		for (i=0; i<nx; i++)
			dpparams[num_func_param+i]=gsl_vector_get(eta_t, i);
		for (i=0; i<nx; i++){
			dpparams[num_func_param+nx+i]=gsl_matrix_get(eta_noise_cov, i, i);
			for (j=i+1; j<nx; j++){
				dpparams[num_func_param+nx+i+j+nx-1] = gsl_matrix_get(eta_noise_cov, i, j);
			}
		}
		
		/*MYPRINT("y_time:\n");
		MYPRINT("%f ",y_time[t-1]);
		MYPRINT("%f",y_time[t]);
		MYPRINT("\n");
		MYPRINT("regime: %lu\n",regime);
		MYPRINT("error_cov_previous:\n");
		print_vector(error_cov_t_vec);
		MYPRINT("\n");
		MYPRINT("parameters:\n");
		print_array(dpparams,num_func_param+nx+(nx+1)*nx/2);
		MYPRINT("\n");*/
		
		
		func_dynam(y_time[t-1], y_time[t], regime, error_cov_t_vec, dpparams, n_dpparams, co_variate, func_dP_dt, Pnewvec);
		
		/*MYPRINT("error_cov_pred:\n");
		print_vector(Pnewvec);
		MYPRINT("\n");*/
		
		for(i=0; i<nx; i++){
			gsl_matrix_set(error_cov_t_plus_1, i, i, gsl_vector_get(Pnewvec, i));
			for (j=i+1; j<nx; j++){
				gsl_matrix_set(error_cov_t_plus_1, i, j, gsl_vector_get(Pnewvec, i+j+nx-1));
				gsl_matrix_set(error_cov_t_plus_1, j, i, gsl_vector_get(Pnewvec, i+j+nx-1));
			}
		}
		
		/*gsl_matrix_set(error_cov_t_plus_1,0,0,gsl_vector_get(Pnewvec,0));
		gsl_matrix_set(error_cov_t_plus_1,0,1,gsl_vector_get(Pnewvec,3));
		gsl_matrix_set(error_cov_t_plus_1,0,2,gsl_vector_get(Pnewvec,4));
		gsl_matrix_set(error_cov_t_plus_1,1,0,gsl_vector_get(Pnewvec,3));
		gsl_matrix_set(error_cov_t_plus_1,1,1,gsl_vector_get(Pnewvec,1));
		gsl_matrix_set(error_cov_t_plus_1,1,2,gsl_vector_get(Pnewvec,5));
		gsl_matrix_set(error_cov_t_plus_1,2,0,gsl_vector_get(Pnewvec,4));
		gsl_matrix_set(error_cov_t_plus_1,2,1,gsl_vector_get(Pnewvec,5));
		gsl_matrix_set(error_cov_t_plus_1,2,2,gsl_vector_get(Pnewvec,2));*/
		
		gsl_vector_free(error_cov_t_vec);
		gsl_vector_free(Pnewvec);
		
	} else if(!isFirstTime & !isContinuousTime) { /* end continuous time "Update P" calculation */
		/*NB if((!isFirstTime || isSmoother) & !isContinuousTime)*/
		
		/*Update P for discrete time*/
		gsl_matrix *jacob_dynam = gsl_matrix_calloc(nx,nx);
		gsl_matrix *p_jacob_dynam = gsl_matrix_calloc(nx, nx);
		/*------------------------------------------------------*\
		* Update P discrete--------error_cov_t_plus_1=eta_noise_cov+jacobdynamic%*%error_cov_t%*%t(jacobdynamic)*
		\*------------------------------------------------------*/
		// error_cov_t_plus_1 = Predicted P for latent variables
		
		func_jacob_dynam(y_time[t-1], y_time[t], regime, eta_t, params, num_func_param, co_variate, func_dF_dx,jacob_dynam);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, error_cov_t, jacob_dynam, 0.0, p_jacob_dynam); /* compute P*jacobdynamic'*/
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, jacob_dynam, p_jacob_dynam, 0.0, error_cov_t_plus_1); /* compute jacobdynamic*P*jacobdynamic'*/
		gsl_matrix_add(error_cov_t_plus_1, eta_noise_cov); /*compute H*P*H'+Q*/
		
		gsl_matrix_free(jacob_dynam);
		gsl_matrix_free(p_jacob_dynam);
	} else if(isFirstTime){
		gsl_matrix_memcpy(error_cov_t_plus_1, error_cov_t);
	}
	/* End "Update P" */
	
	/*------------------------------------------------------*\
	* innovation vector-- will be used to calculate loglikelihood*
	\*------------------------------------------------------*/
	
	/** step 2.1: compute measurement y_hat(t+1|t) **/
	
	/*NB following line is in smoother but not here*/
	//    gsl_matrix_memcpy(error_cov_pred,error_cov_t_plus_1);
	
	if(DEBUG_EKF){
		MYPRINT("eta(%d):",t);
		print_vector(eta_t_plus_1);
		MYPRINT("\n");
	}
	func_measure(t, regime, params, eta_t_plus_1, co_variate, H_t_plus_1, innov_v);/*innov_v <- y_pred*/
	if(DEBUG_EKF){
		MYPRINT("y_hat(%d):", t);
		print_vector(innov_v);
		MYPRINT("\n");
	}
	filter_matrix_rows(H_t_plus_1, y_non_miss_index, H_small);
	filter_vector(innov_v, y_non_miss_index, innov_v_small);
	
	/*------------------------------------------------------*\
	// Predicted Covariance matrix for observed variables
	* innovation variance--------Rek=Rk+H_t_plus_1%*%Pnew%*%t(H_t_plus_1) -- will be used to calculate loglikelihood*
	\*------------------------------------------------------*/
	if(num_non_miss > 0){
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, error_cov_t_plus_1, H_small, 0.0, ph); /* compute P*H'*/
		/*mathfunction_matrix_mul(error_cov_t_plus_1, jacob_measure, false, true, ph);*/
	}
	
	if(DEBUG_EKF){
		MYPRINT("error_cov(%lu):\n", t);
		print_matrix(error_cov_t_plus_1);
		MYPRINT("\n");
	}
	if(DEBUG_EKF){
		MYPRINT("Hk:\n");
		print_matrix(H_t_plus_1);
		MYPRINT("\n");
		MYPRINT("Hk_filtered:\n");
		print_matrix(H_small);
		MYPRINT("\n");
		MYPRINT("ph:\n");
		print_matrix(ph);
		MYPRINT("\n");
	}
	
	if(num_non_miss > 0){
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H_small, ph, 0.0, innov_cov); /* compute H*P*H'*/
	}
	if(DEBUG_EKF){
		/*MYPRINT("jacob_measure:\n");
		print_matrix(jacob_measure);
		MYPRINT("\n");*/
		MYPRINT("ph:\n");
		print_matrix(ph);
		MYPRINT("\n");
	}
	/*mathfunction_matrix_mul(jacob_measure, ph, false, false, innov_cov);*/
	
	if(DEBUG_EKF){
		MYPRINT("H*P(%d)*H':\n", t);
		print_matrix(innov_cov);
		MYPRINT("R(%d):\n", t);
		print_matrix(y_noise_cov);
		MYPRINT("\n");
		MYPRINT("Filtered R(%d):\n", t);
		print_matrix(y_noise_cov_small);
		MYPRINT("\n");
	}
	
	if(num_non_miss > 0){
		gsl_matrix_add(innov_cov, y_noise_cov_small); /*compute H*P*H'+R*/
	}
	if(DEBUG_EKF){
		MYPRINT("y_cov(%d):\n", t);
		print_matrix(innov_cov);
		MYPRINT("\n");
	}
	
	
	/*------------------------------------------------------*\
	* Kalman Gain------Kk=Pnew%*%t(H_t_plus_1)%*%solve(Rek) *
	\*------------------------------------------------------*/
	
	
	if(num_non_miss > 0){
		det = mathfunction_inv_matrix_det(innov_cov, inv_innov_cov_small);
	}
	
	if(DEBUG_EKF){
		MYPRINT("inv_innov_cov filtered:\n");
		print_matrix(inv_innov_cov_small);
		MYPRINT("\n");
	}
	
	if(num_non_miss > 0){
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ph, inv_innov_cov_small, 0.0, kalman_gain); /* compute P*H*S^{-1}*/
	}
	
	if(DEBUG_EKF){
		MYPRINT("ph:\n");
		print_matrix(ph);
		MYPRINT("\n");
		MYPRINT("inverse residu:\n");
		print_matrix(inv_innov_cov_small);
		MYPRINT("\n");
		MYPRINT("kalman_gain:\n");
		print_matrix(kalman_gain);
		MYPRINT("\n");
	}
	
	/*------------------------------------------------------*\
	* Fitered Estimate *
	\*------------------------------------------------------*/
	
	
	/** step 2.2: compute prediction residual y-y_hat **/
	if(num_non_miss > 0){
		gsl_vector_sub(innov_v_small, y_small);
		gsl_vector_scale(innov_v_small, -1); /* now innov_v stores the residual, i.e, v(t+1)*/
		gsl_blas_dgemv(CblasNoTrans, 1.0, kalman_gain, innov_v_small, 1.0, eta_t_plus_1); /* x(k+1|k+1)=x(k+1|k)+W(k+1)*v(k+1)*/
	}
	
	if(DEBUG_EKF){
		MYPRINT("eta_corrected(%lf):", y_time[t]);
		print_vector(eta_t_plus_1);
		MYPRINT("\n");
	}
	
	
	/*------------------------------------------------------*\
	* Filtered Error Cov Matrix *
	\*------------------------------------------------------*/
	/*P_kplus1=Pnew-Kk%*%H_t_plus_1%*%Pnew
	P[,k] = c(P_kplus1[1,1], P_kplus1[2,2], P_kplus1[3,3], P_kplus1[1,2], P_kplus1[1,3], P_kplus1[2,3]) */
	
	
	
	if(num_non_miss > 0){
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, ph, kalman_gain, -1.0, error_cov_t_plus_1); /* W*S*W'-P=P*H'*W'-P=Pnew*H_t_plus_1'*Kk'-Pnew*/
		/*MYPRINT("error_cov(%d):\n", t_plus_1);
		print_matrix(error_cov_t_plus_1);
		MYPRINT("\n");*/
		
		gsl_matrix_scale(error_cov_t_plus_1, -1.0); /* compute P-W*S*W'*/
	}
	
	/*MYPRINT("error_cov_corrected(%lf):", y_time[t]);
	print_matrix(error_cov_t_plus_1);
	MYPRINT("\n");*/
	
	
	/*if(miss_case!=0){
		MYPRINT("innov_v(%lf)", y_time[t]);
		print_vector(innov_v);
		MYPRINT("\n");
		MYPRINT("inv_innov_cov(%lf):\n", y_time[t]);
		print_matrix(inv_innov_cov);
		MYPRINT("\n");
		MYPRINT("eta_corrected(%lf):", y_time[t]);
		print_vector(eta_t_plus_1);
		MYPRINT("\n");
		MYPRINT("error_cov_corrected(%lf):\n", y_time[t]);
	print_matrix(error_cov_t_plus_1);
		MYPRINT("\n");
	}*/
	
	if(DEBUG_EKF){
		MYPRINT("About to compute likelihood for time %lf", y_time[t]);
		MYPRINT("\nInputs 1 through 4");
		MYPRINT("1. Innovation vector");
		print_vector(innov_v_small);
		MYPRINT("2. Inverse innovation covariance");
		print_matrix(inv_innov_cov_small);
		MYPRINT("3. number of non-missing variables");
		MYPRINT("%lu", num_non_miss);
		MYPRINT("4. Determinant of innovation covariance");
		MYPRINT("%lf", det);
	}
	double neg_log_p = mathfunction_negloglike_multivariate_normal_invcov(innov_v_small, inv_innov_cov_small, num_non_miss, det);
	
	if(DEBUG_EKF){
		MYPRINT("neg_log_p(%lf): %lf", y_time[t], neg_log_p);
		MYPRINT("\n");
	}
	
	/** free allocated space **/
	gsl_matrix_free(ph);
	gsl_matrix_free(innov_cov); /*NB not freed in ext_kalmanfilter_smoother because it is passed as an argument rather than allocated */
	gsl_matrix_free(kalman_gain);
	gsl_matrix_free(H_t_plus_1);
	gsl_vector_free(y_non_miss_index);
	gsl_vector_free(y_small);
	gsl_vector_free(innov_v_small);
	gsl_matrix_free(H_small);
	gsl_matrix_free(y_noise_cov_small);
	gsl_matrix_free(inv_innov_cov_small);
	
	return neg_log_p;
}


/**
 * This method finds the missing data and creates a vector that indicates which entries are missing
 * @param y the data. missing ones are represented as NA in R
 * @param non_miss the indication. 1 for not missing and 0 for missing
 *   initialized to all 1, modified on return
 * @return the number of non-missing entries of y
 */
size_t find_miss_data(const gsl_vector *y, gsl_vector *non_miss){
	size_t sum = 0;
	gsl_vector_set_all(non_miss, 1);
	for(size_t col_index=0; col_index < y->size; col_index++){
		if(ISNA(gsl_vector_get(y, col_index))){ // if y[col_index] is missing
			gsl_vector_set(non_miss, col_index, 0);
			sum++;
		}
	}
	return(y->size - sum);
}

// Filter out elements of a vector
void filter_vector(const gsl_vector *y, const gsl_vector *y_ind, gsl_vector *ysmall){
	size_t miss_index = 0;
	size_t col_index;
	for(col_index=0; col_index < y->size; col_index++){
		if(gsl_vector_get(y_ind, col_index) == 1){
			gsl_vector_set(ysmall, miss_index, gsl_vector_get(y, col_index));
			miss_index++;
		}
	}
}

// Filter out some ROWS of a matrix, keeping all columns
void filter_matrix_rows(const gsl_matrix *X, const gsl_vector *y_ind, gsl_matrix *Xsmall){
	size_t miss_index = 0;
	size_t row_index;
	size_t col_index;
	for(row_index=0; row_index < X->size1; row_index++){
		if(gsl_vector_get(y_ind, row_index) == 1){
			for(col_index=0; col_index < X->size2; col_index++){
				gsl_matrix_set(Xsmall, miss_index, col_index, gsl_matrix_get(X, row_index, col_index));
			}
			miss_index++;
		}
	}
}

// Filter out some COLUMNS of a matrix, keeping all rows
void filter_matrix_cols(const gsl_matrix *X, const gsl_vector *y_ind, gsl_matrix *Xsmall){
	size_t miss_index = 0;
	size_t row_index;
	size_t col_index;
	for(col_index=0; col_index < X->size2; col_index++){
		if(gsl_vector_get(y_ind, col_index) == 1){
			for(row_index=0; row_index < X->size1; row_index++){
				gsl_matrix_set(Xsmall, row_index, miss_index, gsl_matrix_get(X, row_index, col_index));
			}
			miss_index++;
		}
	}
}

// Filter both ROWS and COLS of a matrix
void filter_matrix_rows_cols(const gsl_matrix *X, const gsl_vector *y_ind, gsl_matrix *Xsmall){
	size_t miss_index1 = 0;
	size_t miss_index2 = 0;
	size_t row_index;
	size_t col_index;
	for(row_index=0; row_index < X->size1; row_index++){
		if(gsl_vector_get(y_ind, row_index) == 1){
			for(col_index=0; col_index < X->size2; col_index++){
				if(gsl_vector_get(y_ind, col_index) == 1){
					gsl_matrix_set(Xsmall, miss_index1, miss_index2, gsl_matrix_get(X, row_index, col_index));
					miss_index2++;
				}
			}
			miss_index2=0;
			miss_index1++;
		}
	}
}


double ext_kalmanfilter_smoother(size_t t, size_t regime,
        gsl_vector *eta_t,  gsl_matrix *error_cov_t,
		const gsl_vector *y_t_plus_1,const gsl_vector *co_variate, const double *y_time,
		const gsl_matrix *eta_noise_cov, const gsl_matrix *y_noise_cov,
        double *params,size_t num_func_param,
		bool isContinuousTime,
        void (*func_measure)(size_t, size_t, double *, const gsl_vector *, const gsl_vector *, gsl_matrix *, gsl_vector *),
        void (*func_dx_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
        void (*func_dP_dt)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
		void (*func_dF_dx)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
        void (*func_dynam)(const double, const double, size_t, const gsl_vector *,double *, size_t, const gsl_vector *, void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),gsl_vector *),
	    void (*func_jacob_dynam)(const double, const double, size_t, const gsl_vector *,
			double *, size_t,const gsl_vector *,
	        void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
			gsl_matrix *),
        gsl_vector *eta_pred, gsl_matrix *error_cov_pred, gsl_vector *eta_t_plus_1, gsl_matrix *error_cov_t_plus_1, 
		gsl_vector *innov_v, gsl_matrix *inv_innov_cov, gsl_matrix *innov_cov_delete_me, bool isFirstTime){

	int DEBUG_EKF = 0; /*0=false/no; 1=true/yes*/
	if(DEBUG_EKF){
		MYPRINT("Called ext_kalmanfilter_smoother\n");
	}

    double det = 0;
    size_t nx=eta_t->size;

    size_t i,j;

    /** handling missing data **/
	gsl_vector *y_non_miss_index = gsl_vector_alloc(y_t_plus_1->size);
	size_t num_non_miss = find_miss_data(y_t_plus_1, y_non_miss_index);
	gsl_vector *y_small = gsl_vector_alloc(num_non_miss);
	filter_vector(y_t_plus_1, y_non_miss_index, y_small);
	gsl_vector *innov_v_small = gsl_vector_alloc(num_non_miss);

    gsl_matrix *H_t_plus_1=gsl_matrix_calloc(y_t_plus_1->size, nx);
    gsl_matrix *H_small=gsl_matrix_calloc(num_non_miss, nx);
    gsl_matrix *ph=gsl_matrix_calloc(eta_t->size, num_non_miss); /* P*H' - error_cov*jacob'*/
	gsl_matrix *innov_cov=gsl_matrix_calloc(num_non_miss, num_non_miss);
	gsl_matrix *y_noise_cov_small=gsl_matrix_calloc(num_non_miss, num_non_miss);
    gsl_matrix *kalman_gain=gsl_matrix_calloc(eta_t->size, num_non_miss);
	gsl_matrix *inv_innov_cov_small=gsl_matrix_calloc(num_non_miss, num_non_miss);

	filter_matrix_rows_cols(y_noise_cov, y_non_miss_index, y_noise_cov_small);
	
	

    /*------------------------------------------------------*\
    * update xk *
    \*------------------------------------------------------*/
	if(DEBUG_EKF){
		MYPRINT("y_time:\n");
		MYPRINT("%f ",y_time[t-1]);
		MYPRINT("%f",y_time[t]);
		MYPRINT("\n");
		MYPRINT("regime: %lu\n",regime);
		MYPRINT("eta_previous:\n");
		print_vector(eta_t);
		MYPRINT("\n");
		MYPRINT("parameters:\n");
		print_array(params,num_func_param);
		MYPRINT("\n");
	}

    if(isFirstTime){
      gsl_vector_memcpy(eta_t_plus_1,eta_t);
	} else {
	    func_dynam(y_time[t-1], y_time[t], regime, eta_t, params, num_func_param, co_variate, func_dx_dt, eta_t_plus_1);
	}

    /*MYPRINT("eta_pred:\n");
    print_vector(eta_t_plus_1);
    MYPRINT("\n");*/
    gsl_vector_memcpy(eta_pred,eta_t_plus_1);

    /*if(t_plus_1==0 || t_plus_1==1){
        MYPRINT("eta(%d):",t_plus_1-1);
        print_vector(eta_t);
        MYPRINT("\n");
    }*/
    /*if(t_plus_1==0 || t_plus_1==1){
        MYPRINT("eta(%d):", t_plus_1);
        print_vector(eta_t_plus_1);
        MYPRINT("\n");
    }*/
	/*------------------------------------------------------*\
	* update P *
	\*------------------------------------------------------*/
	if (!isFirstTime & isContinuousTime){
		
	    gsl_vector *Pnewvec=gsl_vector_calloc(nx*(nx+1)/2);
	    gsl_vector *error_cov_t_vec=gsl_vector_calloc(nx*(nx+1)/2);
		

		/*print_matrix(error_cov_t);
		MYPRINT("\n");
		exit(0);*/
		for(i=0; i<nx; i++){
			gsl_vector_set(error_cov_t_vec,i,gsl_matrix_get(error_cov_t,i,i));
			for (j=i+1;j<nx;j++){
				gsl_vector_set(error_cov_t_vec,i+j+nx-1,gsl_matrix_get(error_cov_t,i,j));
				/*MYPRINT("%lu",i+j+nx-1);}*/
			}
		}
		/*print_vector(error_cov_t_vec);
		MYPRINT("\n");*/

		/*dpparams include params, eta_t, and eta_noise_cov_vec*/
		size_t n_dpparams=num_func_param+nx+(nx+1)*nx/2;
		double dpparams[n_dpparams];
		for (i=0;i<num_func_param;i++)
			dpparams[i]=params[i];
		for (i=0;i<nx;i++)
			dpparams[num_func_param+i]=gsl_vector_get(eta_t,i);
		for (i=0; i<nx; i++){
			dpparams[num_func_param+nx+i]=gsl_matrix_get(eta_noise_cov,i,i);
			for (j=i+1;j<nx;j++){
				dpparams[num_func_param+nx+i+j+nx-1]=gsl_matrix_get(eta_noise_cov,i,j);
			}
		}

		/*MYPRINT("y_time:\n");
		MYPRINT("%f ",y_time[t-1]);
		MYPRINT("%f",y_time[t]);
		MYPRINT("\n");
		MYPRINT("regime: %lu\n",regime);
		MYPRINT("error_cov_previous:\n");
		print_vector(error_cov_t_vec);
		MYPRINT("\n");
		MYPRINT("parameters:\n");
		print_array(dpparams,num_func_param+nx+(nx+1)*nx/2);
		MYPRINT("\n");*/


		func_dynam(y_time[t-1], y_time[t], regime, error_cov_t_vec, dpparams, n_dpparams, co_variate, func_dP_dt, Pnewvec);

		/*MYPRINT("error_cov_pred:\n");
		print_vector(Pnewvec);
		MYPRINT("\n");*/

		for(i=0; i<nx; i++){
			gsl_matrix_set(error_cov_t_plus_1,i,i,gsl_vector_get(Pnewvec,i));
			for (j=i+1;j<nx;j++){
				gsl_matrix_set(error_cov_t_plus_1,i,j,gsl_vector_get(Pnewvec,i+j+nx-1));
				gsl_matrix_set(error_cov_t_plus_1,j,i,gsl_vector_get(Pnewvec,i+j+nx-1));
			}
		}


		/*gsl_matrix_set(error_cov_t_plus_1,0,0,gsl_vector_get(Pnewvec,0));
		gsl_matrix_set(error_cov_t_plus_1,0,1,gsl_vector_get(Pnewvec,3));
		gsl_matrix_set(error_cov_t_plus_1,0,2,gsl_vector_get(Pnewvec,4));
		gsl_matrix_set(error_cov_t_plus_1,1,0,gsl_vector_get(Pnewvec,3));
		gsl_matrix_set(error_cov_t_plus_1,1,1,gsl_vector_get(Pnewvec,1));
		gsl_matrix_set(error_cov_t_plus_1,1,2,gsl_vector_get(Pnewvec,5));
		gsl_matrix_set(error_cov_t_plus_1,2,0,gsl_vector_get(Pnewvec,4));
		gsl_matrix_set(error_cov_t_plus_1,2,1,gsl_vector_get(Pnewvec,5));
		gsl_matrix_set(error_cov_t_plus_1,2,2,gsl_vector_get(Pnewvec,2));*/
		
	    gsl_vector_free(error_cov_t_vec);
	    gsl_vector_free(Pnewvec);
		
	}else if(!isFirstTime & !isContinuousTime){
		
	    gsl_matrix *jacob_dynam=gsl_matrix_calloc(nx,nx);
	    gsl_matrix *p_jacob_dynam=gsl_matrix_calloc(nx, nx);
	    /*------------------------------------------------------*\
	    * Update P discrete--------error_cov_t_plus_1=eta_noise_cov+jacobdynamic%*%error_cov_t%*%t(jacobdynamic)*
	    \*------------------------------------------------------*/
		
		func_jacob_dynam(y_time[t-1], y_time[t], regime, eta_t, params, num_func_param, co_variate, func_dF_dx,jacob_dynam);
	    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, error_cov_t, jacob_dynam, 0.0, p_jacob_dynam); /* compute P*jacobdynamic'*/
	    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, jacob_dynam, p_jacob_dynam, 0.0, error_cov_t_plus_1); /* compute jacobdynamic*P*jacobdynamic'*/
	    gsl_matrix_add(error_cov_t_plus_1, eta_noise_cov); /*compute H*P*H'+Q*/
		
		gsl_matrix_free(jacob_dynam);
		gsl_matrix_free(p_jacob_dynam);
	} else if(isFirstTime){
      gsl_matrix_memcpy(error_cov_t_plus_1,error_cov_t);
	}


    gsl_matrix_memcpy(error_cov_pred,error_cov_t_plus_1);

    /*gsl_matrix_set(error_cov_t_plus_1,0,0,gsl_vector_get(Pnewvec,0));
    gsl_matrix_set(error_cov_t_plus_1,0,1,gsl_vector_get(Pnewvec,3));
    gsl_matrix_set(error_cov_t_plus_1,0,2,gsl_vector_get(Pnewvec,4));
    gsl_matrix_set(error_cov_t_plus_1,1,0,gsl_vector_get(Pnewvec,3));
    gsl_matrix_set(error_cov_t_plus_1,1,1,gsl_vector_get(Pnewvec,1));
    gsl_matrix_set(error_cov_t_plus_1,1,2,gsl_vector_get(Pnewvec,5));
    gsl_matrix_set(error_cov_t_plus_1,2,0,gsl_vector_get(Pnewvec,4));
    gsl_matrix_set(error_cov_t_plus_1,2,1,gsl_vector_get(Pnewvec,5));
    gsl_matrix_set(error_cov_t_plus_1,2,2,gsl_vector_get(Pnewvec,2));*/

    /*------------------------------------------------------*\
    * innovation vector-- will be used to calculate loglikelihood*
    \*------------------------------------------------------*/

    /** step 2.1: compute measurement y_hat(t+1|t) **/
    
        /*MYPRINT("eta(%d):",t_plus_1);
        print_vector(eta_t_plus_1);
        MYPRINT("\n");*/
	if(DEBUG_EKF){
		MYPRINT("About to call measurement function\n");
	}
        func_measure(t, regime, params, eta_t_plus_1, co_variate, H_t_plus_1, innov_v);
	if(DEBUG_EKF){
		MYPRINT("About to filter measurement\n");
	}
	filter_matrix_rows(H_t_plus_1, y_non_miss_index, H_small);
	filter_vector(innov_v, y_non_miss_index, innov_v_small);
	if(DEBUG_EKF){
		MYPRINT("Just filtered measurement for missing data\n");
	}
	
    /*MYPRINT("y_hat(%d):", t_plus_1);
    print_vector(innov_v);
    MYPRINT("\n");*/

    /*------------------------------------------------------*\
    * innovation variance--------Rek=Rk+H_t_plus_1%*%Pnew%*%t(H_t_plus_1) -- will be used to calculate loglikelihood*
    \*------------------------------------------------------*/

	if(DEBUG_EKF){
		/*MYPRINT("Hk:\n");
		print_matrix(H_t_plus_1);
		MYPRINT("\n");
		MYPRINT("Hk_filtered:\n");
		print_matrix(H_small);
		MYPRINT("\n");
		MYPRINT("ph:\n");
		print_matrix(ph);
		MYPRINT("\n");
		*/
		MYPRINT("Gonna multiply me some P and H\n");
	}

	if(num_non_miss > 0){
	    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, error_cov_t_plus_1, H_small, 0.0, ph); /* compute P*H'*/
	}
    /*mathfunction_matrix_mul(error_cov_t_plus_1, jacob_measure, false, true, ph);*/

	if(DEBUG_EKF){
		/*MYPRINT("jacob_measure:\n");
		print_matrix(jacob_measure);
		MYPRINT("\n");*/
		/*MYPRINT("ph:\n");
		print_matrix(ph);
		MYPRINT("\n");
		*/
		MYPRINT("Just multiplied P and H to make PH\n");
	}

    /*if(t==0){
        MYPRINT("error_cov(%lu):\n", t);
        print_matrix(error_cov_t_plus_1);
        MYPRINT("\n");
        MYPRINT("Hk:\n");
        print_matrix(H_t_plus_1);
        MYPRINT("\n");
        MYPRINT("ph:\n");
        print_matrix(ph);
        MYPRINT("\n");
    }*/

	if(num_non_miss > 0){
	    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H_small, ph, 0.0, innov_cov); /* compute H*P*H'*/
	}
    /*MYPRINT("jacob_measure:\n");
    print_matrix(jacob_measure);
    MYPRINT("\n");
    MYPRINT("ph:\n");
    print_matrix(ph);
    MYPRINT("\n");*/
    /*mathfunction_matrix_mul(jacob_measure, ph, false, false, innov_cov);*/

    /*MYPRINT("y_cov(%d):\n", t_plus_1);
    print_matrix(innov_cov);
    MYPRINT("\n");*/


	if(num_non_miss > 0){
    	gsl_matrix_add(innov_cov, y_noise_cov_small); /*compute H*P*H'+R*/
	}

    /*print_matrix(y_noise_cov);
    print_matrix(innov_cov);
                 MYPRINT("\n");*/


    /*------------------------------------------------------*\
    * Kalman Gain------Kk=Pnew%*%t(H_t_plus_1)%*%solve(Rek) *
    \*------------------------------------------------------*/


	if(num_non_miss > 0){
    	det=mathfunction_inv_matrix_det(innov_cov, inv_innov_cov_small);
	}


        /*MYPRINT("inv_innov_cov:\n");
        print_matrix(inv_innov_cov);
        MYPRINT("\n");
        exit(0);*/

	if(num_non_miss > 0){
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ph, inv_innov_cov_small, 0.0, kalman_gain); /* compute P*H*S^{-1}*/
	}

        /*MYPRINT("ph:\n");
        print_matrix(ph);
        MYPRINT("\n");
        MYPRINT("inverse residu:\n");
        print_matrix(inv_innov_cov);
        MYPRINT("\n");
        MYPRINT("kalman_gain:\n");
        print_matrix(kalman_gain);
        MYPRINT("\n");*/

    /*------------------------------------------------------*\
    * Fitered Estimate *
    \*------------------------------------------------------*/



    /** step 2.2: compute prediction residual y-y_hat **/
	if(num_non_miss > 0){
	    gsl_vector_sub(innov_v_small, y_small);
	    gsl_vector_scale(innov_v_small, -1); /* now innov_v stores the residual, i.e, v(t+1)*/
	}

	if(num_non_miss > 0){
    gsl_blas_dgemv(CblasNoTrans, 1.0, kalman_gain, innov_v_small, 1.0, eta_t_plus_1); /* x(k+1|k+1)=x(k+1|k)+W(k+1)*v(k+1)*/
	}
	// else leave eta_t_plus_1 as is at its predicted value

        /*MYPRINT("eta_corrected(%lf):", y_time[t]);
        print_vector(eta_t_plus_1);
        MYPRINT("\n");*/
    /*------------------------------------------------------*\
    * Filtered Error Cov Matrix *
    \*------------------------------------------------------*/
    /*P_kplus1=Pnew-Kk%*%H_t_plus_1%*%Pnew
    P[,k]=c(P_kplus1[1,1],P_kplus1[2,2],P_kplus1[3,3],P_kplus1[1,2],P_kplus1[1,3],P_kplus1[2,3])*/



	if(num_non_miss > 0){
	    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, ph, kalman_gain, -1.0, error_cov_t_plus_1); /* W*S*W'-P=P*H'*W'-P=Pnew*H_t_plus_1'*Kk'-Pnew*/

    /*MYPRINT("error_cov(%d):\n", t_plus_1);
    print_matrix(error_cov_t_plus_1);
    MYPRINT("\n");*/

	    gsl_matrix_scale(error_cov_t_plus_1, -1.0); /* compute P-W*S*W'*/
	}
	// else leave error_cov_t_plus_1 as is at its predicted value

        /*MYPRINT("error_cov_corrected(%lf):", y_time[t]);
        print_matrix(error_cov_t_plus_1);
        MYPRINT("\n");*/

	if(DEBUG_EKF){
		MYPRINT("About to compute likelihood\n");
	}
	double neg_log_p = mathfunction_negloglike_multivariate_normal_invcov(innov_v_small, inv_innov_cov_small, num_non_miss, det);
	if(DEBUG_EKF){
		MYPRINT("Finished row likelihood computation\n");
	}

	
    /** free allocated space **/
	gsl_matrix_free(ph);
	gsl_matrix_free(kalman_gain);
	gsl_matrix_free(H_t_plus_1);
	gsl_vector_free(y_non_miss_index);
	gsl_vector_free(y_small);
	gsl_vector_free(innov_v_small);
	gsl_matrix_free(H_small);
	gsl_matrix_free(y_noise_cov_small);
	gsl_matrix_free(inv_innov_cov_small);


    return neg_log_p;
}


