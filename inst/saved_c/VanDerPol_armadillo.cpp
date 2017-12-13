#include <math.h>
#include <armadillo>
//#include <armadillo/include/armadillo_bits/mul_gemm.hpp>

using namespace arma ;
#define omega 61.68503

//#include <armadillo/include/armadillo_bits/mul_gemm.hpp>


/*Tested: mathfunction_logistic, mathfunction_softmax, function_dx_dt, function_dF_dx, function_noise_cov*/
/*Converted done: done (compiled successfully)*/

/**
 * This function takes a double and gives back a double
 * It computes the logistic function (i.e. the inverse of the logit link function)
 * @param x, the double value e.g. a normally distributed number
 * @return logistic(x), the double value e.g. a number between 0 and 1
 */
double mathfunction_logistic(const double x){
	double value = 1.0/(1.0 + exp(-x));
	return value;
}

/**
 * This function takes a gsl_vector and modifies its second argument (another gsl_vector)
 * It computes the softmax function (e.g. for multinomial logistic regression)
 * @param x, vector of double values e.g. a vector of normally distributed numbers
 * @param result, softmax(x), e.g. a vector of numbers between 0 and 1 that sum to 1
 */
void mathfunction_softmax(const vec *x, vec *result){
	double scale;
    
    /* Elementwise exponentiation */    
	*result = exp(*x); 
	
	/* Sum for the scaling coeficient */
    scale = accu(*result);
	
	/* Multiply all elements of result by 1/scale */
    *result = 1/scale * (*result);
	
}


void function_dx_dt(double t, size_t regime, const vec *x, double *param, size_t n_param, const vec *co_variate, vec *F_dx_dt){
    (*F_dx_dt)(0) = (*x)(1);
    (*F_dx_dt)(1) = -omega * (*x)(0) + param[0] * (1 - pow((*x)(0), 2)) * (*x)(1);
}



/**
* The dF/dx function
* The partial derivative of the jacobian of the DE function with respect to the variable x
* @param param includes at the end the current state estimates in the same order as the states following the model parameters
*/
void function_dF_dx(double t, size_t regime, double *param, const vec *co_variate, mat *F_dx_dt_dx){
	(*F_dx_dt_dx)(0,0) = 0;
    (*F_dx_dt_dx)(0,1) = 1;
    (*F_dx_dt_dx)(1,0) = -(param[0] * (2 * param[9+0]) * param[9+1] + omega);
    (*F_dx_dt_dx)(1,1) = param[0] * (1 - pow(param[9+0], 2));
}



void function_measurement(size_t t, size_t regime, double *param, const vec *eta, const vec *co_variate, mat *Ht, vec *y){

    //gsl_vector *intVector = gsl_vector_calloc(3);
    vec intVector;
    intVector.set_size(3);
    
	//gsl_matrix_set(Ht, 0, 0, 1);
	(*Ht)(0, 0) = 1;
    //gsl_matrix_set(Ht, 1, 0, param[1]);
	(*Ht)(1, 0) = param[1];
    //gsl_matrix_set(Ht, 2, 0, param[2]);
    (*Ht)(2,0) = param[2];
    
	//gsl_vector_set(intVector, 0, param[3]);
    intVector(0) = param[3];
	//gsl_vector_set(intVector, 1, param[4]);
    intVector(1) = param[4];
	//gsl_vector_set(intVector, 2, param[5]);
    intVector(2) = param[5];
 
	//HJ: the function is a wrapper for GSL to call BLAS functions. Armadiilo also use BLAS, so we need to find/write the wrapper for armadillo.
    //gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);
    //dgemv_arma::apply(y, Ht, eta, 1.0, 0.0);
 
	//gsl_vector_add(y, intVector);
    (*y) = (*y) + intVector;
    
 	//gsl_vector_free(intVector);
 
}



void function_noise_cov(size_t t, size_t regime, double *param, mat *y_noise_cov, mat *eta_noise_cov){


	//gsl_matrix_set(eta_noise_cov, 0, 0, -13.8155105579643);
    (*eta_noise_cov)(0, 0) = -13.8155105579643;
	//gsl_matrix_set(eta_noise_cov, 1, 1, -13.8155105579643);
    (*eta_noise_cov)(1, 1) = -13.8155105579643;
    
	//gsl_matrix_set(y_noise_cov, 0, 0, param[6]);
    (*y_noise_cov)(0, 0) = param[6];
	//gsl_matrix_set(y_noise_cov, 1, 1, param[7]);
    (*y_noise_cov)(1, 1) = param[7];
	//gsl_matrix_sdet(y_noise_cov, 2, 2, param[8]);
    (*y_noise_cov)(2, 2) = param[8];
 
}





void function_initial_condition(double *param, vec **co_variate, vec **pr_0, vec **eta_0, mat **error_cov_0, size_t *index_sbj){
	
	//gsl_vector *Pvector = gsl_vector_calloc(1);
	//gsl_vector *Pintercept = gsl_vector_calloc(1);
	//gsl_vector *Padd = gsl_vector_calloc(1);
	//gsl_vector *Preset = gsl_vector_calloc(1);
    //gsl_vector *eta_local = gsl_vector_calloc(2);
    vec Pvector(1), Pintercept(1), Padd(1), Preset(1), eta_local(2);

    
    //gsl_vector_set(Pvector, 0, 1);
    Pvector(0) = 1;
	//gsl_vector_add(Padd, Pvector);
    Padd = Padd + Pvector;
	//gsl_vector_add(Padd, Pintercept);
    Padd = Padd + Pintercept;
	//gsl_vector_add(Preset, Pvector);
    Preset = Preset + Pvector;
	//gsl_vector_add(Preset, Pintercept);
    Preset = Preset + Pintercept;
	
	//size_t num_regime=pr_0[0]->size;
    size_t num_regime=(*pr_0[0]).n_elem;
	//size_t dim_latent_var=error_cov_0[0]->size1;
    size_t dim_latent_var=(*error_cov_0[0]).n_elem;
	//size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
    size_t num_sbj=((**eta_0).n_elem)/(dim_latent_var);
	size_t i;
	size_t regime;
	for(regime=0; regime < num_regime; regime++){
		for(i=0; i < num_sbj; i++){
			//gsl_vector_set(eta_local, 0, 3);
            eta_local(0) = 3;
			//gsl_vector_set(eta_local, 1, 1);
            eta_local(1) = 1;
			//gsl_vector_set(eta_0[regime], i*dim_latent_var+0, gsl_vector_get(eta_local, 0));
            (*eta_0[regime])(regime, i*dim_latent_var+0) = eta_local(0);
			//gsl_vector_set(eta_0[regime], i*dim_latent_var+1, gsl_vector_get(eta_local, 1));
            (*eta_0[regime])(regime, i*dim_latent_var+1) = eta_local(1);
			//gsl_vector_set_zero(eta_local);
            eta_local.zeros();
		}
        //gsl_matrix_set((error_cov_0)[regime], 0, 0, -4.60517018598809);
        (*error_cov_0[regime])(0, 0) = -4.60517018598809;
        //gsl_matrix_set((error_cov_0)[regime], 1, 1, -4.60517018598809);
        (*error_cov_0[regime])(1, 1) = -4.60517018598809;
	}
	for(i=0; i < num_sbj; i++){
        //HJ: [to be checked] row vector or column vector?
		mathfunction_softmax(&Padd, pr_0[i]);
	}
    
	//gsl_vector_free(Pvector);
	//gsl_vector_free(Pintercept);
	//gsl_vector_free(Padd);
	//gsl_vector_free(Preset);
	//gsl_vector_free(eta_local);
}


/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 * Do not include parameters in noise_cov matrices 
 */
void function_transform(double *param){
}



void function_regime_switch(size_t t, size_t type, double *param, const vec *co_variate, mat *regime_switch_mat){
	//gsl_matrix_set_identity(regime_switch_mat);
    (*regime_switch_mat) = eye<mat>( size(*regime_switch_mat) );
}



/**
 * The dP/dt function: depend on function_dF_dx, needs to be compiled on the user end
 * but user does not need to modify it or care about it.
 */
void mathfunction_mat_to_vec(const mat *mat_, vec *vec_){
	size_t i,j;
	//size_t nx=mat->size1;
    size_t nx=(*mat_).n_rows;
    
	/*convert matrix to vector*/
	for(i=0; i<nx; i++){
		//gsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));
        (*vec_)(i) = (*mat_)(i, i);
		for (j=i+1;j<nx;j++){
			//gsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));
            (*vec_)(i+j+nx-1) = (*mat_)(i, j);
			/*printf("%lu",i+j+nx-1);}*/
		}
	}
}
void mathfunction_vec_to_mat(const vec *vec_, mat *mat_){
	size_t i,j;
	//size_t nx=mat->size1;
    size_t nx=(*mat_).n_rows;
	/*convert vector to matrix*/
	for(i=0; i<nx; i++){
		//gsl_matrix_set(mat,i,i,gsl_vector_get(vec,i));
        (*mat_)(i,i) = (*vec_)(i);
		for (j=i+1;j<nx;j++){
			//gsl_matrix_set(mat,i,j,gsl_vector_get(vec,i+j+nx-1));
            (*mat_)(i,j) = (*vec_)(i+j+nx-1);
			//gsl_matrix_set(mat,j,i,gsl_vector_get(vec,i+j+nx-1));
            (*mat_)(j,i) = (*vec_)(i+j+nx-1);
		}
	}
}

void function_dP_dt(double t, size_t regime, const vec *p, double *param, size_t n_param, const vec *co_variate, vec *F_dP_dt, int flag){
	
	size_t nx;
    //nx = (size_t) floor(sqrt(2*(double) p->size));
	nx = (size_t) floor(sqrt(2*(double) (*p).n_elem));
    //gsl_matrix *P_mat=gsl_matrix_calloc(nx,nx);
    mat P_mat(nx,nx);
    //mathfunction_vec_to_mat(p,&_mat);
	mathfunction_vec_to_mat(p,&P_mat);
    //gsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);
    mat F_dx_dt_dx(nx,nx);
    // function_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);
	function_dF_dx(t, regime, param, co_variate, &F_dx_dt_dx);
	//gsl_matrix *dFP=gsl_matrix_calloc(nx,nx);
    mat dFP(nx,nx);
	//gsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);
    mat dP_dt(nx,nx);
    //HJ: the function is a wrapper for GSL to call BLAS functions. Armadiilo also use BLAS, so we need to find/write the wrapper for armadillo.
	//gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);
    //dgemm_arma::apply(dFP, F_dx_dt_dx, P_mat, 1.0, 0.0);
	//gsl_matrix_transpose_memcpy(dP_dt, dFP);
    dP_dt = dFP.t();
	//gsl_matrix_add(dP_dt, dFP);
    dP_dt = dP_dt + dFP;
	size_t n_Q_vec=(1+nx)*nx/2;
	//gsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);
    vec Q_vec(n_Q_vec);
	size_t i;
	for(i=1;i<=n_Q_vec;i++){
			//gsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);
            Q_vec(n_Q_vec) = param[n_param-i];
	}
	//gsl_matrix *Q_mat=gsl_matrix_calloc(nx,nx);
    mat Q_mat(nx, nx);
    //mathfunction_vec_to_mat(Q_vec,Q_mat);
	mathfunction_vec_to_mat(&Q_vec,&Q_mat);
	//gsl_matrix_add(dP_dt, Q_mat);
    dP_dt = dP_dt + Q_mat;
	//mathfunction_mat_to_vec(dP_dt, F_dP_dt);
    mathfunction_mat_to_vec(&dP_dt, F_dP_dt);
	
    //gsl_matrix_free(P_mat);
	//gsl_matrix_free(F_dx_dt_dx);
	//gsl_matrix_free(dFP);
	//gsl_matrix_free(dP_dt);
	//gsl_vector_free(Q_vec);
	//gsl_matrix_free(Q_mat);
}

int main(){
    /*
    testing function_dx_dt 
    double t = 1;
    size_t regime =0; 
    vec x(2);
    x(0) = 1;
    x(1) = 2;
    
    double param[2] = {2, 2}; 
    size_t n_param = 2;
    vec co_variate(2);
    vec F_dx_dt(2);
    F_dx_dt.zeros();
    
    x.print("x");
    function_dx_dt(t, regime, &x, param, n_param, &co_variate, &F_dx_dt);
    F_dx_dt.print("F_dx_dt");
    */
    
    /*
    //testing F_dx_dt_dx
    double t = 1;
    size_t regime =0; 
    mat F_dx_dt_dx(2,2);
    vec x(2);
    x(0) = 1;
    x(1) = 2;
    
    double param[16] = {1,1,2,3,4,5,6,7,8,9,2}; 
    size_t n_param = 2;
    vec co_variate(2);
    
    
    F_dx_dt_dx.zeros();
    
    x.print("x");
    function_dF_dx( t,  regime,  param,  &co_variate, &F_dx_dt_dx);
    F_dx_dt_dx.print("F_dx_dt_dx");
    */
    
    /*
    //testing function_noise_cov
    double t = 1;
    size_t regime =0; 
    mat F_dx_dt_dx(2,2);
    vec x(2);
    x(0) = 1;
    x(1) = 2;
    
    double param[16] = {1,1,2,3,4,5,6,7,8,9,2}; 
    size_t n_param = 2;
    mat eta_noise_cov(2,2), y_noise_cov(3,3); 
    eta_noise_cov.zeros();
    y_noise_cov.zeros();
    
    function_noise_cov(t, regime, param, &y_noise_cov, &eta_noise_cov);
    eta_noise_cov.print("eta_noise_cov");
    y_noise_cov.print("y_noise_cov");
    */
    
    mat mat_(3,3);
    vec vec_(9);
    
    mat(0, 0) = 1;
    mat(1, 0) = 2;
    mat(2, 0) = 3;
    mat(0, 1) = 4;
    mat(1, 1) = 5;
    mat(2, 1) = 6;
    mat(0, 2) = 7;
    mat(1, 2) = 8;
    mat(2, 2) = 9;
    mathfunction_mat_to_vec(&mat_, &vec_);
    mat_.print("mat_");
    vec_.print("vec_");
    
    vec(8)= 1;
    vec(7)= 2;
    vec(6)= 3;
    vec(5)= 4;
    vec(4)= 5;
    vec(3)= 6;
    vec(2)= 7;
    vec(1)= 8;
    vec(0)= 9;
    mathfunction_vec_to_mat(&vec_, &mat_);
    mat_.print("mat_");
    vec_.print("vec_");
    
    
    return 0; 
}


