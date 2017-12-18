#include <math.h>
#include <armadillo>
#include <stdlib.h>
//#include <armadillo/include/armadillo_bits/mul_gemm.hpp>

using namespace arma ;
#define omega 61.68503

//#include <armadillo/include/armadillo_bits/mul_gemm.hpp>


/*Tested: mathfunction_logistic, mathfunction_softmax, function_dx_dt, function_dF_dx, function_noise_cov mathfunction_mat_to_vec, mathfunction_vec_to_mat,  function_regime_switch, function_measurement, function_dP_dt, function_initial_condition(done)*/
/*Converted done:  done (compiled successfully)*/

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
 

    //gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);
    //interface of dgemv: int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
    // behavior of dgemv y := alpha*A*x + beta*y
    (*y) = 1.0* (*Ht) * (*eta) + 0.0 * (*y);
    
 
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





void function_initial_condition(double *param, vec *co_variate, vec *pr_0, vec *eta_0, mat *error_cov_0, size_t *index_sbj){

//void function_initial_condition(double *param, vec *co_variate, vec *pr_0, vec *eta_0, mat *error_cov_0, size_t *index_sbj){
	
	//gsl_vector *Pvector = gsl_vector_calloc(1);
	//gsl_vector *Pintercept = gsl_vector_calloc(1);
	//gsl_vector *Padd = gsl_vector_calloc(1);
	//gsl_vector *Preset = gsl_vector_calloc(1);
    //gsl_vector *eta_local = gsl_vector_calloc(2);
    vec Pvector(1), Pintercept(1), Padd(1), Preset(1), eta_local(2);
    Pvector.zeros();
    Pintercept.zeros();
    Padd.zeros();
    eta_local.zeros();
    
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
	//Preset.print("Preset");
    
	//size_t num_regime=pr_0[0]->size;
    size_t num_regime=(pr_0[0]).n_elem;
	//size_t dim_latent_var=error_cov_0[0]->size1;
    size_t dim_latent_var=(error_cov_0[0]).n_elem;
	//size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
    size_t num_sbj=((eta_0[0]).n_elem)/(dim_latent_var);
    
    //printf("%d %d %d\n",int(num_regime),int(dim_latent_var), int (num_sbj));
	size_t i;
	size_t regime;
	for(regime=0; regime < num_regime; regime++){
		for(i=0; i < num_sbj; i++){
			//gsl_vector_set(eta_local, 0, 3);
            eta_local(0) = 3;
			//gsl_vector_set(eta_local, 1, 1);
            eta_local(1) = 1;
			//gsl_vector_set(eta_0[regime], i*dim_latent_var+0, gsl_vector_get(eta_local, 0));
            (eta_0[regime])(regime, i*dim_latent_var+0) = eta_local(0);
			//gsl_vector_set(eta_0[regime], i*dim_latent_var+1, gsl_vector_get(eta_local, 1));
            (eta_0[regime])(regime, i*dim_latent_var+1) = eta_local(1);
			//gsl_vector_set_zero(eta_local);
            eta_local.zeros();
		}
        //gsl_matrix_set((error_cov_0)[regime], 0, 0, -4.60517018598809);
        (error_cov_0[regime])(0, 0) = -4.60517018598809;
        //gsl_matrix_set((error_cov_0)[regime], 1, 1, -4.60517018598809);
        (error_cov_0[regime])(1, 1) = -4.60517018598809;
	}
	for(i=0; i < num_sbj; i++){
		mathfunction_softmax(&Padd, &pr_0[i]);
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
    
    // GSL will initialize vectors and matrices as zero ones
    (*vec_).zeros();
    
	/*convert matrix to vector*/
	for(i=0; i<nx; i++){
		//gsl_vector_set(vec,i,gsl_matrix_get(mat,i,i));
        (*vec_)(i) = (*mat_)(i, i);
		for (j=i+1;j<nx;j++){
			//gsl_vector_set(vec,i+j+nx-1,gsl_matrix_get(mat,i,j));
            (*vec_)(i+j+nx-1) = (*mat_)(i, j);
			//printf("%lu",i+j+nx-1);
		}
	}
}
void mathfunction_vec_to_mat(const vec *vec_, mat *mat_){
	size_t i,j;
	//size_t nx=mat->size1;
    size_t nx=(*mat_).n_rows;
	/*convert vector to matrix*/
    
    // GSL will initialize vectors and matrices as zero ones
    (*mat_).zeros();
    
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
    P_mat.zeros();
    //mathfunction_vec_to_mat(p,&_mat);
	mathfunction_vec_to_mat(p,&P_mat);
    //P_mat.print("P_mat");
    
    //gsl_matrix *F_dx_dt_dx=gsl_matrix_calloc(nx,nx);
    mat F_dx_dt_dx(nx,nx);
    F_dx_dt_dx.zeros();
    // function_dF_dx(t, regime, param, co_variate, F_dx_dt_dx);
	function_dF_dx(t, regime, param, co_variate, &F_dx_dt_dx);
    //F_dx_dt_dx.print("F_dx_dt_dx");
	//gsl_matrix *dFP=gsl_matrix_calloc(nx,nx);
    mat dFP(nx,nx);
    dFP.zeros();
	//gsl_matrix *dP_dt=gsl_matrix_calloc(nx,nx);
    mat dP_dt(nx,nx);
    dP_dt.zeros();
    
	//gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, F_dx_dt_dx, P_mat, 0.0, dFP);
    (dFP) = 1.0*(F_dx_dt_dx)*(P_mat)+ 0.0 * (dFP);
    //dFP.print("dFP");
    
	//gsl_matrix_transpose_memcpy(dP_dt, dFP);
    dP_dt = dFP.t();
	//gsl_matrix_add(dP_dt, dFP);
    dP_dt = dP_dt + dFP;
    //dP_dt.print("dP_dt");
    
	size_t n_Q_vec=(1+nx)*nx/2;
	//gsl_vector *Q_vec=gsl_vector_calloc(n_Q_vec);
    //printf("%d\n", (int)n_Q_vec);
    vec Q_vec(n_Q_vec);
    Q_vec.zeros();
	size_t i;
	for(i=1;i<=n_Q_vec;i++){
			//gsl_vector_set(Q_vec,n_Q_vec-i,param[n_param-i]);
            Q_vec(n_Q_vec-i) = param[n_param-i];
            
	}
    //Q_vec.print("Q_vec");
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
    
    /*
    mat mat_(3,3);
    vec vec_(9);
    vec_.zeros();
    
    //mathfunction_mat_to_vec
    mat_(0, 0) = 1;
    mat_(1, 0) = 2;
    mat_(2, 0) = 3;
    mat_(0, 1) = 2;
    mat_(1, 1) = 5;
    mat_(2, 1) = 6;
    mat_(0, 2) = 3;
    mat_(1, 2) = 6;
    mat_(2, 2) = 9;
    mathfunction_mat_to_vec(&mat_, &vec_);
    mat_.print("mat_");
    vec_.print("vec_");
    */
    
    /*
    //mathfunction_vec_to_mat
    mat_.zeros();
    vec_(8)= 1;
    vec_(7)= 2;
    vec_(6)= 3;
    vec_(5)= 4;
    vec_(4)= 5;
    vec_(3)= 6;
    vec_(2)= 7;
    vec_(1)= 8;
    vec_(0)= 9;
    mathfunction_vec_to_mat(&vec_, &mat_);
    vec_.print("vec_");
    mat_.print("mat_");
    
    size_t regime =0;
    double param[2] = {2, 2}; 
    vec co_variate(2);
    mat regime_switch_mat(4,4);
    function_regime_switch(regime, regime, param, &co_variate, &regime_switch_mat);
    regime_switch_mat.print("regime_switch_mat");
    */
    
    //function_initial_condition
    //printf("beginning");
    //mat error_cov_0[10];
    vec *co_variate;
    //co_variate = (vec **)malloc(sizeof(vec*));
    co_variate = (vec *)malloc(5* sizeof(vec));
    co_variate[0].set_size(2,2);
    co_variate[0].ones();
    //co_variate[0].print("co_variate");
    
    vec *pr_0;
    //pr_0 = (vec **)malloc(sizeof(vec*));
    pr_0 = (vec *)malloc(5* sizeof(vec));
    pr_0[0].set_size(2);
    (pr_0[0])(0) = 2;
    (pr_0[0])(1) = 2;
    //pr_0[0].print("pr_0[0]");
    
    vec *eta_0;
    //eta_0 = (vec **)malloc(sizeof(vec*));
    eta_0 = (vec *)malloc(5 * sizeof(vec));
    eta_0[0].set_size(2);
    eta_0[0].zeros();
    eta_0[0](0) = 2;
    (eta_0[0])(1) = 2;
    //eta_0[0].print("(eta_0[0])");
    
    mat *error_cov_0;
    //error_cov_0 = (mat **)malloc(sizeof(mat*));
    error_cov_0 = (mat *)malloc(5* sizeof(mat));
    error_cov_0[0].set_size(2,2);
    error_cov_0[1].set_size(2,2);
    error_cov_0[0].zeros(2,2);
    error_cov_0[1].zeros(2,2);
    //(error_cov_0[0])(0,0) = 2;
    //(error_cov_0[0])(1,0) = 2; 
    //error_cov_0[0].print("error_cov_0");
    
    double param[16] = {1,1,2,3,4,5,6,7,8,9,2}; 
    size_t index_sbj[10] = {1,1,1,1,1,1,1,1,1};

    //printf("here\n");
    function_initial_condition(param, co_variate, pr_0, eta_0, error_cov_0, index_sbj);
    error_cov_0[0].print("error");
    error_cov_0[1].print("error2");
    
    
    /*
    //function_measurement
    size_t t = 1, regime = 0;
    double param[16] = {1,1,2,3,4,5,6,7,8,9,2};
    mat Ht(3,3); 
    vec eta(3), co_variate(3), y(3);
    Ht.zeros();
    eta.ones();
    co_variate.zeros();
    y.zeros();
    function_measurement( t,  regime,  param,  &eta,  &co_variate, &Ht, &y);
    y.print("y");
    co_variate.print("co_variate");
    eta.print("eta");
    Ht.print("Ht");
    */
    
    /*
    //function_dP_dt
    double t = 1;
    size_t regime = 0, n_param = 16;
    double param[16] = {1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    mat Ht(3,3); 
    vec eta(9), co_variate(9), y(9), p(9), F_dP_dt(9);
    Ht.zeros();
    eta.ones();
    co_variate.zeros();
    y.zeros();
    p.zeros();
    p(1) = 1;
    p(0) = 2;
    
    function_dP_dt(t,  regime, &p, param, n_param, &co_variate, &F_dP_dt, 0);
    F_dP_dt.print("F_dP_dt");
    */
    return 0; 
}


