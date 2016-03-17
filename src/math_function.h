#ifndef MATH_FUNCTION_H_INCLUDED
#define MATH_FUNCTION_H_INCLUDED

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdbool.h>

/**
 * This method computes the log-likelihood of a multivariate normal distribution.
 * @param x the variable (column) vector
 * @param cov_matrix the covariance matrix
 */
double mathfunction_multivariate_normal(const gsl_vector *x, const gsl_matrix *cov_matrix);

/**
 * This method computes the determinant of the given matrix.
 * The algorithm will firstly convert the given matrix to a upper triangle one. Then calculate the multiplication of diagonal values.
 * Please refer to http://en.wikipedia.org/wiki/Determinant for details.
 * @matrix the given matrix, should be a square one.
 * @return the determinant of the matrix.
 */
double mathfunction_matrix_determinant(const gsl_matrix *matrix);

/**
 * Given a matrix of log-values, this method normalize each element to v~(0,1) and the sum of them will be equal to 1.
 * For example, given (-1, -2), the resulted value will be (e^{-1}, e^{-2})/(e^{-1}+e^{-2}).
 * @param log_v the target matrix values. After this function, the values in original one will be modified.
 * @return normalizer
 */
double mathfunction_normalize_log(gsl_matrix *log_v);

/**
 * Given a vector of log-values, this method normalize each element to v~(0,1) and the sum of them will be equal to 1.
 * For example, given (-1, -2), the resulted value will be (e^{-1}, e^{-2})/(e^{-1}+e^{-2}).
 * @param log_v the target vector values. After this function, the values in original one will be modified.
 * @return the normalizer
 */
double mathfunction_normalize_log_vector(gsl_vector *log_v);

/**
 * This method normalize the given matrix's values so that the sum of them is equal to 1.
 * @param v the given matrix to be normalized.
 * @return the normalizer
 *
 */
double mathfunction_matrix_normalize(gsl_matrix *v);

/**
 * This method normalizes the given vector's values so that sum of them is equal to 1.
 * @param v the given vector to be normalized.
 * @return the normalizer
 */
double mathfunction_vector_normalize(gsl_vector *v);

/**
 * compute the inverse of a given matrix
 * @param mat the given matrix
 * @param inov_mat the matrix where the inverse one is stored.
 */
void mathfunction_inv_matrix(const gsl_matrix *mat, gsl_matrix *inv_mat);

/**
 * This method computes the trace of the given matrix
 * @param mat the target matrix, make sure the matrix is a square one.
 * @return the trace
 */
double mathfunction_mat_trace(const gsl_matrix *mat);

/**
 * print the given vector's value to the console
 */
void print_vector(const gsl_vector *y);

/**
 * print given array to the console
 * format is [v1, ..., vn]
 */
void print_array(const double *v, int n);

/**
 * print the given matrix;
 */
void print_matrix(const gsl_matrix *mat);

void print_buffer(const void *buffer, size_t len);

/**
 * This method computes C=A*B or C=A'*B or C=A*B' or C=A'*B'
 * @param mat_a the matrix of A
 * @param mat_b the matrix of B
 * @param transpose_a whether A will be transposed
 * @param transpose_b whether B will be transposed
 * @param mat_c the result of multiplication.
 */
void mathfunction_matrix_mul(const gsl_matrix *mat_a, const gsl_matrix *mat_b, bool transpose_a, bool transpose_b, gsl_matrix *mat_c);

/**
 * This function sums the given vector
 */
double mathfunction_sum_vector(const gsl_vector *vec);
double mathfunction_min(const double x,const double y,const double z);
double mathfunction_inv_matrix_det(const gsl_matrix *mat, gsl_matrix *inv_mat); /*via Cholesky decomp*/
double mathfunction_cholesky_det(const gsl_matrix *mat);
double mathfunction_inv_matrix_det_lu(const gsl_matrix *mat, gsl_matrix *inv_mat); /*via LU decomp*/
double mathfunction_negloglike_multivariate_normal_invcov(const gsl_vector *x, const gsl_matrix *inv_cov_matrix, double det);
/**
 * convert a matrix (e.g.,
 * [1 4 5
 *  4 2 6
 *  5 6 3]) to a vector (e.g., [1 2 3 4 5 6])
 * @param mat the given matrix
 */

void mathfunction_mat_to_vec(const gsl_matrix *mat, gsl_vector *vec);
/**
 * convert a vector (e.g., [1 2 3 4 5 6]) to a matrix (e.g.,
 * [1 4 5
 *  4 2 6
 *  5 6 3])
 * @param vec the given vector
 */

void mathfunction_vec_to_mat(const gsl_vector *vec, gsl_matrix *mat);


/**
 * generate a random number ranges from 0 to 1
 */
double drand();

/**
 * generate a random number satisfies standard normal distribution.
 */
double random_std_normal();

/**
 * generate a random number satisfies any normal distribution
 */
double random_normal(double mu, double sigma);

/**
 * generate a vector of random numbers satisfies given normal distribution with zero mean.
 */
void white_noise(const gsl_vector *sigma, gsl_vector *noise);

void random_pos_id_mat(gsl_matrix *mat);
/**
 * This function scales vector A and store the result in B
 * @param vec_a vector A
 * @param vec_b the A*c of scaling
 */
void mathfunction_vec_scale(const gsl_vector *vec_a, const double c, gsl_vector *vec_b);
/**
 * This function scales matrix A and store the result in B
 * @param vec_a vector A
 * @param vec_b the A*c of scaling
 */
void mathfunction_mat_scale(const gsl_matrix *mat_a, const double x, gsl_matrix *mat_b);
/**
 * This function convert a vector A into a matrix B so that the diagonal of matrix B is A*x
 * @param vec_a vector A
 * @param mat_b matrix B
 */
void mathfunction_diagin_scale(const gsl_vector *vec_a, const double x, gsl_matrix *mat_b);
void mathfunction_diagout_scale(const gsl_matrix *mat_a, const double x, gsl_vector *vec_b);


/**
 * This function takes a double and gives back a double
 * It computes the logistic function (i.e. the inverse of the logit link function)
 * @param x, the double value e.g. a normally distributed number
 * @return logistic(x), the double value e.g. a number between 0 and 1
 */
double mathfunction_logistic(const double x);


/**
 * This function takes a gsl_vector and modifies its second argument (another gsl_vector)
 * It computes the softmax function (e.g. for multinomial logistic regression)
 * @param x, the double value e.g. a vector of normally distributed numbers
 * @param result, softmax(x), e.g. a numbers between 0 and 1 that sum to 1
 */
void mathfunction_softmax(const gsl_vector *x, gsl_vector *result);

#endif
