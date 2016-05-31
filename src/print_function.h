#ifndef PRINT_FUNCTION_H_INCLUDED
#define PRINT_FUNCTION_H_INCLUDED

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <R.h>

#define DYNRPRINT(verbose_flag, ...) if(verbose_flag){Rprintf(__VA_ARGS__);}
#define MYPRINT(...) Rprintf(__VA_ARGS__)
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


#endif
