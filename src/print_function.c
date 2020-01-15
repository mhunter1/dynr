#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <R.h>

#define DYNRPRINT(verbose_flag, ...) if(verbose_flag){Rprintf(__VA_ARGS__);}
#define MYPRINT(...) Rprintf(__VA_ARGS__)
/**
 * print the given vector's value to the console
 * format (v1, v2, ... )
 */
void print_vector(const gsl_vector *y){
	if(y==NULL){
		MYPRINT("( NULL )");
		return;
	}
	size_t index;
	if(y->size<=0)
		return;
	MYPRINT("(%.3f",gsl_vector_get(y, 0));
	for(index=1; index<y->size; index++)
		MYPRINT(", %.3f", gsl_vector_get(y, index));
	MYPRINT(")\n");
}

/**
 * print given array to the console
 * format is [v1, ..., vn]
 */
void print_array(const double *v, int n){
	size_t index;
	if(n<=0)
		return;
	MYPRINT("[%.3f", v[0]);
	for(index=1; index<n; index++)
		MYPRINT(", %.3f", v[index]);
	MYPRINT("]\n");
}

/**
 * print the given matrix;
 */
void print_matrix(const gsl_matrix *mat){
	size_t ri, ci;
	if(mat->size1<=0 || mat->size2<=0)
		return;
	for(ri=0; ri<mat->size1; ri++){
		MYPRINT("  %.7f", gsl_matrix_get(mat, ri, 0));
		for(ci=1; ci<mat->size2; ci++)
			MYPRINT(", %.7f", gsl_matrix_get(mat, ri, ci));
		MYPRINT("\n");
	}
}

