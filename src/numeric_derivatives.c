#include "numeric_derivatives.h"
#include "wrappernegloglike.h"
#include "data_structure.h"
#include <string.h>
#include <math.h>/*sqrt(double),pow*/



void forward_diff_grad(double *grad_approx, double ref_fit, const double *x, void * data, double (*func_obj)(const double *, void *))
{
    Data_and_Model data_model=*((Data_and_Model *)data);/*dereference the void pointer*/
    
    double eps = 1e-4;
    double new_point[data_model.pc.num_func_param];
    memcpy(new_point, x, sizeof(new_point));
    int i;
    for(i=0; i < data_model.pc.num_func_param; i++){
        new_point[i] += eps;
        grad_approx[i] = (func_obj(new_point, data) - ref_fit)/eps;
        /*use exp(log(numerator) - log(denominator)) */
        new_point[i] = x[i];
    }
}


double neg_log_like_with_grad(unsigned n, const double *x, double *grad, void *my_func_data)
{
    double fitval = function_neg_log_like(x, my_func_data);
    if (grad) {
        forward_diff_grad(grad, fitval, x, my_func_data, function_neg_log_like);
    }
    return fitval;
}



void hessianR(const double *x,void *data,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian){
    
    Data_and_Model data_model=*((Data_and_Model *)data);/*dereference the void pointer*/
    
    double eps = 1e-4;
    size_t row_index,col_index;
    double xadd[data_model.pc.num_func_param];
    double xsub[data_model.pc.num_func_param];
    double xaa[data_model.pc.num_func_param];
    double xas[data_model.pc.num_func_param];
    double xsa[data_model.pc.num_func_param];
    double xss[data_model.pc.num_func_param];
    
    double fxaa,fxas,fxsa,fxss,hvalue;
    
    
    for(col_index=0; col_index<Hessian->size1; col_index++){
        memcpy(xadd, x, sizeof(xadd));
        memcpy(xsub, x, sizeof(xsub));
        xadd[col_index]=x[col_index] + eps;
        xsub[col_index]=x[col_index] - eps;
        
        gsl_matrix_set(Hessian,col_index,col_index,(func_obj(xadd,data) - 2*fx + func_obj(xsub, data))/(eps*eps));
        
        for(row_index=0; row_index<col_index; row_index++){
            memcpy(xaa, xadd, sizeof(xaa));
            memcpy(xas, xsub, sizeof(xas));
            memcpy(xsa, xadd, sizeof(xsa));
            memcpy(xss, xsub, sizeof(xss));
            
            xaa[row_index] = x[row_index] + eps;
            xas[row_index] = x[row_index] + eps;
            xsa[row_index] = x[row_index] - eps;
            xss[row_index] = x[row_index] - eps;
            
            fxaa=func_obj(xaa,data);
            fxas=func_obj(xas,data);
            fxsa=func_obj(xsa,data);
            fxss=func_obj(xss,data);
            
            hvalue=(fxaa - fxas -fxsa + fxss)/(4 * eps * eps);
            
            gsl_matrix_set(Hessian, row_index, col_index, hvalue);
            gsl_matrix_set(Hessian, col_index, row_index, hvalue);
            
            }
        }
     
     }

void hessianRichardson(const double *x,void *data,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian){
	int i,j;
	for(i=0; i < Hessian->size1; i++){
		hessianOnDiagonal(x, data, func_obj, fx, Hessian, i);
	}
	
	for(j=0; j < Hessian->size1; j++){
		for(i=j+1; i < Hessian->size1; i++){
			hessianOffDiagonal(x, data, func_obj, fx, Hessian, i, j);
		}
	}
}

void hessianOnDiagonal(const double *x,void *data,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian, int index){
	Data_and_Model data_model=*((Data_and_Model *)data);/*dereference the void pointer*/
	
	/*TODO Consider elaborating this to 
	 * double d=0.1; max(fabs(d * x[index]), stepSize)
	 *
	*/
	double stepSize = 1e-4;
	int numIter = 4;
	int k,m;
	static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
	double iOffset = (fabs(stepSize * x[index])) > stepSize ? (fabs(stepSize * x[index])) : stepSize; /* if a > b, then a, else b*/
	
	double xWiggle[data_model.pc.num_func_param];
	memcpy(xWiggle, x, sizeof(xWiggle));
	double Happrox[numIter];
	
	for(k = 0; k < numIter; k++) {
		xWiggle[index] = x[index] + iOffset;
		double f1 = func_obj(xWiggle, data);
		xWiggle[index] = x[index] - iOffset;
		double f2 = func_obj(xWiggle, data);
		Happrox[k] = (f1 - 2.0 * fx + f2) / (iOffset * iOffset);
		xWiggle[index] = x[index];
		iOffset /= v;
	}
	
	for(m = 1; m < numIter; m++) {						// Richardson Step
		for(k = 0; k < (numIter - m); k++) {
			// NumDeriv Hard-wires 4s for r here. Don't remember why.
			Happrox[k] = (Happrox[k+1] * pow(4.0, m) - Happrox[k])/(pow(4.0, m)-1);
		}
	}
	gsl_matrix_set(Hessian, index, index, Happrox[0]);

}

void hessianOffDiagonal(const double *x,void *data,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian, int row_index, int col_index){
	Data_and_Model data_model=*((Data_and_Model *)data);/*dereference the void pointer*/
	
	double stepSize = 1e-4;
	int numIter = 4;
	int k,m;
	static const double v = 2.0; //Note: NumDeriv comments that this could be a parameter, but is hard-coded in the algorithm
	double iOffset = (fabs(stepSize * x[row_index])) > stepSize ? (fabs(stepSize * x[row_index])) : stepSize; /* if a > b, then a, else b*/
	double jOffset = (fabs(stepSize * x[col_index])) > stepSize ? (fabs(stepSize * x[col_index])) : stepSize;
	
	double xWiggle[data_model.pc.num_func_param];
	memcpy(xWiggle, x, sizeof(xWiggle));
	double Happrox[numIter];
	
	for(k = 0; k < numIter; k++) {
		xWiggle[row_index] = x[row_index] + iOffset;
		xWiggle[col_index] = x[col_index] + jOffset;
		double f1 = func_obj(xWiggle, data);
		xWiggle[row_index] = x[row_index] - iOffset;
		xWiggle[col_index] = x[col_index] - jOffset;
		double f2 = func_obj(xWiggle, data);
		Happrox[k] = (f1 - 2.0 * fx + f2 - gsl_matrix_get(Hessian, row_index, row_index)*iOffset*iOffset - gsl_matrix_get(Hessian, col_index, col_index)*jOffset*jOffset ) / (2.0 * iOffset * jOffset);
		xWiggle[row_index] = x[row_index];
		xWiggle[col_index] = x[col_index];
		iOffset /= v;
		jOffset /= v;
	}
	
	for(m = 1; m < numIter; m++) {						// Richardson Step
		for(k = 0; k < (numIter - m); k++) {
			// NumDeriv Hard-wires 4s for r here. Don't remember why.
			Happrox[k] = (Happrox[k+1] * pow(4.0, m) - Happrox[k])/(pow(4.0, m)-1);
		}
	}
	gsl_matrix_set(Hessian, row_index, col_index, Happrox[0]);
	gsl_matrix_set(Hessian, col_index, row_index, Happrox[0]);
}


