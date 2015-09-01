#include "numeric_derivatives.h"
#include "wrappernegloglike.h"
#include <string.h>
#include <math.h>/*sqrt(double),pow*/



void forward_diff_grad(double *grad_approx, double ref_fit, const double *x, void * data, double (*func_obj)(const double *, void *))
{
    double eps = 1e-4;
    double new_point[6];
    memcpy(new_point, x, sizeof(new_point));
    int i;
    for(i=0; i < 6; i++){
        new_point[i] += eps;
        grad_approx[i] = (func_obj(new_point, data) - ref_fit)/eps;
        /*use exp(log(numerator) - log(denominator)) */
        new_point[i] = x[i];
    }
}

void hessian(const double *x,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian){
    double eps = 1e-4;
    size_t row_index,col_index;
    double xadd[6];
    double xsub[6];
    double xaa[6];
    double xas[6];
    double xsa[6];
    double xss[6];
    
    double fxaa,fxas,fxsa,fxss,hvalue;

    
    for(col_index=0; col_index<Hessian->size1; col_index++){
        memcpy(xadd, x, sizeof(xadd));
        memcpy(xsub, x, sizeof(xsub));
        xadd[col_index]=x[col_index] + eps;
        xsub[col_index]=x[col_index] - eps;
        
        gsl_matrix_set(Hessian,col_index,col_index,(func_obj(xadd,NULL) - 2*fx + func_obj(xsub, NULL))/pow(eps,2.0));
        
        for(row_index=0;row_index<col_index;row_index++){
            memcpy(xaa, xadd, sizeof(xaa));
            memcpy(xas, xsub, sizeof(xas));
            memcpy(xsa, xadd, sizeof(xsa));
            memcpy(xss, xsub, sizeof(xss));
            
            xaa[row_index] = x[row_index] + eps;
            xas[row_index] = x[row_index] + eps;
            xsa[row_index] = x[row_index] - eps;
            xss[row_index] = x[row_index] - eps;
            
            fxaa=func_obj(xaa,NULL);
            fxas=func_obj(xas,NULL);
            fxsa=func_obj(xsa,NULL);
            fxss=func_obj(xss,NULL);
            
            hvalue=(fxaa - fxas -fxsa + fxss)/(4 * pow(eps,2.0));
            
            gsl_matrix_set(Hessian, row_index, col_index, hvalue);
            gsl_matrix_set(Hessian, col_index, row_index, hvalue);
            
            }
        }
    
     }


double myfunc_wrapper(unsigned n, const double *x, double *grad, void *my_func_data)
{
    double fitval = function_neg_log_like(x, my_func_data);
    if (grad) {
        forward_diff_grad(grad, fitval, x, my_func_data, function_neg_log_like);
    }
    return fitval;
}



void hessianR(const double *x,void *data,double (*func_obj)(const double *, void *), double fx, gsl_matrix *Hessian){
    double eps = 1e-4;
    size_t row_index,col_index;
    double xadd[6];
    double xsub[6];
    double xaa[6];
    double xas[6];
    double xsa[6];
    double xss[6];
    
    double fxaa,fxas,fxsa,fxss,hvalue;

    
    for(col_index=0; col_index<Hessian->size1; col_index++){
        memcpy(xadd, x, sizeof(xadd));
        memcpy(xsub, x, sizeof(xsub));
        xadd[col_index]=x[col_index] + eps;
        xsub[col_index]=x[col_index] - eps;
        
        gsl_matrix_set(Hessian,col_index,col_index,(func_obj(xadd,data) - 2*fx + func_obj(xsub, data))/pow(eps,2.0));
        
        for(row_index=0;row_index<col_index;row_index++){
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
            
            hvalue=(fxaa - fxas -fxsa + fxss)/(4 * pow(eps,2.0));
            
            gsl_matrix_set(Hessian, row_index, col_index, hvalue);
            gsl_matrix_set(Hessian, col_index, row_index, hvalue);
            
            }
        }
    
     }


