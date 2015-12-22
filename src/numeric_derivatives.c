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


double myfunc_wrapper(unsigned n, const double *x, double *grad, void *my_func_data)
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


