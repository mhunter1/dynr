#install.packages("RcppGSL")
require(RcppGSL)
require(inline)
inctxt="
#include <gsl/gsl_matrix.h>
void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){
    size_t i;
    for (i=0;i<eta_noise_cov->size1;i++){
        gsl_matrix_set(eta_noise_cov,i,i,-20);
    }
    gsl_matrix_set(y_noise_cov,0,0, param[4]);
    gsl_matrix_set(y_noise_cov,1,1, param[5]);

}
"
bdtxt=""
#funct_noise_cov<- cxxfunction(signature(),inc=inctxt,plugin="RcppGSL")
funct_noise_cov<-cfunction(sig=character(), body=bdtxt, includes =inctxt,language="C")
