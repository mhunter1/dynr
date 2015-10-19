#install.packages("RcppGSL")
require(RcppGSL)
require(inline)
inctxt="#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>"
bodytxt="
RcppGSL::matrix<double> M = sM;
int k = M.ncol();
Rcpp::NumericVector n(k);
// create gsl data structures from SEXP
// to store results
for (int j = 0; j < k; j++) {
  RcppGSL::vector_view<double> colview = gsl_matrix_column (M, j);
  n[j] = gsl_blas_dnrm2(colview);
}
M.free() ;
return n;                           // return vector
"
foo <- cxxfunction(signature(sM="numeric"), body=bodytxt, inc=inctxt, plugin="RcppGSL")
## see Section 8.4.13 of the GSL manual: create M as a sum of two outer products
M <- outer(sin(0:9), rep(1,10), "*") + outer(rep(1, 10), cos(0:9), "*") 
print(foo(M))

