// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

//This function prints out "Hello Hui-Ju!" and returm an identity matrix
arma::mat _rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    Rprintf("Hello Hui-Ju!\n");
    return m1;
}
