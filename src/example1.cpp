#define ARMA_USE_SUPERLU 1
// [[Rcpp::depends(RcppArmadillo)]]
//#include <RcppArmadillo.h>
#include <armadillo>
//using namespace Rcpp;


extern "C" void example1(void);
void example1(void) {
    printf("hello world from C++\n");

    arma::mat A(2,2);

    A(0,0) = 1;
    A(0,1) = 0;
    A(1,0) = 0;
    A(1,1) = 4;
    A.print("A");

    A = A % A;
    A.print("A^2");


    using namespace arma ;
    sp_mat A2 = sprandu<sp_mat>(4, 4, 0.5);
    A2.zeros();
    A2(0, 0) = 0.1140;
    A2(1, 0) = 0.5299;
    A2(2, 0) = 0.9758;
    A2(3, 0) = 0.6896;
    A2(1, 1) = 0.9641;
    A2(0, 2) = 0.6484;
    A2(2, 2) = 0.7027;
    A2(3, 3) = 0.7040;

    A2.print("A2");
    vec b = randu<vec>(4);
    b(0) = 0.4924;
    b(1) = 0.8198;
    b(2) = 0.1233;
    b(3) = 0.7975;
    b.t().print("b");
    vec x;
    // = spsolve(A, b);  // solve one system
    //x.print("x");
    
    printf("solve\n");
   // bool status = spsolve(x, A, b);  // use default solver
   // if (status == false)  { printf("no solution"); }
    spsolve(x, A2, b, "lapack");   // use LAPACK  solver
    x.print("x(LAPACK)");
    (A2 * x).print("AX");
    //spsolve(x, A, b, "superlu");  // use SuperLU solver
    //x.t().print("x(SUPERLU)");

    printf("Done.\n");
}
/*
int main(){
    using namespace arma ;
    arma_rng::set_seed_random();
    superLuDemo();
    return 0;
}*/

