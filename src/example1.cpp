#include <stdio.h>


#include <armadillo>
/*
using namespace arma;
*/

extern "C" void example1(void);

void example1(void){
	printf("hello world from C++!\n");
	using namespace arma;
	
	arma::mat A(2,2);

	A(0,0) = 1;
	A(0,1) = 0;
	A(1,0) = 0;
	A(1,1) = 4;
	A.print("A");
	
	A = A % A;
	A.print("A^2");
	
	return;
}

