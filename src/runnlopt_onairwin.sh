#!/bin/bash
rm *.o
gcc -c wrappernegloglike.c estimation_nlopt.c PANAmodel.c functions/*.c -I/usr/local/include -I/usr/include
gcc wrappernegloglike.o estimation_nlopt.o PANAmodel.o adaodesolver.o brekfis.o ekf.o math_function.o -L/usr/lib -lgsl -lgslcblas -L/usr/local/lib -lnlopt -lm -o estimate
./estimate

