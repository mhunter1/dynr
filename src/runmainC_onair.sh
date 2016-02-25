rm *.o
gcc -c -g wrappernegloglike.c mainC.c PANAmodel.c functions/*.c -I/opt/local/include
gcc -g wrappernegloglike.o estimation_nlopt.o PANAmodel.o adaodesolver.o brekfis.o ekf.o math_function.o -L/opt/local/lib -lgsl -lgslcblas -o estimate
./estimate

