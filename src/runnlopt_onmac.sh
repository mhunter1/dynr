export TMPDIR=/private/var/folders/Ej/EjKBBrnaEAKYcCQi+UdvK++++Tc
sudo mkdir /private/var/folders/Ej
sudo mkdir /private/var/folders/Ej/EjKBBrnaEAKYcCQi+UdvK++++Tc
sudo chmod 777 /private/var/folders/Ej/EjKBBrnaEAKYcCQi+UdvK++++Tc
rm *.o
gcc -c -g wrappernegloglike.c numeric_derivatives.c estimation_nlopt.c PANAmodel.c adaodesolver.c brekfis.c ekf.c math_function.c -I/opt/local/include -I/usr/local/include
gcc -g wrappernegloglike.o numeric_derivatives.o estimation_nlopt.o PANAmodel.o adaodesolver.o brekfis.o ekf.o math_function.o -L/opt/local/lib -lgsl -lgslcblas -L/usr/local/lib -lnlopt -lm -o dynr
./dynr

