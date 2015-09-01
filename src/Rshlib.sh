rm *.o
R CMD SHLIB -c -g wrappernegloglike.c numeric_derivatives.c estimation_nloptR.c PANAmodel.c functions/*.c -I/opt/local/include -I/usr/local/include -L/opt/local/lib -lgsl -lgslcblas -L/usr/local/lib -lnlopt -lm



