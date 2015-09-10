rm *.o
R CMD SHLIB -c -g *.c -I/opt/local/include -L/opt/local/lib -lgsl -lgslcblas



