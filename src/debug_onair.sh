rm *.o
gcc -c -g *.c -I/opt/local/include
gcc -g *.o -L/opt/local/lib -lgsl -lgslcblas -o estimate
valgrind --tool=memcheck --leak-check=yes --leak-check=full --show-leak-kinds=all --dsymutil=yes --track-origins=yes  ./estimate

