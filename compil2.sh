gfortran -c *.f
gcc -Wall calElmt.c int.c maille.c impcalel.c alloctab.c freetab.c assemblage.c testTP3.c *.o -lm -lgfortran
