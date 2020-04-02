gfortran -c *.f
gcc -Wall calElmt.c int.c fonc.c maille.c impcalel.c alloctab.c syslin.c freetab.c assemblage.c resol.c Tests.c *.o -lm -lgfortran
