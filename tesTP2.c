#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"
#include "maille.h"

int main(){
  

  int typel ;
  float *w;
  float x[2];
  x[0]=1;
  x[1]=1./2; 
  printf("x=%f %f  \n", x[0], x[1]);
	
  typel=2;
  w=malloc(3*sizeof(float));
  calFbase(typel, x, w);
  printf("pout tout i, wi(x) = %f %f %f %f\n", w[0], w[1], w[2], w[3]);
  
  float **dw=alloctabf(3,2);
  calDerFbase(typel, x, dw);
  printf("i=1,2 , dw1(x)/dxi = %f %f \n", dw[0][0], dw[0][1]);
  printf("i=1,2 , dw2(x)/dxi = %f %f \n", dw[1][0], dw[1][1]);
  printf("i=1,2 , dw3(x)/dxi = %f %f \n", dw[2][0], dw[2][1]);
  
  float **S = alloctabf(3,2);
  S[0][0]=0; S[0][1]=0;
  S[1][0]=1; S[1][1]=0;
  S[2][0]=0; S[2][1]=1;
  
  float *Y = malloc(2*sizeof(float));
  transFk(typel, S, Y, w);
  printf("Y=transFk(x) : Y = (%f,%f)\n", Y[0], Y[1]);
  
  float **Jac;
  Jac = alloctabf(2,2);
  matJacob(typel, Jac, S, dw);
  printf("Matrice jacobienne : %f %f // %f %f \n", Jac[0][0], Jac[0][1], Jac[1][0], Jac[1][1]); 
  
  float **Inv;
  float det;
  Inv = alloctabf(2,2);
  invertM2x2(Jac, &det, Inv);
  printf("Inverse jacobienne : %f %f // %f %f \n", Inv[0][0], Inv[0][1], Inv[1][0], Inv[1][1]); 
  printf("det : %f \n", det);
  
  int aretNum = 1;
  int *ptsNum = malloc(2*sizeof(int));
  numNaret(typel, aretNum, ptsNum);
  printf("ptsNum : %d %d \n", ptsNum[0], ptsNum[1]);
}
