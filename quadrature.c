#include <stdio.h>
#include <stdlib.h>
#include "maille.h"

void ppquad(int t, float *W, float **X){
  switch(t){
    case 1 :
      for(int i=0; i<4; i++){
        W[i]=1./36;
      }
      for(int i=4; i<8; i++){
        W[i]=1./9;
      }
      W[8]=4./9;
      X[0][0]=1; X[0][1]=0;
      X[1][0]=1; X[1][1]=1;
      X[2][0]=0; X[2][1]=1;
      X[3][0]=0; X[3][1]=0;
      X[4][0]=1; X[4][1]=1./2;
      X[5][0]=1./2; X[5][1]=1;
      X[6][0]=0; X[6][1]=1./2;
      X[7][0]=1./2; X[7][1]=0;
      X[8][0]=1./2; X[8][1]=1./2;
    break;
    case 2 :
      for(int i=0; i<2; i++){
        W[i]=1./6;
      }
      X[0][0]=1./2; X[0][1]=1./2;
      X[1][0]=0; X[1][1]=1./2;
      X[2][0]=1./2; X[2][1]=0;
    break;
    case 3 :
      W[0]=1./6; W[1]=1./6; W[2]=2./3; 
      X[0][0] = 1; X[1][0] = 0; X[2][0] = 1./2;
    break;
  }
}

int main(){
  int t=2;
  float *W;
  float **X;
  W = malloc(3*sizeof(float));
  X=alloctabf(3,2);
  ppquad(t,W,X);
  for(int i=0;i<3;i++){
    float*x = X[i];
    printf("%f,%f\n",x[0],x[1]);
  }
  free(W);
  freetab(X);
}
