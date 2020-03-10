#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"

void ppquad(int nbneel, float *W, float **X){
  switch(nbneel){
    case 4 :
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
    case 3 :
      for(int i=0; i<3; i++){
        W[i]=1./6;
      }
      X[0][0]=1./2; X[0][1]=1./2;
      X[1][0]=0; X[1][1]=1./2;
      X[2][0]=1./2; X[2][1]=0;
      break;
    case 2 :
      W[0]=1./6; W[1]=1./6; W[2]=2./3; 
      X[0][0] = 1; X[1][0] = 0; X[2][0] = 1./2;
      break;
  }
}

void invertM2x2(float **M, float *det, float **MInv){
  /* Calcul de l'inverse d'une matrice M et de son déterminant */
  *det=M[0][0]*(M[1][1])-M[0][1]*(M[1][0]);
  float alpha = (1./(*det));
  MInv[0][0]=alpha*M[1][1];
  MInv[0][1]=-alpha*M[0][1];
  MInv[1][0]=-alpha*M[1][0];
  MInv[1][1]=alpha*M[0][0];
}

void matJacob(int nbneel, int dimDom, float **JacFk, float **S, float **Dfoncbase){
  /* Calcul de la matrice Jacobienne JacFk à partir des sommets S
  et des gradients des fonctions de bases fournies par calDerFbase */
  float elmt;
  for(int j=0;j<2;j++){
    for(int k=0;k<dimDom;k++){
      elmt = 0;
      for(int i=0;i<nbneel;i++){
        elmt += Dfoncbase[i][k]*S[i][j];
      }
      JacFk[j][k] = elmt;
    }
  }      
}

void transFk(int nbneel, float **S, float Y[], float *foncbase){
  /* Calcul de Y, l'image par la transformation Fk à partir des sommets S
  et des valeurs foncbase fournies par calFbase */
  for(int j=0;j<2;j++){
    Y[j] = 0;
    for(int i=0; i<nbneel; i++){ 
      Y[j] += foncbase[i]*S[i][j];
    }
  }
}

void calDerFbase(int nbneel, float *x, float **w){
  /* 
  On renvoie un vecteur de taille n qui correspond 
  aux valeurs des gradients fonctions de base en un point x
  et on a : w_k(s_i)=delta_ki pour une sommet si donné
  */
  switch(nbneel){
  case 4 : // quadrangle 
    /*  4  3
        1  2 */
    w[3][0]=-1+x[1]; w[3][1]=-1+x[0]; 
    w[0][0]=1-x[1];  w[0][1]=-x[0];
    w[1][0]=x[1];   w[1][1]=x[0];
    w[2][0]=-x[1];    w[2][1]=1-x[0];
    break;
  case 3 : // triangle 	
    w[2][0]=-1; w[2][1]=-1; 
    w[0][0]=1;  w[0][1]=0;
    w[1][0]=0;  w[1][1]=1;
    break;
  case 2 : // segment
    /*  2  1    */		
    w[0][0]=1;
    w[0][1]=-1;
    break;
  }
}

void calFbase(int nbneel, float *x, float *w){
  /* 
  On renvoie un vecteur de taille n qui correspond 
  aux valeurs des fonctions de base en un point x
  et on a : w_k(s_i)=delta_ki pour un sommet si donné
  */
  switch(nbneel) {
  case 4 : // quadrangle 
  /*  4  3
      1  2      */
    w[3]=1-x[0]-x[1]+x[0]*x[1];
    w[0]=x[0]*(1-x[1]);
    w[1]=x[0]*x[1];
    w[2]=x[1]*(1-x[0]);
    break;
  case 3 : // triangle 
  /*  3  
      1  2      */	
    w[2]=1-x[0]-x[1];
    w[0]=x[0];
    w[1]=x[1];
    break;
  case 2 : // segment
  /*  2  1    */			
    w[0]=x[0];
    w[1]=1-x[0];
    break;
  }
}

void numNaret(int nbneel, int aretNum, int ptsNum[]){
  switch(nbneel){
  case 4 : // quadrangle
    switch(aretNum){
    case 1 :
      ptsNum[0] = 1;
      ptsNum[1] = 2;
      break;
    case 2:
      ptsNum[0] = 2;
      ptsNum[1] = 3;
      break;
    case 3:
      ptsNum[0] = 3;
      ptsNum[1] = 4;
      break;
    case 4:
      ptsNum[0] = 4;
      ptsNum[1] = 1;
      break;
    }
    break;
  case 3 : // triangle
    switch(aretNum){
    case 1:
      ptsNum[0] = 1;
      ptsNum[1] = 2;
      break;
    case 2:
      ptsNum[0] = 2;
      ptsNum[1] = 3;
      break;
    case 3:
      ptsNum[0] = 3;
      ptsNum[1] = 1;
      break;
    }
    break;
  }
}

void selectPts(int nb, int num[], float *coorEns[], float *coorSel[]){
  for(int i=0;i<nb;i++){
    coorSel[i] = coorEns[num[i]-1];
  }
}
