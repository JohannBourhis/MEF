#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"

void invertM2x2(float **M, float *det, float **M_inv){
  /* Calcul de l'inverse d'une matrice M et de son déterminant */
  *det=M[0][0]*(M[1][1])-M[0][1]*(M[1][0]);
  float alpha = (1./(*det));
  M_inv[0][0]=alpha*M[1][1];
  M_inv[0][1]=-alpha*M[0][1];
  M_inv[1][0]=-alpha*M[1][0];
  M_inv[1][1]=alpha*M[0][0];
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

void calDerFbase(int typel, float *x, float **w){
  /* 
  On renvoie un vecteur de taille n qui correspond 
  aux valeurs des gradients fonctions de base en un point x
  et on a : w_k(s_i)=delta_ki pour une sommet si donné
  */
  switch(typel){
  case 1 : // quadrangle 
    /*  4  3
        1  2 */
    w[0][0]=-1+x[1]; w[0][1]=-1+x[0]; 
    w[1][0]=1-x[1];  w[1][1]=-x[0];
    w[2][0]=-x[1];   w[2][1]=1-x[0];
    w[3][0]=x[1];    w[3][1]=x[0];
  break;
  case 2 : // triangle 	
    w[0][0]=-1; w[0][1]=-1; 
    w[1][0]=1;  w[1][1]=0;
    w[2][0]=0;  w[2][1]=1;
  break;
  case 3 : // segment
    /*  2  1    */		
    w[0][0]=1;
    w[0][1]=-1;
  break;
  }
}

void calFbase(int typel, float *x, float *w){
  /* 
  On renvoie un vecteur de taille n qui correspond 
  aux valeurs des fonctions de base en un point x
  et on a : w_k(s_i)=delta_ki pour un sommet si donné
  */
  switch(typel) {
  case 1 : // quadrangle 
  /*  4  3
      1  2      */
    w[0]=1-x[0]-x[1]+x[0]*x[1];
    w[1]=x[0]*(1-x[1]);
    w[2]=x[0]*x[1];
    w[3]=x[1]*(1-x[0]);
  break;
  case 2 : // triangle 
  /*  3  
      1  2      */	
    w[0]=1-x[0]-x[1];
    w[1]=x[0];
    w[2]=x[1];
  break;
  case 3 : // segment
  /*  2  1    */			
    w[0]=x[0];
    w[1]=1-x[0];
  break;
  }
}

void numNaret(int typel, int aretNum, int ptsNum[]){
  switch(typel){
  case 1:
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
      ptsNum[1] = 4;
    break;
    case 4:
      ptsNum[0] = 4;
      ptsNum[1] = 1;
    break;
    }
  break;
  case 2:
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
