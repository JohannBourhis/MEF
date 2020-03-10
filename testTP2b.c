#include <stdlib.h>
#include <stdio.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
 // gcc * -lm   pour utiliser math.h
int main(){
  int typel;
  int nbaret;
  int nbneel;
  int res;
  int nbtng;
  int ntel; 
  float **coord;
  int **ngnel;
  int **nRefAr;
   
  char *ficmai = "car1x1t_1";

  lecfima(ficmai, &typel, &nbtng, &coord, &ntel, &ngnel, &nbneel, &nbaret, &nRefAr);
  
  float **coorEl = alloctabf(nbneel,2);

// Numero des ref des aretes
  int *nRefArEl=malloc(nbneel*sizeof(int));
  int nRefDom=0;
  // Dirichlet
  int nbRefD0=1;
  int numRefD0[nbRefD0];
  numRefD0[0]=1;
  // Dirichlet non-homogène
  int nbRefD1=1;
  int numRefD1[nbRefD1];
  numRefD1[0]=4;
  // Neumann
  int nbRefF1=2;
  int numRefF1[nbRefF1];
  numRefF1[0]=2; numRefF1[1]=3;
  float **MatElem;
  MatElem=alloctabf(nbneel,nbneel);
  float *SMbrElem=malloc(nbneel*sizeof(float));
  int *NuDElem=malloc(nbneel*sizeof(int));
  float *uDElem=malloc(nbneel*sizeof(float));
  
// Boucle sur les elements
  for(int i=0;i<ntel;i++){
    for(int j=0;j<nbneel;j++){
      coorEl[j][0]=coord[ngnel[i][j]-1][0];
      coorEl[j][1]=coord[ngnel[i][j]-1][1];
      nRefArEl[j]=nRefAr[i][j];
    }
  res = cal1Elem(nbneel, nbaret, nRefDom, coorEl, nRefArEl, nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1, MatElem, SMbrElem, NuDElem, uDElem);
  printf("res = %d \n",res);
  
  impCalEl(i+1, typel, nbneel, MatElem, SMbrElem, NuDElem, uDElem);
  }
}
