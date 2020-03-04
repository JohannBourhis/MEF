#include <stdlib.h>
#include <stdio.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
 // gcc * -lm   pour utiliser math.h
int main(){
  int typel=1;
  int nbaret=4;
  int res;
  int nbneel=4;
  int nRefDom=0;
  int nbRefD0=1;
  int numRefD0[nbRefD0];
  numRefD0[0]=1;
  int nbRefD1=1;
  int numRefD1[nbRefD1];
  numRefD1[0]=4;
  int nbRefF1=2;
  int numRefF1[nbRefF1];
  numRefF1[0]=2; numRefF1[1]=3;
  float **MatElem;
  MatElem=alloctabf(nbneel,nbneel);
  float *SMbrElem=malloc(nbneel*sizeof(float));
  int *NuDElem=malloc(nbneel*sizeof(int));
  float *uDElem=malloc(nbneel*sizeof(float));

  // Sommets
  float **coorEl=alloctabf(nbneel,2); 
  coorEl[0][0]=0; coorEl[0][1]=0;
  coorEl[1][0]=0; coorEl[1][1]=1;
  coorEl[2][0]=1; coorEl[2][1]=1;
  coorEl[3][0]=1; coorEl[3][1]=0;
  
  int *nRefArEl=malloc(nbneel*sizeof(int));
  nRefArEl[0]=1;nRefArEl[1]=2;nRefArEl[2]=3;nRefArEl[3]=4;

  res = cal1Elem(nbneel, nbaret, nRefDom, coorEl, nRefArEl, nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1, MatElem, SMbrElem, NuDElem, uDElem);
  
  printf("res = %d \n",res);
  
  impCalEl(1, typel, nbneel, MatElem, SMbrElem, NuDElem, uDElem);
}
