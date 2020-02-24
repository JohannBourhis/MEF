#include <stdio.h>
#include <stdlib.h>
#include "maille.h"

void transFk(int t, float **S, float **Y, float *foncbase){
  /* Calcul de Y, l'image par la transformation Fk à partir des sommets S
  et des valeurs fournies par calFbase */
  float *res;
  res[0] = 0;
  res[1] = 0;
  for(int i=0; i<5-t; i++){
    for(int j=0;j<2;j++){
      res[j] += foncbase[i]*S[i][j];
    }
  }
  *Y = res;
}
