#include <stdio.h>
#include "maille.h"

void etiqAr(int typel, int n1, int n2, int nrefdom, const int *nrefcot, int nbtel, int nbaret, int **nRefAr){
  /* n1 : subdivision en x
  n2 : subdivision en y
  typel = 2 (triangle) ou 1 (quadrangle)
  nrefdom : num ref de l'interieur
  nrefcot : num ref des 4 cotes
  nbtel : nombre d'elements
  nbaret : nombre d'aretes
  nRefAr : tableau des r_i(k)*/
  int N = typel*(n1-1);
  int M = (n2-1);
  int NM = N*M;
  for(int i = 0; i<nbtel; i++){
    for(int j = 0; j<nbaret; j++){
      nRefAr[i][j] = nrefdom;
    }
  } 
  switch(typel){
    case 1 :
      /* lower side */
      for(int i = 0; i < N; i+=typel){
        nRefAr[i][0] = nrefcot[0];
      }
      /* right side */
      for(int i = N-1; i < NM; i+=N){	     
        nRefAr[i][1] = nrefcot[1];
      }
      
      /* upper side */
      for(int i = M*(n1-2)*typel; i < NM; i+=typel){
        nRefAr[i][2] = nrefcot[2];
      }
      /* left side */
      for(int i = 0; i < NM; i+=N){
        nRefAr[i][3] = nrefcot[3];
      }
    break;
    case 2 :
      /* lower side */
      for(int i = 0; i < N; i+=typel){
        nRefAr[i][0] = nrefcot[0];
      }
      /* right side */
      for(int i = N - 1; i < NM; i+=N){	     
        nRefAr[i][2] = nrefcot[1];
      }
      /* upper side */
      for(int i = M*(n1-2)*typel + 1; i < NM; i+=typel){
        nRefAr[i][1] = nrefcot[2];
      }
      /* left side */
      for(int i = 0; i < NM; i+=N){
        nRefAr[i][2] = nrefcot[3];
      }
    break;
  }
}

int lecfima(char *ficmai, int *typel, int *nbtng, float ***pcoord,
  int *nbtel, int ***pngnel, int *nbneel, int *nbaret, int ***pnRefAr){
  
  
  FILE* maillage = NULL;

  if ((maillage = fopen(ficmai,"r+")) == NULL){
    printf("Le fichier ne peut pas etre lu\n");
    return 1;
  }
  else {
    int j;
    
    // remplissage de coord
    fscanf(maillage, "%d", nbtng);      
    float **coord=alloctabf(*nbtng,2);
    if(coord == NULL){return 2;}
    for(j = 0; j < *nbtng; j++){ 
      fscanf(maillage, "%f %f", &coord[j][0], &coord[j][1]);
    }
    
    // remplisage de ngnel et nRefAr
    fscanf(maillage, "%d %d %d %d", nbtel, typel, nbneel, nbaret);  
    int **ngnel=alloctab(*nbtel,*nbneel);
    if(ngnel == NULL){return 2;}
    int **nRefAr=alloctab(*nbtel,*nbaret); 
    if(nRefAr == NULL){return 2;}
    switch(*typel)  {
      case 1 : // quadrangle
      for (j=0; j<*nbtel; j++){
        fscanf(maillage, "%d %d %d %d", &ngnel[j][0], &ngnel[j][1], &ngnel[j][2], &ngnel[j][3]);
	fscanf(maillage, "%d %d %d %d", &nRefAr[j][0], &nRefAr[j][1], &nRefAr[j][2], &nRefAr[j][3]);

	
      }
      break;
      case 2 : // triangle
      for (j=0; j<*nbtel; j++){
        fscanf(maillage, "%d %d %d", &ngnel[j][0], &ngnel[j][1], &ngnel[j][2]);
        fscanf(maillage, "%d %d %d", &nRefAr[j][0], &nRefAr[j][1], &nRefAr[j][2]);
      }
      break;
    }
  *pcoord = coord;
  *pnRefAr = nRefAr;
  *pngnel = ngnel;
  }
  return 0;
}
