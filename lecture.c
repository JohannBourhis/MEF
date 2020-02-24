#include <stdio.h>
#include "maille.h"
  
int main(){
  char *ficmai = "maillage.txt";
  int typel; 
  int nbtng;
  float **pcoord;
  int ntel;
  int **pngnel;
  int nbneel;
  int nbaret;
  int **pnRefAr;

  if(lecfima(ficmai, &typel, &nbtng, &pcoord, &ntel, &pngnel, &nbneel, &nbaret, &pnRefAr)==0){

  printf("%d\n",typel);
  printf("%d\n",nbtng);
  printf("%d\n",ntel);
  printf("%d\n",nbneel);
  printf("%d\n",nbaret);

  // affichage de coord
  for(int j = 0; j < nbtng; j++){ 
    printf("%f %f \n", pcoord[j][0], pcoord[j][1]);
  }
	
// affichage de ngnel et nRefAr
  switch(typel){
    case 1 : // quadrangle
      for (int j=0; j<ntel; j++){
        printf( "%d %d %d %d", pngnel[j][0], pngnel[j][1], pngnel[j][2], pngnel[j][3]);
        printf("%d %d %d %d \n", pnRefAr[j][0], pnRefAr[j][1], pnRefAr[j][2], pnRefAr[j][3]);
      }
      break;
     case 2 : // triangle
       for (int j=0; j<ntel; j++){
         printf("%d %d %d \n", pngnel[j][0], pngnel[j][1], pngnel[j][2]);
         printf("%d %d %d \n", pnRefAr[j][0], pnRefAr[j][1], pnRefAr[j][2]);
       }
       break;
    }
  }
}
