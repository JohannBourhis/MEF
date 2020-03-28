#include <stdlib.h>
#include <stdio.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
#include "assemblage.h"
#include "syslin.h"
#include "forfun.h"

 // gcc * -lm   pour utiliser math.h
int main(){
  printf("start\n");
  int typel;
  int nbaret;
  int nbneel;
  int nbtng;
  int ntel;
  float **coord;
  int **ngnel;
  int **nRefAr;
  //lecture du fichier de maillage
  char *ficmai = "car3x3t_3";
  printf("lecfima\n");
  lecfima(ficmai, &typel, &nbtng, &coord, &ntel, &ngnel,
          &nbneel, &nbaret, &nRefAr);
  // Numero des ref des aretes
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
  printf("declaration des tableaux de la SMD\n");
  //déclaration des tableaux de la SMD
  int dimLmat;  
  dimLmat=(2+2*typel)*nbtng;
  /*switch(typel){
    case 1:
      dimLmat = 4*nbtng;
      break;
    case 2:
      dimLmat = 6*nbtng;
      break;
  } */
  float *Matrice = malloc((dimLmat+nbtng)*sizeof(float)); if(Matrice==NULL){return 1;}
  int *NumCol = malloc(dimLmat*sizeof(int));              if(NumCol==NULL){return 1;}
  int *AdSuccLi = malloc(dimLmat*sizeof(int));            if(AdSuccLi==NULL){return 1;}
  int *AdPrCoefLi = malloc(nbtng*sizeof(int));            if(AdPrCoefLi==NULL){return 1;}
  float *SecMembre = malloc(nbtng*sizeof(float));         if(SecMembre==NULL){return 1;}
  int *NumDLDir = malloc(nbtng*sizeof(int));              if(NumDLDir==NULL){return 1;}
  float *ValDLDir = malloc(nbtng*sizeof(float));          if(ValDLDir==NULL){return 1;}
  printf("assemblage\n");
  assemblage(ntel, typel, nbneel, nbaret, nRefDom, coord, nRefAr, ngnel,
             nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1,
             &nbtng, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
             NumDLDir, ValDLDir);
  printf("EcrSMD\n");
  EcrSMD(&nbtng, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
         NumDLDir, ValDLDir);
  // Réinitialisation (car tout est enregistré dans le fichier SMD)
  free(Matrice); free(AdSuccLi); free(AdPrCoefLi);
  free(NumCol); free(SecMembre); free(NumDLDir); free(ValDLDir);
  Matrice = malloc((dimLmat+nbtng)*sizeof(float));         if(Matrice==NULL){return 1;}
  NumCol = malloc(dimLmat*sizeof(int));                    if(NumCol==NULL){return 1;}
  AdSuccLi = malloc(dimLmat*sizeof(int));                  if(AdSuccLi==NULL){return 1;}
  AdPrCoefLi = malloc(nbtng*sizeof(int));                  if(AdPrCoefLi==NULL){return 1;}
  SecMembre = malloc(nbtng*sizeof(float));                 if(SecMembre==NULL){return 1;}
  NumDLDir = malloc(nbtng*sizeof(int));                    if(NumDLDir==NULL){return 1;}
  ValDLDir = malloc(nbtng*sizeof(float));                  if(ValDLDir==NULL){return 1;}
  float *MatriceO = malloc((dimLmat+nbtng)*sizeof(float)); if(MatriceO==NULL){return 1;}
  int *NumColO = malloc(dimLmat*sizeof(int));              if(NumColO==NULL){return 1;}
  printf("Assemblage de SMO\n");
  dSMDaSMO(&nbtng, MatriceO, NumColO, Matrice, SecMembre,
           AdPrCoefLi, AdSuccLi, ValDLDir, NumDLDir, NumCol);
  printf("LecSMO\n");
  LecSMO(&nbtng, MatriceO, AdPrCoefLi, NumColO, SecMembre);
  printf("End\n");
}
