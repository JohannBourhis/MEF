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

  char *ficmai = "car3x3t_3";
  printf("lecfima\n");
  lecfima(ficmai, &typel, &nbtng, &coord, &ntel, &ngnel, &nbneel, &nbaret, &nRefAr);
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
  int NbLign[1];
  NbLign[0]=nbtng;
  //allocation
  int dimLmat = (*NbLign)*(typel+1);
  float *Matrice=malloc((dimLmat+*NbLign)*sizeof(float)); if(Matrice==NULL){return 1;}
  int *NumCol=malloc(dimLmat*sizeof(int)); if(NumCol==NULL){return 1;}
  int *AdSuccLi=malloc(dimLmat*sizeof(int)); if(AdSuccLi==NULL){return 1;}
  int *AdPrCoefLi=malloc(*NbLign*sizeof(int)); if(AdPrCoefLi==NULL){return 1;}
  float *SecMembre=malloc(*NbLign*sizeof(float)); if(SecMembre==NULL){return 1;}
  int *NumDLDir=malloc(*NbLign*sizeof(int)); if(NumDLDir==NULL){return 1;}
  float *ValDLDir=malloc(*NbLign*sizeof(float)); if(ValDLDir==NULL){return 1;}

  //tests des fonctions
  printf("assemblage\n");
  assemblage(ntel, typel, nbneel, nbaret, nRefDom, coord, nRefAr, ngnel,
             nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1,
             NbLign, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
             NumDLDir, ValDLDir);
  printf("affichage après assemblage (en dehors de LecSMD)\n");
  affsmd_(NbLign, AdPrCoefLi, NumCol, AdSuccLi, Matrice,
          SecMembre, NumDLDir, ValDLDir);
  printf("EcrSMD\n");
  EcrSMD(NbLign, Matrice, AdPrCoefLi, AdSuccLi,
         NumCol, SecMembre, NumDLDir, ValDLDir);

  int NbCoef = AdPrCoefLi[*NbLign-1]-1;
  float *MatriceO = malloc((NbCoef+*NbLign)*sizeof(float));
  int *NumColO = malloc(NbCoef*sizeof(int));

  printf("Assemblage de SMO\n");
  dSMDaSMO(NbLign, SecMembre, AdPrCoefLi, MatriceO, NumColO);
  printf("LecSMO\n");
  int LNbLign, *LAdPrCoefLiO, *LNumColO;
  float *LMatriceO, *LSecMembreO;
  LecSMO(&LNbLign, &LMatriceO, &LAdPrCoefLiO,
         &LNumColO, &LSecMembreO);
}
