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

  char *ficmai = "car1x1t_1";
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
  int NbCoef[1],NbLign[1], *AdPrCoefLi, *AdSuccLi, *NumDLDir, *NumCol;
  float *Matrice, *SecMembre, *ValDLDir;
  NbLign[0]=nbtng;
  //allocation
  int dimDiag = NbLign[0];
  int dimLmat = dimDiag*(typel+1);
  Matrice=malloc((dimLmat+dimDiag)*sizeof(float)); if(Matrice==NULL){return 1;}
  NumCol=malloc(dimLmat*sizeof(int)); if(NumCol==NULL){return 1;}
  AdSuccLi=malloc(dimLmat*sizeof(float)); if(AdSuccLi==NULL){return 1;}
  AdPrCoefLi=malloc(dimDiag*sizeof(float)); if(AdPrCoefLi==NULL){return 1;}
  SecMembre=malloc(dimDiag*sizeof(float)); if(SecMembre==NULL){return 1;}
  NumDLDir=malloc(dimDiag*sizeof(int)); if(NumDLDir==NULL){return 1;}
  ValDLDir=malloc(dimDiag*sizeof(float)); if(ValDLDir==NULL){return 1;}

  //tests des fonctions
  printf("assemblage\n");
  assemblage(ntel, typel, nbneel, nbaret, nRefDom, coord, nRefAr, ngnel,
             nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1,
             NbLign, NbCoef, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
             NumDLDir, ValDLDir);

  printf("EcrSMD\n");
  EcrSMD(NbCoef, NbLign, Matrice, AdPrCoefLi, AdSuccLi,
         NumCol, SecMembre, NumDLDir, ValDLDir);

  float *MatriceO = malloc(*NbCoef*sizeof(float));
  int *NumColO = malloc(*NbCoef*sizeof(int));

  printf("Assemblage de SMO");
  dSMDaSMO(NbLign, SecMembre, AdPrCoefLi, MatriceO, NumColO);
  printf("LecSMO\n");
  int LNbCoef, LNbLign, *LAdPrCoefLiO, *LNumColO;
  float *LMatriceO, *LSecMembreO;
  LecSMO(&LNbCoef, &LNbLign, &LMatriceO, &LAdPrCoefLiO,
         &LNumColO, &LSecMembreO);
}
