#include <stdlib.h>
#include <stdio.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
#include "assemblage.h"
#include "syslin.h"
#include "forfun.h"
#include "resol.h"

int nucas;

int main(){
  nucas = 3;
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
  char *ficmai = "maillage/d1tl_64.txt";
  printf("lecfima\n");
  lecfima(ficmai, &typel, &nbtng, &coord, &ntel, &ngnel,
          &nbneel, &nbaret, &nRefAr);
  // Numero des ref des aretes
  int nRefDom=0;
  // Dirichlet
  int nbRefD0=4;
  int numRefD0[nbRefD0];
  numRefD0[0]=1; numRefD0[1]=2; numRefD0[2]=3; numRefD0[3]=4;
  // Dirichlet non-homogène
  int nbRefD1=4;
  int numRefD1[nbRefD1];
  numRefD1[0]=1; numRefD1[1]=2; numRefD1[2]=3; numRefD1[3]=4;
  //numRefD1[0]=4;
  // Neumann
  int nbRefF1=0;
  int numRefF1[nbRefF1];

  switch(nucas){
    case 3:
      nbRefD0=0;
      break;
    default:
      nbRefD1=0;
      break;
  }

  printf("declaration des tableaux de la SMD\n");
  //déclaration des tableaux de la SMD
  int dimLmat;
  dimLmat=(2+2*typel)*nbtng;

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
             nbtng, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
             NumDLDir, ValDLDir);
  printf("EcrSMD\n");
  EcrSMD(nbtng, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
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
  dSMDaSMO(MatriceO, NumColO, Matrice, SecMembre,
           AdPrCoefLi, AdSuccLi, ValDLDir, NumDLDir, NumCol);

  // Réinitialisation (car tout est enregistré dans le fichier SMO)
  free(Matrice); free(AdSuccLi); free(AdPrCoefLi);
  free(NumCol); free(SecMembre); free(NumDLDir); free(ValDLDir);
  free(MatriceO); free(NumColO);
  MatriceO = malloc((dimLmat+nbtng)*sizeof(float)); if(MatriceO==NULL){return 1;}
  NumColO = malloc(dimLmat*sizeof(int));            if(NumColO==NULL){return 1;}
  AdPrCoefLi = malloc(nbtng*sizeof(int));           if(AdPrCoefLi==NULL){return 1;}
  SecMembre = malloc(nbtng*sizeof(float));          if(SecMembre==NULL){return 1;}
  printf("LecSMO\n");
  LecSMO(&nbtng, MatriceO, AdPrCoefLi, NumColO, SecMembre);
  printf("affichage SMO : (testTP3)\n");
  affsmo_(&nbtng,AdPrCoefLi,NumColO,MatriceO,SecMembre);

  // Détermination de la taille de MatProf
  int indCol;
  int nbcprx=0;
  for(int i=0; i<nbtng-1; i++){
    if( (AdPrCoefLi[i+1]-AdPrCoefLi[i]) !=0 ){
      indCol = NumColO[AdPrCoefLi[i]-1];
      nbcprx += (i+1) - (indCol-1);
    }
  }
  printf("nbcprx = %d\n", nbcprx);

  float* MatProf = malloc((nbtng+nbcprx)*sizeof(float)); if(MatProf==NULL){return 1;}
  int* Profil = malloc(nbtng*sizeof(int));               if(Profil==NULL){return 1;}
  int codret;
  printf("Conversion stockage profile \n");
  dSMOaPR(nbtng, AdPrCoefLi, NumColO, MatriceO, nbcprx,
          Profil, MatProf, &codret);
  if(codret){
    printf("la taille de MatProf est insuffisante (nbcprx)\n");
    return 1;
  }

  dimLmat = Profil[nbtng-1];
  float eps = 1e-6;
  float *ad = &MatProf[0];
  float *al = &MatProf[nbtng];
  float *ld = malloc(nbtng*sizeof(float));             if(ld==NULL){return 1;}
  float *ll = malloc(nbtng*(nbtng-1)/2*sizeof(float)); if(ll==NULL){return 1;} // matrice pleine ?

  printf("Factorisation A=L*Lt/n");
  ltlpr_(&nbtng, Profil, ad, al, &eps, ld, ll);
  printf("Résolution du premier système triangulaire\n");
  float* Y = malloc(nbtng*sizeof(float)); if(Y==NULL){return 1;}
  rsprl_(&nbtng, Profil, ld, ll, SecMembre, Y);
  printf("Résolution du second système triangulaire\n");
  float* U = malloc(nbtng*sizeof(float)); if(U==NULL){return 1;}
  rspru_(&nbtng, Profil, ld, ll, Y, U);
  printf("Calcul de la solution théorique\n");
  float* Uex = malloc(nbtng*sizeof(float));   if(Uex==NULL){return 1;}
  //Calcul de la soltuion exacte
  Calsol(nbtng, coord, Uex);
  printf("Calcul de l'erreur\n");
  //Affichage de la solution
  int impfch=-1;
  affsol_(&nbtng, coord[0], U, Uex, &impfch);
  printf("End\n");
}
