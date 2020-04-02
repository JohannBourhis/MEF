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
  char *FichiersMaillages[6];
  int test;
  printf("Sur quel maillage exécuter le test ?\n");
  printf("1 : d1ql , 2 : d1tl , 3 : d2ql , 4 : d2tl\n");
  scanf("%d", &test);
  printf("Sur quelle fonction exécuter le test ?\n");
  printf("1 : xy(1-x)(1-y) , 2 : sin(pi*x)*sin(pi*y) , 3 : cos(pi*x)*cos(pi*y)\n");
  scanf("%d", &nucas);
  // Gestion des fichiers de maillage
  switch(test){
    case 1 :
      FichiersMaillages[0] = "maillage/d1ql_2.txt";
      FichiersMaillages[1] = "maillage/d1ql_4.txt";
      FichiersMaillages[2] = "maillage/d1ql_8.txt";
      FichiersMaillages[3] = "maillage/d1ql_16.txt";
      FichiersMaillages[4] = "maillage/d1ql_32.txt";
      FichiersMaillages[5] = "maillage/d1ql_64.txt";
      break;
    case 2 :
      FichiersMaillages[0] = "maillage/d1tl_2.txt";
      FichiersMaillages[1] = "maillage/d1tl_4.txt";
      FichiersMaillages[2] = "maillage/d1tl_8.txt";
      FichiersMaillages[3] = "maillage/d1tl_16.txt";
      FichiersMaillages[4] = "maillage/d1tl_32.txt";
      FichiersMaillages[5] = "maillage/d1tl_64.txt";
      break;
    case 3 :
      FichiersMaillages[0] = "maillage/d2ql_2.txt";
      FichiersMaillages[1] = "maillage/d2ql_4.txt";
      FichiersMaillages[2] = "maillage/d2ql_8.txt";
      FichiersMaillages[3] = "maillage/d2ql_16.txt";
      FichiersMaillages[4] = "maillage/d2ql_32.txt";
      FichiersMaillages[5] = "maillage/d2ql_64.txt";
      break;
    case 4 :
      FichiersMaillages[0] = "maillage/d2tl_2.txt";
      FichiersMaillages[1] = "maillage/d2tl_4.txt";
      FichiersMaillages[2] = "maillage/d2tl_8.txt";
      FichiersMaillages[3] = "maillage/d2tl_16.txt";
      FichiersMaillages[4] = "maillage/d2tl_32.txt";
      FichiersMaillages[5] = "maillage/d2tl_64.txt";
      break;
    default :
      printf("Problème non traité\n");
      break;
  }
  // Gestion des conditions aux limites
  int nRefDom = 0;
  int nbRefD0, nbRefD1, nbRefF1;
  int *numRefD0, *numRefD1, * numRefF1;
  switch(nucas){
    case 1 :
      if(test > 2){
        nbRefD0 = 2;
        numRefD0 = malloc(2*sizeof(int));
        numRefD0[0] = 1; numRefD0[1] = 4;
        nbRefD1 = 2;
        numRefD1 = malloc(2*sizeof(int));
        numRefD1[0] = 2; numRefD1[1] = 3;
        nbRefF1 = 0;
      }else{
        nbRefD0 = 4;
        numRefD0 = malloc(4*sizeof(int));
        numRefD0[0] = 1; numRefD0[1] = 2;
        numRefD0[2] = 3; numRefD0[3] = 4;
        nbRefD1 = 0; nbRefF1 = 0;
      }
      break;
    case 2 :
      if(test > 2){
        nbRefD0 = 2;
        numRefD0 = malloc(2*sizeof(int));
        numRefD0[0] = 1; numRefD0[1] = 4;
        nbRefD1 = 2;
        numRefD1 = malloc(2*sizeof(int));
        numRefD1[0] = 2; numRefD1[1] = 3;
        nbRefF1 = 0;
      }else{
        nbRefD0 = 4;
        numRefD0 = malloc(4*sizeof(int));
        numRefD0[0] = 1; numRefD0[1] = 2;
        numRefD0[2] = 3; numRefD0[3] = 4;
        nbRefD1 = 0; nbRefF1 = 0;
      }
      break;
    case 3 :
      if(test > 2){
        nbRefD1 = 2;
        numRefD1 = malloc(2*sizeof(int));
        numRefD1[0] = 1; numRefD1[1] = 4;
        nbRefF1 = 2;
        numRefF1 = malloc(2*sizeof(int));
        numRefF1[0] = 2; numRefF1[1] = 3;
        nbRefD0 = 0;
      }else{
        nbRefF1 = 4;
        numRefF1 = malloc(4*sizeof(int));
        numRefF1[0] = 1; numRefF1[1] = 2;
        numRefF1[2] = 3; numRefF1[3] = 4;
        nbRefD0 = 0; nbRefD1 = 0;
      }
      break;
    default :
      printf("Fonction non définie\n");
      break;
  }
  /*-----------------------------------
  -------------------------------------
  -----------------------------------*/
  // Fichier courant
  char *ficmai;
  // Données du maillage
  int typel;
  int nbaret;
  int nbneel;
  int nbtng;
  int ntel;
  float **coord;
  int **ngnel;
  int **nRefAr;
  // Données pour le stockage matrice
  int dimLmat;
  float *Matrice; int *NumCol;
  int *AdSuccLi; int *AdPrCoefLi;
  float *SecMembre; int *NumDLDir;
  float *ValDLDir;
  float *MatriceO; int *NumColO;
  int codret;
  float *MatProf; int *Profil;
  // Données pour la résolution du sysyème
  float eps = 1e-6;
  float *U; float *Uex; float *Y;
  float *ll; float *ld;
  int aff = -1;
  // Boucle sur les différents maillages testés
  for(int iter=0 ; iter<6 ; iter++){
    ficmai = FichiersMaillages[iter];
    lecfima(ficmai, &typel, &nbtng, &coord, &ntel, &ngnel,
            &nbneel, &nbaret, &nRefAr);
    // Dimensionnement des matrices pour la SMD
    dimLmat=(2+2*typel)*nbtng;
    Matrice = malloc((dimLmat+nbtng)*sizeof(float)); if(Matrice==NULL){return 1;}
    NumCol = malloc(dimLmat*sizeof(int));              if(NumCol==NULL){return 1;}
    AdSuccLi = malloc(dimLmat*sizeof(int));            if(AdSuccLi==NULL){return 1;}
    AdPrCoefLi = malloc(nbtng*sizeof(int));            if(AdPrCoefLi==NULL){return 1;}
    SecMembre = malloc(nbtng*sizeof(float));         if(SecMembre==NULL){return 1;}
    NumDLDir = malloc(nbtng*sizeof(int));              if(NumDLDir==NULL){return 1;}
    ValDLDir = malloc(nbtng*sizeof(float));          if(ValDLDir==NULL){return 1;}
    assemblage(ntel, typel, nbneel, nbaret, nRefDom, coord, nRefAr, ngnel,
               nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1,
               nbtng, Matrice, AdPrCoefLi, AdSuccLi, NumCol, SecMembre,
               NumDLDir, ValDLDir);
    /*affsmd_(&nbtng, AdPrCoefLi, NumCol, AdSuccLi, Matrice,
            SecMembre, NumDLDir, ValDLDir);*/
    // Conversion en SMO
    MatriceO = malloc((dimLmat+nbtng)*sizeof(float)); if(MatriceO==NULL){return 1;}
    NumColO = malloc(dimLmat*sizeof(int));            if(NumColO==NULL){return 1;}
    cdesse_(&nbtng, AdPrCoefLi, NumCol, AdSuccLi, Matrice, SecMembre, NumDLDir, ValDLDir,
            AdPrCoefLi, NumColO, MatriceO, SecMembre);
    /*affsmo_(&nbtng,AdPrCoefLi,NumColO,MatriceO,SecMembre);*/
    // Détermination de la taille de MatProf
    int indCol;
    int nbcprx=0;
    for(int i=0; i<nbtng-1; i++){
      if( (AdPrCoefLi[i+1]-AdPrCoefLi[i]) !=0 ){
        indCol = NumColO[AdPrCoefLi[i]-1];
        nbcprx += (i+1) - (indCol-1);
      }
    }

    MatProf = malloc((nbtng+nbcprx)*sizeof(float)); if(MatProf==NULL){return 1;}
    Profil = malloc(nbtng*sizeof(int));             if(Profil==NULL){return 1;}

    //Conversion stockage profile
    dSMOaPR(nbtng, AdPrCoefLi, NumColO, MatriceO, nbcprx,
            Profil, MatProf, &codret);
    if(codret){
      printf("la taille de MatProf est insuffisante (nbcprx)\n");
      return 1;
    }
    // Factorisation A=L*Lt
    ld = malloc(nbtng*sizeof(float));             if(ld==NULL){return 1;}
    ll = malloc(nbtng*(nbtng-1)/2*sizeof(float)); if(ll==NULL){return 1;} // matrice pleine ?
    ltlpr_(&nbtng, Profil, &MatProf[0],
           &MatProf[nbtng], &eps, ld, ll);
    // Résolution des deux systèmes linéaires
    Y = malloc(nbtng*sizeof(float)); if(Y==NULL){return 1;}
    rsprl_(&nbtng, Profil, ld, ll,
           SecMembre, Y);
    U = malloc(nbtng*sizeof(float)); if(U==NULL){return 1;}
    rspru_(&nbtng, Profil, ld, ll,
           Y, U);
    Uex = malloc(nbtng*sizeof(float)); if(Uex==NULL){return 1;}
    //Calcul de la soltuion exacte
    Calsol(nbtng, coord, Uex);
    //Affichage de la solution
    affsol_(&nbtng, coord[0], U, Uex, &aff);
    //Libération mémoire
    freetab(coord); freetab(ngnel); freetab(nRefAr);
    free(Matrice); free(NumCol); free(AdSuccLi);
    free(AdPrCoefLi); free(SecMembre); free(NumDLDir);
    free(ValDLDir); free(MatriceO); free(NumColO);
    free(MatProf); free(Profil);
    free(U); free(Uex); free(Y);
    free(ll); free(ld);
  }
}
