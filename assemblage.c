#include <stdlib.h>
#include <stdio.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
#include "syslin.h"
#include "forfun.h"
/* --------------------------------------
 * assemblage de la matrice
Arguments d'entrée :
  ntel     : nombre d'elements
  typel    : type d'element
  nbneel   : nombre de noeuds de l'element
  nbaret   : nombre d'aretes de l'element
  coord    : tableau contenant les sommets j de l'élément k (coord[k][j])
  nRefAr   : tableau contenant les num de ref des nbaret j de chaque élément i -> nRefAr[i][j]
  ngnel    : tableau contenant les num de ref globaux des noeuds j de chaque élément
  nRefDom  : etiquette coorespondant aux aretes internes
  nbRefD0  : taille de numRefD0
  numRefDO : numéro locaux des arretes coorespondant a une condition de Dirichlet  u=O
  nbRefD1  : taille de numRefD1
  numRefD1 : numéros locaux des arretes coorespondant a une condition de Dirichlet non homogène   u=ud(x)
  nbRefF1  : taille de numRefF1
  numRefF1 : numéros locaux des arretes coorespond a une condition de Neuman ou Fourier        eq diff
Arguments de sortie :
  AdPrCoefLi : adresses des 1er coeff. des lignes de la partie triang. inf. stricte
  NumCol     : numeros de colonnes des coeff. de la partie triang. inf. stricte
  AdSuccLi   : adresse du coeff. suivant sur la meme ligne
  Matrice    : tableau contenant les éléments diagonaux puis la partie triangulaire inferieure stricte
  NbLign     : nombre de lignes de la matrice A
  NbCoef     : nimbre de coefficients mémorisés de la matrice A
  SecMembre  : tableau contenant le second membre
  NumDLDir   : tableau de description des conditions de Dirichlet
  ValDLDir   : valeur de blocage aux noeuds Dirichlet non homogenes
 ------------------------------------------- */
int assemblage(int ntel, int typel, int nbneel, int nbaret, int nRefDom,
               float **coord, int **nRefAr, int **ngnel, int nbRefD0,
               int *numRefD0, int nbRefD1, int *numRefD1, int nbRefF1,
               int *numRefF1, int *NbLign, float *Matrice,
               int *AdPrCoefLi, int *AdSuccLi, int *NumCol, float *SecMembre,
               int *NumDLDir, float *ValDLDir){
  int i, j, k, res;
  float *DiagMat = &Matrice[0];
  float *LowMat = &Matrice[*NbLign];
  // Variables nécessaires aux calculs élémentaires
  int nextad = 1;
  int I, Il, J;
  int noeuds[nbneel];
  float **coorEl = alloctabf(nbneel, 2); if(coorEl==NULL){return 1;}
  float **MatElem = alloctabf(nbneel, nbneel); if(MatElem==NULL){return 1;}
  float SMbrElem[nbneel];
  int NuDElem[nbneel];
  float uDElem[nbneel];
  // Initialisation des vecteurs
  for(i=0 ; i<*NbLign ; i++){
    DiagMat[i] = 0;
    SecMembre[i] = 0;
    NumDLDir[i] = i+1;
    AdPrCoefLi[i] = 0;
  }
  // Boucle sur les elements
  for(k=0 ; k<ntel ; k++){
    // sélections des noeuds de l'éléments
    for(j=0 ; j<nbneel ; j++){
      noeuds[j] = ngnel[k][j];
    }
    selectPts(nbneel,noeuds,coord,coorEl);
    // calculs élémentaires
    res = cal1Elem(nbneel, nbaret, nRefDom, coorEl, nRefAr[k], nbRefD0,
          numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1,
          MatElem, SMbrElem, NuDElem, uDElem);
    if(res){return res;}
    // assemblage de la matrice
    for(i=0 ; i<nbneel ; i++){
      I = ngnel[k][i];
      Il = I-1;
      DiagMat[Il] += MatElem[i][i];
      SecMembre[Il] += SMbrElem[i];
      if(NumDLDir[Il]>0){
        NumDLDir[Il] *= NuDElem[i];
        ValDLDir[Il] = uDElem[i];
      }
      for(j=0;j<i;j++){
        J = ngnel[k][j];
        if(I<J){  //symétrie : I=max(I,J), J=min(I,J)
          assmat_(&J, &I, &MatElem[i][j], AdPrCoefLi, NumCol, AdSuccLi, LowMat, &nextad);
        }else{
          assmat_(&I, &J, &MatElem[i][j], AdPrCoefLi, NumCol, AdSuccLi, LowMat, &nextad);
        }
      }
    }
  }
  AdPrCoefLi[*NbLign-1] = nextad;
  return 0;
}

/*---------------------------------------------------
 Procédure qui sauvegarde la SMD de A dans un fichier
 *
 *
 *
 ---------------------------------------------------*/
 void EcrSMD(int *NbLign, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
            int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir){
  FILE* SMD;
  if((SMD = fopen("SMD.txt", "w")) != NULL){
    fwrite(NbLign, sizeof(int), 1, SMD);
    fwrite(SecMembre, sizeof(float), *NbLign, SMD);
    fwrite(NumDLDir, sizeof(int), *NbLign, SMD);
    fwrite(ValDLDir, sizeof(float), *NbLign, SMD);
    fwrite(AdPrCoefLi, sizeof(int), *NbLign, SMD);

    printf("AdPrCoefLi :(EcrSMD)\n");
    for(int i=0; i<*NbLign ; i++){
      printf("%d  %d\n", i, AdPrCoefLi[i]);
    }

    int NbCoef = AdPrCoefLi[*NbLign-1]-1;
    fwrite(Matrice, sizeof(float), *NbLign+NbCoef, SMD);
    fwrite(NumCol, sizeof(int), NbCoef, SMD);
    fwrite(AdSuccLi, sizeof(int), NbCoef, SMD);
  }
  fclose(SMD);
}

/*--------------------------------------------------
 Procédure qui lit la SMD de A à partir d'un fichier
 *
 *
 *
 --------------------------------------------------*/
int LecSMD(int *NbLign, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
            int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir){
  FILE* SMD;
  if((SMD = fopen("SMD.txt", "r")) != NULL){
    printf("NbLign : %d (LecSMD)\n", *NbLign);
    fread(NbLign, sizeof(int), 1, SMD);
  	// lecture
  	fread(SecMembre, sizeof(float), *NbLign, SMD);
  	fread(NumDLDir, sizeof(int), *NbLign, SMD);
  	fread(ValDLDir, sizeof(float), *NbLign, SMD);
  	fread(AdPrCoefLi, sizeof(int), *NbLign, SMD);
    printf("AdPrCoefLi :(LecSMD)\n");
    for(int i=0; i<*NbLign ; i++){
      printf("%d  %d\n", i, AdPrCoefLi[i]);
    }
  	int NbCoef=AdPrCoefLi[*NbLign-1]-1;
    printf("NbCoef : %d\n", NbCoef);
  	fread(Matrice, sizeof(float), *NbLign+NbCoef, SMD);
  	fread(NumCol, sizeof(int), NbCoef, SMD);
  	fread(AdSuccLi, sizeof(int), NbCoef, SMD);
  	//transmission des tableaux en dehors de la fonction
    printf("affichage SMD (LecSMD)\n");
    affsmd_(NbLign, AdPrCoefLi, NumCol, AdSuccLi, Matrice,
              SecMembre, NumDLDir, ValDLDir);
    fclose(SMD);
  }
  else {
    printf("Erreur ouverture du fichier SMD\n");
    return 2;
  }
  return 0;
}
