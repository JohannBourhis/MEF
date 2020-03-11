#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"
#include "int.h"
#include "maille.h"
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
  float **coord, int **nRefAr, int **ngnel,
  int nbRefD0, int *numRefD0, int nbRefD1, int *numRefD1, int nbRefF1, int *numRefF1,
  int *NbLign, int *NbCoef, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
  int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir){
    int res;
	int dimDiag = ntel;
	int dimLmat = dimDiag*(typel+1);
	Matrice = malloc((dimLmat+dimDiag)*sizeof(float));
	NumCol=malloc(dimLmat*sizeof(float));
	AdSuccLi=malloc(dimLmat*sizeof(float));
	AdPrCoefLi=malloc((dimDiag-1)*sizeof(float));
	SecMembre=malloc(dimDiag*sizeof(float));
	NumDLDir=malloc(dimDiag*sizeof(int));
	ValDLDir=malloc(dimDiag*sizeof(float));
	float *DiagMat=&Matrice[0];
	float *LowMat=&Matrice[dimDiag];
	// Initialisation des vecteurs
	for(int i=0;i<dimDiag;i++){
	  DiagMat[i]=0;
	  SecMembre[i]=0;
	}
	// Variables nécessaires aux calculs élémentaires
	int nextad=1;
	int I,J,L;
    int num[nbneel];
    float **coorEl = alloctabf(nbneel,2);
    float **MatElem=alloctabf(nbneel,nbneel);
    float SMbrElem[nbneel];
    int NuDElem[nbneel];
    float uDElem[nbneel]; 
    // Boucle sur les elements
    for(int k=0;k<ntel;k++){
	  for(int j=0;j<nbneel;j++){
		num[j]=ngnel[k][j];
	  }
	selectPts(nbneel,num,coord,coorEl);
	res = cal1Elem(nbneel, nbaret, nRefDom, coorEl, nRefAr[k], nbRefD0,
	  numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1,
	  MatElem, SMbrElem, NuDElem, uDElem);
	if(res){return res;}
	//impCalEl(i+1, typel, nbneel, MatElem, SMbrElem, NuDElem, uDElem); // Affichage pour vérification
	for(int i=0;i<nbneel;i++){
	  I=ngnel[k][i];
	  DiagMat[I-1]+=MatElem[i][i];
	  SecMembre[I-1]+=SMbrElem[i];
	  NumDLDir[I-1]=NuDElem[i]*I;
      ValDLDir[I-1]=uDElem[i];
	  for(int j=0;j<i;j++){
		J=ngnel[k][j];
		if(I<J){ L=I; I=J; J=L;}  //symétrie : I=max(I,J), J=min(I,J)
		assmat_(&I, &J, &MatElem[i][j], AdPrCoefLi, NumCol, AdSuccLi, LowMat, &nextad);
	  }
    }
  }
  return 0;
}

/*----------------------------------------------
 Procédure qui sauvegarde la SMD de A
 *
 *
 *
 ----------------------------------------------*/
 void EcrSMD(int *NbCoef, int *NbLign, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
  int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir){
          FILE* SMD;
      if((SMD = fopen("SMD.txt", "w")) != NULL){
                fwrite(NbLign, sizeof(int), 1, SMD);
                fwrite(SecMembre, sizeof(float), *NbLign, SMD);
                fwrite(NumDLDir, sizeof(int), *NbLign, SMD);
                fwrite(ValDLDir, sizeof(float), *NbLign, SMD);
                fwrite(AdPrCoefLi, sizeof(int), *NbLign, SMD);
                fwrite(Matrice, sizeof(float), *NbLign+*NbCoef, SMD);
                fwrite(NumCol, sizeof(float), *NbCoef, SMD);
                fwrite(AdSuccLi, sizeof(int), *NbCoef, SMD);
      }
      fclose(SMD);
 } 

 void LecSMD(int *NbCoef, int *NbLign, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
  int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir){
    FILE* SMD;
    if((SMD = fopen("SMD.txt", "r")) != NULL){
                fread(NbLign, sizeof(int), 1, SMD);
                fread(SecMembre, sizeof(float), *NbLign, SMD);
                fread(NumDLDir, sizeof(int), *NbLign, SMD);
                fread(ValDLDir, sizeof(float), *NbLign, SMD);
                fread(AdPrCoefLi, sizeof(int), *NbLign, SMD);
                *NbCoef=AdPrCoefLi[*NbLign-1]-1;
                fread(Matrice, sizeof(float), *NbLign+*NbCoef, SMD);
                fread(NumCol, sizeof(float), *NbCoef, SMD);
                fread(AdSuccLi, sizeof(int), *NbCoef, SMD);
    } 
    affsmd_(NbLign,AdPrCoefLi,NumCol,AdSuccLi,Matrice,SecMembre,NumDLDir,ValDLDir);
    fclose(SMD);
}

