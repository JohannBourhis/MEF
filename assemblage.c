#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"
#include "int.h"
/* --------------------------------------
Calcul elementaire sur l'element courant K :
Arguments d'entrée :
  N1 : nb de points en x
  N2 : nb de points en y
  typel    : type d'element
  nbneel   : nombre de noeuds de l'element
  nbaret   : nombre d'aretes de l'element
  coorEl   : tableau contenant les sommets de l'element courant K (à initialiser avec la fct selectPts à partir du fichier maillage)
  nRefArEl : tableau contenant les num de ref des nbaret de l'element
             int *nRefAr[nbaret]
  nRefDom  : etiquette coorespondant aux aretes internes
  nbRefD0  : taille de numRefD0
  numRefDO : numéro locaux des arretes coorespondant a une condition de Dirichlet  u=O
  nbRefD1  : taille de numRefD1
  numRefD1 : numéros locaux des arretes coorespondant a une condition de Dirichlet non homogène   u=ud(x)
  nbRefF1  : taille de numRefF1
  numRefF1 : numéros locaux des arretes coorespond a une condition de Neuman ou Fourier        eq diff
Arguments de sortie :
  AdPrCoefLi : adresses des 1er coeff. des lignes de la partie triang. inf. stricte
  NumCol : numeros de colonnes des coeff. de la partie triang. inf. stricte
  AdSuccLi : adresse du coeff. suivant sur la meme ligne
  Matrice : tableau contenant les éléments diagonaux puis la partie triangulaire inferieure stricte
  NbLign : nombre de lignes de la matrice A
  NbCoef : nimbre de coefficients mémorisés de la matrice A 
  SecMembre : tableau contenant le second membre
  NumDlDir : tableau de description des conditions de Dirichlet
  ValDLDir : valeur de blocage aux noeuds Dirichlet non homogenes
 ------------------------------------------- */ 
void assemblage(int N1, int N2, int ntel, int typel, int nbneel, int nbaret, int nRefDom, float **coorEl, int *nRefArEl, 
  int nbRefD0, int *numRefD0, int nbRefD1, int *numRefD1, int nbRefF1, int *numRefF1,
  int *NbLign, int *NbCoef, float *Matrice, int *AdPrCoefLi, int *AdSuccLi, 
  int *NumCol, float *SecMembre, int *NumDlDir, float *ValDLDir){
	  int res;
	  int dimDiag = N1*N2;
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
	  int nextad=1;
      int num[nbneel];
      for(int i=0;i<ntel;i++){
        for(int j=0;j<nbneel;j++){
	      num[j]=ngnel[i][j];
        }
      selectPts(nbneel, num,coord,coorEl);
      res = cal1Elem(nbneel, nbaret, nRefDom, coorEl, nRefAr[i], nbRefD0, numRefD0, nbRefD1, numRefD1, nbRefF1, numRefF1, 
        MatElem, SMbrElem, NuDElem, uDElem);
      if(res){return res;}
      //impCalEl(i+1, typel, nbneel, MatElem, SMbrElem, NuDElem, uDElem);
      }
	  
	  
}
