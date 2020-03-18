#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"
#include "int.h"
#include "maille.h"
#include "forfun.h"
#include "assemblage.h"
/* --------------------------------------
 * Systeme lin�aire
/--> Pas d'arguments de entr�e car on lit un fichier SMD
Arguments du stokage desordon� :
  NbLign     : nombre de lignes de la matrice A
  AdPrCoefLi : adresses des 1er coeff. des lignes de la partie triang. inf. stricte (AdPrCoefLi[NbLign]=NbCoef)
  NumCol     : numeros de colonnes des coeff. de la partie triang. inf. stricte
  AdSuccLi   : adresse du coeff. suivant sur la meme ligne
  Matrice    : tableau contenant les �l�ments diagonaux puis la partie triangulaire inferieure stricte
  SecMembre  : tableau contenant le second membre
  NumDLDir   : tableau de description des conditions de Dirichlet
  ValDLDir   : valeur de blocage aux noeuds Dirichlet non homogenes
/--> Pas d'arguments de sortie car on �crit dans un fichier SMO
Arguments du stokage ordon� :
  AdPrCoefLiO : adresses des 1er coeff. des lignes de la partie triang. inf. stricte
  NumColO     : numeros de colonnes des coeff. de la partie triang. inf. stricte
  AdSuccLiO   : adresse du coeff. suivant sur la meme ligne
  Matrice    : tableau contenant les �l�ments diagonaux puis la partie triangulaire inferieure stricte
  SecMembre  : tableau contenant le second membre
  NumDLDir   : tableau de description des conditions de Dirichlet
  ValDLDir   : valeur de blocage aux noeuds Dirichlet non homogenes
 ------------------------------------------- */
void dSMDaSMO(int *NbLign, float *SecMembreO, int *AdPrCoefLiO,
              float *MatriceO, int *NumColO){
  // SMD
  int i, IAd, N;
  int NbCoef;
  float *Matrice;
  int *AdPrCoefLi;
  int *AdSuccLi;
  int *NumCol;
  float *SecMembre;
  int *NumDLDir;
  float *ValDLDir;
  //Lecture de la SMD
  printf("Lecture SMD\n");
  LecSMD(&NbCoef, NbLign, &Matrice, &AdPrCoefLi, &AdSuccLi,
         &NumCol, &SecMembre, &NumDLDir, &ValDLDir);
  //Appel � la proc�dure fortran pour creer le SMO
  for(i=0 ; i<*NbLign ; i++){
    AdPrCoefLiO[i]=0;
  }
  cdesse_(NbLign, AdPrCoefLi, NumCol, AdSuccLi, Matrice, SecMembre, NumDLDir, ValDLDir,
          AdPrCoefLiO, NumColO, MatriceO, SecMembreO);
  //Rangement des coefficients
  float *LowMat = &MatriceO[*NbLign];
  for(i=0; i<*NbLign; i++){
    IAd = AdPrCoefLiO[i]-1;
    N = AdPrCoefLiO[i]-AdPrCoefLiO[i+1];
    tri_(&N, &NumColO[IAd], &LowMat[IAd]);
  }
  //�criture de la SMO
  FILE* SMO;
  if((SMO = fopen("SMO.txt", "w")) != NULL){
	fwrite(NbLign, sizeof(int), 1, SMO);
    fwrite(SecMembreO, sizeof(float), *NbLign, SMO);
	fwrite(AdPrCoefLiO, sizeof(int), *NbLign, SMO);
	fwrite(MatriceO, sizeof(float), *NbLign+NbCoef, SMO);
	fwrite(NumColO, sizeof(float), NbCoef, SMO);
	fclose(SMO);
  }
  else{
    printf("Erreur ouverture du fichier SMO\n");
  }
  printf("b\n");
}



/*----------------------------------------------
 Proc�dure qui lit la SMO de A
 Tous les arguments sont des arguments de sorties, toutes les informations sont �crites sur le fichier
 ----------------------------------------------*/
int LecSMO(int *NbCoef, int *NbLign, float **MatriceO, int **AdPrCoefLiO,
            int **NumColO, float **SecMembreO){
  FILE* SMO;
  // on peut utiliser une chaine de caract�re pour transmettre le nom du fichier texte
  if((SMO = fopen("SMO.txt", "r")) != NULL){
    fread(NbLign, sizeof(int), 1, SMO);
    // allocation des tableaux temporaires
    float *SecMembre_temp = malloc((*NbLign)*sizeof(float)); if(SecMembre_temp==NULL){return 1;}
    int *AdPrCoefLi_temp = malloc((*NbLign)*sizeof(int)); if(AdPrCoefLi_temp==NULL){return 1;}
    fread(SecMembre_temp, sizeof(float), *NbLign, SMO);
    fread(AdPrCoefLi_temp, sizeof(int), *NbLign, SMO);
    (*NbCoef) = AdPrCoefLi_temp[*NbLign-1]-1;
	float *Matrice_temp = malloc((*NbCoef)*sizeof(float)); if(Matrice_temp==NULL){return 1;}
	int *NumCol_temp = malloc((*NbCoef)*sizeof(float)); if(NumCol_temp==NULL){return 1;}
	fread(Matrice_temp, sizeof(float), (*NbLign)+(*NbCoef), SMO);
	fread(NumCol_temp, sizeof(float), *NbCoef, SMO);
    //transmission des tableaux en dehors de la fonction
    *SecMembreO = SecMembre_temp;
    *AdPrCoefLiO = AdPrCoefLi_temp;
    *MatriceO = Matrice_temp;
    *NumColO = NumCol_temp;
    affsmo_(NbLign,*AdPrCoefLiO,*NumColO,*MatriceO,*SecMembreO);
    fclose(SMO);
  }
  else {
    printf("Erreur ouverture du fichier SMD\n");
    return 2;
  }
  return 0;
}

