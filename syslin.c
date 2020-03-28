#include <stdio.h>
#include <stdlib.h>
#include "calElmt.h"
#include "int.h"
#include "maille.h"
#include "forfun.h"
#include "assemblage.h"
/* --------------------------------------
 * Systeme linéaire
/--> Pas d'arguments de entrée car on lit un fichier SMD
Arguments du stokage desordoné :
  NbLign     : nombre de lignes de la matrice A
  AdPrCoefLi : adresses des 1er coeff. des lignes de la partie triang. inf. stricte (AdPrCoefLi[NbLign]=NbCoef)
  NumCol     : numeros de colonnes des coeff. de la partie triang. inf. stricte
  AdSuccLi   : adresse du coeff. suivant sur la meme ligne
  Matrice    : tableau contenant les éléments diagonaux puis la partie triangulaire inferieure stricte
  SecMembre  : tableau contenant le second membre
  NumDLDir   : tableau de description des conditions de Dirichlet
  ValDLDir   : valeur de blocage aux noeuds Dirichlet non homogenes
/--> Pas d'arguments de sortie car on écrit dans un fichier SMO
Arguments du stokage ordoné :
  AdPrCoefLiO : adresses des 1er coeff. des lignes de la partie triang. inf. stricte
  NumColO     : numeros de colonnes des coeff. de la partie triang. inf. stricte
  AdSuccLiO   : adresse du coeff. suivant sur la meme ligne
  Matrice    : tableau contenant les éléments diagonaux puis la partie triangulaire inferieure stricte
  SecMembre  : tableau contenant le second membre
  NumDLDir   : tableau de description des conditions de Dirichlet
  ValDLDir   : valeur de blocage aux noeuds Dirichlet non homogenes
 ------------------------------------------- */
int dSMDaSMO(int *NbLign, float *MatriceO, int *NumColO, float *Matrice,
             float *SecMembre, int *AdPrCoefLi, int *AdSuccLi, float *ValDLDir,
             int *NumDLDir, int *NumCol){
  //Lecture de la SMD
  printf("Lecture SMD\n");
  LecSMD(NbLign, Matrice, AdPrCoefLi, AdSuccLi,
         NumCol, SecMembre, NumDLDir, ValDLDir);
  //Appel à la procèdure fortran pour creer le SMO
  cdesse_(NbLign, AdPrCoefLi, NumCol, AdSuccLi, Matrice, SecMembre, NumDLDir, ValDLDir,
          AdPrCoefLi, NumColO, MatriceO, SecMembre);
  //écriture de la SMO
  FILE* SMO;
  if((SMO = fopen("SMO.txt", "w")) != NULL){
    fwrite(NbLign, sizeof(int), 1, SMO);
    fwrite(SecMembre, sizeof(float), *NbLign, SMO);
    fwrite(AdPrCoefLi, sizeof(int), *NbLign, SMO);
    int NbCoef = AdPrCoefLi[*NbLign-1]-1;
    fwrite(MatriceO, sizeof(float), *NbLign+NbCoef, SMO);
    fwrite(NumColO, sizeof(int), NbCoef, SMO);
    fclose(SMO);
  }
  else{
    printf("Erreur ouverture du fichier SMO\n");
    return 2;
  }
  return 0;
}



/*----------------------------------------------
 Procédure qui lit la SMO de A
 Tous les arguments sont des arguments de sorties, toutes les informations sont écrites sur le fichier
 ----------------------------------------------*/
int LecSMO(int *NbLign, float *MatriceO, int *AdPrCoefLiO,
            int *NumColO, float *SecMembreO){
  FILE* SMO;
  // on peut utiliser une chaine de caractère pour transmettre le nom du fichier texte
  if((SMO = fopen("SMO.txt", "r")) != NULL){
    fread(NbLign, sizeof(int), 1, SMO);
    // allocation des tableaux temporaires
    fread(SecMembreO, sizeof(float), *NbLign, SMO);
    fread(AdPrCoefLiO, sizeof(int), *NbLign, SMO);
    int NbCoef = AdPrCoefLiO[*NbLign-1]-1;
	  fread(MatriceO, sizeof(float), *NbLign+NbCoef, SMO);
	  fread(NumColO, sizeof(int), NbCoef, SMO);
    //transmission des tableaux en dehors de la fonction
    affsmo_(NbLign,AdPrCoefLiO,NumColO,MatriceO,SecMembreO);
    fclose(SMO);
  }
  else {
    printf("Erreur ouverture du fichier SMO\n");
    return 2;
  }
  return 0;
}

