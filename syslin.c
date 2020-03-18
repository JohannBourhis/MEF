#include <stdio.h>
#include <stdlib.h>
#include "forfun.h"



/*----------------------------------------------
 Procédure qui lit la SMD de A
 Tous les arguments sont des arguments de sorties, toutes les informations sont écrites sur le fichier
 ----------------------------------------------*/
void LecSMD_TP4(int *NbCoef, int *NbLign, float **Matrice, int **AdPrCoefLi, int **AdSuccLi,
  int **NumCol, float **SecMembre, int **NumDLDir, float **ValDLDir){
	FILE* SMD;
    // on peut utiliser une chaine de caractère pour transmettre le nom du fichier texte
    if((SMD = fopen("SMD.txt", "r")) != NULL){
		fread(NbLign, sizeof(int), 1, SMD);
		// allocation des tableaux temporaires
		float *SecMembre_temp;   SecMembre_temp=malloc((*NbLign)*sizeof(float));
		int *NumDLDir_temp;   NumDLDir_temp=malloc((*NbLign)*sizeof(int));
		float *ValDLDir_temp;  ValDLDir_temp=malloc((*NbLign)*sizeof(float));
		int *AdPrCoefLi_temp;  AdPrCoefLi_temp=malloc((*NbLign)*sizeof(int));
		
		fread(SecMembre_temp, sizeof(float), *NbLign, SMD);
		fread(NumDLDir_temp, sizeof(int), *NbLign, SMD);
		fread(ValDLDir_temp, sizeof(float), *NbLign, SMD);
		fread(AdPrCoefLi_temp, sizeof(int), *NbLign, SMD);
		
		(*NbCoef)=AdPrCoefLi_temp[*NbLign-1]-1;
		float *Matrice_temp;  Matrice_temp=malloc((*NbCoef)*sizeof(float));
		int *NumCol_temp;  NumCol_temp=malloc((*NbCoef)*sizeof(int));
		int *AdSuccLi_temp;  AdSuccLi_temp=malloc((*NbCoef)*sizeof(int));
		
		fread(Matrice_temp, sizeof(float), (*NbLign)+(*NbCoef), SMD);
		fread(NumCol_temp, sizeof(int), *NbCoef, SMD);
		fread(AdSuccLi_temp, sizeof(int), *NbCoef, SMD);
	    
	    //transmission des tableaux en dehors de la fonction
	    *SecMembre= SecMembre_temp; 
	    *NumDLDir= NumDLDir_temp;
	    *ValDLDir= ValDLDir_temp;
	    *AdPrCoefLi= AdPrCoefLi_temp;
	    *Matrice= Matrice_temp;
	    *NumCol= NumCol_temp;
	    *AdSuccLi= AdSuccLi_temp;

	    //Liberation de la mémoire
	    free(SecMembre_temp);
	    free(NumDLDir_temp);
	    free(ValDLDir_temp);
	    free(AdPrCoefLi_temp);
	    free(Matrice_temp);
	    free(NumCol_temp);
	    free(AdSuccLi_temp);

	    //affsmd_(NbLign,AdPrCoefLi,NumCol,AdSuccLi,Matrice,SecMembre,NumDLDir,ValDLDir);
	    fclose(SMD);
    }
    else {printf("Erreur ouverture du fichier SMD\n"); }
}




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
void dSMDaSMO(/*char *fic_SMD*/){
  // SMD
  int NbLign, NbCoef;
  float *Matrice;
  int *AdPrCoefLi;
  int *AdSuccLi;
  int *NumCol;
  float *SecMembre;
  int *NumDLDir;
  float *ValDLDir;
  // remplissage des tableaux
  LecSMD_TP4(&NbCoef, &NbLign, &Matrice, &AdPrCoefLi, &AdSuccLi, 
    &NumCol, &SecMembre, &NumDLDir, &ValDLDir);

  // SMO
  float *MatriceO;
  int *AdPrCoefLiO;
  int *NumColO;
  float *SecMembreO;
  MatriceO=malloc(NbCoef*sizeof(float));
  SecMembreO=malloc(NbLign*sizeof(float));
  AdPrCoefLiO=malloc(NbLign*sizeof(int));
  NumColO=malloc(NbCoef*sizeof(float));

  //Appel à la procèdure fortran pour creer le SMO
  cdesse_(&NbLign, AdSuccLi, NumCol, AdSuccLi, Matrice, SecMembre, NumDLDir, ValDLDir,
   AdPrCoefLiO, NumColO, MatriceO, SecMembreO);

  // Création du système linéaire ???

  // écriture de la SMO
  FILE* SMO;
  if((SMO = fopen("SMD.txt", "w")) != NULL){
	fwrite(&NbLign, sizeof(int), 1, SMO);
    fwrite(SecMembreO, sizeof(float), NbLign, SMO);
	fwrite(AdPrCoefLiO, sizeof(int), NbLign, SMO);
	fwrite(MatriceO, sizeof(float), NbLign+NbCoef, SMO);
	fwrite(NumColO, sizeof(float), NbCoef, SMO);
	fclose(SMO);
  }
  else {printf("Erreur ouverture du fichier SMO\n"); }

 
}



/*----------------------------------------------
 Procédure qui lit la SMO de A
 Tous les arguments sont des arguments de sorties, toutes les informations sont écrites sur le fichier
 ----------------------------------------------*/
void LecSMO(int *NbCoef, int *NbLign, float **MatriceO, int **AdPrCoefLiO,
  int **NumColO, float **SecMembreO){
	FILE* SMO;
    // on peut utiliser une chaine de caractère pour transmettre le nom du fichier texte
    if((SMO = fopen("SMO.txt", "r")) != NULL){
		fread(NbLign, sizeof(int), 1, SMO);
		int NbL = (*NbLign);
		// allocation des tableaux temporaires
		float *SecMembre_temp;   SecMembre_temp=malloc((*NbLign)*sizeof(float));
		int *AdPrCoefLi_temp;  AdPrCoefLi_temp=malloc((*NbLign)*sizeof(int));
		
		fread(SecMembre_temp, sizeof(float), *NbLign, SMO);
		fread(AdPrCoefLi_temp, sizeof(int), *NbLign, SMO);
		
		(*NbCoef)=AdPrCoefLi_temp[NbL-1]-1;
		float *Matrice_temp;  Matrice_temp=malloc((*NbCoef)*sizeof(float));
		int *NumCol_temp;  NumCol_temp=malloc((*NbCoef)*sizeof(float));
		
		fread(Matrice_temp, sizeof(float), (*NbLign)+(*NbCoef), SMO);
		fread(NumCol_temp, sizeof(float), *NbCoef, SMO);
	    
	    //transmission des tableaux en dehors de la fonction
	    *SecMembreO= SecMembre_temp;
	    *AdPrCoefLiO= AdPrCoefLi_temp;
	    *MatriceO= Matrice_temp;
	    *NumColO= NumCol_temp;

	    //Liberation de la mémoire
	    free(SecMembre_temp);
	    free(AdPrCoefLi_temp);
	    free(Matrice_temp);
	    free(NumCol_temp);

	    //affsmd_(NbLign,AdPrCoefLi,NumCol,AdSuccLi,Matrice,SecMembre,NumDLDir,ValDLDir);
	    fclose(SMO);
    }
    else {printf("Erreur ouverture du fichier SMD\n"); }
}

