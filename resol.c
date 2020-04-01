#include "forfun.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
--------------------------------------------------------------------------------
   Passage de la structure morse ordonnee (S.M.O.) a la structure profil.

   Arguments d'entree :
  NbLign : nombre de lignes de la matrice
  AdPrCoefLi : adresses (dans lmatri) des premiers coefficients des lignes de
           la partie triangulaire inferieure stricte de la matrice
  NumCol : numeros de colonnes des coeff. de la partie triang. inf. stricte
  Matrice : tableau contenant la diagonale de la matrice (dmatri) puis la
           partie triangulaire inferieure stricte (lmatri) sous forme S.M.0.
  nbcprx : taille max du tableau MatProf

   Arguments de sortie :
  Profil : adresses (dans lmatrip) des premiers coefficients des lignes de
           la partie triangulaire inferieure stricte de la matrice
  MatProf : tableau contenant la diagonale de la matrice (dmatrip) puis la
           partie triangulaire inferieure stricte (lmatrip) sous forme profil
  codret : code de retour egal a 0 si R.A.S., 1 si la taille de MatProf est
           insuffisante (nbcprx)
--------------------------------------------------------------------------------
*/
void dSMOaPR (const int NbLign, const int *AdPrCoefLi, const int *NumCol,
              const float *Matrice, const int nbcprx, int *Profil,
              float *MatProf, int *codret){
  int FirstLiProf=1, counterSMO=0;
  int nbCoefLi, nbCoefProf, Iref, I;
  // Diagonale
  for(int i=0; i<NbLign; i++){
    MatProf[i] = Matrice[i];
  }
  // Initialisation à 0 de la partie inférieure
  for(int i=NbLign ; i<nbcprx ; i++){
    MatProf[i] = 0;
  }
  *codret = 0;
  float *LMatProf = &MatProf[NbLign];
  for(int i=0; i<NbLign-1; i++){
    nbCoefLi = AdPrCoefLi[i+1]-AdPrCoefLi[i]; // nombre d'éléments non nuls sur la ligne i
    Profil[i] = FirstLiProf;
    if(nbCoefLi != 0){
      Iref = NumCol[counterSMO] - 1; //nombre d'éléments nuls avant le premier élément non-nul
      nbCoefProf = (i+1) - Iref;
      if(FirstLiProf-1 + nbCoefProf > nbcprx){ // taille insuffisante
        *codret=1;
        break;
      }
      for(int j=0 ; j<nbCoefLi ; j++){
        I = NumCol[counterSMO+j] - Iref;
        LMatProf[FirstLiProf+I-2] = Matrice[NbLign+counterSMO+j];
      }
      counterSMO += nbCoefLi;
      FirstLiProf += nbCoefProf; // nombre d'éléments par ligne - Iref
    }
  }
  Profil[NbLign-1]=FirstLiProf;
}

/*
--------------------------------------------------------------------------------
    Evaluation de la solution exacte
--------------------------------------------------------------------------------
*/
extern int nucas;

float solex(float *coor){
  const float PI=3.141592;
  float val;

  switch (nucas) {
    case 1 :
  val=16.*coor[0]*coor[1]*(1-coor[0])*(1-coor[1]);
  break;
    case 2 :
  val=sin(PI*coor[0])*sin(PI*coor[1]);
  break;
    case 3 :
  val=cos(PI*coor[0])*cos(PI*coor[1]);
  break;
    default :
  printf("*** SOLEX : exemple non prevu. Abandon.\n");
  exit(1);
  break;
  }
  return(val);
}


/*
--------------------------------------------------------------------------------
  Calcul de solution en chaque neuds contenu dans  coord

   Arguments d'entree :
  NbLign : nombre de noeuds du domaine
  coord : coordonéées de chauqe noeuds
          coord[i][j] --> valeur de la jème composantes du i_ème noeud

   Arguments de sortie :
  Uex : vecteur contenant la valeur de la solution en chaque point de coord
--------------------------------------------------------------------------------
*/
void Calsol( int NbLign, float** coord, float *Uex){
  for(int i=0; i<NbLign; i++){
    Uex[i]=solex(coord[i]);
  }
}
