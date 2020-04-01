#include <math.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
#include <stdlib.h>
#include <stdio.h>
#define EPS 10e-6 //sensibilité pour l'inversibilité de la Jacobienne

/*==============================================
Arguments d'entrée :
	int nbneel : nombre de noeud par élément (induit le type d'élément)
	float **coorEl : les noeuds de l'element K

Arguments de sortie :
        float** MatElem : Les éléments de la matrice élémentaire, pour tout i,j
	float* SmbreElm : Les éléments du second membre, pour tout i
===============================================*/
int intElem(int nbneel, float **coorEl, float *SMbrElm, float **MatElem){
  int i, nQuad;
  float detJac, eltdif;
  float cofvarW, cofvarWW;
  float eps=EPS;
  switch(nbneel){
    case 4 :
      nQuad=9;
      break;
    case 3 :
      nQuad=3;
      break;
    default :
      printf("Erreur dans le type transmis a intElem");
      return 2;
      break;
  }
  float Wq[nQuad];
  float **Xq=alloctabf(nQuad,2); if(Xq==NULL){return 1;};
  float fctbase[nbneel];
  float **Derfctbase=alloctabf(nbneel,2); if(Derfctbase==NULL){return 1;};
  float **cofvarADWDW=alloctabf(2,2); if(cofvarADWDW==NULL){return 1;};
  float Fk[2];
  float **Jac=alloctabf(2,2); if(Jac==NULL){return 1;};
  float **InvJac=alloctabf(2,2); if(InvJac==NULL){return 1;};
  // Calcul des points de quadrature de l'élément de référence
  ppquad(nbneel, Wq, Xq);
  for(i=0; i<nQuad; i++){
    calFbase(nbneel, Xq[i], fctbase);
    calDerFbase(nbneel, Xq[i], Derfctbase);
    transFk(nbneel, coorEl, Fk, fctbase); // Fk : valeur de Fk en xq[i]
    matJacob(nbneel, 2, Jac, coorEl, Derfctbase);
    invertM2x2(Jac, &detJac, InvJac);
    detJac = fabsf(detJac);
    if(detJac < eps){
      printf("Matrice singuliere \n");
      return 3;
    }
    eltdif=detJac*Wq[i];
    cofvarW=FOMEGA(Fk);
    cofvarWW=A00(Fk);
    cofvarADWDW[0][0]=A11(Fk); cofvarADWDW[0][1]=A12(Fk);
    cofvarADWDW[1][0]=A21(Fk); cofvarADWDW[1][1]=A22(Fk);

    // Calcul des integrales
    W(nbneel, fctbase, eltdif, cofvarW, SMbrElm);
    WW(nbneel, fctbase, eltdif, cofvarWW, MatElem);
    ADWDW(nbneel, Derfctbase, InvJac, eltdif, cofvarADWDW, MatElem);
  }
  freetab(Jac); freetab(InvJac);
  freetab(Xq); freetab(cofvarADWDW);
  freetab(Derfctbase);
  return 0;
}


/*==============================================
Arguments d'entrée :
	float **coorAr : coordonnÃ©es des noeuds de l'arrete selectionÃ©e

Arguments de sortie :
	float* SMbrAret : les valeurs de l'intégrale de f*w_i sur K, pour tout i
	float** MatAret : les valeurs de l'intégrale de a*w_i*w_j sur K, pour tout i,j
===============================================*/
int intAret(float *coorAr[], int numNoeuds[], float *SMbrAret, float **MatAret){
  int i, nbneel=2, nQuad=3;
  float L, eltdif;
  float cofvarW, cofvarWW;

  float Fk[2];
  float **Jac=alloctabf(2,1); if(Jac==NULL){return 1;};
  float Wq[nQuad]; // points et poids de quadrature
  float **Xq=alloctabf(nQuad,1); if(Xq==NULL){return 1;}; // sur les aretes les points de quadrature sont des reels
  float fctbase[nbneel];
  float **Derfctbase=alloctabf(nbneel,1); if(Derfctbase==NULL){return 1;};

  // Calcul des points de quadrature
  ppquad(nbneel, Wq, Xq);
  for(i=0; i<nQuad; i++){
    calFbase(nbneel, Xq[i], fctbase);
    calDerFbase(nbneel, Xq[i], Derfctbase);
    matJacob(nbneel, 1, Jac, coorAr, Derfctbase);
    transFk(nbneel, coorAr, Fk, fctbase); // Fk : valeur de Fk en xq[i]
    L=sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);
    eltdif=L*Wq[i];
    cofvarW=FN(Fk);
    cofvarWW=BN(Fk);
    // Calcul des intégrales
    W(nbneel, fctbase, eltdif, cofvarW, SMbrAret);
    WW(nbneel, fctbase, eltdif, cofvarWW, MatAret);
  }
  freetab(Jac); freetab(Xq);
  freetab(Derfctbase);
  return 0;
}

/* --------------------------------------
Calcul elementaire sur l'element courant K :
Arguments d'entrée :
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
  MatElem  : éléments de la matrice élémentaire
  SMbrElem : éléments du second membre
  NuDElem  : noeuds qui ont la condition de dirichlet (0, 1, ou -1)
  uDElem   : valeur de la condition de dirichlet
 ------------------------------------------- */
int cal1Elem(int nbneel, int nbaret, int nRefDom, float **coorEl, int *nRefArEl,
  int nbRefD0, int *numRefD0, int nbRefD1, int *numRefD1, int nbRefF1, int *numRefF1,
    float **MatElem, float *SMbrElem, int *NuDElem, float *uDElem){

  int i, j, k, l, nk, nl, R, condAr;
  int numNoeuds[2]; float *coorAr[2];

  // Mini-matrice pour les aretes
  float *SMbrAret, **MatAret;
  SMbrAret = malloc(2*sizeof(float)); if(SMbrAret == NULL){return 1;}
  MatAret = alloctabf(2,2); if(MatAret == NULL){return 1;}
  //Initialisations de NuDElem, uDElem, SMbrElem et MatElem
  for(i=0; i<nbneel; i++){
    NuDElem[i]=1;
    uDElem[i]=0;
    SMbrElem[i]=0;
    for(j=0; j<nbneel; j++){
      MatElem[i][j]=0;
    }
  }

  //Calcul des intégrales surfaciques
  R = intElem(nbneel, coorEl, SMbrElem, MatElem);
  if(R){return R;};
  //Prise en compte des conditions aux limites
  for (i=0; i<nbaret ; i++){   // (i+1) numéro local de l'arrete
    condAr=nRefArEl[i];
    if (condAr == nRefDom){
      continue; // on passe au i suivant
    }
    // Dirichlet homogène
    for (j=0; j<nbRefD0; j++){
      if (numRefD0[j]==condAr){
        numNaret(nbneel, i+1, numNoeuds);
        NuDElem[numNoeuds[0]-1]=0;
        NuDElem[numNoeuds[1]-1]=0;
        break;
      }
    }
    // Dirichlet non-homogène
    for(j=0; j<nbRefD1; j++){
      if (numRefD1[j]==condAr) {
        numNaret(nbneel, i+1, numNoeuds);
        NuDElem[numNoeuds[0]-1]=-1;
        NuDElem[numNoeuds[1]-1]=-1;
        selectPts(2,numNoeuds, coorEl, coorAr);
        uDElem[numNoeuds[0]-1]=UD(coorAr[0]);
        uDElem[numNoeuds[1]-1]=UD(coorAr[1]);
        break;
      }
    }
    // Neumann ou Fourier     // num-noeuds decale les indices
    for(j=0; j<nbRefF1; j++){
      if (numRefF1[j]==condAr) {
        numNaret(nbneel, i+1, numNoeuds);
        selectPts(2, numNoeuds, coorEl, coorAr);
        // initialisation de MatAret et SmbAret
        for(k=0; k<2; k++){
		  SMbrAret[k]=0;
		  for(l=0; l<2; l++){
			MatAret[k][l]=0;
		  }
		}
        //Calcul des intégrales linéiques
        R = intAret(coorAr, numNoeuds, SMbrAret, MatAret);
        if(R){return R;}
        for(k=0; k<2; k++){  /// changement k<2
          nk = numNoeuds[k]-1;
          SMbrElem[nk] += SMbrAret[k];
          for(l=0; l<2; l++){  /// changement k<2
            nl = numNoeuds[l]-1;
            MatElem[nk][nl]+=MatAret[k][l];
          }
        }
         break;
      }
    }
  }
  free(SMbrAret);
  freetab(MatAret);
  return 0;
}
extern int nucas;
// fonctions
float A12(float *x){
  return 0.0;
}
float A11(float *x){
  return 1.0;
}
float A22(float *x){
  return 1.0;
}
float A21(float *x){
  return 0.0;
}
float A00(float *x){
  return 0.0;
}
float BN(float *x){
  return 0.0;
}
float FOMEGA(float *x){
  float val;
  const float PI=3.141592;
  return -2*PI*PI*sin(PI*x[0])*sin(PI*x[1]);
  switch(nucas){
   case 1:
     val = 32.*( x[1]*(1-x[1]) + x[0]*(1-x[0]) );
     break;
   case 2:
     val = 2.*PI*PI*sin(PI*x[0])*sin(PI*x[1]);
     break;
   case 3:
     val = 2.*PI*PI*cos(PI*x[0])*cos(PI*x[1]);
     break;
   default:
     printf("problème non traité\n");
     exit(1);
     break;
  }
  return val;
}
float FN(float *x){
  return 0.0;
}
float UD(float *x){
  //return 100*x[0]+x[1];
  return 1.0;
}

// intégration

/*--------------------------------------
ON CALCULE LA QUADRATURE POUR UN POINT DE LA QUADARTURE xq
C'EST UN ELEMENT DE LA SOMME
Arguments d'entree :
  nbneel : nombre de noeuds de l'element
  fctbas : valeurs des fonctions de base au point de quadrature courant
  		   dimension utiles : fctbas(nbeel)
  eltdif : element differentiel(jacobienne de FK en xq) multipliÃ© par le poids de quadrature
  cofvar : valeur du coefficient variable (fonction Ã  integrer calculÃ©e en l'image par FK
  			du point de quad courant xq)
  vectelm : vecteur elementaire de masse a  actualiser
  			Dimension utile : vectelm(nbneel)
 Argument de sortie :
  vectelm : vecteur elementaire de masse actualisee
 ------------------------------------------- */
void W(int nbneel, float *fctbas, float eltdif, float cofvar, float *vectelm){
  int i;
  for(i=0 ; i<nbneel ; i++){
    vectelm[i] = vectelm[i] + eltdif*cofvar*fctbas[i];
  }
}

/*--------------------------------------
ON CALCULE UN LA QUADRATURE POUR UN POINT DE LA QUADARTURE xq
C'EST UN ELEMENT DE LA SOMME sur les points de la quadrature
Arguments d'entrée :
  nbneel : nombre de noeuds de l'element
  fctbas : valeurs des fonctions de base au point de quadrature courant
  		   dimension utiles : fctbas(nbeel)
  eltdif : élément différentiel(jacobienne en xq) multiplié par le poids de quadrature
  cofvar : valeur du coefficient variable (fonction à integrer calculée en l'image par FK
  			du point de quad courant xq)
  matelm : matrice élémentaire de masse à actualiser
  			Dimension utile : matelm(nbneel, nbneel)
 Argument de sortie :
  matelm : matrice élémentaire de masse actualisée
 ------------------------------------------- */
void WW(int nbneel, float *fctbas, float eltdif, float cofvar, float **matelm){
  int i,j;
  float coeff;
  for(i=0; i<nbneel; i++){
    coeff=eltdif*cofvar*fctbas[i];
    for(j=0;j<nbneel;j++){
      matelm[i][j]=matelm[i][j] + coeff*fctbas[j];
      // on somme sur les != points de quad à chaque appel de la fonction
    }
  }
}

/*--------------------------------------
ON CALCULE LA QUADRATURE POUR UN POINT DE LA QUADARTURE xq
C'EST UN ELEMENT DE LA SOMME
Arguments d'entrée :
  nbneel : nombre de noeuds de l'element
  fctbas : matrice des valeurs des fonctions de base au point de quadrature courant
  		   dimension utiles : fctbas(2,nbeel)
  InvJac : matrice de jacobienne de Fk^-1
  		   dimension : InvJac(2,2)
  eltdif : élément différentiel(jacobienne en xq) multiplié par le poids de quadrature
  cofvar : valeurs du coefficient variable (fonction à integrer calculée en l'image par FK
  			du point de quad courant xq)
  			Dimension utile : cofvar(4)
  			cofvar[0]=a_11(Fk(xq)) ; cofvar[1]=a_12(Fk(xq))  ; cofvar[2]=a_21(Fk(xq)) ; cofvar[3]=a_22(Fk(xq))
  matelm : matrice élémentaire de masse à actualiser
  			Dimension utile : matelm(nbneel, nbneel)
 Argument de sortie :
  matelm : matrice élémentaire de masse actualisée
 ------------------------------------------- */
void ADWDW(int nbneel, float **Derfctbas, float **InvJac, float eltdif, float **cofvar, float **matelm){
  int i,j,alpha,beta;
  float DWi, DWj, coeff;
  for(i=0; i<nbneel; i++){
    for(alpha=0;alpha<2;alpha++){
		DWi=Derfctbas[i][0]*InvJac[0][alpha]+Derfctbas[i][1]*InvJac[1][alpha];
        for(beta=0;beta<2;beta++){
          coeff=DWi*eltdif*cofvar[alpha][beta];
          for(j=0;j<nbneel;j++){
            DWj=Derfctbas[j][0]*InvJac[0][beta]+Derfctbas[j][1]*InvJac[1][beta];
            matelm[i][j]=matelm[i][j] + DWj*coeff;
            // on somme sur les != points de quad grâce aux appels de la fonction
          }
       }
    }
  }
}

