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
  int i, j, k, nQuad;
  float *Wq, **Xq, *fctbase, **Derfctbase;
  float *Fk, **Jac, detJac, **InvJac, eltdif;
  float cofvarW, cofvarWW, *cofvarADWDW;
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
  cofvarADWDW=malloc(4*sizeof(float));
  Fk=malloc(2*sizeof(float)); if(Fk==NULL){return 1;}
  Jac=alloctabf(2,2); if(Jac==NULL){return 1;}
  InvJac=alloctabf(2,2); if(InvJac==NULL){return 1;}
  Wq=malloc(nQuad*sizeof(float)); if(Wq==NULL){return 1;} // points et poids de quadrature
  Xq=alloctabf(nQuad,2); if(Xq==NULL){return 1;}
  fctbase=malloc(nbneel*sizeof(float)); if(fctbase==NULL){return 1;}
  Derfctbase=alloctabf(nbneel,2); if(Derfctbase==NULL){return 1;}
  // Calcul des points de quadrature de l'élément de référence
  ppquad(nbneel, Wq, Xq);
  for(i=0; i<nQuad; i++){
    calFbase(nbneel, Xq[i], fctbase); 
    calDerFbase(nbneel, Xq[i], Derfctbase); 
    transFk(nbneel, coorEl, fctbase, Fk); // Fk : valeur de Fk en xq[i]
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
    cofvarADWDW[0]=A11(Fk); cofvarADWDW[1]=A12(Fk); 
    cofvarADWDW[2]=A21(Fk); cofvarADWDW[3]=A22(Fk);
    
    // Calcul des integrales
    W(nbneel, fctbase, eltdif, cofvarW, SMbrElm);
    WW(nbneel, fctbase, eltdif, cofvarWW, MatElem);
    ADWDW(nbneel, Derfctbase, InvJac, eltdif, cofvarADWDW, MatElem);
  }
  free(cofvarADWDW); free(Fk);
  freetab(Jac); freetab(InvJac);
  free(Wq); freetab(Xq);
  free(fctbase); freetab(Derfctbase);
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
  float L;
  float *Wq, **Xq;
  float *fctbase, **Derfctbase;
  float *Fk, **Jac, eltdif; 
  float cofvarW, cofvarWW;
  printf("a\n");

  Fk=malloc(2*sizeof(float)); if(Fk == NULL){return 1;}
  Jac=alloctabf(2,1); if(Jac == NULL){return 1;}
  Wq=malloc(nQuad*sizeof(float)); if(Wq == NULL){return 1;} // points et poids de quadrature
  Xq=alloctabf(nQuad,1);  if(Xq == NULL){return 1;} // sur les aretes les points de quadrature sont des reels
  fctbase=malloc(nbneel*sizeof(float)); if(fctbase == NULL){return 1;}
  Derfctbase=alloctabf(nbneel,1); if(Derfctbase == NULL){return 1;}

  // Calcul des points de quadrature
  ppquad(nbneel, Wq, Xq);
  for(i=0; i<nQuad; i++){
	printf("inAret i %d \n",i);  // affichage
    calFbase(nbneel, Xq[i], fctbase); 
    calDerFbase(nbneel, Xq[i], Derfctbase); 
    matJacob(nbneel, 1, Jac, coorAr, Derfctbase);
    Fk[0]=Xq[i][0]*(coorAr[0][0]-coorAr[1][0])+coorAr[0][0];   // coorAr[0][0]-coorAr[1][0]=Jac[0]
    Fk[1]=Xq[i][0]*(coorAr[0][1]-coorAr[1][1])+coorAr[0][1];   // Jac[1]
    L=sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);	
    eltdif=L*Wq[i];
    cofvarW=FN(Fk);
    cofvarWW=BN(Fk);
    // Calcul des intégrales
    W(nbneel, fctbase, eltdif, cofvarW, SMbrAret);
    WW(nbneel, fctbase, eltdif, cofvarWW, MatAret);
  }
  free(Fk); free(Jac);
  free(Wq); free(Xq);
  free(fctbase); freetab(Derfctbase);
  printf("End inAret\n");
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
    float ** MatElem, float *SMbrElem, int *NuDElem, float *uDElem){

  int i, j, k, l, nk, nl, R, numAr, condAr;
  int numNoeuds[2]; float *coorAr[2];
  
  // Mini-matrice pour les aretes
  float *SMbrAret, **MatAret;
  SMbrAret = malloc(2*sizeof(float)); if(SMbrAret == NULL){return 1;}
  MatAret = alloctabf(2,2); if(MatAret == NULL){return 1;}
  
  //Initialisations de NuDElem, uDElem, SMbrElem et MatElem
  for(i=0; i<nbneel; i++){
    NuDElem[i]=0;
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
	printf("i=%d\n",i);
    condAr=nRefArEl[i];
    if (condAr == nRefDom){
      continue; // on passe au i suivant 
    }
    // Dirichlet homogène
    for (j=0; j<nbRefD0; j++){
	  printf("D0 j= %d\n",j);   ///////
      if (numRefD0[j]==condAr){  
		printf(" d\n");
        numNaret(nbneel, i+1, numNoeuds);    //// erreur
        printf(" f\n");
        NuDElem[numNoeuds[0]]=1;
        NuDElem[numNoeuds[1]]=1;        
        break;
      }
    }
    // Dirichlet non-homogène
    for(j=0; j<nbRefD1; j++){
	  printf("D1 j= %d\n",j);
      if (numRefD1[j]==condAr) {  
        numNaret(nbneel, i+1, numNoeuds);  
        NuDElem[numNoeuds[0]]=-1;
        NuDElem[numNoeuds[1]]=-1;
        selectPts(2,numNoeuds, coorEl, coorAr);
        uDElem[numNoeuds[0]]=UD(coorAr[0]);
        uDElem[numNoeuds[1]]=UD(coorAr[1]);
        break;
      }
    }
    
    // Neumann ou Fourier     // num-noeuds decale les indices
    for(j=0; j<nbRefF1; j++){
	  printf("F1 j= %d\n",j);
      if (numRefF1[j]==condAr) {  
        numNaret(nbneel, i+1, numNoeuds); 
        printf(" Neumann 1\n"); 
        selectPts(2, numNoeuds, coorEl, coorAr); // erreur
        printf(" Neumann 2\n");
        //Calcul des intégrales linéiques 
        R = intAret(coorAr, numNoeuds, SMbrAret, MatAret);
        printf("Neumann 3 \n");
        if(R){return R;}
        for(k=0; k<2; k++){  /// changement k<2
		  printf("D0 k %d\n",k);
          nk = numNoeuds[k]-1;
          SMbrElem[nk] += SMbrAret[k];
          for(l=0; l<nbneel; l++){
            nl = numNoeuds[l]-1;
            MatElem[nk][nl]+=MatAret[k][l]; 
            printf("MAtElem %f\n",MatAret[k][l]);
          }
        }
       break;
      }
    }
  }
  free(SMbrAret);
  freetab(MatAret);
  printf("End cal1Elem \n");
  return 0;
}

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
  return 1.0;
}
float BN(float *x){
  return 0.0;
}
float FOMEGA(float *x){
  return 1.0;
}
float FN(float *x){
  return 1.0;
}
float UD(float *x){
  return 100*x[0]+x[1];
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
    // on somme sur les points de quad a  chaque appel de la fonction
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
ON CALCULE UN LA QUADRATURE POUR UN POINT DE LA QUADARTURE xq
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
void ADWDW(int nbneel, float **Derfctbas, float **InvJac, float eltdif, float *cofvar, float **matelm){
  int i,j;
  float coeff0, coeff1, Derfct0, Derfct1;
  for(i=0; i<nbneel; i++){
    Derfct0=Derfctbas[0][i]*InvJac[0][0] + Derfctbas[1][i]*InvJac[0][1]; 
    Derfct1=Derfctbas[0][i]*InvJac[1][0] + Derfctbas[1][i]*InvJac[1][1];
    coeff0=eltdif*(cofvar[0]*Derfct0 + cofvar[2]*Derfct1); // JFk*(a_11*dwi/dx1 + a_21*dwi/dx2) 
    coeff1=eltdif*(cofvar[1]*Derfct0 + cofvar[3]*Derfct1); // JFk*(a_12*dwi/dx1 + a_22*dwi/dx2) 
    for(j=0;j<nbneel;j++){
      matelm[i][j]=matelm[i][j] + coeff0*(Derfctbas[0][j]*InvJac[0][0] + Derfctbas[1][j]*InvJac[0][1]) + 
	      coeff1*(Derfctbas[0][i]*InvJac[1][0] + Derfctbas[1][i]*InvJac[1][1]);
      // on somme sur les != points de quad grâce aux appels de la fonction
    }
  }
}
