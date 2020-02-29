#include <math.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
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
  int i, j, k, nbneel, n_quad;
  float *Wq, **Xq, *fctbase, **Derfctbase;
  float *Fk_xq, **Jac, det_Jac, **Inv_Jac, eltdif;
  float cofvar_W, cofvar_WW, *cofvar_ADWDW;
  float eps=EPS
  switch(nbneel){
    case 4 :
      n_quad=9;
      break;
    case 3 :
      n_quad=3;
      break;
    default :
      printf("Erreur dans le type transmis a intElem");
      return 2;
      break;
  }
  cofvar_ADWDW=malloc(4*sizeof(float));
  Fk_xq=malloc(2*sizeof(float)); if(Fk_xq==NULL){return 1;}
  Jac=alloctabf(2,2); if(Jac==NULL){return 1;}
  Inv_Jac=alloctabf(2,2); if(Inv_Jac==NULL){return 1;}
  Wq=malloc(n_quad*sizeof(float)); if(Wq==NULL){return 1;} // points et poids de quadrature
  Xq=alloctabf(n_quad,2); if(Xq==NULL){return 1;}
  fctbase=malloc(nbneel*sizeof(float)); if(fctbase==NULL){return 1;}
  Derfctbase=alloctabf(nbneel,2); if(Derfctbase==NULL){return 1;}
  // Calcul des points de quadrature de l'élément de référence
  ppquad(nbneel, Wq, Xq);
  for(i=0 ; i<n_quad ; i++){
    calFbase(nbneel, Xq[i], fctbase); 
    calDerFbase(nbneel, Xq[i], Derfctbase); 
    transFK(nbneel, coorEl, fctbase, Fk_xq); // Fk_xq : valeur de Fk en xq[i]
    matJacob(nbneel, 2, Jac, coorEl, Derfctbase);
    invertM2x2(Jac, det_Jac, Inv_Jac);
    det_Jac = fabsf(det_Jac)
    if(det_Jac < 10e-6){
      printf("Matrice singuliere")
      return 3;
    }
    eltdif=det_Jac*Wq[i];
    cofvar_W=FOMEGA(Fk_xq);
    cofvar_WW=A00(Fk_xq);
    cofvar_ADWDW[0]=A11(Fk_xq); cofvar_ADWDW[1]=A12(Fk_xq); 
    cofvar_ADWDW[2]=A21(Fk_xq); cofvar_ADWDW[3]=A22(Fk_xq);
    // Calcul des integrales
    W(nbneel, fctbas, eltdif, cofvar_W, SMbrElm);
    WW(nbneel, fctbas, eltdif, cofvar_WW, MatElem);
    ADWDW(nbneel, Derfctbase, Inv_Jac, eltdif, cofvar_ADWDW, MatElem);
  }
  free(cofvar_ADWDW); free(Fk_xq);
  freetab(Jac); freetab(Inv_Jac);
  free(Wq); freetab(Xq);
  free(fctbase); freetab(Derfctbase);
  return 0;
}


/*==============================================
Arguments d'entrée :
	float **coorAr : coordonnÃ©es des noeuds de l'arrete selectionÃ©e
				
Arguments de sortie :
	float* I_w : les valeurs de l'intÃ©grale de f*w_i sur K, pour tout i
	float** I_ww : les valeurs de l'intÃ©grale de a*w_i*w_j sur K, pour tout i,j
===============================================*/
int intAret(float **coorAR, int *num_noeuds, float* I_w, float **I_ww){
  int i, nbneel=2, n_quad=3;
  float *Wq, *Xq;
  float *fctbase, **Derfctbase;
  float *Fk_xq, **Jac, eltdif; 
  float cofvar_W, cofvar_WW;

  Fk_xq=malloc(2*sizeof(float)); if(Fk_xq == NULL){return 1;}
  Jac=alloctabf(2,1); if(Jac == NULL){return 1;}
  Wq=malloc(n_quad*sizeof(float)); if(Wq == NULL){return 1;} // points et poids de quadrature
  Xq=malloc(n_quad*sizeof(float));  if(Xq == NULL){return 1;} // sur les aretes les points de quadrature sont des reels
  fctbase=malloc(nbneel*sizeof(float)); if(fctbase == NULL){return 1;}
  Derfctbase=allocatbf(nbneel,1); if(Derfctbase == NULL){return 1;}

  // Calcul des points de quadrature
  ppquad(nbneel, Wq, Xq);
  for(i=0; i<n_quad; i++){
    calFbase(nbneel, Xq[i], fctbase); 
    calDerFbase(nbneel, Xq[i], Derfctbase); 
    matJacob(nbneel, 1, Jac, coorAr, Derfctbase);
    xq_Fk[0]=xq[i]*(coorAr[0][0]-coorAr[1][0])+coorAr[0][0];   // coorAr[0][0]-coorAr[1][0]=Jac[0]
    xq_Fk[1]=xq[i]*(coorAr[0][1]-coorAr[1][1])+coorAr[0][1];   // Jac[1]
    L=sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);	
    eltdif=L*Wq[i];
    cofvar_W=FN(Fk_xq);
    cofvar_WW=BN(Fk_xq);
    // Calcul des intégrales
    W(nbneel, fctbase, eltdif, cofvar_W, I_w);
    WW(nbneel, fctbase, eltdif, cofvar_WW, I_ww);
  }
  free(Fk_xq); free(Jac);
  free(Wq); free(Xq);
  free(fctbase); freetab(Derfctbase);
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
             float *nRefAr[nbaret]
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
int cal1Elem(int nbneel, int nbaret, int nRefDom, float **coorEl, float *nRefArEl, 
  int nbRefD0, float *numRefD0, int nbRefD1, float *numRefD1, int nbRefF1, float *numRefF1, 
    float ** MatElem, float *SMbrElem, int *NuDElem, float *uDElem){

  int i, j, k, l, nk, nl, R, numAr, cond_Ar;
  int *num_noeuds, float *coorAr;
  
  // Mini-matrice pour les aretes
  float SMbrAret, MatAret;
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
  if(R){return R};
  
  //Prise en compte des conditions aux limites
  for (i=0; i<nbaret ; i++){   // (i+1) numéro local de l'arrete
    cond_Ar=nRefArEl[i];
    if (cond_Ar == nRefDom){
      continue; // on passe au i suivant 
    }
    // Dirichlet homogène
    for (j=0; j<nbRefD0; j++){
      if (numRefD0[j]==cond_Ar){  
        numNaret(nbneel, i+1, num_noeuds);  
        NuDElem[num_noeuds[0]]=1;
        NuDElem[num_noeuds[1]]=1;        
        break;
      }
    }
    // Dirichlet non-homogène
    for(j=0; j<nbRefD1; j++){
      if (numRefD1[j]==cond_Ar) {  
        numNaret(nbneel, i+1, num_noeuds);  
        NuDElem[num_noeuds[0]]=-1;
        NuDElem[num_noeuds[1]]=-1;
        selectPts(2,num_noeuds, coorEl, coorAr);
        uDElem[num_noeuds[0]]=UD(coorAr[0]);
        uDElem[num_noeuds[1]]=UD(coorAr[1]);
        break;
      }
    }
    // Neumann ou Fourier     // num-noeuds decale les indices
    for(j=0; j<nbRefF1; j++){
      if (numRefF1[j]==cond_Ar) {  
        numNaret(nbneel, i+1, num_noeuds);  
        selectPts(2,num_noeuds, coorEl, coorAr);
        //Calcul des intégrales linéiques 
        R = intAret(coorAr, SMbrAret, MatAret) 
        if(R){return R;}
        for(k=0; k<nbneel; k++){
          nk = num_noeuds[k]-1;
          SMbrElem[nk] += SMbrAret[k];
          for(l=0; l<nbneel; l++){
            nl = num_noeuds[l]-1;
            MatElem[nk][nl]+=MatAret[k][l]; 
          }
        }
       break;
      }
    }
  }
  free(SMbrAret)
  freetab(MatAret)
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
	return 100*(x[0])+x[1];
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
 		Derfct0=Derfctbas[0,i]*InvJac[0,0] + Derfctbas[1,i]*InvJac[0,1]; 
 		Derfct1=Derfctbas[0,i]*InvJac[1,0] + Derfctbas[1,i]*InvJac[1,1];
 		coeff0=eltdif*(cofvar[0]*Derfct0 + cofvar[2]*Derfct1); // JFk*(a_11*dwi/dx1 + a_21*dwi/dx2) 
 		coeff1=eltdif*(cofvar[1]*Derfct0 + cofvar[3]*Derfct1); // JFk*(a_12*dwi/dx1 + a_22*dwi/dx2) 
 		for(j=0;j<nbneel;j++){
 			matelm[i][j]=matelm[i][j] + coeff0*(Derfctbas[0,j]*InvJac[0,0] + Derfctbas[1,j]*InvJac[0,1]) + coeff1*(Derfctbas[0,i]*InvJac[1,0] + Derfctbas[1,i]*InvJac[1,1]);

 			// on somme sur les != points de quad grâce aux appels de la fonction
 		}
 	}
}
