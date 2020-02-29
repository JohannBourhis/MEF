/*==============================================
Arguments d'entrée :
	int nbneel : nombre de noeud par élément (induit le type d'élément)
	float **coorEl : les noeuds de l'element K
	
Arguments de sortie :
        float** MatElem : Les éléments de la matrice élémentaire, pour tout i,j
	float* SmbreElm : Les éléments du second membre, pour tout i
	????
	float* I_w : les valeurs des l'intégrale de f*w_i sur K, pour tout i
	float** I_ww : les valeurs des l'intégrale de a*w_i*w_j sur K, pour tout i,j
	float** I_adwdw : les valeurs des l'intégrale de adwdw sur K, pour tout i,j
	????
===============================================*/
#include <math.h>
#include "maille.h"
#include "calElmt.h"
#include "int.h"
#define EPS 10e-6 //sensibilité pour l'inversibilité de la Jacobienne

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
      printf("Erreur dans le type transmis a intElem");
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
