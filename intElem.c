/*==============================================
Arguments d'entrée :
	int typel : le type d'élément
	float **coorEl : les noeuds de l'element K
	
Arguments de sortie :
	float* I_w : les valeurs des l'intégrale de f*w_i sur K, pour tout i
	float** I_ww : les valeurs des l'intégrale de a*w_i*w_j sur K, pour tout i,j
	float** I_adwdw : les valeurs des l'intégrale de adwdw sur K, pour tout i,j
===============================================*/
#include <math.h>
#include "maille.h"
#include "calElmt.h"

// changer transFK
int intElem(int typel, float **coorEl, float *SMbrElm, float **MatElem){
  int i, j, k, nbneel, n_quad;
  float *Wq, **Xq;
  float *fctbase;
  float **Derfctbase;
  float *Fk_xq;
  float **Jac;
  float det_Jac;
  float **Inv_Jac;
  float eltdif;
  float cofvar_W, cofvar_WW;
  float *cofvar_ADWDW;
  if (typel==1){ // quadrangle
    nbneel=4;
    n_quad=9;
  }
  else if (typel==2){ // triangle
    nbneel=3;
    n_quad=3;
  }
  else{
    printf("Erreur dans le type transmis a�intElem");
    return 2;
  }
  cofvar_ADWDW=malloc(4*sizeof(float));
  Fk_xq=malloc(2*sizeof(float));
  if(Fk_xq==NULL){return 1;}
  Jac=alloctabf(2,2);
  if(Jac==NULL){return 1;}
  Inv_Jac=alloctabf(2,2);
  if(Inv_Jac==NULL){return 1;}
  Wq=malloc(n_quad*sizeof(float)); // points et poids de quadrature
  if(Wq==NULL){return 1;}
  Xq=alloctabf(n_quad,2);
  if(Xq==NULL){return 1;}
  fctbase=malloc(nbneel*sizeof(float));
  if(fctbase==NULL){return 1;}
  Derfctbase=alloctabf(nbneel,2);
  if(Derfctbase==NULL){return 1;}
  // Calcul des points de quadrature de l'element
  ppquad(typel, Wq, Xq);
  for(i=0 ; i<n_quad ; i++){
    calFbase(typel, Xq[i], fctbase); 
    calDerFbase(typel, Xq[i], Derfctbase); 
    transFK(nbneel, coorEl, fctbase, Fk_xq); // Fk_xq : valeur de Fk en xq[i]
    matJacob(typel, 2, Jac, coorEl, Derfctbase);
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
  free(cofvar_ADWDW);
  free(Fk_xq);
  freetab(Jac);
  freetab(Inv_Jac);
  free(Wq);
  freetab(Xq);
  free(fctbase);
  freetab(Derfctbase);
  return 0;
}
