
/*==============================================
Arguments d'entrée :
	float **coorAr : coordonnÃ©es des noeuds de l'arrete selectionÃ©e
				
Arguments de sortie :
	float* I_w : les valeurs de l'intÃ©grale de f*w_i sur K, pour tout i
	float** I_ww : les valeurs de l'intÃ©grale de a*w_i*w_j sur K, pour tout i,j
===============================================*/
#include "maille.h"
#include "calElmt.h"
#include "int.h"

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
