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

