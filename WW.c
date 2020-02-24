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
