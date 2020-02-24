/*--------------------------------------
ON CALCULE LA QUADRATURE POUR UN POINT DE LA QUADARTURE xq
C'EST UN ELEMENT DE LA SOMME
Arguments d'entree :
  nbneel : nombre de noeuds de l'element
  fctbas : valeurs des fonctions de base au point de quadrature courant
  		   dimension utiles : fctbas(nbeel)
  eltdif : element differentiel(jacobienne de FK en xq) multipli√© par le poids de quadrature 
  cofvar : valeur du coefficient variable (fonction √† integrer calcul√©e en l'image par FK 
  			du point de quad courant xq)
  vectelm : vecteur elementaire de masse a† actualiser
  			Dimension utile : vectelm(nbneel)
 Argument de sortie :
  vectelm : vecteur elementaire de masse actualisee
 ------------------------------------------- */ 

void W(int nbneel, float *fctbas, float eltdif, float cofvar, float *vectelm){
  int i;
  for(i=0 ; i<nbneel ; i++){
    vectelm[i] = vectelm[i] + eltdif*cofvar*fctbas[i];
    // on somme sur les points de quad a† chaque appel de la fonction
  }
}
