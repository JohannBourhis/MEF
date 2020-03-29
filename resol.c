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
     const float *Matrice, const int nbcprx, int *Profil, float *MatProf,
     int *codret){
  float *LMat = &Matrice[NbLign];
  int k=0;
  int n, I, Iref;
  MatProf[NbLign-1]=Matrice[NbLign-1];//diag
  for(i=0; i<NbLign-1; i++){ 
    MatProf[i]=Matrice[i]; //diag
    n=AdPrCoefLi[i+1]-AdPrCoefLi[i];
    Profil[i]=k;
    if(n!=0){
      Iref=NumCol[k];
      for(int j=0; j<n; j++){ 
        I=NumCol[k+j]-Iref;
        MatProf[NbLign+k+I]=LMat[k+j];
      }
      k+=i+1-Iref;
      if( (NbLign+k-1) > nbcprx){ codret=1;} //doute sur l'argument du test (+/-1)
    }
  }
  Profil[NbLign-1]=NbLign+k;
}
