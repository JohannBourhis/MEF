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


// Lignes à rajouter dans TP3 :

/*
int nbcprx = NbLign + 3*dimLmat; // ?????
MatProf = malloc(nbcprx*sizeof(float));         if(MatProf==NULL){return 1;}
Profil = malloc(nbtng*sizeof(int));             if(Profil==NULL){return 1;}
int codret=0;

dSMOaPR(nbtng, AdPrCoefLi, NumColO, MatriceO, nbcprx, 
        Profil, MatProf, codret);
if(codret){printf("la taille de MatProf est insuffisante (nbcprx)\n")

dimLmat=Profil[nbtng];
float eps=1e-6;
float* ad = malloc(nbtng*sizeof(int));          if(ad==NULL){return 1;}
ad=&MatProf[0];
float* al = malloc(dimLmat*sizeof(int));        if(al==NULL){return 1;}
al=&MatProf[nbtng];
float* ld = malloc(nbtng*sizeof(int));          if(ld==NULL){return 1;}
float* ll = malloc(dimLmat*sizeof(int));        if(ll==NULL){return 1;}

//factorisation A=LLt
ltlpr(nbtng, profil, ad, al, eps, ld, ll);


*/
