int assemblage(int ntel, int typel, int nbneel, int nbaret, int nRefDom,
  float **coord, int **nRefAr, int **ngnel,
  int nbRefD0, int *numRefD0, int nbRefD1, int *numRefD1, int nbRefF1, int *numRefF1,
  int *NbLign, int *NbCoef, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
  int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir);

void EcrSMD(int *NbCoef, int *NbLign, float *Matrice, int *AdPrCoefLi, int *AdSuccLi,
  int *NumCol, float *SecMembre, int *NumDLDir, float *ValDLDir);

 void LecSMD(int *NbCoef, int *NbLign, float **Matrice, int **AdPrCoefLi, int **AdSuccLi,
  int **NumCol, float **SecMembre, int **NumDLDir, float **ValDLDir);


