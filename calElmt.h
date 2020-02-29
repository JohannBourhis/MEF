void invertM2x2(float **M, float *det, float **M_inv);

void matJacob(int nbneel, int dimDom, float **JacFk, float **S, float **Dfoncbase);

void calDerFbase(int typel, float *x, float **dw);

void calFbase(int typel, float *x, float *w);

void transFk(int nbneel, float **S, float Y[], float *foncbase);

void numNaret(int typel, int aretNum, int ptsNum[]);

void selectPts(int nb, int num[], float *coorEns[], float *coorSel[]);

void ppquad(int t, float *W, float **X);
