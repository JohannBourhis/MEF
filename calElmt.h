void invertM2x2(float **M, float *det, float **M_inv);

void matJacob(int nbneel, int dimDom, float **JacFk, float **S, float **Dfoncbase);

void calDerFbase(int nbneel, float *x, float **dw);

void calFbase(int nbneel, float *x, float *w);

void transFk(int nbneel, float **S, float Y[], float *foncbase);

void numNaret(int nbneel, int aretNum, int ptsNum[]);

void selectPts(int nb, int num[], float *coorEns[], float *coorSel[]);

void ppquad(int nbneel, float *W, float **X);
