void etiqAr(int typel, int n1, int n2, int nrefdom, 
  const int *nrefcot, int nbtel, int nbaret, int **nRefAr);
  
int lecfima(char *ficmai, int *typel, int *nbtng, float ***pcoord,
  int *nbtel, int ***pngnel, int *nbneel, int *nbaret, int ***pnRefAr);
  
int **alloctab(int dim1, int dim2);

float **alloctabf(int dim1, int dim2);

void freetab(void *ptr);
