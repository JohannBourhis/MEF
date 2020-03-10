void etiqAr(int typel, int n1, int n2, int nrefdom, 
  const int *nrefcot, int nbtel, int nbaret, int **nRefAr);
  
void freetab(void *ptr);

float **alloctabf(int dim1, int dim2);

int **alloctab(int dim1, int dim2);

int lecfima(char *ficmai, int *typel, int *nbtng, float ***pcoord,
  int *nbtel, int ***pngnel, int *nbneel, int *nbaret, int ***pnRefAr);

int lecfima(char *ficmai, int *typel, int *nbtng, float ***pcoord,
  int *nbtel, int ***pngnel, int *nbneel, int *nbaret, int ***pnRefAr);
