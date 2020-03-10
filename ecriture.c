#include <stdio.h>
#include "maille.h"


int main(){
  int n1 = 2;
  int n2 = 2;
  int n = n1*n2;
  int t = 2;
  int p = 5-t; /* nombre de noeud par element */
  int q = p;  
  float a = 0;
  float b = 1;
  float c = 0;
  float d = 1;
  int nrefdom = 0;
  int nrefcot[4];
  nrefcot[0]=1;
  nrefcot[1]=2; 
  nrefcot[2]=3; 
  nrefcot[3]=4;
  int k = 0;
  int nbtel = (n1-1)*(n2-1)*t;
  FILE* maillage;
  if((maillage = fopen("maillage.txt", "w")) != NULL){ 
    fprintf(maillage, "%d \n", n);
    float h1 = (b-a)/(n1-1);
    float h2 = (d-c)/(n2-1);
    for(int j = 0; j < n2; j++){ 
      for(int i = 0; i < n1; i++){
        fprintf(maillage, "%f %f \n", a+h1*i, c+h2*j);
      }
    }
    fprintf(maillage, "%d %d %d %d \n", nbtel, t, p, q);
    int **nRefAr = alloctab(nbtel, q);
    etiqAr(t, n1, n2, nrefdom, nrefcot, nbtel, q, nRefAr);
    switch(t){
      case 1 :
        for(int j=0; j<n1*(n2-1); j+=n1){
          for(int i=1; i<n1; i++){
	    fprintf(maillage, "%d %d %d %d ", (i+j+1), (i+j+n1), (i+j+n1+1), (i+j));
	    fprintf(maillage, "%d %d %d %d \n", nRefAr[k][0], nRefAr[k][1], nRefAr[k][2], nRefAr[k][3]);
	    k++;
	  }
        }
      break;
      case 2 :
        for(int j=0; j<n1*(n2-1); j+=n1){
	    for(int i=1; i<n1; i++){
	      fprintf(maillage, "%d %d %d ", (i+j+1), (i+j+n1), (i+j));
	      fprintf(maillage, "%d %d %d \n", nRefAr[k][0], nRefAr[k][1], nRefAr[k][2]);
          fprintf(maillage, "\n");
	      k++;
	      fprintf(maillage, "%d %d %d ", (i+j+n1), (i+j+1), (i+j+n1+1));
	      fprintf(maillage, "%d %d %d \n", nRefAr[k][0], nRefAr[k][1], nRefAr[k][2]);
	      fprintf(maillage, "\n");
	      k++;	    
          }
        }
      break;
    }
    freetab(nRefAr);
  }
}
