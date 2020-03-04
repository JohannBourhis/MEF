#include "maille.h"
#include "calElmt.h"
#include "int.h"

int main(){
  int nbneel=4;
  int nRefDom=0;
  int nbRefD0=1;
  int numReD0[nbRefD0];
  numReD0[0]=1;
  int nbRefD1=1;
  int numRefD1[nbRefD1];
  numRefD1[0]=4;
  int nbRefF1=2;
  int numRefF1[nbRefF1];
  numRefF1[0]=2; numRefF1[1]=3;
  float **MatElem;
  MatElem=alloctabf(nbneel,nbneel);
  float *SMbrElem=malloc(nbneel*sizeof(float));
  int *NuDElem=malloc(nbneel*sizeof(int));
  float *uDElem=malloc(nbneel*sizeof(float));

  // Sommets
  float **coorEl=alloctabf(nbneel,2); 
  coorEl[0][0]=0; coorEl[0][1]=0;
  coorEl[1][0]=0; coorEl[1][1]=1;
  coorEl[2][0]=1; coorEl[2][1]=1;
  coorEl[3][0]=1; coorEl[3][1]=0;
  
  int *nRefArEl=malloc(nbneel*sizeof(int));
  nRefArEl[0]=1;nRefArEl[1]=2;nRefArEl[2]=3;nRefArEl[3]=4;

  int cal1Elem(nbneel, nRefDom, float **coorEl, float *nRefArEl, 
    int nbRefD0, float *numRefD0, int nbRefD1, float *numRefD1, int nbRefF1, float *numRefF1, 
      float ** MatElem, float *SMbrElem, int *NuDElem, float *uDElem);

}
/*
* Numeros de references a donner dans l'ordre :
* 
* nrefdm
* nbred0, (nured0(i),i=1,nbred0)
* nbred1, (nured1(i),i=1,nbred1)
* nbref1, (nuref1(i),i=1,nbref1)
* 
* avec
* nrefdm : numero de reference du domaine
* nbred0 : nombre de numeros de reference Dirichlet homogene
* nured0 : tableau des numeros de reference Dirichlet homogene
* nbred1 : nombre de numeros de reference Dirichlet non homogene
* nured1 : tableau des numeros de reference Dirichlet non homogene
* nbref1 : nombre de numeros de reference Neumann
* nuref1 : tableau des numeros de reference Neumann
*/
