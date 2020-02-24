/* --------------------------------------
Calcul elementaire sur l'element courant K :
Arguments d'entrée :
  tyepl    : type d'element
  nbneel   : nombre de noeuds de l'element
  nbaret   : nombre d'aretes de l'element
  coorEl   : tableau contenant des sommets de l'element courant K (à initialiser avec la fct selectPts à partir du fichier maillage)
  nRefArEl : tableau contenant les num de ref des nbaret de l'element
             float *nRefAr[nbaret]
  nRefDom  : etiquette coorespondant aux aretes internes
  nbRefD0  : taille de numRefD0
  numRefDO : numéro locaux des arretes coorespondant a une condition de Dirichlet  u=O
  nbRefD1  : 
  numRefD1 : numéro locaux des arretes coorespondant a une condition de Dirichlet non homogène   u=ud(x)
  nbRefF1  : 
  numRefF1 : numéro locaux des arretes coorespond a une condition de Neuman ou Fourier        eq diff

Arguments de sortie :
  MatElem  : 
  SMbrElem : 
  NuDElem  : 
  uDElem   : 
 ------------------------------------------- */ 
#include "maille.h"
#include "calElmt.h"

int cal1Elem( int typel, int nbneel, int nbaret, int nRefDom, float **coorEl, float *nRefArEl, 
  int nbRefD0, float *numRefD0, int nbRefD1, float *numRefD1, int nbRefF1, float *numRefF1, 
    float ** MatElem, float *SMbrElem, int *NuDElem, float *uDElem){

  int i, j, k, R, numAr, cond_Ar;
  int *num_noeuds;
  float **coorAr;
  float I_ar_w, I_ar_ww;
  
  I_ar_w = malloc(2*sizeof(float));
  if(I_ar_w == NULL){return 1;}
  I_ar_ww = alloctab(2,2);
  if(I_ar_ww == NULL){return 1;}
  
  R = intElem(typel, coorEl, SMbrElem, MatElem)
  if(R){return R};

  for(i=0; i<nbaret ; i++){   // (i+1) numéro local de l'arrete
    cond_Ar=nRefArEl[i];
    if (cond_Ar == nRefDom){
      continue; // on passe au i suivant 
    }
    // Dirichlet homogène
    for (j=0; j<nbRefD0; j++){
      if (numRefD0[j]==cond_Ar){  
        numNaret(typel, i+1, num_noeuds);  
        NuDElem[num_noeuds[0]]=0;
        NuDElem[num_noeuds[1]]=0;
        break;
      }
    }
    // Dirichlet non-homogène
    for(j=0; j<nbRefD1; j++){
      if (numRefD1[j]==cond_Ar) {  
        numNaret(typel, i+1, num_noeuds);  
        NuDElem[num_noeuds[0]]=-1;
        NuDElem[num_noeuds[1]]=-1;
        selectPts(2,num_noeuds, coorEl, coorAr);
        uDElem[num_noeuds[0]]=UD(coorAr[0]);
        uDElem[num_noeuds[1]]=UD(coorAr[1]);
        break;
      }
    }
    // Neumann ou Fourier     // num-noeuds decale les indices
    for(j=0; j<nbRefF1; j++){
      if (numRefF1[j]==cond_Ar) {  
        numNaret(typel, i+1, num_noeuds);  
        selectPts(2,num_noeuds, coorEl, coorAr);
        R = intAret(coorAr, I_ar_w, I_ar_ww) 
        if(R){return R;}
        for(j=0; j<nbneel; j++){
          SMbrElem[num_noeuds[j]] += I_ar_w[j];
          for(k=0; k<nbneel; k++){
            MatElem[num_noeuds[j]][num_noeuds[k]-1]+=I_ar_ww[j][k]; 
          }
        }
       break;
      }
    }
  }
  free(I_ar_w)
  free(I_ar_ww)
  return 0;
}
