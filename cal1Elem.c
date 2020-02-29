/* --------------------------------------
Calcul elementaire sur l'element courant K :
Arguments d'entrée :
  typel    : type d'element
  nbneel   : nombre de noeuds de l'element
  nbaret   : nombre d'aretes de l'element
  coorEl   : tableau contenant les sommets de l'element courant K (à initialiser avec la fct selectPts à partir du fichier maillage)
  nRefArEl : tableau contenant les num de ref des nbaret de l'element
             float *nRefAr[nbaret]
  nRefDom  : etiquette coorespondant aux aretes internes
  nbRefD0  : taille de numRefD0
  numRefDO : numéro locaux des arretes coorespondant a une condition de Dirichlet  u=O
  nbRefD1  : taille de numRefD1
  numRefD1 : numéros locaux des arretes coorespondant a une condition de Dirichlet non homogène   u=ud(x)
  nbRefF1  : taille de numRefF1
  numRefF1 : numéros locaux des arretes coorespond a une condition de Neuman ou Fourier        eq diff

Arguments de sortie :
  MatElem  : éléments de la matrice élémentaire
  SMbrElem : éléments du second membre
  NuDElem  : noeuds qui ont la condition de dirichlet (0, 1, ou -1)
  uDElem   : valeur de la condition de dirichlet
 ------------------------------------------- */ 
#include "maille.h"
#include "calElmt.h"
#include "int.h"

int cal1Elem(int nbneel, int nbaret, int nRefDom, float **coorEl, float *nRefArEl, 
  int nbRefD0, float *numRefD0, int nbRefD1, float *numRefD1, int nbRefF1, float *numRefF1, 
    float ** MatElem, float *SMbrElem, int *NuDElem, float *uDElem){

  int i, j, k, l, nk, nl, R, numAr, cond_Ar;
  int *num_noeuds, float *coorAr;
  
  // Mini-matrice pour les aretes
  float SMbrAret, MatAret;
  SMbrAret = malloc(2*sizeof(float)); if(SMbrAret == NULL){return 1;}
  MatAret = alloctabf(2,2); if(MatAret == NULL){return 1;}
  
  //Initialisations de NuDElem, uDElem, SMbrElem et MatElem
  for(i=0; i<nbneel; i++){
    NuDElem[i]=0;
    uDElem[i]=0;
    SMbrElem[i]=0;
    for(j=0; j<nbneel; j++){
      MatElem[i][j]=0;
    }
  }
  
  //Calcul des intégrales surfaciques
  R = intElem(nbneel, coorEl, SMbrElem, MatElem);
  if(R){return R};
  
  //Prise en compte des conditions aux limites
  for (i=0; i<nbaret ; i++){   // (i+1) numéro local de l'arrete
    cond_Ar=nRefArEl[i];
    if (cond_Ar == nRefDom){
      continue; // on passe au i suivant 
    }
    // Dirichlet homogène
    for (j=0; j<nbRefD0; j++){
      if (numRefD0[j]==cond_Ar){  
        numNaret(nbneel, i+1, num_noeuds);  
        NuDElem[num_noeuds[0]]=1;
        NuDElem[num_noeuds[1]]=1;        
        break;
      }
    }
    // Dirichlet non-homogène
    for(j=0; j<nbRefD1; j++){
      if (numRefD1[j]==cond_Ar) {  
        numNaret(nbneel, i+1, num_noeuds);  
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
        numNaret(nbneel, i+1, num_noeuds);  
        selectPts(2,num_noeuds, coorEl, coorAr);
        //Calcul des intégrales linéiques 
        R = intAret(coorAr, SMbrAret, MatAret) 
        if(R){return R;}
        for(k=0; k<nbneel; k++){
          nk = num_noeuds[k]-1;
          SMbrElem[nk] += SMbrAret[k];
          for(l=0; l<nbneel; l++){
            nl = num_noeuds[l]-1;
            MatElem[nk][nl]+=MatAret[k][l]; 
          }
        }
       break;
      }
    }
  }
  free(SMbrAret)
  freetab(MatAret)
  return 0;
}
