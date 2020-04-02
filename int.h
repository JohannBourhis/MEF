int intElem(int typel, float **coorEl, float *SMbrElm, float **MatElem);
int intAret(float **coorAr, int numNoeuds[], float *SMbrAret,
            float **MatAret, int condar);
int cal1Elem(int nbneel, int nbaret, int nRefDom, float **coorEl,
             int *nRefArEl, int nbRefD0, int *numRefDO, int nbRefD1,
             int *numRefD1, int nbRefF1, int *numRefF1, float **MatElem,
             float *SMbrElem, int *NuDElem, float *uDElem);

// intégration
void W(int nbneel, float *fctbas, float eltdif,
       float cofvar, float *vectelm);
void WW(int nbneel, float *fctbas, float eltdif,
        float cofvar, float **matelm);
void ADWDW(int nbneel, float **Derfctbas, float **InvJac,
           float eltdif, float **cofvar, float **matelm);

// affichage
void impCalEl(int K, int typEl, int nbneel, float **MatElem,
              float *SMbrElem, int *NuDElem, float *uDElem);
