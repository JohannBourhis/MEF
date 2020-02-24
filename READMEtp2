Remarques concernant le TP2
---------------------------
1. Le fichier NUMREF.Test contient les numeros de reference utilises pour les tests
   dont le resultat est fourni dans les fichiers Tests.*.

2. Canevas de la fonction cal1Elem qui calcule les integrales
   elementaires surfaciques et lineiques relatives a un element.

   IMPORTANT : la prise en compte des conditions aux limites (4e point ci-dessous)
               ne doit etre faite que lorsque la fonction intElem est au point.

========		 
Declaration et allocation des variables et tableaux locaux

Initialisations de NuDElem, uDElem, SMbrElem et MatElem

Calcul des integrales surfaciques (appel a intElem)

Prise en compte des conditions aux limites du domaine
 pour chaque arete de l'element
  si l'arete se situe sur un bord du domaine

    determination des noeuds de l'arete (appel a numNaret, voir paragraphe 3 ci-dessous)
    extraction des coordonnees coorAr de ses points (appel a selectPts)

    -> cas d'un bord de type Neumann
       . allocation et initialisation de MatAret et SMbrAret
       . calcul des integrales lineiques (appel a intAret)
       . mise a jour de la matrice elementaire grace aux valeurs calculees
	 aux noeuds d'aretes (voir paragraphe 3 ci-dessous).

    -> cas d'un bord de type Dirichlet homogene
       . mise a jour de NuDElem

    -> cas d'un bord de type Dirichlet non homogene
	. pour chaque noeud de l'arete (non deja traite) :
            - mise a jour de NuDElem
            - calcul de la valeur donnee en ce noeud (appel a UD et stockage dans uDElem)
========

3. Utilisation de la fonction numNaret

Notons nunear le tableau d'entiers fourni par la fonction numNaret ; ce tableau
contient les numeros locaux des noeuds relatifs a une arete.
On en a besoin pour mettre a jour les tableaux NuDElem et uDElem dans le cas ou il y a
des conditions de Dirichlet et pour calculer les integrales sur cette arete, dans le cas
ou il y a des conditions de Neumann.

Ces integrales sont calculees par intAret, qui est construit sur le modele de intElem.
En effet, dans cal1Elem, on appelle intElem en transmettant la matrice MatElem (a
l'echelle de l'element). De meme, on appelle intAret en transmettant une mini-matrice,
disons MatAret (a l'echelle de l'arete). Il faut ensuite incorporer les valeurs
de MatAret dans MatElem justement par l'intermediaire du tableau nunear, par une
double boucle en i,j ou on aura :
  MatElem[ni][nj] = MatElem[ni][nj] + MatAret[i][j];
avec ni = nunear[i]-1; et nj = nunear[j]-1;

Il faut bien sur faire le meme type de traitement sur le second membre (SMbrElem
pour l'element, SMbrAret pour l'arete).
