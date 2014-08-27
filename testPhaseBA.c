#include "petscksp.h"
#include "phaseBA.h"


/* Test implementation of Giordano et al 2008 viscosity model for silicate melts */

int main(){  
  Composition c;

  c.melt.Si = 62.38;
  c.melt.Ti = 0.41;
  c.melt.Al = 11.79;
  c.melt.Fe = 0.03;
  c.melt.Mn = 0.02;
  c.melt.Mg = 4.80;
  c.melt.Ca = 9.73;
  c.melt.Na = 3.41;
  c.melt.K = 0.59;
  c.melt.P = 0.05;
  c.melt.H = 6.80;
  
  PetscScalar T = 1273;
  PetscScalar eta;
  eta = giordanoViscosity( &c,  T );
  printf("eta=%e\n Pa-s\n",eta);
  PetscScalar phi = 0.1;
  printf("eta at phi=%lf = %e\n",phi,einsteinRoscoe( eta, phi ));

  printf("loading lookup table from disk\n");
  char filename[80] = "BAPhaseDiagrams/Las2HG.txt";
  CompositionLookupTableBA clt;

  clt = loadCompositionLookupTableBAFromFile( filename );
  PetscScalar Ttest = 806;
  printf("looking up composition at T=%e\n",Ttest);
  lookupCompositionBA( &clt, Ttest, &c);
  printf("solid fraction %lf\n",c.phi);

  printComposition( &c );

  return 0;
}


