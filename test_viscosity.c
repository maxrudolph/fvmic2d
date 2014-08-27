#include "petscksp.h"
#include "fdcode.h"
#include "viscosity.h"

inline PetscScalar goldsbyKohlstedt(PetscScalar , PetscScalar , PetscScalar , PetscScalar );


int main(){
  /* test the goldsby-kohlstedt 2001 viscosity implementation */
  PetscScalar pressure = -1.6541e5;
  PetscScalar sii = 2.0740e2;
  PetscScalar d = 1e-3;
  PetscScalar T = 110.68; 
  PetscScalar etaeff;
  etaeff=goldsbyKohlstedt(sii, T, d, pressure);
  printf("Stress (Pa)\t Strain-rate (1/s)\t etaeff(Pa-s)\n");
  printf("%e\t%e\t%e\n",sii,sii/2.0/etaeff,etaeff);

}
