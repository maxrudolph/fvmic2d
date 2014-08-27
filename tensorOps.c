#include "fdcode.h"
#include "tensorOps.h"

PetscScalar IIdev( Tensor33s *T ){
  PetscScalar xx = T->T11;
  PetscScalar yy = T->T22;
  PetscScalar zz = T->T33;
  PetscScalar xy = T->T12;
  PetscScalar xz = T->T13;
  PetscScalar yz = T->T23;
  PetscScalar tr = (xx+yy+zz)/3.0;
  xx -= tr;
  yy -= tr;
  zz -= tr;
  PetscScalar dii = sqrt( 0.5*(xx*xx+yy*yy+zz*zz) + xy*xy + xz*xz + yz*yz);

  return dii;

}
