#include "fdcode.h"

PetscErrorCode updateShearHeating(NodalFields *nodalFields, GridData *grid, Vec Hs){
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  PetscScalar *exx = nodalFields->edotxx;
  PetscScalar *exy = nodalFields->edotxy;
  PetscScalar *sdevxx = nodalFields->sdevxx;
  PetscScalar *sdevxy = nodalFields->sdevxy;
  PetscInt i,j;
  PetscErrorCode ierr;
  //ierr = VecZeroEntries( Hs);CHKERRQ(ierr);
  PetscFunctionBegin;
  for(i=1;i<NY-1;i++){
    for(j=1;j<NX-1;j++){
      PetscInt idxnode = j*NY+i;
      PetscScalar nodalexx = 0.25*(exx[idxnode] + exx[idxnode+1] + exx[idxnode+NY] + exx[idxnode+NY+1]);
      PetscScalar nodalsdevxx = 0.25*(sdevxx[idxnode] + sdevxx[idxnode+1] + sdevxx[idxnode+NY] + sdevxx[idxnode+NY+1]);
      ierr = VecSetValue( Hs, idxnode, 2.0*nodalsdevxx*nodalexx + 2.0*sdevxy[idxnode]*exy[idxnode], INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = VecAssemblyBegin( Hs);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(Hs);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}
