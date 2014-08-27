#include "fdcode.h"

/* this routine takes the markers and projects all fields onto a regular grid for visualization*/
/* NOT PARALLELIZED */

PetscErrorCode gridMarkerField(Markers *markers, GridData *grid, PetscScalar *fieldToGrid, PetscScalar *griddedField, const PetscInt NX, const PetscInt NY){
  PetscErrorCode ierr;

  PetscInt m;
  PetscScalar LX = grid->LX;
  PetscScalar LY = grid->LY;

  PetscScalar dx = LX/((PetscScalar) NX-1.0);
  PetscScalar dy = LY/((PetscScalar) NY-1.0);

  PetscScalar *weights;
  ierr = PetscMalloc( NX*NY*sizeof(PetscScalar), &weights);CHKERRQ(ierr);

  /* begin by zeroing out griddedField*/
  PetscInt i;
  for(i=0;i<NX*NY;i++){
    griddedField[i] = 0.0;
  }

  for(m=0;m<markers->nMark;m++){
    /* find closest node and add value*/
    PetscInt mj = (markers->X[m] + dx/2.0)/dx;
    PetscInt mi = (markers->Y[m] + dy/2.0)/dy;
    if(mi < 0){ mi = 0;} else if(mi > NY-1) {mi = NY-1;}
    if(mj < 0){ mj = 0;} else if(mj > NX-1) {mj = NX-1;}
    PetscScalar mwt = fabs( markers->X[m] - (mj+0.5)*dx) * fabs( markers->Y[m] - (mi+0.5)*dy);
    griddedField[mi + NY*mj] += fieldToGrid[m]*mwt;
    weights[mi + NY*mj] += mwt;
  }
  for(i=0;i<NX*NY;i++){
    if( weights[i] != 0.0){
      griddedField[i]/=weights[i];
    }else{
      griddedField[i] = -1.0;
    }
  }


  ierr=PetscFree(weights);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
