#include "fdcode.h"
#include "reggrid.h"
/* set up a regular grid.*/

PetscErrorCode initializeIrregularGrid(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY){
  PetscErrorCode ierr=0;
  //  PetscScalar dx,dy,dx2,dy2,dxdy;

  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);


  grid->dx=LX/((PetscScalar) NX - 1);
  grid->dy=LY/((PetscScalar) NY - 1);
  grid->dx2 = (grid->dx*grid->dx);
  grid->dy2 = grid->dy*grid->dy;
  grid->dxdy=grid->dx*grid->dy;
  grid->NX = NX;
  grid->NY= NY;
  grid->LX=LX;
  grid->LY=LY;

  PetscInt i,j;
  /* basic node gridlines*/
  for(i=0;i<NY;i++){
    grid->y[i] = 0+grid->dy*((PetscScalar) i);
  }
  for(j=0;j<NX;j++){
    grid->x[j] = 0+grid->dx*((PetscScalar) j);
  }
  /* cell center locations*/
  for(i=0;i<NY;i++){
    grid->yc[i] = grid->dy*(-0.5+((PetscScalar) i));
  }
  for(j=0;j<NX;j++){
    grid->xc[j] = grid->dx*(-0.5+((PetscScalar) j));
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyRegularGrid(GridData *grid){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscFree(grid->x);CHKERRQ(ierr);
  ierr = PetscFree(grid->y);CHKERRQ(ierr);
  ierr = PetscFree(grid->xc);CHKERRQ(ierr);
  ierr = PetscFree(grid->yc);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}

/* routine to set locations of basic nodes in DA for spatial fields*/

