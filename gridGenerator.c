#include "fdcode.h"
#include "gridSpacing.h"
#include "gridGenerator.h"

PetscErrorCode allocateGrid(GridData *, Options *);

PetscErrorCode initializeRegularGrid(GridData *grid, Options *options){
  PetscFunctionBegin;
  PetscInt NX = options->NX;
  PetscInt NY = options->NY;
  PetscScalar LX = options->LX;
  PetscScalar LY = options->LY;
  PetscErrorCode ierr=0;

  ierr = allocateGrid( grid, options);CHKERRQ(ierr);
  if( options->mechBCLeft.type[0] == 2){    
    grid->xperiodic = 1;
    gridSpacingUniform( grid->x, LX, NX+1 );
  } else{
    gridSpacingUniform( grid->x, LX, NX );
  }

  gridSpacingUniform( grid->y, LY, NY );
  getCellCenters( grid->x, grid->xc, LX, NX);
  getCellCenters( grid->y, grid->yc, LY, NY);
  if( grid->xperiodic ){
    makePeriodic( grid->x, grid->xc, LX, NX);
  }  
  PetscFunctionReturn(ierr);
}

PetscErrorCode initializeIrregularGridConstantInnerOuter(GridData *grid, Options *options){
  PetscFunctionBegin;
  PetscInt NX = options->NX;
  PetscInt NY = options->NY;
  PetscScalar LX = options->LX;
  PetscScalar LY = options->LY;
  PetscErrorCode ierr=0;
  PetscScalar hcontrast = options->gridRatio;
  ierr = allocateGrid( grid, options);CHKERRQ(ierr);
  if( options->mechBCLeft.type[0] == 2){    
    grid->xperiodic = 1;
  }  
  if( grid->xperiodic ){
    gridSpacingUniform( grid->x, LX, NX+1 );
  }else{
    gridSpacingConstantInnerOuter( grid->x, LX, NX ,hcontrast);
  }
  gridSpacingConstantInnerOuter( grid->y, LY, NY ,hcontrast);
  getCellCenters( grid->x, grid->xc, LX, NX);
  getCellCenters( grid->y, grid->yc, LY, NY);
  if( grid->xperiodic ){
    makePeriodic( grid->x, grid->xc, LX, NX);
  }  
  PetscFunctionReturn(ierr);
}



PetscErrorCode allocateGrid(GridData *grid, Options *options){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscInt NX = options->NX;
  PetscInt NY = options->NY;
  PetscScalar LX = options->LX;
  PetscScalar LY = options->LY;
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);
  grid->x += 2; /* increment pointer by 2 elements to allow for ghost region */
  grid->xc += 2;
  grid->y += 2;
  grid->yc += 2;
  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;
  grid->xperiodic = 0;
  grid->yperiodic = 0;
  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyGrid(GridData *grid){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  grid->x -= 2; grid->y -=2; grid->xc -= 2; grid->yc -= 2;
  ierr = PetscFree(grid->x);CHKERRQ(ierr);
  ierr = PetscFree(grid->y);CHKERRQ(ierr);
  ierr = PetscFree(grid->xc);CHKERRQ(ierr);
  ierr = PetscFree(grid->yc);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}
