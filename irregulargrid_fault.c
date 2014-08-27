/* This is an irregular grid for studying strike slip faults. 
It has high resolution in the middle and low resolution at the boundaries in the x-direction.
In the y direction, the point of maximum resolution can be specified. */
#include "fdcode.h"
#include "grid.h"
#include "gridSpacing.h"

PetscErrorCode initializeIrregularGridFault(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY, Options *options){

  PetscScalar hcontrast = options->gridRatio;
  PetscErrorCode ierr=0;
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);
  grid->x += 2; /* increment pointer by 2 elements */
  grid->xc += 2;
  grid->y += 2;
  grid->yc += 2;


  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;
  grid->xperiodic = 0;
  grid->yperiodic = 0;

  PetscScalar f=1.2407234;
  PetscInt ix,jy;
  

  /* x gridlines */
  /* check to see if DA is periodic */
  if( options->mechBCLeft.type[0] == 2){
    if( NX%2){
      printf("Need even NX for periodic X boundary conditions\n");
      abort();
    }
    grid->xperiodic = 1;
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX ++;
  }else{
  }

  PetscScalar L = LX/2.0;
  PetscInt NE = (NX-1)/2;
  PetscInt Nr=NE/hcontrast;/* number of elements in refined region */
  PetscInt Nc=NE-Nr-5;
  PetscScalar a = L/(10.0+Nr+Nc*hcontrast);

  grid->x[0] = 0.0;
  for(ix=1;ix<=Nc;ix++){
    grid->x[ix] = grid->x[ix-1] + a*hcontrast;
  }
  for(ix=Nc+1;ix<=Nc+5;ix++){
    grid->x[ix] = grid->x[ix-1] + a*pow(f,6.0-(ix-Nc));
  }
  for(ix=Nc+6;ix<=Nc+Nr+6;ix++){
    grid->x[ix] = grid->x[ix-1] + a;
  }
  
  /* mirror other half of grid */
  grid->x[NE] = LX/2;
  for(ix=(NX+1)/2;ix<NX-1;ix++){
    grid->x[ix] = LX - grid->x[NX-1 - ix];
  }
  grid->x[NX-1] = LX;
  
  /* y-gridlines */
  L = LY/2.0;
  NE = (NY-1)/2;
  Nr=NE/hcontrast;/* number of elements in refined region */
  Nc=NE-Nr-5;
  a = L/(10.0+Nr+Nc*hcontrast);

  /* uniformly spaced y nodes */
  grid->y[0] = 0.0;
  for(jy=1;jy<NY-1;jy++){
    grid->y[jy] = grid->y[jy-1] + LY/(NY-1);
  }
  grid->y[NY-1] = LY;
  
  /* cell center locations*/
  for(jy=1;jy<NY;jy++){
    grid->yc[jy] = (grid->y[jy-1]+grid->y[jy])/2.0;
  }
  grid->yc[0] = - grid->yc[1];
  grid->yc[NY] = 2.0*grid->LY - grid->yc[NY-1];

  for(ix=1;ix<NX;ix++){
    grid->xc[ix] = (grid->x[ix-1]+grid->x[ix])/2.0;
  }
  grid->xc[0] = -grid->xc[1];
  grid->xc[NX] = 2.0*grid->LX - grid->xc[NX-1];


  /* check to see if DA is periodic */
  if( grid->xperiodic ){
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX --;
    grid->x[NX] = LX;
    grid->x[NX+1] = LX+grid->x[1];
    grid->x[-1] = grid->x[NX-1]-LX;
    grid->x[-2] = grid->x[NX-2]-LX;
    grid->xc[NX] = (LX + grid->x[NX-1])/2.0;
    grid->xc[NX+1] = grid->xc[1]+LX;
    grid->xc[0] = 0.0-(LX-grid->xc[NX]);
    grid->xc[-1] = 0.0-(LX-grid->xc[NX-1]);
  }else{
  }  

  PetscFunctionReturn(ierr);
}


/* this mesh has also enhanced resolution near the free-surface */
PetscErrorCode initializeIrregularGridFaultFS(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY, Options *options, PetscScalar fsdepth){

  PetscScalar hcontrast = options->gridRatio;
  PetscErrorCode ierr=0;
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);
  grid->x += 2; /* increment pointer by 2 elements */
  grid->xc += 2;
  grid->y += 2;
  grid->yc += 2;

  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;
  grid->xperiodic = 0;
  grid->yperiodic = 0;

  /* x gridlines */
  /* check to see if DA is periodic */
  if( options->mechBCLeft.type[0] == 2){
    if( NX%2){
      printf("Need even NX for periodic X boundary conditions\n");
      abort();
    }
    grid->xperiodic = 1;
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX ++;
  }else{
  }

  /* x-grid should be uniform */
  gridSpacingUniform( grid->x, LX, NX );

  /* y grid should have constant inner outer spacing */
  gridSpacingRefinedCenter( grid->y, LY, NY, hcontrast );
  /* look for first place where resolution changes and record position */
  PetscInt ref_start=0;
  PetscInt ref_end=NY-1;
  PetscScalar small = 1e-12*grid->LY;
  PetscInt ix = 0;
  for(ix=2;ix<NY;ix++){
    if( fabs((grid->y[ix-1]-grid->y[ix-2])-(grid->y[ix]-grid->y[ix-1])) > small){
      ref_start = ix;
      break;
    }
  }
  for(ix=NY-3;ix>=0;ix--){
    if( fabs((grid->y[ix+2]-grid->y[ix+1])-(grid->y[ix+1]-grid->y[ix])) > small){
      ref_end = ix;
      break;
    }
  }

  /* calculate refined region center */
  PetscScalar ref_center = (2*grid->y[ref_start] + grid->y[ref_end])/3.0;
  ix = 0;
  while(ix < ref_start){
    if( grid->y[ix] > (ref_center - (LY*fsdepth)) ) break;
    ix++;
  }
  circShift( grid->y, LY, NY, -ix );

  getCellCenters( grid->x, grid->xc, LX, NX );
  getCellCenters( grid->y, grid->yc, LY, NY );

  grid->xc[0] = -grid->xc[1];
  grid->yc[0] = -grid->yc[1];
  grid->xc[NX] = 2.0*grid->LX - grid->xc[NX-1];
  grid->yc[NY] = 2.0*grid->LY - grid->yc[NY-1];

  /* check to see if DA is periodic */
  if( grid->xperiodic ){
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX --;
    makePeriodic( grid->x, grid->xc, LX, NX);
  }else{
  }  

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
