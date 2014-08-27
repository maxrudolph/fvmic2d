#include "fdcode.h"
#include "reggrid.h"
#include "gridSpacing.h"
/* set up a regular grid.*/

PetscErrorCode initializeRegularGrid(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY, Options *options){
  PetscErrorCode ierr=0;
  //  PetscScalar dx,dy,dx2,dy2,dxdy;
  PetscFunctionBegin;
  grid->xperiodic = 0;

  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);
  grid->x += 2; /* increment pointer by 2 elements */
  grid->xc += 2;
  grid->y += 2;
  grid->yc += 2;

  grid->NX = NX;
  grid->NY= NY;
  grid->LX=LX;
  grid->LY=LY;

  /* check for periodic grid */
  if( options->mechBCLeft.type[0] == 2){    
    grid->xperiodic = 1;
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX ++;
  }

  PetscScalar dx=LX/((PetscScalar) NX - 1.0);
  PetscScalar dy=LY/((PetscScalar) NY - 1.0);

  PetscInt i,j;
  /* basic node gridlines*/
  for(i=0;i<NY;i++){
    grid->y[i] = 0+dy*((PetscScalar) i);
  }
  for(j=0;j<NX;j++){
    grid->x[j] = 0+dx*((PetscScalar) j);
  }
  /* cell center locations*/
  for(i=1;i<NY;i++){
    grid->yc[i] = (grid->y[i-1]+grid->y[i])/2.0;
  }
  for(j=1;j<NX;j++){
    grid->xc[j] = (grid->x[j-1]+grid->x[j])/2.0;
  }

  grid->xc[0] = -grid->xc[1];
  grid->yc[0] = -grid->yc[1];
  grid->xc[NX] = 2.0*grid->LX - grid->xc[NX-1];
  grid->yc[NY] = 2.0*grid->LY - grid->yc[NY-1];
  
  /* check to see if DA is periodic */
  if( grid->xperiodic){
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
#ifdef DEBUG
    for(j=-2;j<NX+2;j++){
      printf("j=%d x=%e\n",j,grid->x[j]);
    }
#endif
  }else{
  }  

  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyRegularGrid(GridData *grid){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  grid->x -= 2; /* increment pointer by 2 elements */
  grid->xc -= 2;
  grid->y -= 2;
  grid->yc -= 2;
  ierr = PetscFree(grid->x);CHKERRQ(ierr);
  ierr = PetscFree(grid->y);CHKERRQ(ierr);
  ierr = PetscFree(grid->xc);CHKERRQ(ierr);
  ierr = PetscFree(grid->yc);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}

/* routine to set locations of basic nodes in DA for spatial fields*/

