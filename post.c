#include "fdcode.h"
#include "post.h"

/* routines for processing results */
PetscErrorCode nusseltNumber( GridData *grid, Options *options, Vec Tg , PetscScalar *Nu){
  /* use Blankenbach 1989 definition of Nusselt number */
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscScalar dT = options->thermalBCBottom.value[0]-options->thermalBCTop.value[0];
  PetscScalar LX = grid->LX;
  PetscScalar LY = grid->LY;
  
  /* get a local temperature vector */
  Vec Tl;
  PetscScalar **T;
  ierr=DMCreateLocalVector(grid->da, &Tl); CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,Tg,INSERT_VALUES,Tl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,Tg,INSERT_VALUES,Tl);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,Tl,&T);CHKERRQ(ierr);


  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  /* loop over cells */
  if(!grid->xperiodic &&  x+m == grid->NX ) m--;
  PetscInt ix;
  
  PetscScalar integrand = 0.0;
  if( y == 0 ){/* this CPU owns surface */
    for(ix=x;ix<x+m;ix++){
      PetscScalar dT = (T[1][ix]-T[0][ix] + T[1][ix+1] - T[0][ix+1])/2.0;
      PetscScalar dz = grid->y[1]-grid->y[0];
      PetscScalar dx = grid->x[ix+1]-grid->x[ix];
      integrand += dT/dz*dx;
    }   
  }

  /* do an allreduce on the integrand */
  PetscScalar integral=0.0;
  MPI_Allreduce( &integrand, &integral, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD );
  Nu[0] = LY/LX/dT*integral;
  ierr=DMDAVecRestoreArray(grid->da,Tl,&T);CHKERRQ(ierr);
  ierr=VecDestroy(&Tl); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode rmsVelocity( GridData *grid,  NodalFields *nodalFields, PetscScalar *vrms){
  /* get velocity from nodal Fields */
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Vec vxl,vyl;
  PetscScalar **vx, **vy;

  ierr = DMCreateLocalVector(grid->da,&vxl);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(grid->da,&vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr = DMDAVecGetArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr = DMDAVecGetArray( grid->da,vyl,&vy);CHKERRQ(ierr);

  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  /* loop over cells */
  if(!grid->xperiodic &&  x+m == grid->NX ) m--;
  if( y+n == grid->NY ) n--;
  PetscInt ix,jy;
  PetscScalar integrand = 0.0;
  for(jy=y;jy<y+n;jy++){
    for(ix=x;ix<x+m;ix++){
      PetscScalar dx=grid->x[ix+1]-grid->x[ix];
      PetscScalar dy=grid->y[jy+1]-grid->y[jy];
      PetscScalar vx2 = (vx[jy][ix]*vx[jy][ix]+vx[jy][ix+1]*vx[jy][ix+1])/2.0;
      PetscScalar vy2 = (vy[jy][ix]*vy[jy][ix]+vy[jy+1][ix]*vy[jy+1][ix])/2.0;
      integrand += (vx2+vy2)*dx*dy;
    }
  }
  PetscScalar integral=0.0;
  MPI_Allreduce( &integrand, &integral, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD );
  
  integral /= grid->LX*grid->LY;
  integral = sqrt(integral);

  vrms[0] = integral;

  ierr = DMDAVecRestoreArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=VecDestroy(&vxl);CHKERRQ(ierr);
  ierr=VecDestroy(&vyl);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}
