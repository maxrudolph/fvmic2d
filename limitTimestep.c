#include "fdcode.h"
#include "limitTimestep.h"

/*const PetscScalar dispFactor=0.1;*//* maximum fraction of a cell dimension that a marker is allowed to travel in a time step*/

PetscErrorCode limitDisplacementTimestep(GridData *grid, NodalFields *nodalFields, PetscScalar *displacementdt, Options *options){
  PetscErrorCode ierr;
  PetscScalar dispFactor = options->displacementStepLimit;
  
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // PetscScalar lastdt = displacementdt[0];
  PetscInt NX=grid->NX;
  PetscInt NY=grid->NY;
  /* retrieve local portions of vx,vy*/
  PetscInt x,y,m,n;
  PetscFunctionBegin;  
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  Vec vxl, vyl;
  ierr=DMCreateLocalVector(grid->da,&vxl);
  ierr=VecDuplicate(vxl,&vyl);
  //  ierr=VecDuplicate(vxl,&vzl);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);

  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr); */
  PetscScalar **vx, **vy, **vz;  
  ierr=DMDAVecGetArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vyl,&vy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,vzl,&vz);CHKERRQ(ierr); */

  /* get local part of rhodot*/
  PetscScalar **rhodot;
  PetscScalar **rho;
  ierr=DMDAVecGetArray(grid->da, nodalFields->rhodot, &rhodot);
  ierr=DMDAVecGetArray(grid->da, nodalFields->rho, &rho);

  /* loop over local portions of vx and vy*/
  PetscInt ix,jy;
  PetscScalar dx,dy;
  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){
      PetscScalar dx1=grid->x[ix+1]-grid->x[ix];
      PetscScalar dx2=grid->x[ix]-grid->x[ix-1];
      if(ix==0){
	dx=dx1;
      }else if(ix==NX-1){
	dx=dx2;
      }else if(dx1 <= dx2 ){
	dx=dx1;
      }else{
	dx=dx2;
      }
      PetscScalar dt1;
      if( vx[jy][ix] != 0.0){
	dt1=dx/(fabs(vx[jy][ix]))*dispFactor;
      }else{
	dt1=displacementdt[0];
      }
      if( displacementdt[0] > dt1 ) displacementdt[0] = dt1;

      /* vy */
      
      PetscScalar dy1=grid->y[jy+1] - grid->y[jy];
      PetscScalar dy2=grid->y[jy] - grid->y[jy-1];
      if(jy==0){
	dy=dy1;
      }else if(jy==NY-1){
	dy=dy2;
      }else if(dy1 <= dy2){
	dy= dy1;
      }else{
	dy=dy2;
      }
      if( vy[jy][ix] != 0.0){
	dt1=dy/(fabs(vy[jy][ix]))*dispFactor;
      }else{
	dt1=displacementdt[0];
      }
      if( displacementdt[0] > dt1 ) displacementdt[0] = dt1;

      /* vz */
      /*       if( jy > 0 && ix > 0){ */
      /* 	if( vz[jy][ix] != 0.0){ */
      /* 	  dt1 = sqrt(dx*dx+dy*dy)/(fabs(vz[jy][ix]))*dispFactor; */
      /* 	}else{ */
      /* 	  dt1 = displacementdt[0]; */
      /* 	} */
      /* 	if( displacementdt[0] > dt1 ) displacementdt[0] = dt1; */
      /*       } */

      if(rhodot[jy][ix] != 0.0){
	dt1 = fabs(0.01*rho[jy][ix]/rhodot[jy][ix]); /* 0.01 is 1 % density change */
      }
      if( displacementdt[0] > dt1 ) displacementdt[0] = dt1;

     
    }
  }
  /* limit timestep based on density changes*/
  /* calculate maximum nodal density change */
  PetscScalar rhodotmax, rhodotmin;
  ierr = VecMax( nodalFields->rhodot, PETSC_NULL, &rhodotmax );
  ierr = VecMin( nodalFields->rhodot, PETSC_NULL, &rhodotmin );
  rhodotmax = fabs(rhodotmax);
  rhodotmin = fabs(rhodotmin);
  if(rhodotmin>rhodotmax) rhodotmax = rhodotmin;
  if( rhodotmax*displacementdt[0] > 0.01 ) displacementdt[0] = 0.01/rhodotmax;

  PetscScalar dtmin;
  
  ierr=MPI_Allreduce( displacementdt, &dtmin, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);CHKERRQ(ierr); 
  displacementdt[0] = dtmin;
   
  if(!rank) printf("Current displacement timestep %e\n",displacementdt[0]);


  ierr=DMDAVecRestoreArray(grid->da,nodalFields->rhodot,&rhodot);
  ierr=DMDAVecRestoreArray(grid->da,nodalFields->rho,&rho);
  ierr=DMDAVecRestoreArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  //  ierr=DMDAVecRestoreArray( grid->da,vzl,&vz);CHKERRQ(ierr);

  ierr=VecDestroy(&vxl);CHKERRQ(ierr);
  ierr=VecDestroy(&vyl);CHKERRQ(ierr);
  //  ierr=VecDestroy(&vzl);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}


