#include "fdcode.h"
#include "pressureNullSpace.h"
/* This routine sets up a nullspace for the pressure solution */

PetscErrorCode createPressureNullSpace( GridData *grid, MatNullSpace *ns){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscInt nvec = 1; /* size of null space */
  Vec V; /* null space */
  ierr = DMCreateGlobalVector( grid->vda, &V );CHKERRQ(ierr);
  ierr = VecZeroEntries( V ); CHKERRQ(ierr);
  PetscInt x,y,m,n;
  PetscScalar ***v;
  ierr = DMDAVecGetArrayDOF(grid->vda, V, &v );CHKERRQ(ierr);
  PetscInt ix,jy;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL); CHKERRQ(ierr);
  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){
/*       if(0 && ((!grid->xperiodic && ix== 0) || jy == 0) ){ */
/*       }else */ {
	v[jy][ix][DOF_P] = 1.0;
      }
    }
  }
  ierr = DMDAVecRestoreArrayDOF(grid->vda, V, &v);CHKERRQ(ierr);
  ierr = VecNormalize( V, PETSC_NULL );CHKERRQ(ierr);

  ierr = MatNullSpaceCreate( PETSC_COMM_WORLD, PETSC_FALSE, nvec, &V, ns);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyPressureNullSpace( MatNullSpace *ns ){
  PetscErrorCode ierr;
  ierr=MatNullSpaceDestroy(ns); CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}
