#include "fdcode.h"
#include "retrieveSolutions.h"

PetscErrorCode retrieveSolutions( GridData *grid,NodalFields *nodalFields, Vec S, Vec Sz, PetscScalar Kcont, PetscInt *foundnan){
  PetscErrorCode ierr;
  
  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  
  ierr = VecStrideScale(S, DOF_P, Kcont);CHKERRQ(ierr); /* Scale pressure field */

  /* retrieve velocity and pressure from mechanical solutions */
  PetscScalar ***Sa;
  PetscScalar **Sza;
  PetscScalar **vx;
  PetscScalar **vy;
  PetscScalar **vz;
  PetscScalar **p;
  ierr = DMDAVecGetArrayDOF(grid->vda,S,&Sa);CHKERRQ(ierr);
#ifndef TEXTURE
  ierr = DMDAVecGetArray(grid->da,Sz,&Sza);CHKERRQ(ierr);
#endif
  ierr = DMDAVecGetArray(grid->da,nodalFields->vx,&vx);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da,nodalFields->vy,&vy);CHKERRQ(ierr);
#ifndef TEXTURE
  ierr = DMDAVecGetArray(grid->da,nodalFields->vz,&vz);CHKERRQ(ierr);
#endif
  ierr = DMDAVecGetArray(grid->da,nodalFields->p,&p);CHKERRQ(ierr);

  PetscInt i,j;
  PetscInt havenan=0;/* flag if we find a nan in the solution*/
  
  for(j=y;j<y+n;j++){
    for(i=x;i<x+m;i++){
      p[j][i] = Sa[j][i][DOF_P];
      vx[j][i] = Sa[j][i][DOF_U];
      vy[j][i] = Sa[j][i][DOF_V];
#ifndef TEXTURE
      /*       vz[idxnode]=Sza[idxnode]; */
      vz[j][i] = Sza[j][i];
      
      if(isnan(vx[j][i]) || isnan(p[j][i]) || isnan(vy[j][i]) || isnan(vz[j][i])){
	havenan=1;
      }
#else
      if(isnan(vx[j][i]) || isnan(p[j][i]) || isnan(vy[j][i])){
	havenan=1;
      }
#endif
    }
  }
  /* correct the pressure solution to zero top pressure */
  PetscScalar p_mean = 0.0;
  for(j=y;j<y+n;j++){
    if( j == 1 ){
      for(i=x;i<x+m;i++){
	if(i > 0 || grid->xperiodic ){
	  PetscScalar dx = grid->x[i]-grid->x[i-1];
	  p_mean += p[j][i]*dx/grid->LX;	 
	}	
      }
    }
  }

  PetscScalar p_meang;
  ierr = MPI_Allreduce( &p_mean, &p_meang, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

  for(j=y;j<y+n;j++){
    if( j > 0 ){
      for(i=x;i<x+m;i++){
	if(i > 0 || grid->xperiodic ){
	  p[j][i] -= p_meang;
	}	
      }
    }
  }

  /* mpi allreduce on foundnan*/
  ierr = MPI_Allreduce( &havenan, foundnan, 1, MPI_INT, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
  
  ierr = DMDAVecRestoreArrayDOF(grid->vda,S,&Sa);CHKERRQ(ierr);
#ifndef TEXTURE
  ierr = DMDAVecRestoreArray(grid->da,Sz,&Sza);CHKERRQ(ierr);
#endif
  ierr = DMDAVecRestoreArray(grid->da,nodalFields->vx,&vx);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,nodalFields->vy,&vy);CHKERRQ(ierr);
#ifndef TEXTURE
  ierr = DMDAVecRestoreArray(grid->da,nodalFields->vz,&vz);CHKERRQ(ierr);
#endif
  ierr = DMDAVecRestoreArray(grid->da,nodalFields->p,&p);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
