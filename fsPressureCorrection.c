#include "fdcode.h"
#include "fsPressureCorrection.h"
#include "markerProjection.h"

/* This file contains routines to detect location of free surface and correct the pressure field accordingly */

PetscErrorCode freeSurfacePressureCorrection( GridData *grid, MarkerSet *markerset, NodalFields *nodalFields, Materials *materials, PetscInt airmat ){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  /* pull out density field */
  Vec rhol, pl;
  PetscScalar **rho, **p;
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


  ierr = DMGetLocalVector( grid->da, &rhol);CHKERRQ(ierr);
  ierr = DMGetLocalVector( grid->da, &pl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->etaS, INSERT_VALUES, rhol);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->etaS, INSERT_VALUES, rhol);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->p, INSERT_VALUES, pl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->p, INSERT_VALUES, pl);CHKERRQ(ierr);
  ierr = DMDAVecGetArray( grid->da, rhol, &rho);CHKERRQ(ierr);
  ierr = DMDAVecGetArray( grid->da, pl, &p);CHKERRQ(ierr);

  Vec mat;
  /* get current viscosity of sticky air material (arithmetic average over sticky air markers) */
  PetscScalar etam_mark;
  {
    PetscScalar etam = 0.0;
    PetscInt nAir = 0;
    PetscInt m;
    Marker *markers = markerset->markers;
    /* count number of locally owned air markers */
    for( m=0;m<markerset->nMark;m++){
      if( markers[m].cellX != -1 && markers[m].Mat == airmat ){
	nAir ++;
      }
    }
    /* if there are any locally owned air particles, continue averaging viscosity */
    for( m=0;m<markerset->nMark;m++){
       if( markers[m].cellX != -1 && markers[m].Mat == airmat ){
	 etam += markers[m].eta/((PetscScalar) nAir);
       }
    }
    /* allreduce both the average and total number */
    PetscScalar nAirLocal = (PetscScalar) nAir;
    PetscScalar nAirGlobal;
    MPI_Allreduce( &nAirLocal , &nAirGlobal, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD );
    /* weight local average by nLocal/nGlobal */
    etam *= nAirLocal/nAirGlobal;    
    
    MPI_Allreduce( &etam , &etam_mark, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD );    
    //printf("Current sticky air viscosity is %le\n",etam_mark);
  }



  /* identify point where density increases 2X from sticky air density */
  PetscScalar rhoair = etam_mark;
  const PetscScalar rhoc = 10.0*rhoair;
  /* look at locally owned columns of density field */
  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL); CHKERRQ(ierr);
  if( y+n > grid->NY-1 ){
    n --;
  }
  PetscInt ix, jy;
  
  /* accumulators for mean pressure at free surface */
  PetscScalar yfs = 0.0; /* free surface depth */
  PetscScalar pfs = 0.0; /* free surface pressure */
  PetscScalar nLocal = 0;
  for(ix=x;ix<x+m;ix++){ /* loop over columns */
    jy = y-1;
    if( jy < 0) jy = 0;
    
    if( ix == 0 || rho[y][ix] >= rhoc || rho[y+n-1][ix] < rhoc ){  
      /* free surface is not locally owned */
    }else{
      nLocal += 1.0;
      while( jy < y+n && rho[jy+1][ix] < rhoc ){ jy++;}
      
      /* free surface should now be between jy and jy+1 */
      /* assume that rho changes linearly between the two nodes */
      PetscScalar rho0 = rho[jy][ix];
      PetscScalar rho1 = rho[jy+1][ix];
      PetscScalar dy = grid->y[jy+1] - grid->y[jy];
      PetscScalar yfsl = (rhoc-rho0)*dy / (rho1-rho0)+grid->y[jy];
      
      if( rho[jy+1][ix] < rhoc ) yfsl = 0.0;
      if( yfsl > grid->yc[jy+1] ) jy++;

      //PetscScalar p0 = p[jy][ix];
      //PetscScalar p1 = p[jy+1][ix];

      //PetscScalar pfsl = p0 + (p1-p0) * (yfsl - grid->yc[jy])/(grid->yc[jy+1]-grid->yc[jy]);      
     
      PetscScalar pfsl = p[jy+2][ix];

      //printf("found free surface idx %d, location %le, pressure %le\n",jy,yfsl,pfsl);
      pfs += pfsl;
      yfs += yfsl;
    }    
  }
  
  { /* allreduce on nLocal and on yfs */
    PetscScalar nTot;
    PetscScalar yfs_sum;
    PetscScalar p_sum;
    MPI_Allreduce( &nLocal , &nTot, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD ) ;
    MPI_Allreduce( &yfs , &yfs_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD ) ;
    MPI_Allreduce( &pfs , &p_sum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD ) ;
    if( nTot > 0.0 ){
      yfs = yfs_sum/nTot;
      pfs = p_sum/nTot;
    }else{
      yfs = 0.0;
      pfs = 0.0;
    }
  }
  if(!rank){
    printf("Free surface pressure correction %e\n, depth %e\n",pfs,yfs);
  }
  PetscScalar pshift = 1.0e5 - pfs;
  
  /* shift the pressure vector */
  ierr = VecShift( nodalFields->p, pshift);CHKERRQ(ierr);
    
  ierr = DMDAVecRestoreArray(grid->da, rhol, &rho);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da, pl, &p);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector( grid->da, &rhol); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector( grid->da, &pl); CHKERRQ(ierr);
 
  PetscFunctionReturn(ierr);
}
