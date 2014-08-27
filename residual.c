/* This file contains subroutines related to tracking a residual during newton iteration */

#include "fdcode.h"
#include "viscosity.h"
#include "markers.h"
#include "residual.h"
#include "markerProjection.h"

#define NRES 100 /* maximum number of residuals to store */
PetscScalar residuals[NRES];

PetscErrorCode updateGlobalStrainRateResidual( GridData *grid, MarkerSet *markerset, PetscInt iter ){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  Marker *markers = markerset->markers;
  /* calculate the current global strain rate residual */
  Vec srrg;
  ierr = DMCreateGlobalVector(grid->da, &srrg);

  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &(markers[0].strainRateResidual), srrg, -1,  -1, ARITHMETIC);
  /* calculate cell contribution to global area = center_value*dx*dy */
  PetscInt x,y,m,n,ix,jy;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  PetscScalar **srr;
  ierr = DMDAVecGetArray( grid->da, srrg, &srr );CHKERRQ(ierr);
  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){
      if( ix > 0 && jy > 0){
	PetscScalar dx = grid->x[ix]-grid->x[ix-1];
	PetscScalar dy = grid->y[jy]-grid->y[jy-1];
	srr[jy][ix] *= dx*dy;
      }else{
	srr[jy][ix] = 0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray( grid->da, srrg, &srr );CHKERRQ(ierr);
  /* calculate global residual by summation of cell residuals */
  PetscScalar srrsum = 0.0;
  ierr = VecSum( srrg, &srrsum);CHKERRQ(ierr);

  /* normalize global residual by area (LX*LY) */
  srrsum /= (grid->LX*grid->LY);

  /* store residual  */
  setGlobalStrainRateResidual( iter, srrsum );

  /* clean up */
  ierr=VecDestroy(&srrg);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}

PetscInt checkConvergence( PetscScalar atol, PetscScalar rtol, PetscInt iter){
  /* this function checks to see whether absolute or relative error tolerances have been met */
  /* check absolute tolerance */
  if( getGlobalStrainRateResidual( iter ) < atol ){
    return( 1 );
  }
  /* check relative tolerance */
  if( iter > 0 ){/* only for nonzero iter */
    PetscScalar thisres = getGlobalStrainRateResidual( iter   );
    PetscScalar lastres = getGlobalStrainRateResidual( iter-1 );
    if( fabs((thisres-lastres)/lastres) < rtol ){
      return(2);
    }
  }
  /* No convergence - return 0 */
  return (0);
}

PetscScalar getGlobalStrainRateResidual( PetscInt iter ){
  return( residuals[iter % NRES] );
}

void setGlobalStrainRateResidual( PetscInt iter, PetscScalar res){
  residuals[ iter % NRES ] = res;
}
