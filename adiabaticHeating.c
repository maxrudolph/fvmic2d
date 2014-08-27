#include "fdcode.h"
#include "adiabaticHeating.h"

/* this is a subroutine to compute adiabatic heating */
/* it actually just computes vy*g*T*alpha*rho on the markers and projects this to the nodes*/

PetscErrorCode adiabaticHeating(GridData *grid, MarkerSet *markerset, NodalFields *nodalFields, Materials *materials, Options *options){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  Vec wtng,wtnl,hal;
  /* allocate global vector for weights */
  ierr=DMCreateGlobalVector(grid->da, &wtng);CHKERRQ(ierr);
  ierr=DMCreateLocalVector(grid->da, &wtnl);CHKERRQ(ierr);
  /* create local vector for adiabatic heating */
  ierr=DMCreateLocalVector(grid->da, &hal);CHKERRQ(ierr);
  ierr=VecZeroEntries(hal);CHKERRQ(ierr);
  ierr=VecZeroEntries(wtnl);
  ierr=VecZeroEntries(nodalFields->ha);CHKERRQ(ierr);

  PetscScalar **wtn,**ha;
  ierr=DMDAVecGetArray(grid->da,hal,&ha);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);

  /* loop over markers */
  Marker *markers = markerset->markers;
  PetscInt m;
  for(m=0;m<markerset->nMark;m++){
    if(markers[m].cellX != -1){
      PetscInt cellI = markers[m].cellX;
      PetscInt cellJ = markers[m].cellY;
      PetscScalar mdx = (markers[m].X - grid->x[cellI])/(grid->x[cellI+1]-grid->x[cellI]);
      PetscScalar mdy = (markers[m].Y - grid->y[cellJ])/(grid->y[cellJ+1]-grid->y[cellJ]);
      /* calculate marker adiabatic heating */
      PetscScalar mrho = materials->materialRho[(PetscInt) markers[m].Mat];
      PetscScalar mha = materials->materialAlpha[(PetscInt) markers[m].Mat] * markers[m].T * markers[m].VY * mrho * options->gy;
      /* contribution of heating and weights to surrounding nodes*/
      ha[cellJ][cellI]      += mha*(1-mdx)*(1-mdy);
      ha[cellJ][cellI+1]    += mha*(mdx)*(1-mdy);
      ha[cellJ+1][cellI]    += mha*(1-mdx)*(mdy);
      ha[cellJ+1][cellI+1]  += mha*(mdx)*(mdy);
      wtn[cellJ][cellI]     +=     (1-mdx)*(1-mdy);
      wtn[cellJ][cellI+1]   +=     (mdx)*(1-mdy);
      wtn[cellJ+1][cellI]   +=     (1-mdx)*(mdy);
      wtn[cellJ+1][cellI+1] +=     (mdx)*(mdy); 
    }
  }
  /* restore arrays */
  ierr = DMDAVecRestoreArray(grid->da,hal,&ha);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);

  /* scatter to global vectors */
  ierr = DMLocalToGlobalBegin(grid->da,wtnl,ADD_VALUES,wtng);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,wtnl,ADD_VALUES,wtng);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,hal,ADD_VALUES,nodalFields->ha);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,wtnl,ADD_VALUES,nodalFields->ha);CHKERRQ(ierr);

  /* normalize nodal ha */
  ierr = VecPointwiseDivide(nodalFields->ha,nodalFields->ha,wtng);CHKERRQ(ierr);

  ierr = VecDestroy(&wtng);CHKERRQ(ierr);
  ierr = VecDestroy(&wtnl);CHKERRQ(ierr);
  ierr = VecDestroy(&hal);CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}
