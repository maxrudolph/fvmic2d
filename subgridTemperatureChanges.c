#include "fdcode.h"
#include "subgridTemperatureChanges.h"
#include "markers.h"

PetscErrorCode subgridTemperatureChanges(Vec thisT, GridData *grid, NodalFields *nodalFields, MarkerSet *markerset, Materials *materials, PetscScalar dt, Options *options){
  PetscErrorCode ierr;
  PetscScalar dsubgrid = options->subgridTemperatureDiffusivity;

  Marker *markers = markerset -> markers;

  /* subtract off last solution from Sa*/
  /*   PetscScalar *deltaT; */
  Vec deltaT;
  ierr = VecDuplicate(thisT, &deltaT);
  ierr = VecCopy( thisT, deltaT);
  ierr = VecAXPY( deltaT, -1.0, nodalFields->lastT);CHKERRQ(ierr);/* deltaT = thisT - lastT */

  /* get T, rho, Cp, k for ghost nodes */
  Vec Tl, rhol, Cpl, kThermall;
  ierr=DMCreateLocalVector(grid->da,&Tl); CHKERRQ(ierr);
  ierr=VecDuplicate(Tl, &rhol); CHKERRQ(ierr);
  ierr=VecDuplicate(Tl, &Cpl); CHKERRQ(ierr);
  ierr=VecDuplicate(Tl, &kThermall); CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->lastT,INSERT_VALUES,Tl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->lastT,INSERT_VALUES,Tl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->Cp,INSERT_VALUES,Cpl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->Cp,INSERT_VALUES,Cpl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->kThermal,INSERT_VALUES,kThermall);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->kThermal,INSERT_VALUES,kThermall);CHKERRQ(ierr);
  /* get arrays with ghost nodes*/
  PetscScalar **T,**rho,**Cp,**kThermal;
  ierr = DMDAVecGetArray(grid->da,Tl,&T);
  ierr = DMDAVecGetArray(grid->da,rhol,&rho);
  ierr = DMDAVecGetArray(grid->da,Cpl,&Cp);
  ierr = DMDAVecGetArray(grid->da,kThermall,&kThermal);

  /* allocate vector for subgrid temperature changes*/
  PetscScalar *deltaTmSubgrid;
  ierr = PetscMalloc( markerset->nMark*sizeof(PetscScalar), &deltaTmSubgrid);CHKERRQ(ierr);
  
  /* begin looping over markers*/
  PetscInt m;
  for(m=0;m<markerset->nMark;m++){
    if( markers[m].cellX != -1){    /* Check to see if marker is in bounds */
      PetscInt cellY = markers[m].cellY;
      PetscInt cellX = markers[m].cellX;
      /* restrict to INTERNAL cells */
      //if(cellX<1) cellX=1;
      //if(cellX>NX-3) cellX=NX-3;
      //if(cellY<1) cellY=1;
      //if(cellY>NY-3) cellY=NY-3;
      /* compute marker mdx, mdy (position in cell) */
      PetscScalar dx = grid->x[cellX+1]-grid->x[cellX];
      PetscScalar dy = grid->y[cellY+1]-grid->y[cellY];
      PetscScalar mdx=(markers[m].X - grid->x[cellX])/dx;
      PetscScalar mdy=(markers[m].Y - grid->y[cellY])/dy;

      PetscScalar Tmnodeslast = T[cellY][cellX]*(1.0-mdx)*(1.0-mdy) + T[cellY][cellX+1]*mdx*(1.0-mdy) + T[cellY+1][cellX]*(1.0-mdx)*(mdy) + T[cellY+1][cellX+1]*mdx*mdy;
      PetscScalar rhonodeslast = rho[cellY][cellX]*(1.0-mdx)*(1.0-mdy) + rho[cellY][cellX+1]*(mdx)*(1.0-mdy) + rho[cellY+1][cellX]*(1.0-mdx)*(mdy) + rho[cellY+1][cellX+1]*mdx*mdy;
      PetscScalar Cpnodeslast =  Cp[cellY][cellX]*(1.0-mdx)*(1.0-mdy) + Cp[cellY][cellX+1]*(mdx)*(1.0-mdy) + Cp[cellY+1][cellX]*(1.0-mdx)*(mdy) + Cp[cellY+1][cellX+1]*mdx*mdy;
      PetscScalar knodeslast= kThermal[cellY][cellX]*(1.0-mdx)*(1.0-mdy) + kThermal[cellY][cellX+1]*(mdx)*(1.0-mdy) + kThermal[cellY+1][cellX]*(1.0-mdx)*(mdy) + kThermal[cellY+1][cellX+1]*mdx*mdy;
      /* Gerya Book Equation 10.16 */
      //PetscScalar dTdiff = Cpnodeslast*rhonodeslast/knodeslast/(2.0/grid.dx2 + 2.0/grid.dy2);
      PetscScalar dTdiff = Cpnodeslast*rhonodeslast/knodeslast/(2.0/(dx*dx) + 2.0/(dy*dy));
      deltaTmSubgrid[m] = (Tmnodeslast - markers[m].T)*(1.0-exp( -dsubgrid*dt/dTdiff));
    }/* end if in domain*/ 
  }/* end loop over markers*/

  /* project deltaTmSubgrid onto the nodes */

  ierr=projectMarkersNodesFromScalar(markerset, grid, deltaTmSubgrid, deltaT);CHKERRQ(ierr);

  /* deltaT[i] = Sa[i]-nodalFields.lastT[i]- deltaT[i]; */
  /* deltaT = +1.0*thisT - 1.0*nodalFields.lastT -1.0*deltaT */
  ierr = VecAXPBYPCZ(deltaT,1.0,-1.0,-1.0,thisT,nodalFields->lastT);/*temperature correction*/
  /* get deltaT as an array */
  Vec deltaTl;
  ierr=DMCreateLocalVector(grid->da,&deltaTl); CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,deltaT,INSERT_VALUES,deltaTl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,deltaT,INSERT_VALUES,deltaTl);CHKERRQ(ierr);
  PetscScalar **dT;
  ierr =  DMDAVecGetArray(grid->da,deltaTl,&dT);

  /* now transfer thermal solution onto the markers */
  /* begin looping over markers*/
  for(m=0;m<markerset->nMark;m++){
    if( markers[m].cellX != -1){    /* Check to see if marker is in bounds */
      PetscInt cellY = markers[m].cellY;
      PetscInt cellX = markers[m].cellX;
      /* compute marker mdx, mdy (position in cell) */
      PetscScalar dx = grid->x[cellX+1]-grid->x[cellX];
      PetscScalar dy = grid->y[cellY+1]-grid->y[cellY];
      PetscScalar mdx=(markers[m].X - grid->x[cellX])/dx;
      PetscScalar mdy=(markers[m].Y - grid->y[cellY])/dy;
      PetscScalar dtm = (dT[cellY][cellX]*(1.0-mdy)*(1.0-mdx)+dT[cellY+1][cellX]*mdy*(1.0-mdx)+dT[cellY][cellX+1]*(1.0-mdy)*mdx+dT[cellY+1][cellX+1]*mdy*mdx) + deltaTmSubgrid[m];
      markers[m].T += dtm;
      markers[m].Tdot = dtm/dt;
      /* UPDATE OF DENSITY, INCLUDING THERMAL EXPANSION GOES HERE - one option. This is currently in main instead to account for damage*/
    }/* end if in bounds*/
  }/* end loop over markers*/

  ierr = DMDAVecRestoreArray(grid->da,deltaTl,&dT);CHKERRQ(ierr);
  ierr = VecDestroy(&deltaTl);CHKERRQ(ierr);

  /* clean up!!! */
  ierr = PetscFree(deltaTmSubgrid); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,Tl,&T);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,rhol,&rho);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,Cpl,&Cp);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,kThermall,&kThermal);CHKERRQ(ierr);

  ierr = VecDestroy(&deltaT);CHKERRQ(ierr);

  ierr = VecDestroy(&Tl);CHKERRQ(ierr);
  ierr = VecDestroy(&rhol);CHKERRQ(ierr);
  ierr = VecDestroy(&Cpl);CHKERRQ(ierr);
  ierr = VecDestroy(&kThermall);CHKERRQ(ierr);


  PetscFunctionReturn(ierr);
}/* end function*/
