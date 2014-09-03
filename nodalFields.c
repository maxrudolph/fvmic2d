#include "petscdmda.h"
#include "fdcode.h"
#include "nodalFields.h"

PetscErrorCode initializeNodalFields( NodalFields *nodalFields, GridData *grid, Options *options){

  PetscErrorCode ierr;
  /* initialize nodal quantity arrays*/
  PetscFunctionBegin;

  DMBoundaryType xwrap;
  if( options->mechBCLeft.type[0] == 2 || options->mechBCRight.type[0] == 2 ){
    /* periodicity exists in x-direction */
    xwrap = DM_BOUNDARY_PERIODIC;
    printf("Using periodic X boundary conditions\n");
  } else {
    xwrap = DM_BOUNDARY_GHOSTED;
  }
  ierr = DMDACreate2d( PETSC_COMM_WORLD, xwrap, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, grid->NX, grid->NY, PETSC_DECIDE,PETSC_DECIDE, 1, 2, PETSC_NULL, PETSC_NULL, &grid->da ); CHKERRQ(ierr);
  ierr = DMDACreate2d( PETSC_COMM_WORLD, xwrap, DM_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, grid->NX, grid->NY, PETSC_DECIDE,PETSC_DECIDE, 3, 2, PETSC_NULL, PETSC_NULL, &grid->vda ); CHKERRQ(ierr);

  /* put coordinate information into DA - from Petsc ex4.c */
  {
    Vec coords, global;
    DM coordsda;
    DMDACoor2d **X;/* coordinates */
    ierr = DMDASetUniformCoordinates(grid->da,0.0,1.0,0.0,1.0,0.0,0.0);
    ierr = DMGetCoordinateDM( grid->da, &coordsda ); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(grid->da, &coords );CHKERRQ(ierr);
    ierr = DMDAVecGetArray(coordsda, coords, &X); CHKERRQ(ierr);
    PetscInt ix,jy;
    PetscInt x,y,m,n;
    ierr=DMDAGetCorners(coordsda,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
    for(jy=y;jy<y+n;jy++){
      for(ix=x;ix<x+m;ix++){
	X[jy][ix].x = grid->x[ix];
	X[jy][ix].y = grid->y[jy];
	//printf("X[%d][%d].x=%e\n",jy,ix,X[jy][ix].x);
      }
    }    
    ierr=DMDAVecRestoreArray(coordsda, coords, &X);CHKERRQ(ierr);
    ierr=DMGetCoordinates( grid->da, &global);CHKERRQ(ierr);
    ierr=DMLocalToGlobalBegin(coordsda,coords,INSERT_VALUES,global);
    ierr=DMLocalToGlobalEnd(coordsda,coords,INSERT_VALUES,global);
  }
  ierr=DMDASetFieldName(grid->da,0,"");
  ierr=PetscObjectSetName( (PetscObject) grid->da,"Scalar DMDA");
  
  /* if texture is enabled, allocate nodal fields for viscoplasticity tensor */
#ifdef TEXTURE
  ierr = DMDACreate2d( PETSC_COMM_WORLD, xwrap, DMDA_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, grid->NX, grid->NY, PETSC_DECIDE,PETSC_DECIDE, 4, 4, PETSC_NULL, PETSC_NULL, &grid->tda ); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(grid->tda, &nodalFields->VPTensorC);CHKERRQ(ierr);/*cell-centered viscoplasticity tensor*/
  ierr = DMCreateGlobalVector(grid->tda, &nodalFields->VPTensorB);CHKERRQ(ierr);/*same but at basic nodes instead of cell centers*/
  ierr = PetscObjectSetName((PetscObject) nodalFields->VPTensorC, "Nc");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->VPTensorB, "Nb");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(grid->da, &nodalFields->strainRateResidual);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->strainRateResidual,"strainRateResidual");CHKERRQ(ierr);
#endif

  ierr = DMCreateGlobalVector(grid->da, &nodalFields->lastT); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->lastT, "lastT");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->thisT);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->thisT, "thisT");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->kThermal);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->kThermal, "kThermal");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->Cp);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->Cp, "Cp");CHKERRQ(ierr);

  /*initial density distribution*/

/*   ierr = DACreate2d( PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_BOX, grid->NY, grid->NX, PETSC_DECIDE,PETSC_DECIDE, 1, 1, PETSC_NULL, PETSC_NULL, &nodalFields->rhodot ); CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->rhodot);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->rhodot, "rhodot");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->rho);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->rho, "rho");CHKERRQ(ierr);
  /* arrays to hold nodal and mid-cell viscosity */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->etaS);CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->etaN);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->etaS,"etaS");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->etaN,"etaN");CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->etavx);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->etavx,"etavx");CHKERRQ(ierr); */
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->etavy);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->etavy,"etavy");CHKERRQ(ierr); */


  /* arrays to hold nodal strain rates*/

  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->edotxx);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->edotxx, "edotxx");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->edotyy);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->edotyy, "edotyy");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->edotzz);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->edotzz, "rho");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->edotxy);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->edotxy, "edotxy");CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->edotxz);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->edotxz, "edotxz");CHKERRQ(ierr); */
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->edotyz);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->edotyz, "edotyz");CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->eii);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->eii, "eii");CHKERRQ(ierr);
  /* same for nodal stresses*/
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->dsxx);CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->dsyy);CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->dszz);CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->dsxy);CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->dsxz);CHKERRQ(ierr); */
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->dsyz);CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->sii);CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->ha);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->ha,"ha");CHKERRQ(ierr);
  /* same for nodal muPETSC_DECIDE, muS*/
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->muN);CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->muS);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->muN, "muN");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->muS, "muS");CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->muvx);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->muvx,"muvx");CHKERRQ(ierr); */
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->muvy);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->muvy,"muvy");CHKERRQ(ierr); */

  /* same for last deviatoric stresses*/
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->soxx);CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->soyy);CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->sozz);CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->soxy);CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->soxz);CHKERRQ(ierr); */
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->soyz);CHKERRQ(ierr); */
  ierr = PetscObjectSetName((PetscObject) nodalFields->soxx, "soxx");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->soyy, "soyy");CHKERRQ(ierr);
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->sozz, "sozz");CHKERRQ(ierr); */
  ierr = PetscObjectSetName((PetscObject) nodalFields->soxy, "soxy");CHKERRQ(ierr);
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->soxz, "soxz");CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->soyz, "soyz");CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->vx);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->vx, "vx");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->vy);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->vy, "vy");CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->vz);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->vz, "vz");CHKERRQ(ierr); */
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->p);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->p, "p");CHKERRQ(ierr);
  ierr = VecDuplicate( nodalFields->lastT, &nodalFields->wxy);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) nodalFields->wxy, "wxy");CHKERRQ(ierr);
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->wxz);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->wxz, "wxz");CHKERRQ(ierr); */
/*   ierr = VecDuplicate( nodalFields->lastT, &nodalFields->wyz);CHKERRQ(ierr); */
/*   ierr = PetscObjectSetName((PetscObject) nodalFields->wyz, "wyz");CHKERRQ(ierr); */

  PetscFunctionReturn(ierr);
}


PetscErrorCode resetNodalFields( NodalFields *nodalFields, GridData *grid, Options *options){
  PetscErrorCode ierr;
  ierr = VecZeroEntries(nodalFields->thisT);CHKERRQ(ierr);
  ierr = VecZeroEntries(nodalFields->lastT);CHKERRQ(ierr);

#ifdef TEXTURE
    ierr=VecZeroEntries(nodalFields->VPTensorB);CHKERRQ(ierr);
    ierr=VecZeroEntries(nodalFields->VPTensorC);CHKERRQ(ierr);
    ierr=VecZeroEntries(nodalFields->strainRateResidual);CHKERRQ(ierr);
#endif

/*     nodalFields->thisT[i] = 0; */
/*     nodalFields->lastT[i] = 0; */
/*     nodalFields->kThermal[i] = 0; */
  ierr = VecZeroEntries(nodalFields->kThermal);CHKERRQ(ierr);
/*     nodalFields->Cp[i] = 0; */
  ierr = VecZeroEntries(nodalFields->Cp);CHKERRQ(ierr);
/*     nodalFields->rho[i] = 0; */
  ierr = VecZeroEntries(nodalFields->rho);CHKERRQ(ierr);
/*     nodalFields->rhodot[i] = 0; */
  ierr = VecZeroEntries(nodalFields->rhodot);CHKERRQ(ierr);
/*     nodalFields->etaS[i] = 0; */
  ierr = VecZeroEntries(nodalFields->etaS);CHKERRQ(ierr);
/*     nodalFields->etaN[i] = 0; */
  ierr = VecZeroEntries(nodalFields->etaN);CHKERRQ(ierr);
/*   ierr = VecZeroEntries(nodalFields->etavx);CHKERRQ(ierr); */
/*   ierr = VecZeroEntries(nodalFields->etavy);CHKERRQ(ierr); */
/*     nodalFields->edotxx[i] = 0; */
  ierr = VecZeroEntries(nodalFields->edotxx);CHKERRQ(ierr);
/*     nodalFields->edotyy[i] = 0; */
  ierr = VecZeroEntries(nodalFields->edotyy);CHKERRQ(ierr);
/*     nodalFields->edotzz[i] = 0; */
  ierr = VecZeroEntries(nodalFields->edotzz);CHKERRQ(ierr);
/*     nodalFields->edotxy[i] = 0; */
  ierr = VecZeroEntries(nodalFields->edotxy);CHKERRQ(ierr);
/*     nodalFields->edotxz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->edotxz);CHKERRQ(ierr); */
/*     nodalFields->edotyz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->edotyz);CHKERRQ(ierr); */
/*     nodalFields->dsxx[i] = 0; */
  ierr = VecZeroEntries(nodalFields->dsxx);CHKERRQ(ierr);
/*     nodalFields->dsyy[i] = 0; */
  ierr = VecZeroEntries(nodalFields->dsyy);CHKERRQ(ierr);
/*     nodalFields->dszz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->dszz);CHKERRQ(ierr); */
/*     nodalFields->dsxy[i] = 0; */
  ierr = VecZeroEntries(nodalFields->dsxy);CHKERRQ(ierr);
/*     nodalFields->dsxz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->dsxz);CHKERRQ(ierr); */
/*     nodalFields->dsyz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->dsyz);CHKERRQ(ierr); */
/*     nodalFields->eii[i] = 0; */
  ierr = VecZeroEntries(nodalFields->eii);CHKERRQ(ierr);
/*     nodalFields->sii[i] = 0; */
  ierr = VecZeroEntries(nodalFields->sii);CHKERRQ(ierr);
/*     nodalFields->ha[i] = 0; */
  ierr = VecZeroEntries(nodalFields->ha);CHKERRQ(ierr);
  //ierr = VecZeroEntries(nodalFields->nodalHeating);CHKERRQ(ierr);
/*     nodalFields->muN[i] = 0; */
  ierr = VecZeroEntries(nodalFields->muN);CHKERRQ(ierr);
/*   ierr = VecZeroEntries(nodalFields->muvx);CHKERRQ(ierr); */
/*   ierr = VecZeroEntries(nodalFields->muvy);CHKERRQ(ierr); */
/*     nodalFields->muS[i] = 0; */
  ierr = VecZeroEntries(nodalFields->muS);CHKERRQ(ierr);
/*     nodalFields->soxx[i] = 0; */
  ierr = VecZeroEntries(nodalFields->soxx);CHKERRQ(ierr);
/*     nodalFields->soyy[i] = 0; */
  ierr = VecZeroEntries(nodalFields->soyy);CHKERRQ(ierr);
/*     nodalFields->sozz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->sozz);CHKERRQ(ierr); */
/*     nodalFields->soxy[i] = 0; */
  ierr = VecZeroEntries(nodalFields->soxy);CHKERRQ(ierr);
/*     nodalFields->soxz[i] = 0; */
/*   ierr = VecZeroEntries(nodalFields->soxz);CHKERRQ(ierr); */
/*   ierr = VecZeroEntries(nodalFields->soyz);CHKERRQ(ierr); */
  ierr = VecZeroEntries(nodalFields->vx);CHKERRQ(ierr);
  ierr = VecZeroEntries(nodalFields->vy);CHKERRQ(ierr);
/*   ierr = VecZeroEntries(nodalFields->vz);CHKERRQ(ierr); */
  ierr = VecZeroEntries(nodalFields->p);CHKERRQ(ierr);
  ierr = VecZeroEntries(nodalFields->wxy);CHKERRQ(ierr);
/*   ierr = VecZeroEntries(nodalFields->wxz);CHKERRQ(ierr); */
/*   ierr = VecZeroEntries(nodalFields->wyz);CHKERRQ(ierr); */

  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyNodalFields( NodalFields *nodalFields, GridData *grid){
  PetscErrorCode ierr;
  ierr = VecDestroy(&nodalFields->thisT);CHKERRQ(ierr);
  ierr = VecDestroy(&nodalFields->lastT);CHKERRQ(ierr);
  
/*     nodalFields->thisT[i] = 0; */
/*     nodalFields->lastT[i] = 0; */
/*     nodalFields->kThermal[i] = 0; */
  ierr = VecDestroy(&nodalFields->kThermal);CHKERRQ(ierr);
/*     nodalFields->Cp[i] = 0; */
  ierr = VecDestroy(&nodalFields->Cp);CHKERRQ(ierr);
/*     nodalFields->rho[i] = 0; */
  ierr = VecDestroy(&nodalFields->rho);CHKERRQ(ierr);
/*     nodalFields->rhodot[i] = 0; */
  ierr = VecDestroy(&nodalFields->rhodot);CHKERRQ(ierr);
/*     nodalFields->etaS[i] = 0; */
  ierr = VecDestroy(&nodalFields->etaS);CHKERRQ(ierr);
/*     nodalFields->etaN[i] = 0; */
  ierr = VecDestroy(&nodalFields->etaN);CHKERRQ(ierr);
/*   ierr = VecDestroy(&nodalFields->etavx);CHKERRQ(ierr); */
/*   ierr = VecDestroy(&nodalFields->etavy);CHKERRQ(ierr); */

/*     nodalFields->edotxx[i] = 0; */
  ierr = VecDestroy(&nodalFields->edotxx);CHKERRQ(ierr);
/*     nodalFields->edotyy[i] = 0; */
  ierr = VecDestroy(&nodalFields->edotyy);CHKERRQ(ierr);
/*     nodalFields->edotzz[i] = 0; */
  ierr = VecDestroy(&nodalFields->edotzz);CHKERRQ(ierr);
/*     nodalFields->edotxy[i] = 0; */
  ierr = VecDestroy(&nodalFields->edotxy);CHKERRQ(ierr);
/*     nodalFields->edotxz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->edotxz);CHKERRQ(ierr); */
/*     nodalFields->edotyz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->edotyz);CHKERRQ(ierr); */
/*     nodalFields->dsxx[i] = 0; */
  ierr = VecDestroy(&nodalFields->dsxx);CHKERRQ(ierr);
/*     nodalFields->dsyy[i] = 0; */
  ierr = VecDestroy(&nodalFields->dsyy);CHKERRQ(ierr);
/*     nodalFields->dszz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->dszz);CHKERRQ(ierr); */
/*     nodalFields->dsxy[i] = 0; */
  ierr = VecDestroy(&nodalFields->dsxy);CHKERRQ(ierr);
/*     nodalFields->dsxz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->dsxz);CHKERRQ(ierr); */
/*     nodalFields->dsyz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->dsyz);CHKERRQ(ierr); */
/*     nodalFields->eii[i] = 0; */
  ierr = VecDestroy(&nodalFields->eii);CHKERRQ(ierr);
/*     nodalFields->sii[i] = 0; */
  ierr = VecDestroy(&nodalFields->sii);CHKERRQ(ierr);
/*     nodalFields->ha[i] = 0; */
  ierr = VecDestroy(&nodalFields->ha);CHKERRQ(ierr);
/*     nodalFields->muN[i] = 0; */
  ierr = VecDestroy(&nodalFields->muN);CHKERRQ(ierr);
/*     nodalFields->muS[i] = 0; */
  ierr = VecDestroy(&nodalFields->muS);CHKERRQ(ierr);
/*   ierr = VecDestroy(&nodalFields->muvx);CHKERRQ(ierr); */
/*   ierr = VecDestroy(&nodalFields->muvy);CHKERRQ(ierr); */

/*     nodalFields->soxx[i] = 0; */
  ierr = VecDestroy(&nodalFields->soxx);CHKERRQ(ierr);
/*     nodalFields->soyy[i] = 0; */
  ierr = VecDestroy(&nodalFields->soyy);CHKERRQ(ierr);
/*     nodalFields->sozz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->sozz);CHKERRQ(ierr); */
/*     nodalFields->soxy[i] = 0; */
  ierr = VecDestroy(&nodalFields->soxy);CHKERRQ(ierr);
/*     nodalFields->soxz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->soxz);CHKERRQ(ierr); */
/*     nodalFields->soyz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->soyz);CHKERRQ(ierr); */
/*     nodalFields->vx[i] = 0; */
  ierr = VecDestroy(&nodalFields->vx);CHKERRQ(ierr);
/*     nodalFields->vy[i] = 0; */
  ierr = VecDestroy(&nodalFields->vy);CHKERRQ(ierr);
/*     nodalFields->vz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->vz);CHKERRQ(ierr); */
/*     nodalFields->p[i] = 0; */
  ierr = VecDestroy(&nodalFields->p);CHKERRQ(ierr);
/*     nodalFields->wxy[i] = 0; */
  ierr = VecDestroy(&nodalFields->wxy);CHKERRQ(ierr);
/*     nodalFields->wxz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->wxz);CHKERRQ(ierr); */
/*     nodalFields->wyz[i] = 0; */
/*   ierr = VecDestroy(&nodalFields->wyz);CHKERRQ(ierr); */

#ifdef TEXTURE
  ierr = VecDestroy(&nodalFields->VPTensorB);CHKERRQ(ierr);
  ierr = VecDestroy(&nodalFields->VPTensorC);CHKERRQ(ierr);
  ierr = VecDestroy(&nodalFields->strainRateResidual);CHKERRQ(ierr);
  ierr = DMDestroy(&grid->tda);CHKERRQ(ierr);
#endif

  ierr=DMDestroy(&grid->da);CHKERRQ(ierr);
  ierr=DMDestroy(&grid->vda);CHKERRQ(ierr);



  PetscFunctionReturn(ierr);
}


