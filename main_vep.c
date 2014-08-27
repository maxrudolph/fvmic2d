/* exercise 10.1 from Gerya's boook */
/* this version checks against gerya's solution to exercise 10.2 */
//#include "../include/fem.h"
//
/*
 *Main routine of FEM code. This copy is specifically for ridge damage problem
 *Do not put any compiler macros here. They should be in fem.h.
 */
static char help[] = "Help statement goes here\n\
 petscmpiexec -n 2 ./Fem ./input/fault.1 output/test -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps >log1\n\
It can be one line or many \n";

#include<stdio.h>
#include<stdlib.h>
//#include<malloc.h>
#include "petscksp.h"
#include "mpi.h"
#include "fdcode.h"
#include "kspLinearSolve.h"
#include "markers.h"
#include "reggrid.h"
#include "nodalFields.h"
#include "markers.h"
#include "materialProperties.h"
#include "thermalSystem.h"
//#include "stokes_varvisc.h"
#include "vep_system.h"
#include "nodalStressStrain.h"/* includes shear heating too*/
//#include "shearHeating.h"
#include "updateMarkerStrainPressure.h"
#include "subgridStresschanges.h"
#include "updateMarkerStress.h"
#include "io.h"

//Include file for PETSc linear solve commands

//extern struct ELEMENT linearHexahedralElementCreate();
//extern double evaluateNodalBasisFunction();

//extern struct ELEMENT quadTriElementCreate();
//extern double evaluateNodalBasisFunction();

#define NADD 12 
//maximum number of entries to add at once.

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **args){
  PetscMPIInt rank,size;
  PetscErrorCode ierr;
  KSP	ksp; /* linear solver context */
  PetscViewer viewer;
  PetscLogStage stages[6];  
  /*end debugging stuff*/
  
  PetscInt NX,NY;/* number of nodes in x,y,directions*/
  //  PetscScalar eta;

  /*   PetscScalar *etaS,*etaN; */ /* center and nodal viscosity */
  /* thermal problem */
  /*   PetscScalar *thisT, *lastT;*//*current and former timestep temperature*/
  /*   PetscScalar *Cp,*kThermal; */
  
  GridData grid;
  NodalFields nodalFields;
  
  //PetscScalar *rho; /*nodal density dist.*/

  PetscScalar LX,LY;/* model dimensions (x,y)*/
  /*   PetscScalar dx,dy; */ /*x,y grid spacing*/
  /*   PetscScalar dx2,dy2,dxdy;  *//*dx*dx, dy*dy*/
  PetscScalar *x,*y,*xc,*yc;/* grid line x and y values, cell center x and y values.*/
  Mat LHS, thermalLHS;
  Vec RHS, thermalRHS, nodalHeating;
  
  /* markers*/
  /*   PetscScalar *markerX, *markerY; */  /* marker nodal positions*/
  /*   PetscScalar *markerVX, *markerVY; */ /* marker velocities*/
  /*   PetscScalar *markerT; */
  Markers markers; /* data structure that holds all of the marker info*/
  Materials materials;
  /*   char *markerMat; *//*marker material*/
  
  //INITIALIZE PETSC
  PetscInitialize(&argc,&args,(char *)0,help); //copied from petsc/src/ksp/ksp/examles/tutorials
  //Set MPI rank
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  /* grid parameters */
  NX=31;
  NY=41;
  PetscInt NMX=10;/* number of markers per cell, x*/
  PetscInt NMY=10;/* markers per cell y*/
  PetscInt ndof = NX*NY*3;
  PetscInt nTime = 100;
  PetscInt maxNumPlasticity=100;

  /* Model Parameters */
  LX=1000000.0;
  LY=1000000.0;
  //  eta=1.0e21;/* viscosity, pa-s*/
  /*   PetscScalar rholeft = 3200; */
  /*   PetscScalar rhoright = 3300; */
  /*   PetscScalar etaleft=1.0e20; */
  /*   PetscScalar etaright=1.0e22; */
  PetscScalar boundaryT = 1000;/* boundary value temperature*/
  nodalFields.boundaryT = boundaryT;


  PetscScalar bumpT = 1300;/* temperature of thermal "bump"*/
  PetscScalar gy = 10;
  /*   PetscScalar materialEta[2] = {1.0e20, 1.0e22}; */
  /*   PetscScalar materialRho[2] = {3200.0, 3300.0}; */
  /*   PetscScalar materialkThermal[2] = {3.0,10.0}; */
  /*   PetscScalar materialCp[2] = {1000.0,1100.0}; */

  PetscScalar etamin = 1.0e20;
  PetscScalar etamax = 1.0e23;
  /* grid setup*/
  /*   PetscInt nMark = NX*NMX*NY*NMY; */
  PetscScalar dtMax = 3.16e13;/*maximum timestep in seconds*/
  dtMax = 1e+2*(365.25*24*3600);/* taras*/
  PetscScalar dt = 7.466666666666666e+13;
  PetscScalar displacementdt;/* this is the displacement timestep*/
  PetscScalar elapsedTime=0.0;

  PetscLogStageRegister("Main",&stages[0]);
  PetscLogStageRegister("Mechanical Assembly",&stages[1]);
  PetscLogStageRegister("Check Plastic Yield",&stages[2]);
  PetscLogStageRegister("Mechanical Solve",&stages[3]);
  PetscLogStagePush(stages[0]);

  /* initialize material properties*/
  setMaterialProperties( &materials);
  printf("made it here\n"); 
  /* initialize grid*/
  ierr = initializeRegularGrid( &grid, LX, LY, NX, NY);CHKERRQ(ierr);

  /* scaling parameters*/
  PetscScalar Kbond = 4*etamin/((grid.dx+grid.dy)*(grid.dx+grid.dy));
  PetscScalar Kcont = 2*etamin/(grid.dx+grid.dy);

  /* allocate the nodal fields*/
  ierr = initializeNodalFields( &nodalFields, &grid);CHKERRQ(ierr);

  /* initialize Markers*/
  ierr=allocateMarkers( NX*NMX*NY*NMY, &markers);CHKERRQ(ierr);
  /* distribute markers*/
  ierr=distributeMarkers( &markers, NX, NY, NMX, NMY, LX, LY);CHKERRQ(ierr);
  /* make sure that markers are inside domain*/
  ierr = checkMarkersOutsideDomain( &markers, LX, LY); CHKERRQ(ierr);

  PetscInt m;
  for(m=0;m<markers.nMark;m++){
    if( fabs(markers.Y[m] - LY/2.0) < 0.2*LY && fabs(markers.X[m] - LX/2.0) < 0.2*LX){/* inside perturbed region*/
      markers.Mat[m] = 1;
      markers.T[m] = bumpT;

    } else{
      markers.Mat[m] = 0;
      markers.T[m] = boundaryT;
    }
    markers.rho[m] = materials.materialRho[(PetscInt) markers.Mat[m]];
    markers.eta[m] = materials.materialEta[(PetscInt) markers.Mat[m]];
  }    

  //  printMarkersAllASCIIMatlab( &markers, -1);
  ierr=saveMarkersBinary( &markers, -1);
  PetscInt iTime;
  /* initialize matrix for left hand side and vector for right hand side of THERMAL PROBLEM*/
  PetscInt nnode = NX*NY;
  ierr = MatCreate( PETSC_COMM_WORLD, &thermalLHS);CHKERRQ(ierr);
  ierr = MatSetSizes( thermalLHS, PETSC_DECIDE,PETSC_DECIDE,nnode,nnode); CHKERRQ(ierr);
  ierr = MatSetType( thermalLHS, MATAIJ); CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) thermalLHS,"CTLHS");
  ierr = VecCreate( PETSC_COMM_WORLD, &thermalRHS); CHKERRQ(ierr);
  ierr = VecSetSizes( thermalRHS, PETSC_DECIDE, nnode);CHKERRQ(ierr);
  ierr = VecSetType( thermalRHS, VECMPI);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) thermalRHS,"CTRHS");
  ierr = VecDuplicate(thermalRHS, &nodalHeating);CHKERRQ(ierr);

  
  /* initialize matrix for left hand side and vector for right hand side*/
  ierr = MatCreate( PETSC_COMM_WORLD, &LHS);CHKERRQ(ierr);
  ierr = MatSetSizes( LHS, PETSC_DECIDE,PETSC_DECIDE,ndof,ndof); CHKERRQ(ierr);
  ierr = MatSetType( LHS, MATAIJ); CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) LHS,"CmechLHS");
  ierr = VecCreate( PETSC_COMM_WORLD, &RHS); CHKERRQ(ierr);
  ierr = VecSetSizes( RHS, PETSC_DECIDE, ndof);CHKERRQ(ierr);
  ierr = VecSetType( RHS, VECMPI);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) RHS,"CmechRHS");

  PetscInt isYielding;
  PetscInt iPlastic;
  PetscInt i,j;
  displacementdt=0.0;/* initially*/
  for(iTime=0;iTime<nTime;iTime++){
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Beginning Timestep %d\n",iTime);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);    
    isYielding = 1;
    iPlastic=0;
    dt=dtMax;/* initial timestep*/
    while(isYielding && iPlastic < maxNumPlasticity){/*begin plastic iteration loop*/
      /* check for plastic yielding */
      isYielding=0;


      ierr = checkPlasticYielding(&grid, &markers, &materials, displacementdt, &isYielding,etamin,etamax);
      /* uncomment these to monitor plastic iterations*/
      /*      saveNodalFieldsASCIIMatlab( &nodalFields, &grid, 1000+100*(iTime)+iPlastic); */
      /*       printMarkersAllASCIIMatlab( &markers, 1000+100*iTime+iPlastic); */

      if(isYielding) printf("iTime=%d,iPlastic = %d Yielding detected\n",iTime,iPlastic);
      if(!isYielding) printf("no more yielding\n");
      /* project marker quantities onto nodes: kThermal, Cp, rho, T*/
      ierr = projectMarkersNodesAll(&markers, &grid, &nodalFields, &materials);

      /* form mechanical problem LHS*/
      /* zero out LHS, RHS*/
      ierr = VecZeroEntries( RHS);CHKERRQ(ierr);
      //    ierr = formStokesVarViscSystem( &nodalFields, &grid, LHS, RHS, Kbond, Kcont, gy);CHKERRQ(ierr);
      ierr=formVEPSystem( &nodalFields, &grid, LHS, RHS, &Kbond, &Kcont, gy, dt);
      
      Vec mechanicalS;/* solution to mechanical problem*/
      ierr = VecDuplicate(RHS, &mechanicalS);
      ierr = kspLinearSolve(LHS, RHS, mechanicalS);

      /* dump solution*/
      /* ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./loadMechLHS.m",&viewer);CHKERRQ(ierr);   */
/*       ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr); */
/*       ierr = MatView( LHS  ,viewer);CHKERRQ(ierr);   */
/*       ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr); */
/*       ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./loadMechRHS.m",&viewer);CHKERRQ(ierr);   */
/*       ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr); */
/*       ierr = VecView( RHS  ,viewer);CHKERRQ(ierr);   */
/*       ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr); */
      
      /* retrieve velocity and pressure from mechanicalS */
      PetscScalar *Sa;
      ierr = VecGetArray(mechanicalS,&Sa); CHKERRQ(ierr);  
      /* check for nans in solution*/
      {PetscInt fn=0;
	for(i=0;i<3*NX*NY;i++){
	  if(isnan(Sa[i])) fn=1;
	}
	if(fn) {
	  printf("found nan in solution\n");
	  goto abort;
	}
      }
      for(i=0;i<NY;i++){
	for(j=0;j<NX;j++){
	  PetscInt idxnode = j*NY+i;
	  PetscInt pdof = 3*idxnode;
	  PetscInt vxdof = 3*idxnode+1;
	  PetscInt vydof = 3*idxnode+2;
	  nodalFields.p[idxnode] = Sa[pdof]*Kcont;
	  nodalFields.vx[idxnode]=Sa[vxdof];
	  nodalFields.vy[idxnode]=Sa[vydof];
	}
      }
      ierr = VecRestoreArray(mechanicalS,&Sa);CHKERRQ(ierr); 
      ierr = VecDestroy(mechanicalS);CHKERRQ(ierr);

      displacementdt = dt;
      /* determine displacement timestep*/
      {
	PetscScalar vxmax=0.0;
	PetscScalar vymax=0.0;
	for(i=0;i<grid.NX*grid.NY;i++){/* compute maximum vx, vy in absolute sense*/
	  if( fabs(nodalFields.vx[i]) > vxmax) vxmax=fabs(nodalFields.vx[i]);
	  if( fabs(nodalFields.vy[i]) > vymax) vymax=fabs(nodalFields.vy[i]);
	}
	//dt = grid.dx/2.0/vxmax;
	if( displacementdt > grid.dy/2.0/vymax) displacementdt=grid.dy/2.0/vymax;
	if( displacementdt > grid.dx/2.0/vxmax) displacementdt=grid.dx/2.0/vxmax;
	if( displacementdt > dtMax) displacementdt=dtMax;
	printf("Using Displacement Timestep %e years\n",displacementdt/3.1557e7);
      }      
      /* compute normal quantities: exx, sxx */
      /* compute shear quantities: exy, sxy, vorticity, adiabatic heating*/
      ierr = VecZeroEntries(nodalHeating);CHKERRQ(ierr);
      ierr=nodalStressStrain(&grid, &nodalFields, displacementdt,nodalHeating, gy);CHKERRQ(ierr);/* no nans here in first timestep*/

      /* compute marker strain and pressure*/
      ierr=updateMarkerStrainPressure( &grid,&nodalFields, &markers, &materials, displacementdt);

      /* update marker effective viscosity to match marker strain reate*/

      iPlastic++;
    }/* end plasticity loop*/
    
    printf("Done with plasticity loops.\n");

    /* do sub-grid stress diffusion*/
    PetscScalar dsubgridstress = 1.0;
    if(dsubgridstress>0.0) ierr=subgridStressChanges(&grid, &nodalFields, &markers, &materials, displacementdt, dsubgridstress);

    /* update marker stress*/
    ierr= updateMarkerStress( &grid, &nodalFields,&markers, &materials);

    /* begin solve thermal problem*/

    /* determine correct time step for thermal problem */

    /* set thermal problem timestep to displacement dt*/
    dt = displacementdt;

    /* form the thermal system*/
    ierr = formThermalSystem( &grid, &nodalFields, thermalLHS, thermalRHS, dt);
    /* add shear heating to RHS*/
    //    ierr = VecAXPY( thermalRHS, 1.0, nodalHeating);
    /* solve linear system */
    //  ierr 
    Vec thermalS;
    ierr = VecDuplicate( thermalRHS, &thermalS);
    ierr = PetscObjectSetName( (PetscObject) thermalS,"CthermalS");
    ierr = kspLinearSolve( thermalLHS, thermalRHS, thermalS);CHKERRQ(ierr);
    /* thermal solution is dT/dt*/

    /*get solution as an array */
    PetscScalar *Sa;
    ierr = VecGetArray(thermalS,&Sa); CHKERRQ(ierr);
    /* copy thermal solution to lastT */
/*     for(i=0;i<nnode;i++){ */
/*       lastT[i] = Sa[i];  */
/*     } */
    /* subtract off last solution from Sa*/
    PetscScalar *deltaT;
    ierr = PetscMalloc( nnode*sizeof(PetscScalar), &deltaT);CHKERRQ(ierr);

    /* do sub-grid temperature correction*/
    PetscScalar *deltaTmSubgrid;
    ierr = PetscMalloc( markers.nMark*sizeof(PetscScalar), &deltaTmSubgrid);CHKERRQ(ierr);
    for(m=0;m<markers.nMark;m++){
      /* calculate temperature that would have been assigned to this marker at last timestep from nodes*/
      /* find my cell */
      PetscInt I,J;
      I = floor( markers.Y[m]/grid.dy);
      J = floor( markers.X[m]/grid.dx);
      PetscScalar dx = grid.dx;
      PetscScalar dy = grid.dy;
      PetscInt idxnode = J*NY+I;
      PetscScalar deltaxm=markers.X[m] - grid.dx*((PetscScalar) J);
      PetscScalar deltaym=markers.Y[m] - grid.dy*(((PetscScalar) I));
      PetscScalar *T = nodalFields.lastT;
      PetscScalar Tmnodeslast = (T[idxnode]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + T[idxnode+NY]*(deltaxm/dx)*(1.0-deltaym/dy) + T[idxnode+1]*(1.0-deltaxm/dx)*(deltaym/dy) + T[idxnode+NY+1]*(deltaxm*deltaym)/dx/dy);
      PetscScalar rhonodeslast = (nodalFields.rho[idxnode]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + nodalFields.rho[idxnode+NY]*(deltaxm/dx)*(1.0-deltaym/dy) + nodalFields.rho[idxnode+1]*(1.0-deltaxm/dx)*(deltaym/dy) + nodalFields.rho[idxnode+NY+1]*(deltaxm*deltaym)/dx/dy);
      PetscScalar Cpnodeslast = (nodalFields.Cp[idxnode]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + nodalFields.Cp[idxnode+NY]*(deltaxm/dx)*(1.0-deltaym/dy) + nodalFields.Cp[idxnode+1]*(1.0-deltaxm/dx)*(deltaym/dy) + nodalFields.Cp[idxnode+NY+1]*(deltaxm*deltaym)/dx/dy);
      PetscScalar knodeslast=(nodalFields.kThermal[idxnode]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + nodalFields.kThermal[idxnode+NY]*(deltaxm/dx)*(1.0-deltaym/dy) + nodalFields.kThermal[idxnode+1]*(1.0-deltaxm/dx)*(deltaym/dy) + nodalFields.kThermal[idxnode+NY+1]*(deltaxm*deltaym)/dx/dy);
      PetscScalar subgridDiffusionConstant=1.0;
      /* gerya equation 10.16 (2) :*/
      PetscScalar dTdiff = Cpnodeslast*rhonodeslast/knodeslast/(2.0/grid.dx2 + 2.0/grid.dy2);
      deltaTmSubgrid[m] = (Tmnodeslast - markers.T[m])*(1.0-exp( -subgridDiffusionConstant*dt/dTdiff));
    }

    /* now take deltaTmSubgrid and project it back onto the nodes */
    
    ierr=projectMarkersNodesFromScalar(&markers, &grid, deltaTmSubgrid, deltaT);
    /* dTremaining = dTnodes-dTsubgrid*/
    for(i=0;i<nnode;i++){
      deltaT[i] = Sa[i]-nodalFields.lastT[i]- deltaT[i];
    }

    /* transfer thermal solution onto markers */
    for(m=0;m<markers.nMark;m++){
      PetscInt I,J;
      /* find cell that this marker belongs to*/
      I = floor( markers.Y[m]/grid.dy);
      J = floor( markers.X[m]/grid.dx);
      PetscScalar dx = grid.dx;
      PetscScalar dy = grid.dy;
      PetscInt idxnode = J*NY+I;
      PetscScalar deltaxm=markers.X[m] - grid.dx*((PetscScalar) J);
      PetscScalar deltaym=markers.Y[m] - grid.dy*(((PetscScalar) I));
      //      if( deltaxm<0.0 || deltaym<0.0 || deltaxm > dx || deltaym > dy) printf("BADBADBAD!!!\n");

      /*             j      j+1   */ 
      /*        i    x       x    */
      /*                m        */
      /* 	i+1  x       x    */
      
      //markers.T[m] += (deltaT[idxnode]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + deltaT[idxnode+NY]*(deltaxm/dx)*(1.0-deltaym/dy) + deltaT[idxnode+1]*(1.0-deltaxm/dx)*(deltaym/dy) + deltaT[idxnode+NY+1]*(deltaxm*deltaym)/dx/dy);
      PetscScalar *dT = deltaT;
      PetscScalar dtm =(dT[idxnode]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + dT[idxnode+NY]*(deltaxm/dx)*(1.0-deltaym/dy) + dT[idxnode+1]*(1.0-deltaxm/dx)*(deltaym/dy) + dT[idxnode+NY+1]*(deltaxm*deltaym)/dx/dy) + deltaTmSubgrid[m]; 
      markers.T[m] += dtm;
      markers.rho[m] += dtm*materials.materialAlpha[(PetscInt) markers.Mat[m]];
      /* use dT to update marker density*/
      
      
    }
    
    ierr = VecRestoreArray(thermalS,&Sa);CHKERRQ(ierr); 
    ierr = VecDestroy( thermalS);CHKERRQ(ierr);
    ierr = PetscFree( deltaT );CHKERRQ(ierr);
    ierr = PetscFree( deltaTmSubgrid); CHKERRQ(ierr);

    saveNodalFieldsASCIIMatlab( &nodalFields, &grid, iTime);

    /* dump shear heating */
/*     FILE *of; */
/*     char name4[80]; */
/*     sprintf(name4,"output/loadT%d.m",iTime); */
/*     of=fopen(name4,"w"); */
/*     fprintf(of,"T=["); */
/*     for(i=0;i<NY;i++){ */
/*       for(j=0;j<NX;j++){ */
/* 	fprintf(of,"%e,",nodalFields.thisT[j*NY+i]); */
/*       } */
/*       fprintf(of,";\n"); */
/*     } */
/*     fprintf(of,"];\n"); */
/*     fclose(of); */

 
    
    /* move the markers*/
    /* 1. project velocity field onto markers*/
    ierr= projectVelocityToMarkers(&markers, &grid, &nodalFields );
    //    printMarkersAllASCIIMatlab( &markers, iTime);
    ierr=saveMarkersBinary( &markers, -1);
    /* 2. rotate marker stresses, total strains*/

    /* 3. advect markers */
    for( m=0;m<markers.nMark;m++){ 
      markers.X[m] +=markers.VX[m]*dt;/* first order*/
      markers.Y[m] +=markers.VY[m]*dt; 

      /*  rotate stresses*/
      PetscScalar espm=markers.w[m] * dt;
      PetscScalar msxx0 = markers.sxx[m];
      PetscScalar msxy0 = markers.sxy[m];
      markers.sxy[m] =msxx0*sin(2.0*espm) + msxy0*cos(2.0*espm);
      markers.sxx[m] =msxx0*(cos(espm)*cos(espm) - sin(espm)*sin(espm)) - msxy0*sin(2.0*espm);

      /* 	if(markers.X[m] > grid.LX) markers.X[m] -= grid.LX; */
      /* 	if(markers.Y[m] > grid.LY) markers.Y[m] -= grid.LY; */
    }	
    elapsedTime+=dt;    
  }/* end time step*/
  
  
  /* clean up*/
  ierr = MatDestroy( thermalLHS);CHKERRQ(ierr);
  ierr = VecDestroy(thermalRHS);CHKERRQ(ierr);
  ierr = VecDestroy(nodalHeating);CHKERRQ(ierr);
  ierr = finalizeNodalFields( &nodalFields);CHKERRQ(ierr);
  // ierr = PetscFree(rho);CHKERRQ(ierr);
  //ierr = PetscFree(etaS);CHKERRQ(ierr);
  //ierr = PetscFree(etaN);CHKERRQ(ierr);
  /* clean up thermal problem */
  //ierr = PetscFree(Cp); CHKERRQ(ierr);
  //ierr = PetscFree(kThermal); CHKERRQ(ierr);
  //ierr = PetscFree(thisT);CHKERRQ(ierr);
  //ierr = PetscFree(lastT); CHKERRQ(ierr);
  /* clean up markers*/
/*   ierr = PetscFree(markerX);CHKERRQ(ierr); */
/*   ierr = PetscFree(markerY);CHKERRQ(ierr); */
/*   ierr = PetscFree(markerVX);CHKERRQ(ierr); */
/*   ierr = PetscFree(markerVY);CHKERRQ(ierr); */
/*  ierr = PetscFree(markerMat);CHKERRQ(ierr);*/
  
  
abort:
  
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}




