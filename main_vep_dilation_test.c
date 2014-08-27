
//#include "../include/fem.h"
//
/*
 *Main routine of FEM code. This copy is specifically for ridge damage problem
 *Do not put any compiler macros here. They should be in fem.h.
 */
static char help[] = "Help statement goes here\n";
/*petscmpiexec -n 2 ./Fem ./input/fault.1 output/test -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps >log1\n \ */
   //It can be one line or many \n";

#include<stdio.h>
#include<stdlib.h>
   //#include<malloc.h>
#include "petscksp.h"
#include "petsctime.h"
#include "mpi.h"
#include "fdcode.h"
#include "kspLinearSolve.h"
#include "markers.h"
#include "grid.h"
#include "nodalFields.h"
#include "markers.h"
#include "materialProperties.h"
#include "thermalSystem.h"
   //#include "stokes_varvisc.h"
#include "vep_system.h"
#include "nodalStressStrain.h"/* includes shear heating too*/
   //#include "shearHeating.h"
#include "updateMarkerStrainPressure.h"
#include "subgridStressChanges.h"
#include "subgridTemperatureChanges.h"
#include "updateMarkerStress.h"
#include "updateDamageViscosity.h"
#include "io.h"
#include "options.h"
#include "randomProperties.h"
#include "retrieveSolutions.h"
#include "limitTimestep.h"
#include "updateBCs.h"
#include "adiabaticHeating.h"
#include "memuse.h"
#include "version.h"
#include "benchmarkInitialConditions.h"
#include "markerProjection.h"
#include "fsPressureCorrection.h"
#include "profile.h"
#include "pressureNullSpace.h"

//Include file for PETSc linear solve commands

//extern struct ELEMENT linearHexahedralElementCreate();
//extern double evaluateNodalBasisFunction();

#define TBREAK 273.15 /*273.15*/ /* stop if temperature exceeds this number */
#define INITIALDT options.dtMax
/*3600.0*/ /* take ten really small timesteps (10 s here) before increasing timestep later */
//extern struct ELEMENT quadTriElementCreate();
//extern double evaluateNodalBasisFunction();

#define NADD 12 
//maximum number of entries to add at once.

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **args){
  PetscMPIInt rank,size;
  PetscErrorCode ierr;

  /* logging, debugging infrastructure */
  PetscLogStage *stages;  

  /* Grid, field variables */
  GridData grid;
  NodalFields nodalFields;
  
  /* Linear Systems */
  Mat LHS, thermalLHS;
  Mat LHSz;
  Vec RHS, thermalRHS, nodalHeating;
  Vec RHSz;

  /* markers */
  MarkerSet markerset; /* data structure that holds all of the marker info*/
  Marker *markers;     /* pointer to the markers themselves */
  /* material properties, boundary values and options */
  BoundaryValues boundaryValues; /* used to store time-varying boundary condition information*/
  Materials materials;
  PetscRandom r;
  Options options;

  /* BEGIN PROGRAM */

  PetscInitialize(&argc,&args,PETSC_NULL,help);  //INITIALIZE PETSC
  MPI_Comm_size(PETSC_COMM_WORLD,&size);  //Get MPI rank
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /* Print welcome message */
  if(!rank) printf("Version Information : %s\n",MARKERCODE_HG_VERSION);
  stages = initializeLogging();

  /* read options from input file and distribute to all nodes */
  {
    PetscBool set = PETSC_FALSE;
    char filename[80];
    ierr = PetscOptionsGetString(PETSC_NULL,"-input_file",filename,sizeof(filename),&set);CHKERRQ(ierr);
    if(set){
      ierr=csvOptions(filename , &options,&materials);
    }else{
      printf("reading from default input file\n");
      ierr=csvOptions((char *)PETSC_NULL, &options,&materials);
    }
  }

  /* initialize material properties*/
  //setMaterialProperties( &materials);

  ierr=  PetscRandomCreate(PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  ierr = PetscRandomSetType(r,PETSCRAND48);CHKERRQ(ierr);
  {
    /* seed random number generator with cpu time */
    PetscLogDouble time;
    unsigned long seed;
    if(!rank) ierr = PetscGetTime(&time); CHKERRQ(ierr);
    ierr = MPI_Bcast( &time, sizeof(PetscLogDouble),MPI_BYTE , 0, PETSC_COMM_WORLD);
    seed = (unsigned long) time;
    printf("seeding random number generator with %ld\n",seed);fflush(stdout);
    ierr = PetscRandomSetSeed(r,seed);CHKERRQ(ierr);
    ierr = PetscRandomSeed( r );CHKERRQ(ierr);/* seed the generator*/
  }
  /* initialize the grid. This version does not require a regular grid*/
  /* note that all processes know the locations of all gridlines and all cell centers */
  //ierr = initializeRegularGrid( &grid, options.LX, options.LY, options.NX, options.NY);CHKERRQ(ierr);
  if( options.gridRatio == 1.0 ){
    ierr=initializeRegularGrid( &grid, options.LX, options.LY, options.NX, options.NY, &options);
  }else{
     ierr=initializeIrregularGridFaultFS( &grid, options.LX, options.LY, options.NX, options.NY, &options, 1.0/5.0 ); 
     //    ierr=initializeIrregularGridFault( &grid, options.LX, options.LY, options.NX, options.NY, &options);
  }
  /* save the grid */
  ierr = saveGrid( &grid );

  /* initialize scaling parameters*/
  PetscInt ndof = grid.NX*grid.NY*3;  
  PetscScalar Kbond = 0.0;
  PetscScalar Kcont = 0.0;
  /* allocate the nodal fields. The nodal fields are all global PETSc Vecs and nodalFields->da is the associated DA*/
  ierr = initializeNodalFields( &nodalFields, &grid, &options);CHKERRQ(ierr);
  {
    /*allocate the markers. First, calculate the right number of markers for each processor.*/
    int x,y,z,m,n,p;
    ierr=DMDAGetCorners(grid.da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
    ierr=allocateMarkers( (PetscInt)(options.maxMarkFactor*m*options.NMX*n*options.NMY), &markerset, &options);CHKERRQ(ierr); 
    markers = markerset.markers;
    markerset.nMark = n*m*options.NMX*options.NMY;
  }
  /* set up local to global mapping*/
  ISLocalToGlobalMapping ltogm;
  ierr=DMGetLocalToGlobalMapping(grid.da, &ltogm);CHKERRQ(ierr);
  
  /* initialize LHS,LHSz matrix and RHS,RHSz matrix for mechanical problem */
  ierr=DMGetMatrix(grid.vda,MATAIJ,&LHS);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) LHS,"CmechLHS");

  MatNullSpace ns;
  ierr = createPressureNullSpace( &grid, &ns );CHKERRQ(ierr); 

  ierr=DMGetMatrix(grid.da,MATAIJ,&LHSz);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) LHSz,"CmechLHSz");
  {
    PetscInt m,n; 
    ierr=MatGetSize(LHS,&m,&n);CHKERRQ(ierr); 
    ierr=MatGetOwnershipRange(LHS,&m,&n);CHKERRQ(ierr);
  }

  ierr=DMGetGlobalVector(grid.vda,&RHS);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) RHS,"CmechRHS"); 
  ierr=DMGetGlobalVector(grid.da,&RHSz);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) RHSz,"CmechRHSz"); 
  /* initialize LHS matrix, RHS Vec and nodalHeating Vec for thermal problem*/
  ierr=DMGetGlobalVector(grid.da,&thermalRHS);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) thermalRHS,"CthermalRHS");
  ierr=DMGetMatrix(grid.da,MATAIJ,&thermalLHS);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) thermalLHS,"CthermalLHS");
  ierr=VecDuplicate(thermalRHS,&nodalHeating);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) nodalHeating,"nodalHeating");

  PetscInt iMonte;
  for(iMonte=0;iMonte<options.nMonte;iMonte++){
    PetscScalar displacementdt;/* this is the displacement timestep*/
    PetscScalar elapsedTime=0.0;
    if(!rank) printf("beginning MonteCarlo run %d\n",iMonte);
    /* randomize parameters*/
    if(options.doMonte){ ierr = getRandomProperties( &materials, &options, &r); CHKERRQ(ierr);}
    /*if(!rank)*/ ierr = saveRunInfo( &options, &materials, iMonte);

    /* reset (zero out) all marker properties*/

    if( iMonte > 0 ){/* destroy old markers and make new ones */
      ierr=  destroyMarkers( &markerset,&options );CHKERRQ(ierr);
      int x,y,z,m,n,p;
      ierr=DMDAGetCorners(grid.da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
      ierr=allocateMarkers( (PetscInt)(options.maxMarkFactor*m*options.NMX*n*options.NMY), &markerset, &options);CHKERRQ(ierr);
      markers = markerset.markers;
      markerset.nMark = n*m*options.NMX*options.NMY;
    }

    resetMarkers( &markerset, &options);
    resetNodalFields( &nodalFields, &grid, &options);

    ierr=distributeMarkersUniformInDomain( &markerset, &options, &grid);CHKERRQ(ierr);
    
    PetscInt m;
    /* INITIAL CONDITIONS */
    ierr= initialConditionsDilationTest2( &markerset, &options, &materials, &grid);
    //ierr = initialConditionsInclusion( &markerset, &options, &materials, &grid);

    updateDamageViscosity( &grid, &markerset, &materials, &options, 0.0);        
    /* calculate cell to which each marker belongs*/
    findMarkerCells( &markerset, &grid);

    /* project initial condition (thermal) from markers to nodes */
    //    ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials, &options);
    ierr = projectMarkersNodesAll2(&markerset, &grid, &nodalFields, &materials, &options);
    ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);
    saveNodalFields( &nodalFields, &grid, -iMonte,-6,0.0,0.0,0);
    projectNodalFieldToMarkersS(&nodalFields,nodalFields.lastT,&markerset,&(markerset.markers[0].T),&grid);


    /* exchange markers among cpus*/

    ierr = exchangeMarkers( &markerset, &grid, &options);

    findMarkerCells( &markerset, &grid);
    
    /* save markers for debugging*/
    saveNodalFields( &nodalFields, &grid, -iMonte,-5,0.0,0.0,0);
    ierr=saveMarkersBinary( &markerset, -iMonte,-5,0.0);

    PetscInt iTime;
    PetscScalar dt;
    PetscInt isYielding;
    PetscInt iPlastic;
    
    displacementdt=0.0;/* initially*/
    for(iTime=0;iTime<options.nTime;iTime++){
      /*if(!rank)*/ printf("Beginning Timestep %d, elapsedTime %e years\n",iTime,elapsedTime/(365.25*24*3600));
      /* get memory use */
      getMemoryUse();

      isYielding = 1;
      iPlastic=0;
      dt=options.dtMax;/* initial timestep*/
      findMarkerCells( &markerset, &grid);
      ierr = exchangeMarkers( &markerset, &grid,&options);
      findMarkerCells( &markerset, &grid);
      /* check marker density, add markers if necessary.*/
      ierr = checkMarkerDensity(&markerset, &grid, &options, r);
      findMarkerCells( &markerset, &grid);
  
      /* project rhodot to nodes*/
      findMarkerCells( &markerset, &grid);

      /* project all marker fields onto nodes */
      ierr = projectMarkersNodesAll2(&markerset, &grid, &nodalFields, &materials,&options);

      updateDamageViscosity( &grid, &markerset, &materials, &options, dt);
      /* check for overheating (melting). abort if Tmelt exceeded*/
      { /* check for exceedence of maximum temperature */
	PetscInt m;
	PetscScalar maxT = 0.0;
	for(m=0;m<markerset.nMark;m++){
	  /* make sticky air isothermal and cold */
	  if(markers[m].Mat == 2) markers[m].T = options.thermalBCTop.value[0];
	  if(markers[m].T > maxT) maxT = markers[m].T;
	}
	PetscScalar maxTg;
	MPI_Allreduce( &maxT, &maxTg, 1, MPI_DOUBLE, MPI_MAX,PETSC_COMM_WORLD);
	if( maxTg > TBREAK){
	  if(!rank) printf("Maximum temperature exceeded. Ending Montecarlo step\n");
	  saveNodalFields( &nodalFields, &grid, iMonte,iTime,dt,elapsedTime,0);
	  //ierr=saveMarkersBinary( &markerset, iMonte,iTime,elapsedTime);
	  //saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	  goto nextMonte;
	}
	for(m=0;m<markerset.nMark;m++){
	  if( markers[m].Mat == 1){
	    markers[m].rhodot = -3.17e-10;
	  }
	}       
      }
      /* update boundary conditions, which can vary in time */
      updateBoundaryConditions( &options, &boundaryValues, elapsedTime );
      if( iTime < 10){
	displacementdt = INITIALDT;
      }else{
	displacementdt = options.dtMax;/* try for maximum displacement timestep */
      }
      while(isYielding && iPlastic < options.maxNumPlasticity){/*begin plastic iteration loop*/
	/* check for plastic yielding */
	isYielding=0;

	//ierr = projectMarkersNodesAll2(&markerset, &grid, &nodalFields, &materials,&options);

	ierr = checkPlasticYielding(&grid, &markerset, &materials, displacementdt, &isYielding, &options, elapsedTime);

	ierr = projectMarkersNodesAll2(&markerset, &grid, &nodalFields, &materials,&options);

	if(!rank) printf("iTime=%d,iPlastic = %d\n",iTime,iPlastic);
	if(!isYielding && !rank) printf("no yielding\n");
	if( isYielding && !rank) printf("yielding detected\n");
	if( iTime < options.plasticDelay && !rank) {isYielding=0; printf("yielding disabled\n");}
		
	/* form mechanical problem LHS*/
	/* zero out LHS, RHS*/
	ierr = VecZeroEntries( RHS);CHKERRQ(ierr);	
	ierr = VecZeroEntries( RHSz); CHKERRQ(ierr);

	ierr=formVEPSystem( &nodalFields, &grid, LHS,PETSC_NULL, RHS,LHSz,RHSz, &Kbond, &Kcont, options.gy, dt, &options, &boundaryValues);

	Vec mechanicalS;/* solution to mechanical problem*/
	ierr = VecDuplicate(RHS, &mechanicalS);

	ierr = kspLinearSolve(LHS,PETSC_NULL, RHS, mechanicalS,"stokes_",ns);

	Vec Sz;/* solution to out-of-plane problem*/
	ierr = VecDuplicate(RHSz, &Sz);
	ierr = VecZeroEntries( Sz );CHKERRQ(ierr);
	ierr = kspLinearSolve(LHSz, PETSC_NULL, RHSz, Sz,"zstokes_",PETSC_NULL);CHKERRQ(ierr);
	/* move solutions from solution vectors to nodalFields*/
	PetscInt fn=0;/* flag for nan in solution*/
	ierr= retrieveSolutions( &grid, &nodalFields, mechanicalS, Sz, Kcont, &fn);
	
	ierr= freeSurfacePressureCorrection( &grid, &markerset, &nodalFields, &materials, 2 );

	ierr = VecDestroy(&mechanicalS); CHKERRQ(ierr);
	ierr = VecDestroy(&Sz); CHKERRQ(ierr);
	if(fn){
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Found nan in solution\n");CHKERRQ(ierr);
	  saveNodalFields( &nodalFields, &grid, iMonte,iTime,dt,elapsedTime,0);
	  ierr=saveMarkersBinary( &markerset, iMonte,iTime,elapsedTime);
	  if(iMonte < options.nMonte-1){
	    goto nextMonte;
	  }else{
	    goto abort;
	  }
	} /* end deal with nan in solution*/
	
	//ierr = VecZeroEntries(nodalHeating);CHKERRQ(ierr);
	ierr=PetscLogStagePush(stages[6]);
	
	displacementdt = dt;/* displacement dt needs to be <= dt */
	ierr=limitDisplacementTimestep( &grid, &nodalFields, &displacementdt, &options); CHKERRQ(ierr);	

	/* nodalStressStrain also forecasts damage rate using calculated strain rate and pressure */
	ierr=nodalStressStrain(&grid, &nodalFields,&options, &boundaryValues, displacementdt, nodalHeating, options.gy);CHKERRQ(ierr); ierr=PetscLogStagePop();
	
	/* compute marker strain and pressure*/
	PetscLogStagePush(stages[12]); ierr=updateMarkerStrainPressure( &grid,&nodalFields, &markerset, &materials, &options, displacementdt); PetscLogStagePop();
	
	/* if(lastdt != displacementdt)  */
	/* 	isYielding = 1; */ /* timestep is still adapting to density changes - stay in nonlinear loop */
	
	/* calculate new effective viscosity without updating damage variable */
	{PetscInt m; /* rhodot and Ddot are stored on marker - use them to calculate a new damage value and effective viscosity for the next plastic iteration */
	  for(m=0;m<markerset.nMark;m++){
	    //markerDamageViscosityNoUpdate( &markers[m], &materials, &options, displacementdt);
	  }
	}
	/* update marker effective viscosity to match marker strain reate*/
	iPlastic++;
      }/* end plasticity loop*/
      if(!rank) printf("Done with plasticity loops.\n");
      /* do sub-grid stress diffusion*/
      PetscLogStagePush(stages[7]);
      if(options.subgridStressDiffusivity>0.0)	ierr=subgridStressChanges(&grid, &nodalFields, &markerset, &materials, displacementdt,options.subgridStressDiffusivity);
      PetscLogStagePop();
      /* update marker stress*/
      ierr= updateMarkerStress( &grid, &nodalFields,&markerset, &materials);CHKERRQ(ierr);
      /* set thermal problem timestep to displacement dt*/
      dt = displacementdt;
     
      /* Form thermal system (diffusion, variable coefficient, non-uniform grid */
      
      if( !options.shearHeating ){
	ierr=VecZeroEntries(nodalHeating);CHKERRQ(ierr); /* this line disables shear heating */
      }
      if( options.adiabaticHeating ){
	adiabaticHeating( &grid, &markerset, &nodalFields, &materials, &options);
	ierr=VecAXPY(nodalHeating,1.0,nodalFields.ha);CHKERRQ(ierr);/* add adiabatic heating to shear heating */
      }
      ierr = VecCopy(nodalHeating,thermalRHS);CHKERRQ(ierr);
      ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);

      PetscLogStagePush( stages[8] );
      ierr = formThermalSystem( &grid, &nodalFields, thermalLHS, thermalRHS, dt, &options);CHKERRQ(ierr);
      PetscLogStagePop();

      /* make the solution vector */
      Vec thermalS;
      ierr = VecDuplicate( thermalRHS, &thermalS);CHKERRQ(ierr);
      /* Solve the thermal system: */
      PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(thermalLHS,PETSC_NULL, thermalRHS, thermalS,"energy_",PETSC_NULL);CHKERRQ(ierr); PetscLogStagePop();
      /* check maximum temperature change */
      Vec dT;
      ierr = VecDuplicate( thermalS, &dT);CHKERRQ(ierr);
      ierr = VecCopy(thermalS, dT);CHKERRQ(ierr);
      ierr = VecAXPY(dT,-1.0,nodalFields.lastT);CHKERRQ(ierr);
      PetscScalar dTmax;/* max abs change in temperature*/
      ierr = VecNorm(dT,NORM_INFINITY,&dTmax);CHKERRQ(ierr);

      if(dTmax > options.maxTempChange){
	if(!rank) printf("[%d] temperature change too large, limiting timestep from %e to %e\n",rank,dt,dt*options.maxTempChange/dTmax);
	dt = dt*options.maxTempChange/dTmax;
	/* repeat temperature solution */
	PetscLogStagePush( stages[8] );
	ierr = VecCopy(nodalHeating,thermalRHS);CHKERRQ(ierr);
	ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);
	ierr = formThermalSystem( &grid, &nodalFields, thermalLHS, thermalRHS, dt, &options);CHKERRQ(ierr); PetscLogStagePop();
	PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(thermalLHS,PETSC_NULL, thermalRHS, thermalS,"energy_",PETSC_NULL);CHKERRQ(ierr); PetscLogStagePop();
	/* end repeat temperature solution */
      }
      ierr = VecDestroy(&dT);CHKERRQ(ierr);

      /* Apply sub-grid temperature corrections and project corrected temperature onto the markers (all integrated into one routine*/
      PetscLogStagePush( stages[10] ); ierr=subgridTemperatureChanges(thermalS, &grid, &nodalFields, &markerset, &materials,dt, &options);      PetscLogStagePop();

      ierr = VecDestroy(&thermalS); /* free the thermal solution vector */

      //updateDamageRate( &grid, &markerset, &materials, &options);

      /* save solution */
      if(!(iTime % (options.saveInterval))) {
	ierr=saveNodalFields( &nodalFields, &grid,iMonte,iTime, dt, elapsedTime,0);
	ierr=saveMarkersBinary( &markerset, iMonte,iTime,elapsedTime);
      }
      if(!(iTime % options.saveInterval)) {
	//saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime,elapsedTime);
      }
      
      elapsedTime+=dt;
      if(elapsedTime > options.totalTime || iTime==options.nTime-1){
	PetscLogStagePush(stages[11]);
	//saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	PetscLogStagePop();
	ierr=saveMarkersBinary( &markerset,iMonte,iTime,elapsedTime);
	goto nextMonte;
      }

      if(dt < 1.0) {	
	if(!rank) printf("timestep less than 1s, aborting run\n");fflush(stdout);
	PetscLogStagePush(stages[11]);
	ierr=saveMarkersBinary( &markerset, iMonte,iTime,elapsedTime);
	saveNodalFields( &nodalFields, &grid, iMonte,iTime,dt,elapsedTime,0);
	PetscLogStagePop();
	goto nextMonte;
      }


      /* advect markers */

      ierr=advectMarkersRK(&markerset,&nodalFields,&grid,&options,&boundaryValues,dt);

      /* update marker density and damage and viscosity using stored rates*/
      {
	PetscInt m;
	for(m=0;m<markerset.nMark;m++){
	  if(markers[m].cellX != -1){
	    if(markers[m].rho > materials.materialRho[(PetscInt) markers[m].Mat]) markers[m].rho = materials.materialRho[(PetscInt) markers[m].Mat];
	    if(materials.hasDilation[ markers[m].Mat ] == 0){
	      //markers[m].rho = materials.materialRho[(PetscInt) markers[m].Mat]*(1-materials.materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
	    }
	    markers[m].rho += dt*markers[m].rhodot;
	    /* update damage variable and viscosity */
	    //updateMarkerDamageViscosity( &markers[m], &materials, &options, dt);
	  }
	}
      }


      
    }/* end time stepping*/
  nextMonte:;


  }/* end montecarlo loop */
 abort:;

  //ierr = VecDestroy(thermalRHS);CHKERRQ(ierr);
  //ierr = VecDestroy(RHSz);CHKERRQ(ierr);
  //ierr = VecDestroy(RHS);CHKERRQ(ierr);
  ierr = VecDestroy(&nodalHeating);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(grid.vda,&RHS);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(grid.da,&RHSz);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(grid.da,&thermalRHS);CHKERRQ(ierr);
  ierr = destroyPressureNullSpace( &ns ); CHKERRQ(ierr);
  ierr = MatDestroy( &thermalLHS );CHKERRQ(ierr);
  ierr = MatDestroy( &LHSz );CHKERRQ(ierr);
  ierr = MatDestroy( &LHS );CHKERRQ(ierr);

  ierr = destroyGrid(&grid);CHKERRQ(ierr);
  ierr = destroyNodalFields(&nodalFields, &grid);CHKERRQ(ierr);
  ierr=  destroyMarkers( &markerset,&options );CHKERRQ(ierr);
  ierr = PetscRandomDestroy( &r);CHKERRQ(ierr);
  finalizeLogging();
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
 
