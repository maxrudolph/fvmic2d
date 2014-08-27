static char help[] = "Help statement goes here\n";

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
#include "benchmarkInitialConditions.h"
#include "restart.h"
#include "residual.h"
#include "phaseBA.h"

#define TBREAK 1e99 /*273.15*/ /* stop if temperature exceeds this number */

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **args){
  PetscMPIInt rank,size;
  PetscErrorCode ierr;

  /* logging, debugging infrastructure */
  PetscLogStage stages[13];  

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

  PetscInitialize(&argc,&args,(char *)0,help);  //INITIALIZE PETSC
  MPI_Comm_size(PETSC_COMM_WORLD,&size);  //Get MPI rank
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /* log stages */
  PetscLogStageRegister("Main",&stages[0]);
  PetscLogStageRegister("Mech Assem",&stages[1]);
  PetscLogStageRegister("Check Yield",&stages[2]);
  PetscLogStageRegister("Mech Solve",&stages[3]);
  PetscLogStageRegister("Exch Mark",&stages[4]);
  PetscLogStageRegister("Proj M to N",&stages[5]);
  PetscLogStageRegister("Node StrStra",&stages[6]);
  PetscLogStageRegister("Sub, Stres",&stages[7]);
  PetscLogStageRegister("Therm Assem",&stages[8]);
  PetscLogStageRegister("Therm Solve",&stages[9]);
  PetscLogStageRegister("Sub. Temp.",&stages[10]); 
  PetscLogStageRegister("I/O",&stages[11]); 
  PetscLogStageRegister("Mark Str Press",&stages[12]);  
  PetscLogStagePush(stages[0]);

  /* read options from input file and distribute to all nodes */
  printf("argc = %d\n",argc);
  printf("%s\n",args[1]);
  /* check to see if first argument looks like a filename*/
  if( !strncmp( "-",&args[1][0],1)){ //not a file name. pass NULL to csvOptions
    printf("reading from default input file\n");
    ierr=csvOptions((char *)PETSC_NULL, &options,&materials);
  } else {
    ierr=csvOptions( (char *) args[1] , &options,&materials);
  }

  /* begin stuff specific to magma mixing problem */
  char ltf0[80] = "BAPhaseDiagrams/Las2IG.txt";
  char ltf1[80] = "BAPhaseDiagrams/Las2HG.txt";
  materials.compositionLookupTable[0] = loadCompositionLookupTableBAFromFile( ltf0 );
  materials.compositionLookupTable[1] = loadCompositionLookupTableBAFromFile( ltf1 );
  /* end stuff specific to magma mixing problem */
  

  ierr=  PetscRandomCreate(PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  ierr = PetscRandomSetType(r,PETSCRAND48);CHKERRQ(ierr);
  {
    /* seed random number generator with cpu time */
    PetscLogDouble time;
    unsigned long seed;
    if(!rank) ierr = PetscGetTime(&time); CHKERRQ(ierr);
    ierr = MPI_Bcast( &time, sizeof(PetscLogDouble),MPI_BYTE , 0, PETSC_COMM_WORLD);
    seed = (unsigned long) time;
    //    seed = 1;
    printf("seeding random number generator with %ld\n",seed);fflush(stdout);
    ierr = PetscRandomSetSeed(r,seed);CHKERRQ(ierr);
    ierr = PetscRandomSeed(r);CHKERRQ(ierr);/* seed the generator*/
  }
  /* initialize the grid. This version does not require a regular grid*/
  /* note that all processes know the locations of all gridlines and all cell centers */
  if( options.gridRatio == 1.0 ){
    ierr = initializeRegularGrid( &grid, options.LX, options.LY, options.NX, options.NY, &options);CHKERRQ(ierr);
  }else{
    /*     ierr=initializeIrregularGridFault( &grid, options.LX, options.LY, options.NX, options.NY, &options); */
    ierr=initializeIrregularGridConstantInnerOuter( &grid, options.LX, options.LY, options.NX, options.NY, &options);

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
    //printf("[%d], %d %d %d %d %d %d\n",rank,x,y,z,m,n,p);CHKERRQ(ierr);
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
  ierr=DMGetMatrix(grid.da,MATAIJ,&LHSz);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) LHSz,"CmechLHSz");
  {
    PetscInt m,n; 
    ierr=MatGetSize(LHS,&m,&n);CHKERRQ(ierr); 
    printf("LHS is %d by %d, ndof=%d\n",m,n,ndof); 
    ierr=MatGetOwnershipRange(LHS,&m,&n);CHKERRQ(ierr);
    printf("[%d] has rows %d-%d\n",rank,m,n);CHKERRQ(ierr);
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
    printf("beginning MonteCarlo run %d\n",iMonte);
    /* randomize parameters*/
    if(options.nMonte > 1){ierr = getRandomProperties( &materials, &options, &r);CHKERRQ(ierr);}
/*     if(!rank) */ ierr = saveRunInfo( &options, &materials, iMonte);
    /* reset (zero out) all marker properties*/
    resetMarkers( &markerset, &options);
    resetNodalFields( &nodalFields, &grid, &options);

    ierr=distributeMarkersUniformInDomain( &markerset, &options, &grid); CHKERRQ(ierr);

    PetscInt m;
    /* INITIAL CONDITIONS */
 
#ifdef SANDBOX
    ierr = initialConditionsSandbox( &markerset, &options, &materials, &grid); CHKERRQ(ierr);
#else
    //ierr=initialConditionsBending( &markerset, &options, &materials, &grid); CHKERRQ(ierr);
    //ierr= initialConditionsOOPCouette( &markerset, &options, &materials, &grid); CHKERRQ(ierr);
    ierr = initialConditionsMagmaMixingBA( &markerset, &options, &materials, &grid); CHKERRQ(ierr);
    //ierr = initialConditionsConvectionTuran( &markerset, &options, &materials, &grid); CHKERRQ(ierr);
#endif

    PetscInt iTime0=0;/* initial timestep */
    if( options.restartStep ){
      /* destroy current markers */
      destroyMarkers( &markerset ,&options);
      restartFromMarkers( &markerset, &options, iMonte, options.restartStep, &elapsedTime);
      markers = &(markerset.markers[0]);
      iTime0 = options.restartStep + 1;
    }
        
    /* calculate cell to which each marker belongs*/
    findMarkerCells( &markerset, &grid);

    /* project initial condition (thermal) from markers to nodes */
    if(options.restartStep == 0){
      ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials, &options);
      ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);
      saveNodalFieldsMatlab( &nodalFields, &grid, -5,-6,0.0,0.0,0);
      projectNodalFieldToMarkersS(&nodalFields,nodalFields.lastT,&markerset,&(markerset.markers[0].T),&grid);
    }

    /* exchange markers among cpus*/
    PetscLogStagePush(stages[4]);
    ierr = exchangeMarkers( &markerset, &grid, &options);
    PetscLogStagePop();
    findMarkerCells( &markerset, &grid);
    
    /* save markers for debugging*/
    PetscLogStagePush(stages[11]);
    ierr=saveMarkersBinary( &markerset, -5,-5,0.0);
    saveNodalFieldsMatlab( &nodalFields, &grid, -5,-5,0.0,0.0,0);
    PetscLogStagePop();

    PetscInt iTime;
    PetscScalar dt;
    PetscInt isYielding;
    PetscInt iPlastic;
    
    displacementdt=0.0;/* initially*/
    for(iTime=iTime0;iTime<options.nTime;iTime++){/* begin time stepping */
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Beginning Timestep %d, elapsedTime %e years\n",iTime,elapsedTime/(365.25*24*3600));
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
      {/* Check global memory usage */
	PetscLogDouble memuse;
	PetscScalar memuses;
	PetscScalar memuset;
	PetscMemoryGetCurrentUsage( &memuse );
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] memory use %e \n",rank,memuse);
	PetscSynchronizedFlush(PETSC_COMM_WORLD);
	memuses = (PetscScalar) memuse;
	MPI_Allreduce( &memuses, &memuset, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	if(!rank) PetscPrintf(PETSC_COMM_SELF,"Global memory use %e\n",memuset); 
      }
      isYielding = 1;
      iPlastic=0;
      dt=options.dtMax;/* initial timestep*/
      findMarkerCells( &markerset, &grid);
      ierr = exchangeMarkers( &markerset, &grid,&options);
      findMarkerCells( &markerset, &grid);
      /* check marker density, add markers if necessary.*/
      ierr = checkMarkerDensity(&markerset, &grid, &options, r);
      findMarkerCells( &markerset, &grid);

      findMarkerCells( &markerset, &grid);

      /* project all marker fields onto nodes */
      PetscLogStagePush(stages[5]);      ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials,&options);      PetscLogStagePop();


      
      /* update damage- or temperature-dependent viscosity */
      updateDamageViscosity( &grid, &markerset, &materials, &options, dt);

      /* check for overheating (melting). abort if Tmelt exceeded*/
      { /* check for exceedence of maximum temperature */
	PetscInt m;
	PetscScalar maxT = 0.0;
	for(m=0;m<markerset.nMark;m++){
	  if(markers[m].T > maxT) maxT = markers[m].T;
	}
	PetscScalar maxTg;
	MPI_Allreduce( &maxT, &maxTg, 1, MPI_DOUBLE, MPI_MAX,PETSC_COMM_WORLD);
	if( maxTg > TBREAK){
	  printf("Maximum temperature exceeded. Ending Montecarlo step\n");
	  saveNodalFieldsMatlab( &nodalFields, &grid, iMonte,iTime,dt,elapsedTime,0);
	  //saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	  goto nextMonte;
	}
      }
      /* update boundary conditions, which can vary in time */
      updateBoundaryConditions( &options, &boundaryValues, elapsedTime );
      PetscInt converged=0;
      while((isYielding || !converged) && iPlastic < options.maxNumPlasticity){/*begin plastic iteration loop*/
	/* check for plastic yielding */
	isYielding=0;
	
	ierr = checkPlasticYielding(&grid, &markerset, &materials, dt, &isYielding, &options, elapsedTime);
	{
	  /* do an allreduce on isYielding - if any cpu is still yielding, all should still think that they're yielding*/
	  PetscInt tmp1=isYielding;
	  ierr=MPI_Allreduce( &tmp1,&isYielding,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
	}
	if(isYielding) printf("iTime=%d,iPlastic = %d Yielding detected\n",iTime,iPlastic);
	if(!isYielding && !rank) printf("no more yielding\n");
	if( iTime < options.plasticDelay && !rank) {isYielding=0; printf("yielding disabled\n");}
		
	/* project marker quantities onto nodes: kThermal, Cp, rho, T*/
	PetscLogStagePush(stages[5]);
	ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials,&options);
	PetscLogStagePop();
	/* form mechanical problem LHS*/
	/* zero out LHS, RHS*/
	ierr = VecZeroEntries( RHS);CHKERRQ(ierr);	
	ierr = VecZeroEntries( RHSz); CHKERRQ(ierr);
	PetscLogStagePush(stages[1]);
	ierr=formVEPSystem( &nodalFields, &grid, LHS, RHS,LHSz,RHSz, &Kbond, &Kcont, options.gy, dt, &options, &boundaryValues);
	PetscLogStagePop();
	Vec mechanicalS;/* solution to mechanical problem*/
	ierr = VecDuplicate(RHS, &mechanicalS);
	ierr = VecZeroEntries(mechanicalS);
	PetscLogStagePush(stages[3]);
	ierr = kspLinearSolve(LHS, RHS, mechanicalS);
	PetscLogStagePop();
	Vec Sz;/* solution to out-of-plane problem*/
	ierr = VecDuplicate(RHSz, &Sz);
	ierr = VecZeroEntries( Sz );CHKERRQ(ierr);
	//ierr = kspLinearSolve(LHSz, RHSz, Sz);CHKERRQ(ierr);
	if(!rank) printf("out=of-plane solution disabled\n");
	/* move solutions from solution vectors to nodalFields*/
	PetscInt fn=0;/* flag for nan in solution*/
	ierr= retrieveSolutions( &grid, &nodalFields, mechanicalS, Sz, Kcont, &fn);
	ierr = VecDestroy(&mechanicalS); CHKERRQ(ierr);
	ierr = VecDestroy(&Sz); CHKERRQ(ierr);
	if(fn) {
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Found nan in solution\n");CHKERRQ(ierr);
	  PetscLogStagePush(stages[11]);
	  saveNodalFieldsMatlab( &nodalFields, &grid, iMonte,iTime,dt,elapsedTime,0);
	  ierr=saveMarkersBinary( &markerset, iMonte,iTime,elapsedTime);
	  PetscLogStagePush(stages[11]);

	  if(iMonte < options.nMonte-1){
	    goto nextMonte;
	  }else{
	    goto abort;
	  }
	}/* end deal with nan in solution*/

	/* calculate global strain rate residual */
	updateGlobalStrainRateResidual( &grid, &markerset, iPlastic );
	converged = checkConvergence( 1e-16, 1e-3, iPlastic);

	displacementdt = dt;
	/* determine displacement timestep*/
	ierr=limitDisplacementTimestep( &grid, &nodalFields, &displacementdt,&options);CHKERRQ(ierr);

	ierr=PetscLogStagePush(stages[6]);
	ierr=nodalStressStrain(&grid, &nodalFields,&options,&boundaryValues, displacementdt,nodalHeating, options.gy); CHKERRQ(ierr);
	ierr=PetscLogStagePop();
	/* save nodalFields for debugging*/
	PetscLogStagePush(stages[11]);
	/*ierr=saveNodalFieldsASCIIMatlab( &nodalFields, &grid,iMonte,iTime, dt, elapsedTime);*/ /* for debugging*/
	PetscLogStagePop();
	/* compute marker strain and pressure*/
	PetscLogStagePush(stages[12]); ierr=updateMarkerStrainPressure( &grid,&nodalFields, &markerset, &materials, &options, displacementdt); PetscLogStagePop();
	if(!rank) printf("Timestep %d iteration %d, displacementdt = %e, yielding %d, residual %e, converged = %d\n",iTime,iPlastic,displacementdt,isYielding,getGlobalStrainRateResidual( iPlastic ),converged);
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
      PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(thermalLHS, thermalRHS, thermalS);CHKERRQ(ierr); PetscLogStagePop();
      /* check maximum temperature change */
      Vec dT;
      ierr = VecDuplicate( thermalS, &dT);CHKERRQ(ierr);
      ierr = VecCopy(thermalS, dT);CHKERRQ(ierr);
      ierr = VecAXPY(dT,-1.0,nodalFields.lastT);CHKERRQ(ierr);
      PetscScalar dTmax;/* max abs change in temperature*/
      ierr = VecNorm(dT,NORM_INFINITY,&dTmax);CHKERRQ(ierr);

      if(dTmax > options.maxTempChange){
	printf("[%d] temperature change too large, limiting timestep from %e to %e\n",rank,dt,dt*options.maxTempChange/dTmax);
	dt = dt*options.maxTempChange/dTmax;
	/* repeat temperature solution */
	PetscLogStagePush( stages[8] );
	ierr = VecCopy(nodalHeating,thermalRHS);CHKERRQ(ierr);
	ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);
	ierr = formThermalSystem( &grid, &nodalFields, thermalLHS, thermalRHS, dt, &options);CHKERRQ(ierr); PetscLogStagePop();
	PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(thermalLHS, thermalRHS, thermalS);CHKERRQ(ierr); PetscLogStagePop();
	/* end repeat temperature solution */
      }
      ierr = VecDestroy(&dT);CHKERRQ(ierr);

      /* Apply sub-grid temperature corrections and project corrected temperature onto the markers (all integrated into one routine*/
      PetscLogStagePush( stages[10] ); ierr=subgridTemperatureChanges(thermalS, &grid, &nodalFields, &markerset, &materials,dt, &options);      PetscLogStagePop();

      ierr = VecDestroy(&thermalS); /* free the thermal solution vector */
      
      /* 1. project velocity field onto markers - only needed if using 1st order advection*/
      //ierr= projectVelocityToMarkers(&markerset, &grid, &nodalFields );
      
      /* update damage rate on markers */
      updateDamageRate( &grid, &markerset, &materials, &options);
      
      /* advect markers */
      ierr=advectMarkersRK(&markerset,&nodalFields,&grid,&options,&boundaryValues,dt);
#ifdef SANDBOX
      /*       ierr= sandboxMobileWall( &markerset, &grid, &materials); */
      ierr= sandboxMobileWallAndTimesteps( &markerset, &grid, &materials, &options, iTime, elapsedTime);
#endif
      /* save solution */
      if(!(iTime % (options.saveInterval))) {
	ierr=saveNodalFieldsMatlab( &nodalFields, &grid,iMonte,iTime, dt, elapsedTime,0);
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
      

    }/* end time stepping*/
  nextMonte:;


  }/* end montecarlo loop */
 abort:;

  //ierr = VecDestroy(&thermalRHS);CHKERRQ(ierr);
  //ierr = VecDestroy(&RHSz);CHKERRQ(ierr);
  //ierr = VecDestroy(&RHS);CHKERRQ(ierr);
  ierr = VecDestroy(&nodalHeating);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(grid.vda,&RHS);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(grid.da,&RHSz);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(grid.da,&thermalRHS);CHKERRQ(ierr);

  ierr = MatDestroy( &thermalLHS );CHKERRQ(ierr);
  ierr = MatDestroy( &LHSz );CHKERRQ(ierr);
  ierr = MatDestroy( &LHS );CHKERRQ(ierr);

  ierr = destroyGrid(&grid);CHKERRQ(ierr);
  ierr = destroyNodalFields(&nodalFields, &grid);CHKERRQ(ierr);
  ierr=  destroyMarkers( &markerset,&options );CHKERRQ(ierr);
  ierr = PetscRandomDestroy( &r);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
 
