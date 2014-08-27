static char help[] = "Help statement goes here\n";
/*petscmpiexec -n 2 ./Fem ./input/fault.1 output/test -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps >log1\n \ */


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
#include "profile.h"
#include "markerProjection.h"


#define TBREAK 273.15 /*273.15*/ /* stop if temperature exceeds this number */
#define INITIALDT problem.options.dtMax
/*3600.0*/ /* take ten really small timesteps (10 s here) before increasing timestep later */


#define NADD 12 
//maximum number of entries to add at once.

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **args){
  PetscMPIInt rank,size;
  PetscErrorCode ierr;

  /* logging, debugging infrastructure */
  PetscLogStage *stages;  

  /* problem context */
  Problem problem;
  problem.mech_system.pc = PETSC_NULL;

  /* material properties, boundary values, options */
  BoundaryValues boundaryValues; /* used to store time-varying boundary condition information*/
  PetscRandom r;

  /* Grid, field variables */
  
  /* BEGIN PROGRAM */

  PetscInitialize(&argc,&args,(char *)0,help);  //INITIALIZE PETSC
  MPI_Comm_size(PETSC_COMM_WORLD,&size);  //Get MPI rank
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /* Print welcome message */
  if(!rank) printf("Version Information : %s\n",MARKERCODE_HG_VERSION);
  /* read options from input file and distribute to all nodes */
  if(!rank) printf("argc = %d\n",argc);
  if(!rank) printf("%s\n",args[1]);

  /* log stages */
  stages = initializeLogging();
  {
    PetscBool set = PETSC_FALSE;
    char filename[80];
    ierr = PetscOptionsGetString(PETSC_NULL,"-input_file",filename,sizeof(filename),&set);CHKERRQ(ierr);
    if(set){
      ierr=csvOptions(filename , &problem.options,&problem.materials);
    }else{
      printf("reading from default input file\n");
      ierr=csvOptions((char *)PETSC_NULL, &problem.options,&problem.materials);
    }
  }

  /* initialize material properties*/
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
    ierr = PetscRandomSeed(r);CHKERRQ(ierr);/* seed the generator*/
  }
  /* initialize the grid. This version does not require a regular grid*/
  /* note that all processes know the locations of all gridlines and all cell centers */
  //ierr = initializeRegularGrid( &grid, options.LX, options.LY, options.NX, options.NY);CHKERRQ(ierr);
  if( problem.options.gridRatio == 1.0 ){
    ierr=initializeRegularGrid( &problem.grid, problem.options.LX, problem.options.LY, problem.options.NX, problem.options.NY, &problem.options);
  }else{
    ierr=initializeIrregularGridFaultFS( &problem.grid, problem.options.LX, problem.options.LY, problem.options.NX, problem.options.NY, &problem.options, 1.0/5.0);
  }
  /* save the grid */
  ierr = saveGrid( &problem.grid );

  /* initialize scaling parameters*/
  PetscInt ndof = problem.grid.NX*problem.grid.NY*3;  
  PetscScalar Kbond = 0.0;
  PetscScalar Kcont = 0.0;
  /* allocate the nodal fields. The nodal fields are all global PETSc Vecs and nodalFields->da is the associated DA*/
  ierr = initializeNodalFields( &problem.nodal_fields, &problem.grid, &problem.options);CHKERRQ(ierr);
  {
    /*allocate the markers. First, calculate the right number of markers for each processor.*/
    int x,y,z,m,n,p;
    ierr=DMDAGetCorners(problem.grid.da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
    ierr=allocateMarkers( (PetscInt)(problem.options.maxMarkFactor*m*problem.options.NMX*n*problem.options.NMY), &problem.markerset, &problem.options);CHKERRQ(ierr); 
    problem.markerset.nMark = n*m*problem.options.NMX*problem.options.NMY;
  }
  /* set up local to global mapping*/
  ISLocalToGlobalMapping ltogm;
  ierr=DMGetLocalToGlobalMapping(problem.grid.da, &ltogm);CHKERRQ(ierr);
  
  /* initialize LHS,LHSz matrix and RHS,RHSz matrix for mechanical problem */
  ierr=DMGetMatrix(problem.grid.vda,MATAIJ,&problem.mech_system.lhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem.mech_system.lhs,"CmechLHS");
  ierr=DMGetMatrix(problem.grid.da,MATAIJ,&problem.z_system.lhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem.z_system.lhs,"CmechLHSz");
  {
    PetscInt m,n; 
    ierr=MatGetSize(problem.mech_system.lhs,&m,&n);CHKERRQ(ierr); 
    ierr=MatGetOwnershipRange(problem.mech_system.lhs,&m,&n);CHKERRQ(ierr);
  }

  ierr=DMGetGlobalVector(problem.grid.vda,&problem.mech_system.rhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem.mech_system.rhs,"CmechRHS"); 
  ierr=DMGetGlobalVector(problem.grid.da,&problem.z_system.rhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem.z_system.rhs,"CmechRHSz"); 
  /* initialize LHS matrix, RHS Vec and nodalHeating Vec for thermal problem*/
  ierr=DMGetGlobalVector(problem.grid.da,&problem.thermal_system.rhs);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem.thermal_system.rhs,"CthermalRHS");
  ierr=DMGetMatrix(problem.grid.da,MATAIJ,&problem.thermal_system.lhs);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem.thermal_system.lhs,"CthermalLHS");
  ierr=VecDuplicate(problem.thermal_system.rhs,&problem.nodal_heating);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem.nodal_heating,"nodalHeating");

  PetscInt iMonte;
  for(iMonte=0;iMonte<problem.options.nMonte;iMonte++){
    PetscScalar displacementdt;/* this is the displacement timestep*/
    PetscScalar elapsedTime=0.0;
    if(!rank) printf("beginning MonteCarlo run %d\n",iMonte);
    /* randomize parameters*/
    if(problem.options.doMonte){ ierr = getRandomProperties( &problem.materials, &problem.options, &r); CHKERRQ(ierr);}
    /*if(!rank)*/ ierr = saveRunInfo( &problem.options, &problem.materials, iMonte);

    /* reset (zero out) all marker properties*/

    if( iMonte > 0 ){/* destroy old markers and make new ones */
      ierr=  destroyMarkers( &problem.markerset,&problem.options );CHKERRQ(ierr);
      int x,y,z,m,n,p;
      ierr=DMDAGetCorners(problem.grid.da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
      ierr=allocateMarkers( (PetscInt)(problem.options.maxMarkFactor*m*problem.options.NMX*n*problem.options.NMY), &problem.markerset, &problem.options);CHKERRQ(ierr);
      problem.markerset.nMark = n*m*problem.options.NMX*problem.options.NMY;
    }

    resetMarkers( &problem.markerset, &problem.options);
    resetNodalFields( &problem.nodal_fields, &problem.grid, &problem.options);

    ierr=distributeMarkersUniformInDomain( &problem.markerset, &problem.options, &problem.grid);CHKERRQ(ierr);
    
    PetscInt m;
    /* INITIAL CONDITIONS */
    for(m=0;m<problem.markerset.nMark;m++){
      /* initial condition for ridges */
      Marker *markers = problem.markerset.markers;
      if(1){
	const PetscScalar initialDamageWidth = 100.0/2;
	const PetscScalar initialDamage = 0.5;

	if( markers[m].Y <= problem.grid.LY/5.0){/* this is the fraction of the domain occupied by sticky air*/
	  /* weak layer */
	  markers[m].Mat = 1;
	  markers[m].T = problem.options.thermalBCTop.value[0];
	}else if( fabs(markers[m].X - problem.grid.LX/2.0) < 0.01*problem.grid.LX && markers[m].Y < 0.9*problem.grid.LY){/* inside perturbed region*/
	  markers[m].Mat = 0;
	  PetscScalar TTop = problem.options.thermalBCTop.value[0];
	  PetscScalar TBtm = problem.options.thermalBCBottom.value[0];
	  
	  markers[m].T = TTop+(TBtm-TTop)/problem.grid.LY*(markers[m].Y-0.1*problem.grid.LY) + problem.options.Tperturb;
	  if( fabs(markers[m].X - problem.grid.LX/2.0) < initialDamageWidth ){
	    markers[m].D = initialDamage;
	  }	
	} else{
	  markers[m].Mat = 0;
	  PetscScalar TTop = problem.options.thermalBCTop.value[0];
	  PetscScalar TBtm = problem.options.thermalBCBottom.value[0];
	  markers[m].T = TTop+(TBtm-TTop)/problem.grid.LY*(markers[m].Y-0.1*problem.grid.LY); /* conductive profile */
	  //markers[m].T = problem.options.TBtm;
	  if( fabs(markers[m].X - problem.grid.LX/2.0) < initialDamageWidth){
	    markers[m].D = initialDamage;
	  }
	}

      }else if(0){/* initial condition for testing out-of-plane couette flow*/
	markers[m].Mat = 0;
	markers[m].T = (markers[m].X/problem.grid.LX)*(problem.options.thermalBCRight.value[0]-problem.options.thermalBCLeft.value[0]) + problem.options.thermalBCLeft.value[0];
      }
      markers[m].rhodot = 0.0;
      markers[m].rho = getDamageDensity( &markers[m], &problem.materials, &problem.options );

      /* calculate density that is self consistent with damage */
      
      markers[m].eta = problem.materials.materialEta[(PetscInt) markers[m].Mat];
      markers[m].mu = problem.materials.materialMu[(PetscInt) markers[m].Mat];
    } 

    /* calculate cell to which each marker belongs*/
    findMarkerCells( &problem.markerset, &problem.grid);

    /* project initial condition (thermal) from markers to nodes */
    ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials, &problem.options);
    ierr = enforceThermalBCs1( &problem.grid, &problem.options, &problem.nodal_fields);CHKERRQ(ierr);
    saveNodalFields( &problem.nodal_fields, &problem.grid, -iMonte,-6,0.0,0.0,0);
    projectNodalFieldToMarkersS(&problem.nodal_fields,problem.nodal_fields.lastT,&problem.markerset,&(problem.markerset.markers[0].T),&problem.grid);


    /* exchange markers among cpus*/
    PetscLogStagePush(stages[4]);
    ierr = exchangeMarkers( &problem.markerset, &problem.grid, &problem.options);
    PetscLogStagePop();
    findMarkerCells( &problem.markerset, &problem.grid);
    
    /* save markers for debugging*/
    PetscLogStagePush(stages[11]);
    ierr=saveMarkersBinary( &problem.markerset, -iMonte,-5,0.0);
    saveNodalFields( &problem.nodal_fields, &problem.grid, -iMonte,-5,0.0,0.0,0);
    PetscLogStagePop();

    PetscInt iTime;
    PetscScalar dt;
    PetscInt isYielding;
    PetscInt iPlastic;
    
    displacementdt=0.0;/* initially*/
    for(iTime=0;iTime<problem.options.nTime;iTime++){
      if(!rank) printf("Beginning Timestep %d, elapsedTime %e years\n",iTime,elapsedTime/(365.25*24*3600));
      /* get memory use */
      getMemoryUse();

      isYielding = 1;
      iPlastic=0;
      dt=problem.options.dtMax;/* initial timestep*/
      findMarkerCells( &problem.markerset, &problem.grid);
      ierr = exchangeMarkers( &problem.markerset, &problem.grid,&problem.options);
      findMarkerCells( &problem.markerset, &problem.grid);
      /* check marker density, add markers if necessary.*/
      ierr = checkMarkerDensity(&problem.markerset, &problem.grid, &problem.options, r);
      findMarkerCells( &problem.markerset, &problem.grid);
  
      /* project rhodot to nodes*/
      findMarkerCells( &problem.markerset, &problem.grid);
      PetscLogStagePush(stages[5]);
      /* project all marker fields onto nodes */
      ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials,&problem.options);      PetscLogStagePop();

      /* check for overheating (melting). abort if Tmelt exceeded*/
      { /* check for exceedence of maximum temperature */
	PetscInt m;
	PetscScalar maxT = 0.0;
	Marker *markers = problem.markerset.markers;
	for(m=0;m<problem.markerset.nMark;m++){
	  /* make sticky air isothermal and cold */
	  if(markers[m].Mat == 1) markers[m].T = problem.options.thermalBCTop.value[0];
	  if(markers[m].T > maxT) maxT = markers[m].T;
	}
	PetscScalar maxTg;
	MPI_Allreduce( &maxT, &maxTg, 1, MPI_DOUBLE, MPI_MAX,PETSC_COMM_WORLD);
	if( maxTg > TBREAK){
	  if(!rank) printf("Maximum temperature exceeded. Ending Montecarlo step\n");
	  saveNodalFields( &problem.nodal_fields, &problem.grid, iMonte,iTime,dt,elapsedTime,0);
	  ierr=saveMarkersBinary( &problem.markerset, iMonte,iTime,elapsedTime);
	  //saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	  goto nextMonte;
	}
      }
      /* update boundary conditions, which can vary in time */
      updateBoundaryConditions( &problem.options, &boundaryValues, elapsedTime );
      if( iTime < 10){
	displacementdt = INITIALDT;
      }else{
	displacementdt = problem.options.dtMax;/* try for maximum displacement timestep */
      }
      while(isYielding && iPlastic < problem.options.maxNumPlasticity){/*begin plastic iteration loop*/
	/* check for plastic yielding */
	isYielding=0;

	ierr = checkPlasticYielding(&problem.grid, &problem.markerset, &problem.materials, displacementdt, &isYielding, &problem.options, elapsedTime);

	ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials,&problem.options);
	
	if(!rank) printf("iTime=%d,iPlastic = %d\n",iTime,iPlastic);
	if(!isYielding && !rank) printf("no more yielding\n");
	if( iTime < problem.options.plasticDelay && !rank) {isYielding=0; printf("yielding disabled\n");}
	  
	ierr=formVEPSystem( &problem.nodal_fields, &problem.grid, problem.mech_system.lhs,PETSC_NULL, problem.mech_system.rhs,problem.z_system.lhs,problem.z_system.rhs, &Kbond, &Kcont, problem.options.gy, dt, &problem.options, &boundaryValues);

	Vec mechanicalS;/* solution to mechanical problem*/
	ierr = VecDuplicate(problem.mech_system.rhs, &mechanicalS);

	ierr = kspLinearSolve(problem.mech_system.lhs,PETSC_NULL, problem.mech_system.rhs, mechanicalS,"stokes_",PETSC_NULL);

	Vec Sz;/* solution to out-of-plane problem*/
	ierr = VecDuplicate(problem.z_system.rhs, &Sz);
	ierr = VecZeroEntries( Sz );CHKERRQ(ierr);
	ierr = kspLinearSolve(problem.z_system.lhs,PETSC_NULL, problem.z_system.rhs, Sz,"zstokes_",PETSC_NULL);CHKERRQ(ierr);
	/* move solutions from solution vectors to nodalFields*/
	PetscInt fn=0;/* flag for nan in solution*/
	ierr= retrieveSolutions( &problem.grid, &problem.nodal_fields, mechanicalS, Sz, Kcont, &fn);
	ierr = VecDestroy(&mechanicalS); CHKERRQ(ierr);
	ierr = VecDestroy(&Sz); CHKERRQ(ierr);
	if(fn){
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Found nan in solution\n");CHKERRQ(ierr);
	  saveNodalFields( &problem.nodal_fields, &problem.grid, iMonte,iTime,dt,elapsedTime,0);
	  ierr=saveMarkersBinary( &problem.markerset, iMonte,iTime,elapsedTime);
	  if(iMonte < problem.options.nMonte-1){
	    goto nextMonte;
	  }else{
	    goto abort;
	  }
	}/* end deal with nan in solution*/
	ierr=limitDisplacementTimestep( &problem.grid, &problem.nodal_fields, &displacementdt, &problem.options); CHKERRQ(ierr);
	//ierr = VecZeroEntries(nodalHeating);CHKERRQ(ierr);

	ierr=nodalStressStrain(&problem.grid, &problem.nodal_fields,&problem.options, &boundaryValues, displacementdt,problem.nodal_heating, problem.options.gy);CHKERRQ(ierr);
	/* compute marker strain and pressure*/
	ierr=updateMarkerStrainPressure( &problem.grid,&problem.nodal_fields, &problem.markerset, &problem.materials, &problem.options, displacementdt);

	/* update marker effective viscosity to match marker strain reate*/
	iPlastic++;
      }/* end plasticity loop*/
      if(!rank) printf("Done with plasticity loops.\n");
      /* do sub-grid stress diffusion*/
      PetscLogStagePush(stages[7]);
      if(problem.options.subgridStressDiffusivity>0.0)	ierr=subgridStressChanges(&problem.grid, &problem.nodal_fields, &problem.markerset, &problem.materials, displacementdt,problem.options.subgridStressDiffusivity);
      PetscLogStagePop();
      /* update marker stress*/
      ierr= updateMarkerStress( &problem.grid, &problem.nodal_fields,&problem.markerset, &problem.materials);CHKERRQ(ierr);
      /* set thermal problem timestep to displacement dt*/
      dt = displacementdt;
     
      /* Form thermal system (diffusion, variable coefficient, non-uniform grid */
      
      if( !problem.options.shearHeating ){
	ierr=VecZeroEntries(problem.nodal_heating);CHKERRQ(ierr); /* this line disables shear heating */
      }
      if( problem.options.adiabaticHeating ){
	adiabaticHeating( &problem.grid, &problem.markerset, &problem.nodal_fields, &problem.materials, &problem.options);
	ierr=VecAXPY(problem.nodal_heating,1.0,problem.nodal_fields.ha);CHKERRQ(ierr);/* add adiabatic heating to shear heating */
      }
      ierr = VecCopy(problem.nodal_heating,problem.thermal_system.rhs);CHKERRQ(ierr);
      ierr = enforceThermalBCs1( &problem.grid, &problem.options, &problem.nodal_fields);CHKERRQ(ierr);

      PetscLogStagePush( stages[8] );
      ierr = formThermalSystem( &problem.grid, &problem.nodal_fields, problem.thermal_system.lhs, problem.thermal_system.rhs, dt, &problem.options);CHKERRQ(ierr);
      PetscLogStagePop();

      /* make the solution vector */
      Vec thermalS;
      ierr = VecDuplicate( problem.thermal_system.rhs, &thermalS);CHKERRQ(ierr);
      /* Solve the thermal system: */
      PetscLogStagePush( stages[9] );  
      ierr = kspLinearSolve(problem.thermal_system.lhs, PETSC_NULL, problem.thermal_system.rhs, thermalS, "energy_",PETSC_NULL);CHKERRQ(ierr); PetscLogStagePop();
      /* check maximum temperature change */
      Vec dT;
      ierr = VecDuplicate( thermalS, &dT);CHKERRQ(ierr);
      ierr = VecCopy(thermalS, dT);CHKERRQ(ierr);
      ierr = VecAXPY(dT,-1.0,problem.nodal_fields.lastT);CHKERRQ(ierr);
      PetscScalar dTmax;/* max abs change in temperature*/
      ierr = VecNorm(dT,NORM_INFINITY,&dTmax);CHKERRQ(ierr);

      if(dTmax > problem.options.maxTempChange){
	if(!rank) printf("[%d] temperature change too large, limiting timestep from %e to %e\n",rank,dt,dt*problem.options.maxTempChange/dTmax);
	dt = dt*problem.options.maxTempChange/dTmax;
	/* repeat temperature solution */
	PetscLogStagePush( stages[8] );
	ierr = VecCopy(problem.nodal_heating,problem.thermal_system.rhs);CHKERRQ(ierr);
	ierr = enforceThermalBCs1( &problem.grid, &problem.options, &problem.nodal_fields);CHKERRQ(ierr);
	ierr = formThermalSystem( &problem.grid, &problem.nodal_fields, problem.thermal_system.lhs, problem.thermal_system.rhs, dt, &problem.options);CHKERRQ(ierr); PetscLogStagePop();
	PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(problem.thermal_system.lhs, PETSC_NULL, problem.thermal_system.rhs, thermalS, "energy_",PETSC_NULL);CHKERRQ(ierr); PetscLogStagePop();
	/* end repeat temperature solution */
      }
      ierr = VecDestroy(&dT);CHKERRQ(ierr);

      /* Apply sub-grid temperature corrections and project corrected temperature onto the markers (all integrated into one routine*/
      ierr=subgridTemperatureChanges(thermalS, &problem.grid, &problem.nodal_fields, &problem.markerset, &problem.materials,dt, &problem.options); 
      ierr = VecDestroy(&thermalS); /* free the thermal solution vector */
                  
      /* advect markers */
      ierr=advectMarkersRK(&problem.markerset,&problem.nodal_fields,&problem.grid,&problem.options,&boundaryValues,dt);
      /* update marker density and damage and viscosity using stored rates*/
      {
	PetscInt m;
	Marker *markers = problem.markerset.markers;
	for(m=0;m<problem.markerset.nMark;m++){	
	  if(markers[m].cellX != -1){
	    if(markers[m].rho > problem.materials.materialRho[(PetscInt) markers[m].Mat]) markers[m].rho = problem.materials.materialRho[(PetscInt) markers[m].Mat];
	    if(problem.materials.hasDilation[ markers[m].Mat ] == 0){
	      markers[m].rho = problem.materials.materialRho[(PetscInt) markers[m].Mat]*(1-problem.materials.materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
	    }
	    markers[m].rho += dt*markers[m].rhodot;
	    /* update damage variable */
	    updateMarkerDamage( &markers[m], &problem.materials, &problem.options, dt);
	  }
	}
      }


      /* save solution */
      if(!(iTime % (problem.options.saveInterval))) {
	ierr=saveNodalFields( &problem.nodal_fields, &problem.grid,iMonte,iTime, dt, elapsedTime,0);
	ierr=saveMarkersBinary( &problem.markerset, iMonte,iTime,elapsedTime);
      }
      if(!(iTime % problem.options.saveInterval)) {
	//saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime,elapsedTime);
      }
      
      elapsedTime+=dt;
      if(elapsedTime > problem.options.totalTime || iTime==problem.options.nTime-1){
	PetscLogStagePush(stages[11]);
	//saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	PetscLogStagePop();
	ierr=saveMarkersBinary( &problem.markerset,iMonte,iTime,elapsedTime);
	goto nextMonte;
      }

      if(dt < 1.0) {	
	if(!rank) printf("timestep less than 1s, aborting run\n");fflush(stdout);
	PetscLogStagePush(stages[11]);
	ierr=saveMarkersBinary( &problem.markerset, iMonte,iTime,elapsedTime);
	saveNodalFields( &problem.nodal_fields, &problem.grid, iMonte,iTime,dt,elapsedTime,0);
	PetscLogStagePop();
	goto nextMonte;
      }
      
    }/* end time stepping*/
  nextMonte:;


  }/* end montecarlo loop */
 abort:;

  //ierr = VecDestroy(problem.thermal_system.rhs);CHKERRQ(ierr);
  //ierr = VecDestroy(problem.z_system.rhs);CHKERRQ(ierr);
  //ierr = VecDestroy(RHS);CHKERRQ(ierr);
  ierr = VecDestroy(&problem.nodal_heating);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(problem.grid.vda,&problem.mech_system.rhs);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(problem.grid.da,&problem.z_system.rhs);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(problem.grid.da,&problem.thermal_system.rhs);CHKERRQ(ierr);

  ierr = MatDestroy( &problem.thermal_system.lhs );CHKERRQ(ierr);
  ierr = MatDestroy( &problem.z_system.lhs );CHKERRQ(ierr);
  ierr = MatDestroy( &problem.mech_system.lhs );CHKERRQ(ierr);

  ierr = destroyGrid(&problem.grid);CHKERRQ(ierr);
  ierr = destroyNodalFields(&problem.nodal_fields, &problem.grid);CHKERRQ(ierr);
  ierr=  destroyMarkers( &problem.markerset,&problem.options );CHKERRQ(ierr);
  //  ierr = PetscRandomDestroy( r);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
 
