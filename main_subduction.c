static const char help[] = "Help statement goes here\n";

#include<stdio.h>
#include<stdlib.h>
//#include<malloc.h>
#include "petscksp.h"
#include "petsctime.h"
#include "mpi.h"
#include "fdcode.h"
#include "kspLinearSolve.h"
#include "markers.h"
#include "gridGenerator.h"
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
#include "viscosity.h"
#include "restart.h"
#include "benchmarkInitialConditions.h"
#include "residual.h"
#include "version.h"
#include "profile.h"
#include "markerProjection.h"
#include "initialPressureGuess.h"
#include "pressureNullSpace.h"
#include "post.h"

#define TBREAK 9999 /* stop if temperature exceeds this number */

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **args){
  PetscMPIInt rank,size;
  PetscErrorCode ierr;

  /* Problem context */
  Problem problem;
  problem.mech_system.pc = PETSC_NULL;/* nullify preconditioner */

  /* material properties, boundary values and options */
  BoundaryValues boundaryValues; /* used to store time-varying boundary condition information*/
  PetscRandom r;

  /* BEGIN PROGRAM */

  PetscInitialize(&argc,&args,NULL,help);  //INITIALIZE PETSC
  MPI_Comm_size(PETSC_COMM_WORLD,&size);  //Get MPI rank
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /* Print welcome message */
  if(!rank) printf("Version Information : %s\n",MARKERCODE_HG_VERSION);
  initializeLogging();
  /* read options from input file and distribute to all nodes */
  printf("argc = %d\n",argc);
  printf("%s\n",args[1]);
  /* Get input filename */
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
  if(!rank) ierr = saveRunInfo( &problem.options, &problem.materials, 0);

  /* initialize material properties*/
  //setMaterialProperties( &materials);

  ierr=  PetscRandomCreate(PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  ierr = PetscRandomSetType(r,PETSCRAND48);CHKERRQ(ierr);
  {
    /* seed random number generator with cpu time */
    PetscLogDouble time;
    unsigned long seed;
    if(!rank) ierr = PetscTime(&time); CHKERRQ(ierr);
    ierr = MPI_Bcast( &time, sizeof(PetscLogDouble),MPI_BYTE , 0, PETSC_COMM_WORLD);
    seed = (unsigned long) time;
    seed = 1;
    printf("seeding random number generator with %ld\n",seed);fflush(stdout);
    ierr = PetscRandomSetSeed(r,seed);CHKERRQ(ierr);
    ierr = PetscRandomSeed(r);CHKERRQ(ierr);/* seed the generator*/
  }
  /* note that all processes know the locations of all gridlines and all cell centers */
  if( problem.options.gridRatio == 1.0 ){
    //ierr = initializeRegularGrid( &problem.grid,  &problem.options);CHKERRQ(ierr);
    ierr=initializeSubductionGrid( &problem.grid, &problem.options ); CHKERRQ(ierr);
  }else{
    //    ierr = initializeIrregularGridConstantInnerOuter( &problem.grid, problem.options.LX, problem.options.LY, problem.options.NX, problem.options.NY, &problem.options);CHKERRQ(ierr);
    initializeIrregularGridConstantInnerOuter( &problem.grid, &problem.options );CHKERRQ(ierr);
  }
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
    //printf("[%d], %d %d %d %d %d %d\n",rank,x,y,z,m,n,p);CHKERRQ(ierr);
    ierr=allocateMarkers( (PetscInt)(problem.options.maxMarkFactor*m*problem.options.NMX*n*problem.options.NMY), &problem.markerset, &problem.options);CHKERRQ(ierr); 
    problem.markerset.nMark = n*m*problem.options.NMX*problem.options.NMY;
  }
  /* set up local to global mapping*/
  ISLocalToGlobalMapping ltogm;
  ierr=DMGetLocalToGlobalMapping(problem.grid.da,&ltogm);CHKERRQ(ierr);
  
  /* initialize LHS,LHSz matrix and RHS,RHSz matrix for mechanical problem */
  ierr=DMCreateMatrix(problem.grid.vda,&problem.mech_system.lhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem.mech_system.lhs,"CmechLHS");
  ierr = createPressureNullSpace( &problem.grid, &problem.mech_system.ns );CHKERRQ(ierr); 

  {
    PetscInt m,n; 
    ierr=MatGetSize(problem.mech_system.lhs,&m,&n);CHKERRQ(ierr); 
    if(!rank) printf("LHS is %d by %d, ndof=%d\n",m,n,ndof); 
    ierr=MatGetOwnershipRange(problem.mech_system.lhs,&m,&n);CHKERRQ(ierr);
    printf("[%d] has rows %d-%d\n",rank,m,n);CHKERRQ(ierr);
  }

  ierr=DMGetGlobalVector(problem.grid.vda,&problem.mech_system.rhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem.mech_system.rhs,"CmechRHS"); 

  /* initialize LHS matrix, RHS Vec and nodalHeating Vec for thermal problem*/
  ierr=DMGetGlobalVector(problem.grid.da,&problem.thermal_system.rhs);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem.thermal_system.rhs,"CthermalRHS");
  ierr=DMCreateMatrix(problem.grid.da,&problem.thermal_system.lhs);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem.thermal_system.lhs,"CthermalLHS");
  ierr=VecDuplicate(problem.thermal_system.rhs,&problem.nodal_heating);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem.nodal_heating,"nodalHeating");
  {
    PetscBool dop = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-stokes_ksp_build_pc_diag",&dop,0);
    if( dop ){
      printf("[%d] Allocating Preconditioner\n",rank);
      /*       ierr=MatDuplicate( LHS, MAT_DO_NOT_COPY_VALUES ,&P ); CHKERRQ(ierr); */
      ierr = DMCreateMatrix(problem.grid.vda, &problem.mech_system.pc); CHKERRQ(ierr);
      ierr = PetscObjectSetName( (PetscObject) problem.mech_system.pc, "MechanicalPreconditioner");
    }
  }

  ierr = VecDuplicate(problem.mech_system.rhs, &problem.mech_system.solution);  
  ierr = VecZeroEntries(problem.mech_system.solution);CHKERRQ(ierr);

  PetscInt iMonte=0;
  
    PetscScalar displacementdt;/* this is the displacement timestep*/
    PetscScalar elapsedTime=0.0;
    printf("beginning MonteCarlo run %d\n",iMonte);
    /* randomize parameters*/

    /* reset (zero out) all marker properties*/
    resetMarkers( &problem.markerset, &problem.options);
    resetNodalFields( &problem.nodal_fields, &problem.grid, &problem.options);

    ierr=distributeMarkersUniformInDomain( &problem.markerset, &problem.options, &problem.grid); CHKERRQ(ierr);

    PetscInt m;
    /* INITIAL CONDITIONS */
    ierr = initialConditionsVanKeken( &problem.markerset, &problem.options, &problem.materials, &problem.grid);CHKERRQ(ierr);
   
    /* check to see if we are restarting from a saved state */
    PetscInt iTime0=0;/* initial timestep */
    if( problem.options.restartStep ){
      /* destroy current markers */
      destroyMarkers( &problem.markerset ,&problem.options);
      restartFromMarkers( &problem.markerset, &problem.grid, &problem.materials, &problem.options, iMonte, problem.options.restartStep, &elapsedTime);
      //      markers = &(problem.markerset.markers[0]);
      iTime0 = problem.options.restartStep + 1;
    }
    
    /* calculate cell to which each marker belongs*/
    findMarkerCells( &problem.markerset, &problem.grid);
    /* project initial condition to nodes, apply bcs, push back to markers*/
    ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials,&problem.options);

    if( problem.options.restartStep == 0 ){/* if we are not restarting, project markers to nodes and back */
      ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials,&problem.options);
      ierr = enforceThermalBCs1( &problem.grid, &problem.options, &problem.nodal_fields);CHKERRQ(ierr);
      projectNodalFieldToMarkersS(&problem.nodal_fields,problem.nodal_fields.lastT,&problem.markerset,&(problem.markerset.markers[0].T),&problem.grid);
    }

    /* exchange markers among cpus*/
    ierr = exchangeMarkers( &problem.markerset, &problem.grid, &problem.options);
 
    findMarkerCells( &problem.markerset, &problem.grid);
    
    saveNodalFields( &problem.nodal_fields, &problem.grid, iMonte, -5,0.0,0.0,0);    
    /* save markers for debugging*/
    ierr=saveMarkersBinary( &problem.markerset, -5,-5,0.0);

    /* make an initial guess at the pressure solution to help iterative solver */
    ierr = initialPressureGuess( &problem.grid, &problem.nodal_fields, &problem.options, problem.mech_system.solution ); CHKERRQ(ierr);

    PetscInt iTime;
    PetscScalar dt;
    PetscInt isYielding;
    PetscInt iPlastic;
    
    displacementdt=problem.options.dtMax;/* initially*/

    Vec steadySolution;
    ierr = VecDuplicate( problem.mech_system.solution, &steadySolution); CHKERRQ(ierr);    

    for(iTime=iTime0;iTime<problem.options.nTime;iTime++){
      {
	PetscLogDouble memuse;
	PetscScalar memuses;
	PetscScalar memuset;
	PetscMemoryGetCurrentUsage( &memuse );
	//PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] memory use %e \n",rank,memuse);
	PetscSynchronizedFlush(PETSC_COMM_WORLD,stdout);
	memuses = (PetscScalar) memuse;
	MPI_Allreduce( &memuses, &memuset, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	//if(!rank) PetscPrintf(PETSC_COMM_SELF,"Global memory use %e\n",memuset); 
	
      }
      isYielding = 1;
      iPlastic=0;
      dt=problem.options.dtMax;/* try to take a large initial timestep*/
      findMarkerCells( &problem.markerset, &problem.grid);

      ierr = exchangeMarkers( &problem.markerset, &problem.grid,&problem.options);

      findMarkerCells( &problem.markerset, &problem.grid);
      /* check marker density, add markers if necessary.*/

      ierr = checkMarkerDensity(&problem, &problem.markerset, &problem.grid, &problem.options, r);

      findMarkerCells( &problem.markerset, &problem.grid);
      /* calculate rhodot for markers*/
      //PetscScalar Vmax = 1.0/problem.options.maxPorosity;
      {
	Marker *markers = problem.markerset.markers;
	for(m=0;m<problem.markerset.nMark;m++){
	  if(markers[m].cellX != -1 ){
	    //markers[m].rho = problem.materials.materialRho[(PetscInt) markers[m].Mat]*(1-problem.materials.materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
	    markers[m].rho = -problem.materials.materialRho[(PetscInt) markers[m].Mat]*(problem.materials.materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
	    markers[m].rhodot = 0.0;
	  }       
	}/* end loop over markers*/
      }
      //updateViscosity( &problem.markerset, &problem.options, &problem.materials );
      
      /* project all marker fields onto nodes */
      ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials,&problem.options);

      /* check for overheating (melting). abort if Tmelt exceeded*/
      { /* check for exceedence of maximum temperature */
	PetscScalar maxT = 0.0;
	Marker *markers = problem.markerset.markers;
	for(m=0;m<problem.markerset.nMark;m++){
	  if(markers[m].cellX != -1){
	    if(markers[m].T > maxT) maxT = markers[m].T;
	    if(markers[m].T > TBREAK){
	      printf("Maximum temperature exceeded. Ending Montecarlo step. Marker position (X,Y) = %e, %e\n", markers[m].X, markers[m].Y);
	    }
	  }
	}
	PetscScalar maxTg;
	MPI_Allreduce( &maxT, &maxTg, 1, MPI_DOUBLE, MPI_MAX,PETSC_COMM_WORLD);
	if( maxTg > TBREAK){
	  printf("Maximum temperature exceeded. Ending Montecarlo step. Marker position (X,Y) = %e, %e\n", markers[m].X, markers[m].Y);
	  //saveGriddedMarkersBinary( &markers, &problem.grid, 5*problem.grid.NX, 5*problem.grid.NY,iMonte,iTime);
	  goto nextMonte;
	}
      }
      /* update boundary conditions, which can vary in time */
      updateBoundaryConditions( &problem.options, &boundaryValues, elapsedTime );
      displacementdt = dt;
      PetscInt converged=0;
      while( (iTime==0 && iPlastic < 1) || ((isYielding || !converged) && iPlastic < problem.options.maxNumPlasticity)){/*begin plastic iteration loop*/
	/* check for plastic yielding */
	//ierr=saveMarkersBinary( &problem.markerset, -6,-iPlastic,0.0);
	isYielding = 0;
	converged = 0;
	ierr = checkPlasticYielding(&problem.grid, &problem.markerset, &problem.materials, displacementdt, &isYielding, &problem.options, elapsedTime);
	limitViscosity( &problem.markerset, &problem.options, &problem.materials);
	
	/* project marker quantities onto nodes: kThermal, Cp, rho, T*/
	ierr = projectMarkersNodesAll2(&problem.markerset, &problem.grid, &problem.nodal_fields, &problem.materials, &problem.options);

	/* form mechanical problem LHS*/
	/* zero out LHS, RHS*/
	if( !problem.options.staticVelocity || iTime == iTime0 ){	
	  ierr = VecZeroEntries( problem.mech_system.rhs);CHKERRQ(ierr);	
	  ierr=formVEPSystem( &problem.nodal_fields, &problem.grid, problem.mech_system.lhs, problem.mech_system.pc, problem.mech_system.rhs, &Kbond, &Kcont, problem.options.gy, displacementdt, &problem.options, &boundaryValues);
	  
	  /* scale the pressure guess */
	  ierr = VecStrideScale( problem.mech_system.solution, DOF_P, 1.0/Kcont );CHKERRQ(ierr);	  
	  ierr = kspLinearSolve(problem.mech_system.lhs, problem.mech_system.pc,problem.mech_system.rhs, problem.mech_system.solution,"stokes_",problem.mech_system.ns);	  
	  ierr = VecCopy( problem.mech_system.solution, steadySolution ); CHKERRQ(ierr);
	}else{
	  ierr = VecCopy( steadySolution, problem.mech_system.solution ); CHKERRQ(ierr);
	}

	PetscInt fn=0;/* flag for nan in solution*/
	ierr= retrieveSolutions( &problem.grid, &problem.nodal_fields, problem.mech_system.solution, Kcont, &fn);


	if(fn) {
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Found nan in solution\n");CHKERRQ(ierr);

	  saveNodalFields( &problem.nodal_fields, &problem.grid, iMonte,iTime,displacementdt,elapsedTime,0);

	  /* 	  saveGriddedMarkersBinary( &markers, &problem.grid, 5*problem.grid.NX, 5*problem.grid.NY,iMonte,iTime); */
	  ierr=saveMarkersBinary( &problem.markerset, iMonte,-iTime,elapsedTime); 

	  goto abort;

	}/* end deal with nan in solution*/

	ierr=nodalStressStrain(&problem.grid, &problem.nodal_fields,&problem.options,&boundaryValues, displacementdt, problem.nodal_heating, problem.options.gy);CHKERRQ(ierr);
	/* save nodalFields for debugging*/

	/*ierr=saveNodalFieldsASCIIMatlab( &problem.nodal_fields, &grid,iMonte,iTime, dt, elapsedTime);*/ /* for debugging*/

	/* compute marker strain and pressure*/
	ierr=updateMarkerStrainPressure( &problem.grid,&problem.nodal_fields, &problem.markerset, &problem.materials, &problem.options, displacementdt);
	/* determine displacement timestep*/
	
	/* calculate global strain rate residual */
	updateGlobalStrainRateResidual( &problem.grid, &problem.markerset, iPlastic );
	converged = checkConvergence( 1e-16, 1e-1, iPlastic);
	{
	  PetscScalar lastdt = displacementdt;
	  ierr=limitDisplacementTimestep( &problem.grid, &problem.nodal_fields, &displacementdt, &problem.options);CHKERRQ(ierr);
	  if( displacementdt != lastdt){
	    converged = 0;
	  }
	}	
	if(!rank) printf("Timestep %d iteration %d, displacementdt = %e, yielding %d, residual %e, converged = %d\n",iTime,iPlastic,displacementdt,isYielding,getGlobalStrainRateResidual( iPlastic ),converged);
	
	/* update marker effective viscosity to match marker strain reate*/
	iPlastic++;
      }/* end plasticity loop*/
      if(!rank) printf("Done with plasticity loops.\n");
      /* set thermal problem timestep to displacement dt*/
      dt = displacementdt;

      /* do sub-grid stress diffusion*/

      if(problem.options.subgridStressDiffusivity>0.0)	ierr=subgridStressChanges(&problem.grid, &problem.nodal_fields, &problem.markerset, &problem.materials, displacementdt,problem.options.subgridStressDiffusivity);

      /* update marker stress*/
      ierr= updateMarkerStress( &problem.grid, &problem.nodal_fields,&problem.markerset, &problem.materials);CHKERRQ(ierr);
      /* 1. project velocity field onto markers*/

      /* calculate adiabatic heating */
      if( !problem.options.shearHeating ){
	ierr=VecZeroEntries(problem.nodal_heating);CHKERRQ(ierr); /* this line disables shear heating */
      }
      if( problem.options.adiabaticHeating) {
	adiabaticHeating( &problem.grid, &problem.markerset, &problem.nodal_fields, &problem.materials, &problem.options);
	ierr=VecAXPY(problem.nodal_heating,1.0,problem.nodal_fields.ha);CHKERRQ(ierr);/* add adiabatic heating to shear heating */
      }

      ierr=VecCopy(problem.nodal_heating,problem.thermal_system.rhs);CHKERRQ(ierr);
      /* Form thermal system (diffusion, variable coefficient, non-uniform grid */
      ierr = enforceThermalBCs1( &problem.grid, &problem.options, &problem.nodal_fields);CHKERRQ(ierr);

      ierr = formThermalSystem(&problem, &problem.grid, &problem.nodal_fields, problem.thermal_system.lhs, problem.thermal_system.rhs, dt, &problem.options, 0);CHKERRQ(ierr);

      /* make the solution vector */
      Vec thermalS;
      ierr = VecDuplicate( problem.thermal_system.rhs, &thermalS);CHKERRQ(ierr);
      ierr = VecCopy( problem.nodal_fields.lastT, thermalS); CHKERRQ(ierr);
      /* Solve the thermal system: */
      ierr = kspLinearSolve(problem.thermal_system.lhs, PETSC_NULL, problem.thermal_system.rhs, thermalS,"energy_",PETSC_NULL);CHKERRQ(ierr); 
      /* check maximum temperature change */
      Vec dT;
      ierr = VecDuplicate( thermalS, &dT);CHKERRQ(ierr);
      ierr = VecCopy(thermalS, dT);CHKERRQ(ierr);
      ierr = VecAXPY(dT,-1.0,problem.nodal_fields.lastT);CHKERRQ(ierr);
      PetscScalar dTmax;/* max abs change in temperature*/
      ierr = VecNorm(dT,NORM_INFINITY,&dTmax);CHKERRQ(ierr);
      if(!rank) printf("Maximum temperature change %le\n",dTmax);
      if(dTmax > problem.options.maxTempChange){
	printf("[%d] temperature change too large, limiting timestep from %e to %e\n",rank,dt,dt*problem.options.maxTempChange/dTmax);
	dt = dt*problem.options.maxTempChange/dTmax;
	/* repeat temperature solution */

	ierr = VecCopy(problem.nodal_heating,problem.thermal_system.rhs);CHKERRQ(ierr);
	ierr = enforceThermalBCs1( &problem.grid, &problem.options, &problem.nodal_fields);CHKERRQ(ierr);
	ierr = formThermalSystem(&problem, &problem.grid, &problem.nodal_fields, problem.thermal_system.lhs, problem.thermal_system.rhs, dt, &problem.options, 0);CHKERRQ(ierr);
	ierr = kspLinearSolve(problem.thermal_system.lhs, PETSC_NULL, problem.thermal_system.rhs, thermalS,"energy_",PETSC_NULL);CHKERRQ(ierr);
	/* end repeat temperature solution */
      }
      ierr = VecDestroy(&dT);CHKERRQ(ierr);
      
      /* Apply sub-grid temperature corrections and project corrected temperature onto the markers (all integrated into one routine*/

      ierr=subgridTemperatureChanges(thermalS, &problem.grid, &problem.nodal_fields, &problem.markerset, &problem.materials,dt, &problem.options);


      PetscScalar Nu;
      PetscScalar vrms;
      ierr = nusseltNumber( &(problem.grid), &(problem.options), thermalS, &Nu );CHKERRQ(ierr);
      ierr = rmsVelocity( &problem.grid, &problem.nodal_fields, &vrms );CHKERRQ(ierr);
      if(!rank) printf("Nusselt Number %e, vrms=%e\n",Nu,vrms);
      
      ierr = VecDestroy(&thermalS); /* free the thermal solution vector */
               
      /* save solution */
      if(!(iTime % (problem.options.saveInterval))){
	ierr=saveNodalFields( &problem.nodal_fields, &problem.grid,iMonte,iTime, dt, elapsedTime,0);
      }
      if( !(iTime % (problem.options.markerSaveInterval))) {
	ierr=saveMarkersBinary( &problem.markerset, iMonte,iTime,elapsedTime);
      }
    
      /* advect markers */
      ierr=advectMarkersRK(&problem.markerset,&problem.nodal_fields, &problem.grid,&problem.options,&boundaryValues,dt);

      {
	PetscInt mm;
	PetscScalar sxxmax=0.0;
	Marker *markers = problem.markerset.markers;
	for(mm=0;mm<problem.markerset.nMark;mm++){
	  if( problem.markerset.markers[mm].cellX != -1 ){
	    if( fabs(markers[mm].s.T11) > sxxmax) sxxmax = fabs(markers[mm].s.T11);
	  }
	}
	printf("[%d] sxx max = %e\n",rank,sxxmax);

      }




      elapsedTime+=dt;
      if(elapsedTime > problem.options.totalTime || iTime==problem.options.nTime-1){

	//saveGriddedMarkersBinary( &markers, &problem.grid, 5*problem.grid.NX, 5*problem.grid.NY,iMonte,iTime);

	ierr=saveMarkersBinary( &problem.markerset,iMonte,iTime,elapsedTime);
	goto nextMonte;
      }
      
      
    }/* end time stepping*/
 nextMonte:;            
 abort:;
  
  //ierr = VecDestroy(&problem.thermal_system.rhs);CHKERRQ(ierr);
  //ierr = VecDestroy(&RHSz);CHKERRQ(ierr);
  //ierr = VecDestroy(&RHS);CHKERRQ(ierr);

  ierr = VecDestroy(&problem.nodal_heating);CHKERRQ(ierr);
  ierr=DMRestoreGlobalVector(problem.grid.vda,&problem.mech_system.rhs);CHKERRQ(ierr);

  ierr=DMRestoreGlobalVector(problem.grid.da,&problem.thermal_system.rhs);CHKERRQ(ierr);
  
  ierr = MatDestroy(&problem.thermal_system.lhs );CHKERRQ(ierr);
  // ierr = MatDestroy(&problem.z_system.lhs );CHKERRQ(ierr);
  ierr = destroyPressureNullSpace( &problem.mech_system.ns ); CHKERRQ(ierr);
  ierr = MatDestroy(&problem.mech_system.lhs );CHKERRQ(ierr);
  if( problem.mech_system.pc != PETSC_NULL){
    ierr = MatDestroy(&problem.mech_system.pc);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&problem.mech_system.solution); CHKERRQ(ierr);
  ierr = destroyGrid(&problem.grid);CHKERRQ(ierr);
  ierr = destroyNodalFields(&problem.nodal_fields, &problem.grid);CHKERRQ(ierr);
  ierr=  destroyMarkers( &problem.markerset,&problem.options );CHKERRQ(ierr);

  ierr = PetscRandomDestroy(& r );CHKERRQ(ierr);
  finalizeLogging();
  ierr = PetscFinalize();
}
