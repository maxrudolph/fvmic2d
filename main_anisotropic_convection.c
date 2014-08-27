
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
#include "anisotropic_system.h"
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
#include "texture.h"
#include "adiabaticHeating.h"
#include "restart.h"
#include "benchmarkInitialConditions.h"
#include "version.h"
#include "profile.h"

//Include file for PETSc linear solve commands

//extern struct ELEMENT linearHexahedralElementCreate();
//extern double evaluateNodalBasisFunction();

//#define TBREAK 273.15 
/* stop if temperature exceeds this number */
#define TBREAK 10000

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
  Marker *markers;

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
  printf("argc = %d\n",argc);
  printf("%s\n",args[1]);
  /* check to see if first argument looks like a filename*/
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
    //    seed = 1;
    printf("seeding random number generator with %ld\n",seed);fflush(stdout);
    ierr = PetscRandomSetSeed(r,seed);CHKERRQ(ierr);
    ierr = PetscRandomSeed(r);CHKERRQ(ierr);/* seed the generator*/
  }
  /* initialize the grid. This version does not require a regular grid*/
  /* note that all processes know the locations of all gridlines and all cell centers */
  //ierr = initializeIrregularGridConstantIncrease( &grid, options.LX, options.LY, options.NX, options.NY,1.02,1.02);CHKERRQ(ierr);
  //  ierr = initializeIrregularGridConstantIncreaseFixedRatio( &grid, options.LX, options.LY, options.NX, options.NY, options.gridRatio);CHKERRQ(ierr);
  if( options.gridRatio == 1.0 ){
    ierr = initializeRegularGrid( &grid, options.LX, options.LY, options.NX, options.NY, &options);CHKERRQ(ierr);
  }else{
    ierr = initializeIrregularGridConstantInnerOuter( &grid, options.LX, options.LY, options.NX, options.NY, &options);CHKERRQ(ierr);
  }
  //  ierr = initializeRegularGrid(  &grid, options.LX, options.LY, options.NX, options.NY);CHKERRQ(ierr);

  ierr = saveGrid( &grid );
  /* initialize scaling parameters*/
  PetscInt ndof = grid.NX*grid.NY*3;  
  PetscScalar Kbond = 0.0;
  PetscScalar Kcont = 0.0;
  /* allocate the nodal fields. The nodal fields are all global PETSc Vecs and nodalFields->da is the associated DA*/
  ierr = initializeNodalFields( &nodalFields, &grid,&options);CHKERRQ(ierr);
  {
    /*allocate the markers. First, calculate the right number of markers for each processor.*/
    int x,y,z,m,n,p;
    ierr=DMDAGetCorners(grid.da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
    //printf("[%d], %d %d %d %d %d %d\n",rank,x,y,z,m,n,p);CHKERRQ(ierr);
    ierr=allocateMarkers( (PetscInt)(options.maxMarkFactor*m*options.NMX*n*options.NMY), &(markerset), &options);CHKERRQ(ierr); 
    markers = &(markerset.markers[0]);
    markerset.nMark = n*m*options.NMX*options.NMY;
  }
  /* set up local to global mapping*/
  ISLocalToGlobalMapping ltogm;
  ierr=DMGetLocalToGlobalMapping(grid.da, &ltogm);CHKERRQ(ierr);
  
  /* initialize LHS,LHSz matrix and RHS,RHSz matrix for mechanical problem */
  ierr=DMCreateMatrix(grid.vda,MATAIJ,&LHS);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) LHS,"CmechLHS");
  ierr=DMCreateMatrix(grid.da,MATAIJ,&LHSz);CHKERRQ(ierr); 
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
  ierr=DMCreateMatrix(grid.da,MATAIJ,&thermalLHS);CHKERRQ(ierr);
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
    /*     options.grainSize = 1e-3; */ /* grain size in meters*/
    if(!rank) ierr = saveRunInfo( &options, &materials, iMonte);
    /* reset (zero out) all marker properties*/
    resetMarkers( &markerset, &options);
    resetNodalFields( &nodalFields, &grid, &options);
    ierr=distributeMarkersUniformInDomain( &markerset, &options, &grid);CHKERRQ(ierr);
    
    PetscInt m;
    /* INITIAL CONDITIONS */

    ierr = initialConditionsConvectionBarr( &markerset, &options, &materials, &grid);
    
    /* check to see if we are restarting from a saved state */
    PetscInt iTime0=0;/* initial timestep */
    if( options.restartStep != 0 ){
      /* destroy current markers */
      destroyMarkers( &markerset ,&options);
      if( options.restartStep > 0){
	restartFromMarkers( &markerset, &grid, &materials, &options, iMonte, options.restartStep, &elapsedTime);
      }else{
	options.restartStep = -options.restartStep;
	restartFromMarkersNoTexture( &markerset, &grid, &materials, &options, iMonte, options.restartStep, &elapsedTime);

      }

      markers = &(markerset.markers[0]);
      iTime0 = options.restartStep + 1;
      /*       ierr=saveMarkersBinary( &markerset, -4,-4,0.0); */
    }
    
    /* calculate cell to which each marker belongs*/
    findMarkerCells( &markerset, &grid);
    ierr = exchangeMarkers( &markerset, &grid, &options);

    /* project initial condition (thermal) from markers to nodes */
    ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials, &options);
    /* calculate cell to which each marker belongs*/
    findMarkerCells( &markerset, &grid);
    /* project initial condition to nodes, apply bcs, push back to markers*/
    if( options.restartStep == 0 ){/* if we are not restarting, project markers to nodes and back */
      ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials,&options);
      ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);
      projectNodalFieldToMarkersS(&nodalFields,nodalFields.lastT,&markerset,&(markerset.markers[0].T),&grid);
    }
    /* exchange markers among cpus*/
    PetscLogStagePush(stages[4]);
    ierr = exchangeMarkers( &markerset, &grid, &options);
    PetscLogStagePop();
    findMarkerCells( &markerset, &grid);

    {
      /* debugging texture reading */
      PetscScalar thetamean=0.0;      
      PetscScalar thetamax = 0.0;
      PetscScalar thetamin = 0.0;
      PetscInt nn  =0;
      PetscScalar NT3 = (PetscScalar) NT*NT*NT;
      for( m=0;m<markerset.nMark;m++){
	if( markerset.markers[m].cellX != -1 ){
	  PetscInt i,j,k;
	  for(i=0;i<NT;i++){
	    for(j=0;j<NT;j++){
	      for(k=0;k<NT;k++){
		thetamean += markerset.markers[m].texture.ctheta[i][j][k]/NT3;
		if( markerset.markers[m].texture.ctheta[i][j][k] > thetamax  ) thetamax =  markerset.markers[m].texture.ctheta[i][j][k];
		if( markerset.markers[m].texture.ctheta[i][j][k] < thetamin  ) thetamin =  markerset.markers[m].texture.ctheta[i][j][k];
	      }
	    }
	  }
	  nn=0;
	}
      }
      thetamean /= (PetscScalar)nn;
      printf("[%d]: thetamean = %e, max=%e, min=%e\n",rank,thetamean,thetamax,thetamin);
    }

    /* initialize the texture */
    if( options.restartStep == 0){
      initializeIsotropicFabric( &markerset );
    }

    /* save markers for debugging*/
    PetscLogStagePush(stages[11]);
    ierr=saveMarkersBinary( &markerset, -5,-5,0.0);
    saveTextureBinary(&markerset,&materials, -5, -5);
    saveNodalFields( &nodalFields, &grid, -5,-5,0.0,0.0,0);

    PetscLogStagePop();

    formViscoplasticM( &markerset, &materials, &options, 0);

    PetscInt iTime;
    PetscScalar dt;
    PetscInt isYielding;
    PetscInt iPlastic;
    
    displacementdt=0.0;/* initially*/
    for(iTime=iTime0;iTime<options.nTime;iTime++){ /* BEGIN TIMESTEP LOOP */
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Beginning Timestep %d, elapsedTime %e years\n",iTime,elapsedTime/(365.25*24*3600));
      PetscSynchronizedFlush(PETSC_COMM_WORLD);    
      isYielding = 1;
      iPlastic=0;
      dt=options.dtMax;/* initial timestep*/
      findMarkerCells( &markerset, &grid);
      ierr = exchangeMarkers( &markerset, &grid, &options);
      findMarkerCells( &markerset, &grid);
      ierr = checkMarkerDensity(&markerset, &grid, &options, r);
      findMarkerCells( &markerset, &grid);

      /* now update marker density using preferred timestep*/
      {
	PetscInt m;
	for(m=0;m<markerset.nMark;m++){
	  if(markers[m].cellX != -1){
	    markers[m].rho = materials.materialRho[(PetscInt) markers[m].Mat]*(1-materials.materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
	    markers[m].rhodot = 0.0;
	  }	      
	}
      }
      
      PetscLogStagePush(stages[5]);
      /* project all marker fields onto nodes */

      PetscLogStagePop();

      //updateDamageViscosity( &grid, &markerset, &materials, &options, dt);
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
	  //saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	  goto nextMonte;
	}
      }
      /* update boundary conditions, which can vary in time */
      updateBoundaryConditions( &options, &boundaryValues, elapsedTime );

      while(isYielding && iPlastic < options.maxNumPlasticity){/*begin plastic iteration loop*/
	/* check for plastic yielding */
	isYielding=1; /* for anisotropy, always do maxNumPlasticity plastic iterations */
	
/* 	ierr = checkPlasticYielding(&grid, &markerset, &materials, displacementdt, &isYielding, &options); */

	{ /* do an allreduce on isYielding - if any cpu is still yielding, all should still think that they're yielding*/
	  PetscInt tmp1=isYielding;
	  ierr=MPI_Allreduce( &tmp1,&isYielding,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);     
	}
	
	if(isYielding) printf("iTime=%d,iPlastic = %d Yielding detected\n",iTime,iPlastic);
	if(!isYielding && !rank) printf("no more yielding\n");
	if( iTime < options.plasticDelay && !rank) {isYielding=0; printf("yielding disabled\n");}
		
	/* Calculate correct viscoplastic response for current assumed marker stress state*/
	PetscLogStagePush(stages[12]);
	if(!rank) printf("Finding correct out-of-plane stresses\n");
	if( iPlastic==0 ){ formViscoplasticMNewton( &markerset, &materials, &options);}
	PetscLogStagePop();
	/* project marker quantities onto nodes: kThermal, Cp, rho, T*/
	ierr = limitAnisotropicViscosity(&markerset,&materials,&options);
	ierr = projectMarkersNodesAll(&markerset, &grid, &nodalFields, &materials, &options);
	
	/* form mechanical problem LHS*/
	/* zero out LHS, RHS*/
	ierr = VecZeroEntries( RHS);CHKERRQ(ierr);	
	ierr = VecZeroEntries( RHSz); CHKERRQ(ierr);
	PetscLogStagePush(stages[1]);
	//ierr=formVEPSystem( &nodalFields, &grid, LHS, RHS,LHSz,RHSz, &Kbond, &Kcont, options.gy, dt, boundaryValues.vbx,boundaryValues.vbz);
	ierr=formAnisotropicSystem( &nodalFields, &grid, LHS, RHS,LHSz,RHSz, &Kbond, &Kcont, options.gy, dt, &boundaryValues ,&options);
	PetscLogStagePop();
	Vec mechanicalS;/* solution to mechanical problem*/
	ierr = VecDuplicate(RHS, &mechanicalS);
	PetscLogStagePush(stages[3]);
	ierr = kspLinearSolve(LHS,PETSC_NULL, RHS, mechanicalS, "stokes_",PETSC_NULL);
	PetscLogStagePop();
	
	/* move solutions from solution vectors to nodalFields*/
	PetscInt fn=0;/* flag for nan in solution*/
	ierr= retrieveSolutions( &grid, &nodalFields, mechanicalS, PETSC_NULL, Kcont, &fn);
	ierr = VecDestroy(&mechanicalS); CHKERRQ(ierr);

	if(fn) {
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Found nan in solution\n");CHKERRQ(ierr);
	  PetscLogStagePush(stages[11]);
	  saveNodalFields( &nodalFields, &grid, iMonte,iTime,dt,elapsedTime,0);
	  PetscLogStagePush(stages[11]);
	  if(iMonte < options.nMonte-1){
	    goto nextMonte;
	  }else{
	    goto abort;
	  }
	}/* end deal with nan in solution*/
	displacementdt = dt;
	/* determine displacement timestep*/
	ierr=limitDisplacementTimestep( &grid, &nodalFields, &displacementdt, &options);CHKERRQ(ierr);
	/* check fabric timestep */
	{
/* 	  if( options.textureDevelopment){ */
/* 	    PetscScalar fabricdt = getFabricDT(&markerset); */
/* 	    if( fabricdt < displacementdt){  */
/* 	      printf("adjusting dt from %e to %e for fabric\n",displacementdt, fabricdt); */
/* 	      displacementdt = fabricdt; */
/* 	    } */
/* 	  } */
	}


	//ierr = VecZeroEntries(nodalHeating);CHKERRQ(ierr);
	ierr=PetscLogStagePush(stages[6]);
	ierr=nodalStressStrain(&grid, &nodalFields,&options, &boundaryValues, displacementdt,nodalHeating, options.gy);CHKERRQ(ierr);
	ierr=PetscLogStagePop();
	

	/* save nodalFields for debugging*/
	PetscLogStagePush(stages[11]);
	/*ierr=saveNodalFieldsMatlab( &nodalFields, &grid,iMonte,iTime, dt, elapsedTime,0);*/ /* for debugging*/
	PetscLogStagePop();
	/* compute marker strain and pressure*/
	ierr=updateMarkerStrainPressure( &grid,&nodalFields, &markerset, &materials, &options, displacementdt);

	/* project residual onto nodes */
	{
	  PetscScalar *srr;
	  ierr = PetscMalloc( markerset.nMark*sizeof(PetscScalar), &srr ); CHKERRQ(ierr);
	  PetscInt mm;
	  for(mm=0;mm<markerset.nMark;mm++){
	    srr[mm] = markerset.markers[mm].strainRateResidual;
	  }
	  ierr= projectMarkersNodesFromScalar(&markerset, &grid, srr, nodalFields.strainRateResidual);
	  ierr = PetscFree(srr);CHKERRQ(ierr);
	}

	/* reduce norm and print value */
	{
	  PetscScalar resnorm;
	  ierr = VecNorm(nodalFields.strainRateResidual,NORM_2,&resnorm);  
	  if(!rank) printf("Strain Rate Resdiual Norm curently %e\n",resnorm);      
	}
#ifdef debug1
	/* dump marker information at this plastic step */
	ierr=saveNodalFields( &nodalFields, &grid,iPlastic,-iTime, dt, elapsedTime,0);
	ierr=saveMarkersBinary( &markerset, iPlastic,-iTime,elapsedTime);
	//	saveTextureBinary(&markerset,&materials, iPlastic, -iTime);	
#endif
	/* update marker effective viscosity to match marker strain reate*/
	iPlastic++;
      }/* end plasticity/non-newtonian loop*/
      
      /* update marker stress - do this OUTSIDE of plasticity loop*/
      ierr= updateMarkerStress( &grid, &nodalFields,&markerset, &materials);CHKERRQ(ierr);

      if(!rank) printf("Done with plasticity loops.\n");
      /* do sub-grid stress diffusion*/
      //PetscLogStagePush(stages[7]);
      /* subgrid stress diffusion for viscoelastic problems would go here */
      //PetscLogStagePop();

      /* project velocity field onto markers*/
      //ierr= projectVelocityToMarkers(&markerset, &grid, &nodalFields );
      /* calculate adiabatic heating */
      if( !options.shearHeating ){
	ierr=VecZeroEntries(nodalHeating);CHKERRQ(ierr); /* this line disables shear heating */
      }
      if( options.adiabaticHeating ){
	adiabaticHeating( &grid, &markerset, &nodalFields, &materials, &options);
	ierr=VecAXPY(nodalHeating,1.0,nodalFields.ha);CHKERRQ(ierr);/* add adiabatic heating to shear heating */
      }
      /* set thermal problem timestep to displacement dt*/
      dt = displacementdt;
      
      /* Form thermal system (diffusion, variable coefficient, non-uniform grid */
      PetscLogStagePush( stages[8] );
      ierr = VecCopy(nodalHeating,thermalRHS);CHKERRQ(ierr);
      ierr = enforceThermalBCs1( &grid, &options, &nodalFields);CHKERRQ(ierr);
      ierr = formThermalSystem( &grid, &nodalFields, thermalLHS, thermalRHS, dt, &options);CHKERRQ(ierr);
      PetscLogStagePop();
      /* make the solution vector */
      Vec thermalS;
      ierr = VecDuplicate( thermalRHS, &thermalS);CHKERRQ(ierr);
      /* Solve the thermal system: */
      PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(thermalLHS,PETSC_NULL, thermalRHS, thermalS, "energy_",PETSC_NULL);CHKERRQ(ierr); PetscLogStagePop();
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
	PetscLogStagePush( stages[9] );  ierr = kspLinearSolve(thermalLHS, PETSC_NULL, thermalRHS, thermalS,"energy_",PETSC_NULL);CHKERRQ(ierr); PetscLogStagePop();
	/* end repeat temperature solution */
      }
      ierr = VecDestroy(&dT);CHKERRQ(ierr);
      /* Apply sub-grid temperature corrections and project corrected temperature onto the markers (all integrated into one routine*/
      PetscLogStagePush( stages[10] );      ierr=subgridTemperatureChanges(thermalS, &grid, &nodalFields, &markerset, &materials,dt, &options);      PetscLogStagePop();
      ierr = VecDestroy(&thermalS); /* free the thermal solution vector */
            
      /* update damage rate on markers */
      //updateDamageRate( &grid, &markerset, &materials, &options, dt);
        
      /* save solution */
      if(!(iTime % (options.saveInterval))) {
	ierr=saveNodalFields( &nodalFields, &grid,iMonte,iTime, dt, elapsedTime,0);
	ierr=saveMarkersBinary( &markerset, iMonte,iTime,elapsedTime);
	saveTextureBinary(&markerset,&materials, iMonte, iTime);
      }
      if(!(iTime % options.saveInterval)) {
	//saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime,elapsedTime);
      }
      /* advect markers */
      //ierr=advectMarkers(&markerset,&grid,dt);
      ierr= advectMarkersRK(&markerset, &nodalFields, &grid, &options, &boundaryValues, dt);
      /* if we are updating fabric, make it so */
      if( options . textureDevelopment ){
	updateFabric( &markerset, &materials, dt);
      }

      elapsedTime+=dt;
      if(elapsedTime > options.totalTime || iTime==options.nTime-1){
	PetscLogStagePush(stages[11]);
	//saveGriddedMarkersBinary( &markers, &grid, 5*grid.NX, 5*grid.NY,iMonte,iTime);
	PetscLogStagePop();
	ierr=saveMarkersBinary( &markerset,iMonte,iTime,elapsedTime);
	saveTextureBinary(&markerset,&materials, iMonte, iTime);
	goto nextMonte;
      }
      

    }/* end time stepping*/
  nextMonte:;


  }/* end montecarlo loop */
 abort:;
  ierr = MatDestroy(& thermalLHS);CHKERRQ(ierr);
  ierr = VecDestroy(&thermalRHS);CHKERRQ(ierr);
  ierr = VecDestroy(&nodalHeating);CHKERRQ(ierr);
  finalizeLogging();
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
