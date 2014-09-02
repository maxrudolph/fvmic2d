#include "fdcode.h"
#include "options.h"
#include <float.h>

#ifndef LINEWIDTH
#define LINEWIDTH 128
#endif

PetscErrorCode csvOptions(char *csvFileName, Options *options, Materials *materials){
  FILE *csvFile;
  //  PetscErrorCode ierr;
  char line[LINEWIDTH];
  PetscMPIInt rank,size;
  char fn[12] = "options.csv";
  char *ifn;
  setOptions (options);

  PetscFunctionBegin;
  /* for testing: */
  if( csvFileName == NULL ){
    ifn = &fn[0];
  } else {
    ifn = csvFileName;
  }

  PetscFunctionBegin;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  if(!rank) printf("reading from %s\n",ifn);

  if(!rank){ 
    csvFile = fopen( ifn, "r");
    /* get a line*/
    while( fgets(line, sizeof(line), csvFile) != NULL ){
      /* parse the line */
      PetscInt idxComma=0;
      /* find the location of the comma */
      if( !strncmp( "#",line,1) ){
	
      } else {
	while( line[++idxComma] != ',');
	//      printf("idxComma=%d\n", idxComma);
	//printf("comma = %c \n",line[idxComma]);
	/* parse stuff before the comma */
	
	/*ignore this line */
	if( !strncmp( "LX", line, 2) ){
	  //	  printf("reading LX\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> LX)); printf("LX = %le\n",options->LX);
	  
	} else if( !strncmp( "LY", line, 2) ){
	  //	  printf("reading LY\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> LY)); printf("LY = %le\n",options->LY);
	} else if( !strncmp( "NX", line, 2) ){
	  //printf("reading NX\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> NX)); printf("NX = %d\n",options->NX);
	  
	} else if( !strncmp( "NY", line, 2) ){
	  //printf("reading NY\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> NY)); printf("NY = %d\n",options->NY);
	
	} else if( !strncmp( "gridRatio", line, 9) ){
	  //printf("reading NY\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> gridRatio)); printf("gridRatio = %e\n",options->gridRatio);
	
	
	/*options->NMX=10;*//* number of markers per cell, x*/
	} else if( !strncmp( "NMX", line, 3) ){
	  //printf("reading NMX\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> NMX)); printf("NMX = %d\n",options->NMX);
	  
	  /* options->NMY=10; */ /* markers per cell y*/
	} else if( !strncmp( "NMY", line, 2) ){
	  //printf("reading NMY\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> NMY)); printf("NMY = %d\n",options->NMY);
	 
	/*   options->maxMarkFactor = 1.5; */ /*memory to allocate for markers is this number*original number of markers*/
	} else if( !strncmp( "maxMarkFactor", line,13 ) ){
	  //printf("reading maxMarkFactor\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> maxMarkFactor)); printf("maxMarkFactor = %le\n",options->maxMarkFactor);
	  
	  /* options->minMarkers = 60; */ /* munimum number of markers per cell */
	} else if( !strncmp( "minMarkers", line, 10) ){
	  //printf("reading minMarkers\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> minMarkers)); printf("minMarkers = %d\n",options->minMarkers);
	} else if( !strncmp( "maxMarkers", line, 10) ){
	  //printf("reading minMarkers\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> maxMarkers)); printf("maxMarkers = %d\n",options->maxMarkers);
	

	} else if( !strncmp( "slabAngle", line, 9) ){
	  sscanf( &line[idxComma+1],"%le\n", &(options -> slabAngle)); printf("slabAngle = %le\n",options->slabAngle);
	  options->slabAngle *= M_PI/180.0;
	  printf("slabAngle in radians %le\n",options->slabAngle);
	} else if( !strncmp( "slabVelocity", line, 12) ){
	  sscanf( &line[idxComma+1],"%le\n", &(options -> slabVelocity)); printf("slabVelocity = %le\n",options->slabVelocity);
	  
	  /* time step, plastic iteration*/

	  /*   options->dtMax = 1e0*(365.25*24*3600); */
	} else if( !strncmp( "dtMax", line, 5) ){
	  //printf("reading dtMax\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> dtMax)); printf("dtMax = %le\n",options->dtMax);
	} else if( !strncmp( "displacementStepLimit", line, 21) ){
	  sscanf( &line[idxComma+1],"%le\n", &(options -> displacementStepLimit)); printf("displacementStepLimit = %le\n",options->displacementStepLimit);
	
	} else if( !strncmp( "maxTempChange", line, 13) ){
	  //printf("reading dtMax\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> maxTempChange)); printf("maxTempChange = %le\n",options->maxTempChange);

	  /*   options->nTime = 2001; */ /* max number of timesteps*/
	} else if( !strncmp( "nTime", line, 5) ){
	  //printf("reading nTime\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> nTime)); printf("nTime = %d\n",options->nTime);
	} else if( !strncmp( "restartStep", line, 11) ){
	  //printf("reading nTime\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> restartStep)); printf("restartStep = %d\n",options->restartStep);
	
	  /* options->totalTime = 4e6*(365.25*24*3600); */ /* maximum model time to run*/
	} else if( !strncmp( "totalTime", line, 9) ){
	  //printf("reading totalTime\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> totalTime)); printf("totalTime = %le\n",options->totalTime);
	
	  /*   options->maxNumPlasticity=1; */
	} else if( !strncmp( "maxNumPlasticity", line, 16) ){
	  //printf("reading maxNumPlasticity\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> maxNumPlasticity)); printf("maxNumPlasticity = %d\n",options->maxNumPlasticity);
	
/*   options->plasticDelay = 1; */
	} else if( !strncmp( "plasticDelay", line,12) ){
	  //printf("reading plasticDelay\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> plasticDelay)); printf("plasticDelay = %d\n",options->plasticDelay);
	
/*   options->etamin = 1.0e12; */ /* minimum viscosity allowed in plastic yielding*/
	} else if( !strncmp( "etamin", line, 6) ){
	  //printf("reading etamin\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> etamin)); printf("etamin = %le\n",options->etamin);

	} else if(  !strncmp( "fractionalEtamin", line, 16) ){
	  //printf("reading etamin\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> fractionalEtamin)); printf("fractionalEtamin = %le\n",options->fractionalEtamin);
	
/*   options->etamax = 1.0e36; */
	} else if( !strncmp( "etamax", line, 6) ){
	  //printf("reading etamax\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> etamax)); printf("etamax = %le\n",options->etamax);
	
	  /*texture - for ice viscoplasticity*/
/* 	} else if(  !strncmp( "doTexture", line, 9) ){ */
/* 	  //printf("reading doMonte\n"); */
/* 	  sscanf( &line[idxComma+1],"%d\n", &(options -> doTexture)); printf("doTexture = %d\n",options->doTexture); */

       

	} else if(  !strncmp( "grainSize", line, 9) ){
	  //printf("reading doMonte\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> grainSize)); printf("grainSize = %e\n",options->grainSize);
	
	  /* Monte-carlo search */
	
	
  /* damage stuff */
	
	  
	  /* model setup/ICs */
	  /*   options->Tperturb = 1; */
 	} else if( !strncmp( "Tperturb", line, 8) ){
	  //printf("reading Tperturb\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> Tperturb)); printf("Tperturb = %le\n",options->Tperturb);
	  
	  
	} else if( !strncmp( "shearHeating", line, 12) ){

	  sscanf( &line[idxComma+1],"%d\n", &(options -> shearHeating)); printf("shearHeating = %d\n",options->shearHeating);

	} else if( !strncmp( "adiabaticHeating", line, 16) ){

	  sscanf( &line[idxComma+1],"%d\n", &(options -> adiabaticHeating)); printf("adiabaticHeating = %d\n",options->adiabaticHeating);

	} else if( !strncmp( "mechBCLeft", line, 10) ){
	  sscanf( &line[idxComma+1],"%d,%d,%d,%le,%le,%le\n", &(options->mechBCLeft.type[0]), &(options->mechBCLeft.type[1]), &(options->mechBCLeft.type[2]), &(options->mechBCLeft.value[0]), &(options->mechBCLeft.value[1]), &(options->mechBCLeft.value[2]));
	} else if( !strncmp( "mechBCRight", line, 10) ){
	  sscanf( &line[idxComma+1],"%d,%d,%d,%le,%le,%le\n", &(options->mechBCRight.type[0]), &(options->mechBCRight.type[1]), &(options->mechBCRight.type[2]), &(options->mechBCRight.value[0]), &(options->mechBCRight.value[1]), &(options->mechBCRight.value[2]));
	} else if( !strncmp( "mechBCTop", line, 9) ){
	  sscanf( &line[idxComma+1],"%d,%d,%d,%le,%le,%le\n", &(options->mechBCTop.type[0]), &(options->mechBCTop.type[1]), &(options->mechBCTop.type[2]), &(options->mechBCTop.value[0]), &(options->mechBCTop.value[1]), &(options->mechBCTop.value[2]));
	} else if( !strncmp( "mechBCBottom", line, 12) ){
	  sscanf( &line[idxComma+1],"%d,%d,%d,%le,%le,%le\n", &(options->mechBCBottom.type[0]), &(options->mechBCBottom.type[1]), &(options->mechBCBottom.type[2]), &(options->mechBCBottom.value[0]), &(options->mechBCBottom.value[1]), &(options->mechBCBottom.value[2]));

	  /* Thermal BCs */
	} else if( !strncmp( "thermalBCBottom", line, 15) ){
	  sscanf( &line[idxComma+1],"%d,%le\n", &(options->thermalBCBottom.type[0]), &(options->thermalBCBottom.value[0]));
	} else if( !strncmp( "thermalBCTop", line, 12) ){
	  sscanf( &line[idxComma+1],"%d,%le\n", &(options->thermalBCTop.type[0]), &(options->thermalBCTop.value[0]));
	} else if( !strncmp( "thermalBCLeft", line, 13) ){
	  sscanf( &line[idxComma+1],"%d,%le\n", &(options->thermalBCLeft.type[0]), &(options->thermalBCLeft.value[0]));
	} else if( !strncmp( "thermalBCRight", line, 14) ){
	  sscanf( &line[idxComma+1],"%d,%le\n", &(options->thermalBCRight.type[0]), &(options->thermalBCRight.value[0]));
	  
	  
	  /* body forces*/
	  /*   options->gy = 1.3; */
	} else if( !strncmp( "gy", line, 2) ){
	  //printf("reading gy\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> gy)); printf("gy = %le\n",options->gy);
	  /*   options->gx = 0; */
	} else if( !strncmp( "gx", line, 2) ){
	  //printf("reading gx\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> gx)); printf("gx = %le\n",options->gx);
	  
	  /* I/O etc*/
	  /*   options->saveInterval = 25; */
	} else if( !strncmp( "saveInterval", line, 12) ){
	  //printf("reading saveInterval\n");
	  sscanf( &line[idxComma+1],"%d\n", &(options -> saveInterval)); printf("saveInterval = %d\n",options->saveInterval);
	  
  /* subgrid diffusion constants*/
/*   options-> subgridStressDiffusivity = 1.0; */
	} else if( !strncmp( "subgridStressDiffusivity", line, 24) ){
	  //printf("reading subgridStressDiffusivity\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> subgridStressDiffusivity)); printf("subgridStressDiffusivity = %le\n",options->subgridStressDiffusivity);
/*   options-> subgridTemperatureDiffusivity = 1.0; */
	} else if( !strncmp( "subgridTemperatureDiffusivity", line, 29) ){
	  //printf("reading subgridTemperatureDiffusivity\n");
	  sscanf( &line[idxComma+1],"%le\n", &(options -> subgridTemperatureDiffusivity)); printf("subgridTemperatureDiffusivity = %le\n",options->subgridTemperatureDiffusivity);
	} else if( !strncmp( "grainSize", line, 9) ){

	  sscanf( &line[idxComma+1],"%le\n", &(options -> grainSize)); printf("grainSize = %le\n",options->grainSize);

	  /* Material Properties */
	} else if( !strncmp( "nMaterials", line, 10) ){
	  //printf("reading nMaterials\n");
	  sscanf( &line[idxComma+1],"%d\n", &(materials -> nMaterials)); printf("nMaterials = %d\n",materials->nMaterials);
	  
	} else if( !strncmp( "materialEta", line,11 ) ){
	  /* first read material number */
	  //printf("reading materialEta\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialEta[thisMat] = val1;
	  printf("materialEta[%d] = %le\n",thisMat,materials->materialEta[thisMat]);


/*   materials-> materialRho[0]=1000.0; */
	} else if( !strncmp( "materialRho", line,11 ) ){
	  /* first read material number */
	  //printf("reading materialRho\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialRho[thisMat] = val1;
	  printf("materialRho[%d] = %le\n",thisMat,materials->materialRho[thisMat]);

	  /*   materials-> materialkThermal[2]=100.0; */
	} else if( !strncmp( "materialkThermal", line,16 ) ){
	  /* first read material number */
	  //printf("reading materialkThermal\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialkThermal[thisMat] = val1;
	  printf("materialkThermal[%d] = %le\n",thisMat,materials->materialkThermal[thisMat]);
  
	  /*   materials->materialCp[2]=1100.0; */
	} else if( !strncmp( "materialCp", line,10 ) ){
	  /* first read material number */
	  //printf("reading materialCp\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialCp[thisMat] = val1;
	  printf("materialCp[%d] = %le\n",thisMat,materials->materialCp[thisMat]);


	  /*   materials->materialAlpha[2] = 0.0; */
	} else if( !strncmp( "materialAlpha", line,13 ) ){
	  /* first read material number */
	  //printf("reading materialAlpha\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialAlpha[thisMat] = val1;
	  printf("materialAlpha[%d] = %le\n",thisMat,materials->materialAlpha[thisMat]);

	  /*   materials->materialMu[2] = 1e99; */
	} else if( !strncmp( "materialMu", line,10 ) ){
	  /* first read material number */
	  //printf("reading materialMu\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialMu[thisMat] = val1;
	  printf("materialMu[%d] = %le\n",thisMat,materials->materialMu[thisMat]);
	  
	  /*   materials->materialCohesion[2] = 1e99; */
	} else if( !strncmp( "materialCohesion", line,16 ) ){
	  /* first read material number */
	  //printf("reading materialCohesion\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialCohesion[thisMat] = val1;
	  printf("materialCohesion[%d] = %le\n",thisMat,materials->materialCohesion[thisMat]);


	  /*   materials->materialFriction[2] = 0.6; */
	} else if( !strncmp( "materialFriction", line,16 ) ){
	  /* first read material number */
	  //printf("reading materialFriction\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->materialFriction[thisMat] = val1;
	  printf("materialFriction[%d] = %le\n",thisMat,materials->materialFriction[thisMat]);

	} else if( !strncmp( "hasPlasticity", line,13 ) ){
	  PetscInt thisMat;
	  PetscInt val1;
	  sscanf( &line[idxComma+1],"%d,%d\n",&thisMat,&val1); /* read material number, idx, value */
	  materials->hasPlasticity[thisMat] = val1;
	  printf("hasPlasticity[%d] = %d\n",thisMat,materials->hasPlasticity[thisMat]);

	} else if( !strncmp( "materialC", line,9 ) ){
	  PetscInt thisMat;
	  PetscInt idx;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%d,%le\n", &thisMat,&idx,&val1); /* read material number, idx, value */
	  materials->C[thisMat][idx] = val1;
	  printf("materialC[%d][%d] = %le\n",thisMat,idx,materials->C[thisMat][idx]);

	} else if( !strncmp( "materialF", line,9 ) ){
	  PetscInt thisMat;
	  PetscInt idx;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%d,%le\n",&thisMat,&idx,&val1); /* read material number, idx, value */
	  materials->F[thisMat][idx] = val1;
	  printf("materialF[%d][%d] = %le\n",thisMat,idx,materials->F[thisMat][idx]);

	} else if(  !strncmp( "binghamYieldStress", line,18 ) ){
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n",&thisMat,&val1); /* read material number, value */
	  materials->binghamYieldStress[thisMat] = val1;
	  printf("binghamYieldStress[%d] = %le\n",thisMat,materials->binghamYieldStress[thisMat]);


	} else if( !strncmp( "materialGamma", line,13 ) ){
	  PetscInt thisMat;
	  PetscInt idx;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%d,%le\n",&thisMat,&idx,&val1); /* read material number, idx, value */
	  materials->gamma[thisMat][idx] = val1;
	  printf("materialGamma[%d][%d] = %le\n",thisMat,idx,materials->gamma[thisMat][idx]);


	} else if( !strncmp( "hasDamage", line,9 ) ){
	  /* first read material number */
	  //printf("reading hasDamage\n");
	  PetscInt thisMat;
	  PetscInt val1;
	  sscanf( &line[idxComma+1],"%d,%d\n", &thisMat,&val1); /* read material number, value */
	  materials->hasDamage[thisMat] = val1;
	  printf("hasDamage[%d] = %d\n",thisMat,materials->hasDamage[thisMat]);

	 

	} else if( !strncmp( "hasDilation", line,11 ) ){
	  /* first read material number */
	  //printf("reading hasDamage\n");
	  PetscInt thisMat;
	  PetscInt val1;
	  sscanf( &line[idxComma+1],"%d,%d\n", &thisMat,&val1); /* read material number, value */
	  materials->hasDilation[thisMat] = val1;
	  printf("hasDilation[%d] = %d\n",thisMat,materials->hasDilation[thisMat]);

/*   materials->hayhurstAlpha[0] = 0.3; */
	} else if( !strncmp( "hayhurstAlpha", line,13 ) ){
	  /* first read material number */
	  //printf("reading hayhurstAlpha\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->hayhurstAlpha[thisMat] = val1;
	  printf("hayhurstAlpha[%d] = %le\n",thisMat,materials->hayhurstAlpha[thisMat]);

/*   materials->hayhurstBeta[0] = 0.6; */
	} else if( !strncmp( "hayhurstBeta", line,12 ) ){
	  /* first read material number */
	  //printf("reading hayhurstBeta\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->hayhurstBeta[thisMat] = val1;
	  printf("hayhurstBeta[%d] = %le\n",thisMat,materials->hayhurstBeta[thisMat]);

/*   materials->damagek[0] = 0.0; */
	} else if( !strncmp( "damagek", line,7 ) ){
	  /* first read material number */
	  //printf("reading damagek\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->damagek[thisMat] = val1;
	  printf("damagek[%d] = %le\n",thisMat,materials->damagek[thisMat]);

/*   materials->damageAlpha3[0] = 5.56806e+01; */
	} else if( !strncmp( "damageAlpha3", line,12 ) ){
	  /* first read material number */
	  //printf("reading damageAlpha3\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->damageAlpha3[thisMat] = val1;
	  printf("damageAlpha3[%d] = %le\n",thisMat,materials->damageAlpha3[thisMat]);

/*   materials->damager[0] = 0.43; */
	} else if( !strncmp( "damager", line,7 ) ){
	  /* first read material number */
	  //printf("reading damager\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->damager[thisMat] = val1;
	  printf("damager[%d] = %le\n",thisMat,materials->damager[thisMat]);

/*   materials->damageB[0] = 1.673757e-27; *//* 1e-30 produced damage at a reasonable rate*/
	} else if( !strncmp( "damageB", line,7 ) ){
	  /* first read material number */
	  //printf("reading damageB\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->damageB[thisMat] = val1;
	  printf("damageB[%d] = %le\n",thisMat,materials->damageB[thisMat]);

/*   materials->damagem[0] = 2.595666e-01; */
	} else if( !strncmp( "damagem", line,7 ) ){
	  /* first read material number */
	  //printf("reading damagem\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->damagem[thisMat] = val1;
	  printf("damagem[%d] = %le\n",thisMat,materials->damagem[thisMat]);
  
  /* temperature dependent viscosity*/
/*   materials -> hasEtaT[0] = 1; */
	} else if( !strncmp( "hasEtaT", line,7 ) ){
	  /* first read material number */
	  //printf("reading hasEtaT\n");
	  PetscInt thisMat;
	  PetscInt val1;
	  sscanf( &line[idxComma+1],"%d,%d\n", &thisMat,&val1); /* read material number, value */
	  materials->hasEtaT[thisMat] = val1;
	  printf("hasEtaT[%d] = %d\n",thisMat,materials->hasEtaT[thisMat]);

/*   materials -> hasEtaT[1] = 0; */
/*   materials -> hasEtaT[2] = 0; */

/*   materials -> QonR[0] = 6000.0; */
	} else if( !strncmp( "QonR", line,4 ) ){
	  /* first read material number */
	  //printf("reading QonR\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->QonR[thisMat] = val1;
	  printf("QonR[%d] = %le\n",thisMat,materials->QonR[thisMat]);

/*   materials -> Tref[0] = 273.15; */ /* kelvin*/
	} else if( !strncmp( "Tref", line,4 ) ){
	  /* first read material number */
	  //printf("reading Tref\n");
	  PetscInt thisMat;
	  PetscScalar val1;
	  sscanf( &line[idxComma+1],"%d,%le\n", &thisMat,&val1); /* read material number, value */
	  materials->Tref[thisMat] = val1;
	  printf("Tref[%d] = %le\n",thisMat,materials->Tref[thisMat]);	  

	} else {
	  
	  printf("WARNING: Unknown input parameter %s !!!\n",line);
	}	
      }
    }
      
    fclose( csvFile );
  }
  /* send parameters to all nodes */
  MPI_Bcast( options, sizeof( Options ), MPI_BYTE, 0, PETSC_COMM_WORLD);
  MPI_Bcast( materials, sizeof( Materials ), MPI_BYTE, 0, PETSC_COMM_WORLD);
  /*   printf("[%d] has mu=%f , rho0 = %f after scatteinrg Parameters\n",rank,params->mu0,params->rho0) */;
  
  PetscFunctionReturn(0);
}

void setOptions( Options *options){
  /* set options to default (often meaningless) values */
  options->NX = 0;
  options->NY = 0;
  options->NMX = 0;
  options->NMY = 0;
  options->gridRatio = 0.0;
  options->nTime = 0;
  options->restartStep = 0;
  options->totalTime = 0.0;
  options->maxNumPlasticity = 0;
  options->plasticDelay = 0;
  options->Tperturb = 0.0;
  options->shearHeating = 0;
  options->adiabaticHeating = 0;
  options->gx = 0;
  options->gy = 0;
  options->gz = 0;
  //BC stuff here
  options->fractionalEtamin = 1.0;
  options->etamin = 0.0;
  options->etamax = DBL_MAX;
  options->dtMax = DBL_MAX;
  options->displacementStepLimit = 1.0;
  options->maxTempChange = DBL_MAX;
  options->saveInterval = 1;
  options->subgridStressDiffusivity = 0.0;
  options->subgridTemperatureDiffusivity = 0.0;
  options->maxMarkFactor = 0.0;
  options->minMarkers = 0;
  options->maxMarkers = INT_MAX;
  options->grainSize = 0.0;
  options->slabAngle = 0.0;
  options->slabVelocity = 0.0;






}
