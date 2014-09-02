#include "fdcode.h"
#include "gridMarkers.h"
#include "io.h"
#include "version.h"
#include "profile.h"
#include "petscbag.h"

PetscErrorCode saveRunInfo( Options *options, Materials *materials, PetscInt iMonte){
  PetscErrorCode ierr=0;
  char fn[80];
  FILE *of;
  PetscMPIInt rank;
  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  setLogStage( LOG_IO );
  sprintf(fn,"output/runInfo.%d_%d.csv",iMonte,rank);
  of=fopen(fn,"w");
  /* spit out the version of the code that we're running */
  fprintf(of,"# VERSION: ");
  fprintf(of,MARKERCODE_HG_VERSION);
  fprintf(of,"\n");
  /* save the material properties and run options as a csv file*/
  fprintf(of,"#Options:\n");
  fprintf(of,"LX,%e\n",options->LX);
  fprintf(of,"LY,%e\n",options->LY);
  fprintf(of,"NX,%d\n",options->NX);
  fprintf(of,"NY,%d\n",options->NY);
  fprintf(of,"gy,%e\n",  options->gy);
  fprintf(of,"gx,%e\n",  options->gx);
  fprintf(of,"gridRatio,%e\n",options->gridRatio); 
  fprintf(of,"NMX,%d\n",options->NMX);
  fprintf(of,"NMY,%d\n",options->NMY); 
  fprintf(of,"maxMarkFactor,%le\n",options->maxMarkFactor);
  fprintf(of,"minMarkers,%d\n",options->minMarkers);
  fprintf(of,"maxMarkers,%d\n",options->maxMarkers);
  fprintf(of,"saveInterval,%d\n",options->saveInterval);

  /* time stepping, plastic iterations */
  fprintf(of,"dtMax,%e\n",options->dtMax);
  fprintf(of,"nTime,%d\n",options->nTime);
  fprintf(of,"restartStep,%d\n",options->restartStep);
  fprintf(of,"totalTime,%e\n",options->totalTime);
  fprintf(of,"maxNumPlasticity,%d\n",options->maxNumPlasticity);
  fprintf(of,"plasticDelay,%d\n",options->plasticDelay);
  fprintf(of,"displacementStepLimit,%e\n",options->displacementStepLimit);
  fprintf(of,"maxTempChange,%e\n",options->maxTempChange);

  /* initial and boundary conditions */
  fprintf(of,"Tperturb,%e\n",options->Tperturb);
  {
    BCInfo *thisBC = &options->mechBCLeft;
    fprintf(of,"mechBCLeft,%d,%d,%d,%e,%e,%e\n",thisBC->type[0],thisBC->type[1],thisBC->type[2],thisBC->value[0],thisBC->value[1],thisBC->value[2]);
    thisBC = &options->mechBCRight;
    fprintf(of,"mechBCRight,%d,%d,%d,%e,%e,%e\n",thisBC->type[0],thisBC->type[1],thisBC->type[2],thisBC->value[0],thisBC->value[1],thisBC->value[2]);
    thisBC = &options->mechBCTop;
    fprintf(of,"mechBCTop,%d,%d,%d,%e,%e,%e\n",thisBC->type[0],thisBC->type[1],thisBC->type[2],thisBC->value[0],thisBC->value[1],thisBC->value[2]);
    thisBC = &options->mechBCBottom;
    fprintf(of,"mechBCBottom,%d,%d,%d,%e,%e,%e\n",thisBC->type[0],thisBC->type[1],thisBC->type[2],thisBC->value[0],thisBC->value[1],thisBC->value[2]);

    thisBC = &options->thermalBCLeft;
    fprintf(of,"thermalBCLeft,%d,%e\n",thisBC->type[0],thisBC->value[0]);
    thisBC = &options->thermalBCRight;
    fprintf(of,"thermalBCRight,%d,%e\n",thisBC->type[0],thisBC->value[0]);
    thisBC = &options->thermalBCTop;
    fprintf(of,"thermalBCTop,%d,%e\n",thisBC->type[0],thisBC->value[0]);
    thisBC = &options->thermalBCBottom;
    fprintf(of,"thermalBCBottom,%d,%e\n",thisBC->type[0],thisBC->value[0]);
  }

  /* physics */
  fprintf(of,"adiabaticHeating,%d\n", options->adiabaticHeating);
  fprintf(of,"shearHeating,%d\n", options->shearHeating);
  fprintf(of,"subgridStressDiffusivity,%e\n", options-> subgridStressDiffusivity);
  fprintf(of,"subgridTemperatureDiffusivity,%e\n", options-> subgridTemperatureDiffusivity);
  /* limits on material properties */
  fprintf(of,"fractionalEtamin,%e\n",options->fractionalEtamin);
  fprintf(of,"etamin,%e\n",options->etamin);
  fprintf(of,"etamax,%e\n",options->etamax);

  /* material properties*/
  fprintf(of,"#Material Properties:\n");

  fprintf(of,"grainSize,%e\n",options->grainSize);
  fprintf(of,"nMaterials,%d\n",materials->nMaterials);

  PetscInt i;
  for(i=0;i<materials->nMaterials;i++){
    fprintf(of,"materialEta,%d,%e\n",i,materials->materialEta[i]);
    fprintf(of,"materialRho,%d,%e\n",i,materials-> materialRho[i]);

    fprintf(of,"materialkThermal,%d,%e\n",i,materials-> materialkThermal[i]);
    fprintf(of,"materialCp,%d,%e\n",i,materials->materialCp[i]);
    fprintf(of,"materialAlpha,%d,%e\n",i,materials->materialAlpha[i]);
    fprintf(of,"materialMu,%d,%e\n",i,materials->materialMu[i]);
    fprintf(of,"hasPlasticity,%d,%d\n",i,materials->hasPlasticity[i]);
    fprintf(of,"materialCohesion,%d,%e\n",i,materials->materialCohesion[i]);
    fprintf(of,"materialFriction,%d,%e\n",i,materials->materialFriction[i]);
    fprintf(of,"materialGamma,%d,%e,%e\n",i,materials->gamma[i][0], materials->gamma[i][1]);
    fprintf(of,"materialF,%d,%e,%e\n",i,materials->F[i][0], materials->F[i][1]);
    fprintf(of,"materialC,%d,%e,%e\n",i,materials->C[i][0], materials->C[i][1]);
    fprintf(of,"binghamYieldStress,%d,%e\n",i,materials->binghamYieldStress[i]);
    fprintf(of,"hasDamage,%d,%d\n",i,materials->hasDamage[i]);
    fprintf(of,"hasDilation,%d,%d\n",i,materials->hasDilation[i]);
#ifdef TEXTURE
    fprintf(of,"hasTexture,%d,%d\n",i,materials->hasTexture[i]);
#endif
/*     if( materials->hasDamage[i]){ */
/*       fprintf(of,"hayhurstAlpha,%d,%e\n",i,materials->hayhurstAlpha[i]); */
/*       fprintf(of,"hayhurstBeta,%d,%e\n",i,materials->hayhurstBeta[i]); */
/*       fprintf(of,"damager,%d,%e\n",i, materials->damager[i]); */
/*       fprintf(of,"damageB,%d,%e\n",i, materials->damageB[i]); */
/*       fprintf(of,"damagem,%d,%e\n",i,  materials->damagem[i]); */
/*       fprintf(of,"damagek,%d,%e\n",i, materials->damagek[i]); */
/*       fprintf(of,"damageAlpha3,%d,%e\n",i,materials->damageAlpha3[i]); */
/*     } */
    fprintf(of,"hasEtaT,%d,%d\n",i,  materials -> hasEtaT[i]);
    if(materials->hasEtaT[i]){
      fprintf(of,"QonR,%d,%e\n",i, materials -> QonR[i]);
      fprintf(of,"Tref,%d,%e\n",i,  materials -> Tref[i]);
    }
  }

  fclose(of);
  PetscLogStagePop();
  PetscFunctionReturn(ierr);
  
}

PetscErrorCode saveGrid( GridData *grid){
  PetscMPIInt rank;
  setLogStage( LOG_IO );
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if(!rank){
    FILE *of;
    of=fopen("output/loadgrid.m","w");
    fprintf(of,"grid.x=[");
    PetscInt i;
    for(i=0;i<grid->NX;i++){
      fprintf(of,"%e\n",grid->x[i]);
    }
    fprintf(of,"];\n");
    fprintf(of,"grid.y=[");
    for(i=0;i<grid->NY;i++){
      fprintf(of,"%e\n",grid->y[i]);
    }
    fprintf(of,"];\n");
    fprintf(of,"grid.xc=[");
    for(i=0;i<grid->NX;i++){
      fprintf(of,"%e\n",grid->xc[i]);
    }
    fprintf(of,"];\n");
    fprintf(of,"grid.yc=[");
    for(i=0;i<grid->NY;i++){
      fprintf(of,"%e\n",grid->yc[i]);
    }
    fprintf(of,"];\n");

    fclose(of);
  }
  PetscLogStagePop();
  PetscFunctionReturn(0);
}

PetscErrorCode saveNodalFields( NodalFields *nodalFields, GridData *grid,PetscInt iMonte,PetscInt iTime, PetscScalar dt, PetscScalar elapsedTime, PetscInt format){
  /* save nodal fields in format specified by format */
  /* check for command line arguments to override format */
  {
    PetscBool set=PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,"-output_ascii_matlab",&set,PETSC_NULL);
    if( set ) format = OUTPUT_NODAL_FIELDS_ASCII_MATLAB;
    set = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,"-output_ascii_vtk",&set,PETSC_NULL);
    if( set ){ format = OUTPUT_NODAL_FIELDS_ASCII_VTK; printf("vtk flag set\n");}
    set = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,"-output_binary_vtk",&set,PETSC_NULL);
    if( set ) format = OUTPUT_NODAL_FIELDS_BINARY_VTK;
    set = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULL,"-output_binary",&set,PETSC_NULL);
    if( set ) format = OUTPUT_NODAL_FIELDS_BINARY;
  }
  printf("output format %d\n",format);
  PetscMPIInt rank;
  setLogStage( LOG_IO );
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscErrorCode ierr;
  FILE *of;
  char name4[80];
  PetscViewer viewer;
  if(format == OUTPUT_NODAL_FIELDS_ASCII_MATLAB ){
    sprintf(name4,"output/loadNodalFields_%d_%d.m",iMonte,iTime);
    ierr=  PetscViewerASCIIOpen(PETSC_COMM_WORLD, &name4[0],&viewer);CHKERRQ(ierr);
    ierr=  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  }else if(format == OUTPUT_NODAL_FIELDS_BINARY){
    sprintf(name4,"output/loadNodalFields_%d_%d.petscbin",iMonte,iTime);
    ierr=  PetscViewerBinaryOpen(PETSC_COMM_WORLD,&name4[0],FILE_MODE_WRITE,&viewer);
  }else if(format == OUTPUT_NODAL_FIELDS_BINARY_MATLAB){
#ifdef PETSC_HAVE_MATLAB_ENGINE
    sprintf(name4,"output/loadNodalFields_%d_%d.mat",iMonte,iTime);
    ierr=  PetscViewerMatlabOpen(PETSC_COMM_WORLD,&name4[0],FILE_MODE_WRITE,&viewer);
#else
    printf("ERROR: matlab binary output is not available\n"); abort();
#endif
  }else if(format == OUTPUT_NODAL_FIELDS_ASCII_VTK){
    sprintf(name4,"output/loadNodalFields_%d_%d.vtk",iMonte,iTime);
    ierr=  PetscViewerASCIIOpen(PETSC_COMM_WORLD, &name4[0],&viewer);CHKERRQ(ierr);    
    ierr=  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr); 
    ierr=  DMView( grid->da, viewer);CHKERRQ(ierr); 
  }else if(format == OUTPUT_NODAL_FIELDS_BINARY_VTK){
    sprintf(name4,"output/loadNodalFields_%d_%d.vts",iMonte,iTime);
    ierr= PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    //ierr= PetscViewerSetType(viewer, PETSCVIEWERVTK );CHKERRQ(ierr);
    printf("VTK support disabled\n"); abort();
    //ierr= PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);CHKERRQ(ierr);
    ierr= PetscViewerFileSetName(viewer, &name4[0] );CHKERRQ(ierr);
    //ierr= DMView( grid->da, viewer);CHKERRQ(ierr);
  }else{
    fprintf(stderr,"Invalid output format selected\n");
    abort();
  }


  ierr=  VecView( nodalFields->vx,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->vy,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->vz,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->p,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->rho,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->rhodot,viewer );CHKERRQ(ierr);
  //ierr=  VecView( nodalFields->D,viewer );CHKERRQ(ierr);
  //ierr=  VecView( nodalFields->Ddot,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->etaN,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->etaS,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->etavx,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->etavy,viewer );CHKERRQ(ierr);

  ierr=  VecView( nodalFields->kThermal,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->lastT,viewer );CHKERRQ(ierr);
  ierr=  VecView( nodalFields->Cp,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->muN,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->muS,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->muvx,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->muvy,viewer);CHKERRQ(ierr);

#ifdef TEXTURE
  ierr=  VecView( nodalFields->VPTensorB,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->VPTensorC,viewer);CHKERRQ(ierr);
#endif
  ierr=  VecView( nodalFields->soxx,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->soyy,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->sozz,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->soxy,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->soxz,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->soyz,viewer);CHKERRQ(ierr);
  /* strain rate components */
  ierr=  VecView( nodalFields->edotxx,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->edotyy,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->edotxy,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->edotxz,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->edotyz,viewer);CHKERRQ(ierr);
  /* rotation tensor */
  ierr=  VecView( nodalFields->wxy,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->wxz,viewer);CHKERRQ(ierr);
  ierr=  VecView( nodalFields->wyz,viewer);CHKERRQ(ierr);


  ierr=  VecView( nodalFields->ha,viewer);CHKERRQ(ierr);
#ifdef TEXTURE
  ierr = VecView( nodalFields->strainRateResidual,viewer);CHKERRQ(ierr);
#endif
  if(format == OUTPUT_NODAL_FIELDS_ASCII_MATLAB){
    ierr=  PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    if(!rank){
      of = fopen(name4,"a");
      fprintf(of,"elapsedTime = %e\n",elapsedTime);
      fclose(of);
    }
  }else if(format == OUTPUT_NODAL_FIELDS_BINARY){      
    /* dump grid structure here */
    Vec X;
    ierr=DMGetCoordinates( grid->da, &X );CHKERRQ(ierr);
    ierr=PetscObjectSetName( (PetscObject) X, "coordinates"); CHKERRQ(ierr);
    ierr=VecView( X ,viewer); CHKERRQ(ierr);
    printf("Preparing bag...\n");
/*     if(!rank){ */
    {
      PetscScalar *info;
      PetscBag bag;
      ierr = PetscBagCreate(PETSC_COMM_WORLD,3*sizeof(PetscScalar),&bag);CHKERRQ(ierr);
      ierr = PetscBagGetData(bag,(void **)&info);CHKERRQ(ierr);
      ierr = PetscBagRegisterScalar(bag,&info[0],(PetscScalar) grid->NX,"NX","NX");CHKERRQ(ierr);
      ierr = PetscBagRegisterScalar(bag,&info[1],(PetscScalar) grid->NY,"NY","NY");CHKERRQ(ierr);
      ierr = PetscBagRegisterScalar(bag,&info[2],elapsedTime,"elapsedTime","Elapsed Time in seconds");CHKERRQ(ierr);
      ierr = PetscBagView(bag,viewer);CHKERRQ(ierr);
      ierr = PetscBagDestroy(&bag);CHKERRQ(ierr);
    }
    ierr=PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    
  }else if(format == OUTPUT_NODAL_FIELDS_BINARY_MATLAB){
#ifdef PETSC_HAVE_MATLAB_ENGINE
    if(!rank) ierr=PetscViewerMatlabPutArray(viewer,1,1,&elapsedTime,"elapsedTime");
#else

#endif
    ierr=  PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }else{
    ierr=  PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
    
  PetscLogStagePop();
  PetscFunctionReturn(ierr);
}


PetscErrorCode saveMarkersBinary( MarkerSet *markerset, PetscInt iMonte, PetscInt iTime, PetscScalar elapsedTime){
  PetscMPIInt rank,size;
  setLogStage( LOG_IO );
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscErrorCode ierr=0;
  char fn[80];
  FILE *of=PETSC_NULL;
  PetscFunctionBegin;
  if(!rank) sprintf(fn,"output/Markers.%d.%d",iMonte,iTime);
  if(!rank) of=fopen(fn,"wb");
  
  PetscInt *nMarks;
  Marker *markers = markerset -> markers;
  PetscInt *nOuts;
  ierr = PetscMalloc( size*sizeof(PetscInt),&nOuts); CHKERRQ(ierr);
  ierr = PetscMalloc( size*sizeof(PetscInt),&nMarks);CHKERRQ(ierr);
  /* count number of markers on each processor*/
  ierr=MPI_Allgather( &(markerset->nMark), 1, MPI_INT, &nMarks[0], 1, MPI_INT, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr=MPI_Allgather( &(markerset->nOut), 1, MPI_INT, &nOuts[0], 1, MPI_INT, PETSC_COMM_WORLD);CHKERRQ(ierr);
  /* modify nMark by nOut */
  {  PetscInt i;
    for(i=0;i<size;i++) nMarks[i] -= nOuts[i];/* nMarks is number of markers to actually save */
  }
  PetscInt nMarkG, nOutG;
  ierr=MPI_Allreduce( &(markerset->nMark), &nMarkG, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr=MPI_Allreduce( &(markerset->nOut), &nOutG, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
#ifdef DEBUG  
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] writing %d markers, omitting %d out of bounds\n",rank,nMarks[rank],nOuts[rank]);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
  nMarkG-=nOutG;
  /*   fwrite( &markers->nMark, sizeof(PetscScalar),1,of); */ /* write nMark*/
  if(!rank) fwrite( &nMarkG, sizeof(PetscInt),1,of);/* number of markers*/
  if(!rank) fwrite( &elapsedTime, sizeof(PetscScalar),1,of); /*elapsed time */
#ifdef DEBUG
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d made it here\n",rank); PetscSynchronizedFlush(PETSC_COMM_WORLD);
#endif
  saveMarkerFieldI(markerset,& markers[0].cpu, &nMarks[0], of);
  saveMarkerFieldI(markerset,& markers[0].cellX, &nMarks[0], of);
  saveMarkerFieldI(markerset,& markers[0].cellY, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].X, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].Y, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].Z, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].VX, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].VY, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].VZ, &nMarks[0], of);
  saveMarkerFieldI(markerset,& markers[0].Mat, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].T, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].Tdot, &nMarks[0], of);
  //eta
  /*   fwrite( markers->eta, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].eta, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].mu, &nMarks[0], of);
#ifdef TEXTURE
  saveMarkerFieldS(markerset,& markers[0].texture.N[0][0],&nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].texture.N[0][1],&nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].texture.N[1][0],&nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].texture.N[1][1],&nMarks[0], of);
#endif
  //Damage
  /*   fwrite( markers->D, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].D, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].Ddot, &nMarks[0], of);
  //exx
  /*   fwrite( markers->exx, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].e.T11, &nMarks[0], of);
  //exy
  /*   fwrite( markers->exy, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].e.T12, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].Eii,   &nMarks[0], of);

  //sxx
  /*   fwrite( markers->sxx, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].s.T11, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].s.T22, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].s.T33, &nMarks[0], of);
  //sxy
  /*   fwrite( markers->sxy, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].s.T13, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].s.T23, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].s.T12, &nMarks[0], of);

  //p
  /*   fwrite( markers->p, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].p, &nMarks[0], of);
  //w
  /*   fwrite( markers->wxy, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->wxz, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->wyz, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->rho, sizeof(PetscScalar), markers->nMark, of); */
  saveMarkerFieldS(markerset,& markers[0].rho, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].rhodot, &nMarks[0], of);
#ifdef TEXTURE
  saveMarkerFieldS(markerset,& markers[0].residual, &nMarks[0], of);
  saveMarkerFieldS(markerset,& markers[0].strainRateResidual, &nMarks[0], of);
#endif
  if(!rank) fclose(of);
  PetscFree(nMarks);CHKERRQ(ierr);
  PetscFree(nOuts);CHKERRQ(ierr);
  PetscLogStagePop();
  PetscFunctionReturn(ierr);
} 

PetscErrorCode saveMarkerFieldS( MarkerSet *markerset, PetscScalar *fieldptr, PetscInt *nMarks, FILE *of){
  /* save a scalar-valued marker field*/
  setLogStage( LOG_IO );
  PetscErrorCode ierr;
  PetscMPIInt rank,size;
  
  const PetscInt stride = sizeof( Marker );/* spacing between adjacent elements in bytes */
  
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscFunctionBegin;
  /* calculate largest number of markers on any node */
  PetscInt nmax=0;
  PetscInt i;
  if(!rank){
    for(i=0;i<size;i++){ 
      if( nMarks[i]>nmax) nmax=nMarks[i];
    }
  }

  PetscScalar *tmpField1, *tmpField;
  /* allocate space for temporary field */
  ierr = PetscMalloc( nMarks[rank]*sizeof(PetscScalar), &tmpField);CHKERRQ(ierr);
  if(!rank){
    ierr = PetscMalloc( nmax*sizeof(PetscScalar), &tmpField1);CHKERRQ(ierr);
  }
  /* copy locally owned markers to tmpField */
  PetscInt m;
  PetscInt mm=0;
  for(m=0;m<markerset->nMark;m++){
    if( markerset->markers[m].cellX != -1 ){
      /* cast fieldptr to char * so that we can seek forward by a number of bytes*/
      /* then cast back to PetscScalar * so that we can dereference to a PetscScalar */
      tmpField[mm] = *((PetscScalar *)((unsigned char *)fieldptr + stride*m));
      mm++;
    }
  }
  if( mm != nMarks[rank] ){ printf("ERROR: number of in-bounds markers is inconsistent with nMarks[rank]\n"); abort();}

  MPI_Status status;
  if(!rank){
    /* write my portion of the field*/
    fwrite( tmpField, sizeof(PetscScalar), nMarks[0], of);/* nMark is number of markers on processor 0*/
    /* receive marker field from each other process*/
    for(i=1;i<size;i++){/* for better performance, use non-blocking MPI calls here to overlap disk access and communication*/
      MPI_Recv( tmpField1, nMarks[i], MPI_DOUBLE, i, i, PETSC_COMM_WORLD, &status);
      fwrite( tmpField1, sizeof(PetscScalar), nMarks[i],of); 
    }
  }else{
    /* rank is not zero*/
    MPI_Send( tmpField, nMarks[rank], MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD);
  }
  if(!rank){  ierr=  PetscFree(tmpField1);CHKERRQ(ierr);}
  ierr=  PetscFree(tmpField);CHKERRQ(ierr);
  PetscFunctionReturn( ierr );
}
PetscErrorCode saveMarkerFieldI(MarkerSet *markerset, PetscInt *fieldptr, PetscInt *nMarks, FILE *of){
  /* save a scalar-valued marker field*/
  PetscErrorCode ierr;
  PetscMPIInt rank,size;
  
  const PetscInt stride = sizeof( Marker );/* spacing between adjacent elements in bytes */
  
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscFunctionBegin;
  /* calculate largest number of markers on any node */
  PetscInt nmax=0;
  PetscInt i;
  if(!rank){
    for(i=0;i<size;i++){ 
      if( nMarks[i]>nmax) nmax=nMarks[i];
    }
  }

  PetscInt *tmpField1, *tmpField;
  /* allocate space for temporary field */
  ierr = PetscMalloc( nMarks[rank]*sizeof(PetscInt), &tmpField);CHKERRQ(ierr);
  if(!rank){
    ierr = PetscMalloc( nmax*sizeof(PetscInt), &tmpField1);CHKERRQ(ierr);
  }
  /* copy locally owned markers to tmpField */
  PetscInt m;
  PetscInt mm=0;
  for(m=0;m<markerset->nMark;m++){
    if( markerset->markers[m].cellX != -1 ){
      /* cast fieldptr to char * so that we can seek forward by a number of bytes*/
      /* then cast back to PetscScalar * so that we can dereference to a PetscScalar */
      tmpField[mm] = *((PetscInt *)((unsigned char *)fieldptr + stride*m));
      mm++;
    }
  }
  if( mm != nMarks[rank] ){ printf("ERROR: number of in-bounds markers is inconsistent with nMarks[rank]\n"); abort();}

  MPI_Status status;
  if(!rank){
    /* write my portion of the field*/
    fwrite( tmpField, sizeof(PetscInt), nMarks[0], of);/* nMark is number of markers on processor 0*/
    /* receive marker field from each other process*/
    for(i=1;i<size;i++){/* for better performance, use non-blocking MPI calls here to overlap disk access and communication*/
      MPI_Recv( tmpField1, nMarks[i], MPI_INT, i, i, PETSC_COMM_WORLD, &status);
      fwrite( tmpField1, sizeof(PetscInt), nMarks[i],of); 
    }
  }else{
  /* rank is not zero*/
    MPI_Send( tmpField, nMarks[rank], MPI_INT, 0, rank, PETSC_COMM_WORLD);
  }
  if(!rank){  ierr=  PetscFree(tmpField1);CHKERRQ(ierr);}
  ierr=  PetscFree(tmpField);CHKERRQ(ierr);
  PetscFunctionReturn( ierr );
}


#ifdef TEXTURE
/* write the c-axis orientations (theta,phi) */
PetscErrorCode saveTextureBinary(MarkerSet *markerset, Materials *materials, PetscInt iMonte, PetscInt iTime){

  PetscMPIInt rank, size;
  PetscErrorCode ierr=0;
  PetscInt m;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  char fnt[160];
  FILE *fh;
  sprintf(fnt,"output/Texture.%d.%d.%d",iMonte,iTime,rank);
  fh = fopen(fnt,"w");
  /*   if(!rank){ */
  fwrite(&size,sizeof(PetscInt),1,fh);
  /*   } */

  PetscInt nMarkTex =0; /* number of markers with texture */
  for(m=0;m<markerset->nMark;m++){
    if(markerset->markers[m].cellX != -1 && materials->hasTexture[(PetscInt) markerset->markers[m].Mat]){
      nMarkTex++;
    }
  }

  //  fwrite(&markerset->nMark,sizeof(markerset->nMark),1,fh);
  fwrite(&nMarkTex,sizeof(PetscInt),1,fh);

  PetscInt nt = NT;
  fwrite(&nt,sizeof(PetscInt),1,fh);

  for(m=0;m<markerset->nMark;m++){
    if(markerset->markers[m].cellX != -1 && materials->hasTexture[(PetscInt) markerset->markers[m].Mat]){
      fwrite(&markerset->markers[m].texture.ctheta[0][0][0],sizeof(PetscScalar),NT*NT*NT,fh);
      fwrite(&markerset->markers[m].texture.cphi[0][0][0],sizeof(PetscScalar),NT*NT*NT,fh);
    }
  }
  fclose(fh);


/*   MPI_File fh; */
/*   MPI_Info info; */
/*   printf("opening MPI file %s\n",fnt); */
/*   ierr=MPI_File_open(PETSC_COMM_WORLD,,MPI_MODE_CREATE,info,&fh);   CHKERRQ(ierr); */
  

/*   PetscInt nWrite = NT*NT*NT; */
/*   float temp = 134; */
/*   MPI_Status status; */

/*   MPI_File_write(fh,&temp,1,MPI_FLOAT,&status); */
/*   PetscInt count; */
/*   MPI_Get_count( &status, MPI_INT, &count ); */
/*   printf("wrote %d\n",count); */
/*   MPI_File_close(&fh); */
  PetscFunctionReturn(ierr);
}
#endif
