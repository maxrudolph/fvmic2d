
#include "fdcode.h"
#include "restart.h"
#include "markers.h"

/* This is a subroutine to restart the simulation from a saved state */

PetscErrorCode restartFromMarkers( Problem *problem, PetscInt iMonte, PetscInt iTime, PetscScalar *elapsedTime){
  MarkerSet *markerset = &(problem->markerset);
  GridData *grid = &(problem->grid);
  Materials *materials = &(problem->materials);
  Options *options = &(problem->options);
  
  /* restart the simulation by reading markers from some other timestep. Distribute to correct processors */

  PetscMPIInt rank,size;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscErrorCode ierr=0;
  char fn[80];
  FILE *of=PETSC_NULL;
  PetscFunctionBegin;
  if(!rank) printf("restart: attempting to read markers from output/Markers.%d.%d\n",iMonte,iTime);
  if(!rank) sprintf(fn,"output/Markers.%d.%d",iMonte,iTime);
  if(!rank) of=fopen(fn,"rb");
  
  PetscInt *nMarks;

  
  ierr = PetscMalloc( size*sizeof(PetscInt),&nMarks); CHKERRQ(ierr);
  PetscInt nMarkG;

  if(!rank) fread( &nMarkG, sizeof(PetscInt),1,of);/* number of markers*/
  if(!rank) fread( elapsedTime, sizeof(PetscScalar),1,of); /*elapsed time */

  if(!rank) printf("restart: reading %d global markers\n",nMarkG);

  /* read marker cpu info */
  if(!rank){
    PetscInt m;
    nMarks[0] = nMarkG/size + nMarkG%size; /* put remainder of nMarks on cpu 0 */
    for(m=1;m<size;m++){
      nMarks[ m ] = nMarkG/size;     
    }
  }
  ierr=MPI_Bcast ( nMarks, size, MPI_INT, 0,PETSC_COMM_WORLD ); /* distribute the number of markers on each cpu to every cpu */
  ierr=MPI_Bcast ( elapsedTime, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD );/* let everyone know what time it is */
  printf("restart: cpu %d restoring %d markers\n",rank,nMarks[rank]);

  ierr=allocateMarkers( problem );/* allocate markers with extra space specified by maxMarkFactor */

  Marker *markers = markerset -> markers;

  loadMarkerFieldI(& markers[0].cpu, &nMarks[0], of);

  loadMarkerFieldI(& markers[0].cellX, &nMarks[0], of);
  loadMarkerFieldI(& markers[0].cellY, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].X, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Y, &nMarks[0], of);
/*   loadMarkerFieldS(& markers[0].Z, &nMarks[0], of); */
  loadMarkerFieldS(& markers[0].VX, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].VY, &nMarks[0], of);
/*   loadMarkerFieldS(& markers[0].VZ, &nMarks[0], of); */
  loadMarkerFieldI(& markers[0].Mat, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].T, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Tdot, &nMarks[0], of);
  //eta
  /*   fwrite( markers->eta, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].eta, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].mu, &nMarks[0], of);
#ifdef TEXTURE
  loadMarkerFieldS(& markers[0].texture.N[0][0],&nMarks[0], of);
  loadMarkerFieldS(& markers[0].texture.N[0][1],&nMarks[0], of);
  loadMarkerFieldS(& markers[0].texture.N[1][0],&nMarks[0], of);
  loadMarkerFieldS(& markers[0].texture.N[1][1],&nMarks[0], of);
#endif
  //Damage
  /*   fwrite( markers->D, sizeof(PetscScalar), markers->nMark, of); */
/*   loadMarkerFieldS(& markers[0].D, &nMarks[0], of); */
/*   loadMarkerFieldS(& markers[0].Ddot, &nMarks[0], of); */
  //exx
  /*   fwrite( markers->exx, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].e.T11, &nMarks[0], of);
  //exy
  /*   fwrite( markers->exy, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].e.T12, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Eii, &nMarks[0], of);
  //sxx
  /*   fwrite( markers->sxx, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].s.T11, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].s.T22, &nMarks[0], of);
  //loadMarkerFieldS(& markers[0].s.T33, &nMarks[0], of);
  //sxy
  /*   fwrite( markers->sxy, sizeof(PetscScalar), markers->nMark, of); */
  //loadMarkerFieldS(& markers[0].s.T13, &nMarks[0], of);
  //loadMarkerFieldS(& markers[0].s.T23, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].s.T12, &nMarks[0], of);

  //p
  /*   fwrite( markers->p, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].p, &nMarks[0], of);
  //w
  /*   fwrite( markers->wxy, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->wxz, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->wyz, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->rho, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].rho, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].rhodot, &nMarks[0], of);
#ifdef TEXTURE
  loadMarkerFieldS(& markers[0].residual, &nMarks[0], of);
#endif
  
  if(!rank) fclose(of);
  /* set markerset parameters */
  markerset->nMark = nMarks[rank];
  /* nOut can be set with a call to findMarkerCells */
  if(!rank) printf("restart: assigning cells and beginning exchange\n");

#ifdef TEXTURE
  ierr = loadTextureBinary( markerset, materials, nMarks, iMonte, options->restartStep);
#endif
  findMarkerCells( markerset, grid );

  exchangeMarkers( markerset, grid, options);

  PetscFree(nMarks);CHKERRQ(ierr);
  if(!rank) printf("restart: done\n");
  PetscFunctionReturn(ierr);
}

#ifdef TEXTURE
PetscErrorCode restartFromMarkersNoTexture( MarkerSet *markerset, GridData *grid, Materials *materials, Options *options, PetscInt iMonte, PetscInt iTime, PetscScalar *elapsedTime){
  /* restart the simulation by reading markers from some other timestep. Distribute to correct processors */
  /* restart by reading markers from the last timestep (without texture) */
  PetscMPIInt rank,size;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  PetscErrorCode ierr=0;
  char fn[80];
  FILE *of=PETSC_NULL;
  PetscFunctionBegin;
  if(!rank) printf("restart: restarting from markers produced by isotropic code\n");
  if(!rank) printf("attempting to read markers from output/Markers.%d.%d\n",iMonte,iTime);
  if(!rank) sprintf(fn,"output/Markers.%d.%d",iMonte,iTime);
  if(!rank) of=fopen(fn,"rb");
  
  PetscInt *nMarks;

  
  ierr = PetscMalloc( size*sizeof(PetscInt),&nMarks); CHKERRQ(ierr);
  PetscInt nMarkG;

  if(!rank) fread( &nMarkG, sizeof(PetscInt),1,of);/* number of markers*/
  if(!rank) fread( elapsedTime, sizeof(PetscScalar),1,of); /*elapsed time */

  if(!rank) printf("reading %d global markers\n",nMarkG);

  if(!rank){
    PetscInt m;
    nMarks[0] = nMarkG/size + nMarkG%size; /* put remainder of nMarks on cpu 0 */
    for(m=1;m<size;m++){
      nMarks[ m ] = nMarkG/size;     
    }
  }
  ierr=MPI_Bcast ( nMarks, size, MPI_INT, 0,PETSC_COMM_WORLD ); /* distribute the number of markers on each cpu to every cpu */
  ierr=MPI_Bcast ( elapsedTime, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD );/* let everyone know what time it is */
  printf("cpu %d restoring %d markers\n",rank,nMarks[rank]);

  ierr=allocateMarkers( (PetscInt)(((PetscScalar) nMarks[rank])*options->maxMarkFactor) , markerset, PETSC_NULL);/* allocate markers with extra space specified by maxMarkFactor */

  Marker *markers = markerset -> markers;

  loadMarkerFieldI(& markers[0].cpu, &nMarks[0], of);

  loadMarkerFieldI(& markers[0].cellX, &nMarks[0], of);
  loadMarkerFieldI(& markers[0].cellY, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].X, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Y, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Z, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].VX, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].VY, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].VZ, &nMarks[0], of);
  loadMarkerFieldI(& markers[0].Mat, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].T, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Tdot, &nMarks[0], of);
  //eta
  /*   fwrite( markers->eta, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].eta, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].mu, &nMarks[0], of);
  /* #ifdef TEXTURE */
  /*   loadMarkerFieldS(& markers[0].texture.N[0][0],&nMarks[0], of); */
  /*   loadMarkerFieldS(& markers[0].texture.N[0][1],&nMarks[0], of); */
  /*   loadMarkerFieldS(& markers[0].texture.N[1][0],&nMarks[0], of); */
  /*   loadMarkerFieldS(& markers[0].texture.N[1][1],&nMarks[0], of); */
  /* #endif */
  //Damage
  /*   fwrite( markers->D, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].D, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Ddot, &nMarks[0], of);
  //exx
  /*   fwrite( markers->exx, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].e.T11, &nMarks[0], of);
  //exy
  /*   fwrite( markers->exy, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].e.T12, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].Eii, &nMarks[0], of);
  //sxx
  /*   fwrite( markers->sxx, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].s.T11, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].s.T22, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].s.T33, &nMarks[0], of);
  //sxy
  /*   fwrite( markers->sxy, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].s.T13, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].s.T23, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].s.T12, &nMarks[0], of);

  //p
  /*   fwrite( markers->p, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].p, &nMarks[0], of);
  //w
  /*   fwrite( markers->wxy, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->wxz, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->wyz, sizeof(PetscScalar), markers->nMark, of); */
  /*   fwrite( markers->rho, sizeof(PetscScalar), markers->nMark, of); */
  loadMarkerFieldS(& markers[0].rho, &nMarks[0], of);
  loadMarkerFieldS(& markers[0].rhodot, &nMarks[0], of);
  /* #ifdef TEXTURE */
  /*   loadMarkerFieldS(& markers[0].residual, &nMarks[0], of); */
  /* #endif */
  
  if(!rank) fclose(of);
  /* set markerset parameters */
  markerset->nMark = nMarks[rank];
  /* nOut can be set with a call to findMarkerCells */
  
  /* make the viscoplasticity tensor */
  PetscInt m;
  for(m=0;m < markerset->nMark; m++){
    markers[m].texture.N[0][0] = 2.0*markers[m].eta;
    markers[m].texture.N[1][1] = 2.0*markers[m].eta;
    markers[m].texture.N[0][1] = 0.0;
    markers[m].texture.N[1][0] = 0.0;
  }

  if(!rank) printf("restart: loading texture\n");
  ierr = loadTextureBinary( markerset, materials, nMarks, iMonte, options->restartStep);
  if(!rank) printf("restart: assigning cells and beginning exchange\n"); fflush(stdout);
  findMarkerCells( markerset, grid );
  exchangeMarkers( markerset, grid, options);
  findMarkerCells( markerset, grid );
  PetscFree(nMarks);CHKERRQ(ierr);
  if(!rank) printf("restart: done\n");
  PetscFunctionReturn(ierr);
}
#endif


#ifdef TEXTURE
/* read the c-axis orientations (theta,phi) */
PetscErrorCode loadTextureBinary(MarkerSet *markerset, Materials *materials, PetscInt *nMarks, PetscInt iMonte, PetscInt iTime){
  PetscMPIInt rank, size;
  PetscErrorCode ierr=0;
  PetscInt m;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  char fnt[160];
  FILE *fh;
  sprintf(fnt,"output/Texture.%d.%d",iMonte,iTime);  /* requires that consolidateTexture be run on output from all cpus */
  if(!rank) printf("Reading Texture from %s\n",fnt);
  
  if(!rank) fh = fopen(fnt,"r");

  /* for i=1:ncpu */
  PetscInt icpu;
  PetscInt nt = NT;
  const PetscInt NT3 = nt*nt*nt;
  PetscInt nMark = 0; /* number of markers in this file */
  PetscInt nMarkTot;
  PetscInt oldsize;/* size of previous run */
  if(!rank){
    fread(&oldsize,sizeof(PetscInt),1,fh);
    fread(&nMarkTot,sizeof(PetscInt),1,fh);
    fread(&nt,sizeof(PetscInt),1,fh);
  }
  if(!rank) printf("restoring %d textures\n",nMarkTot);
  if( nt != NT) printf("ERROR: NT=%d in texture file not consistent with compiled-in value %d\n",nt,NT);

  /* count number of local markers with texture */
  PetscInt nMarkTex = 0;
  for( m=0;m<markerset->nMark;m++){
    if(materials->hasTexture[(PetscInt) markerset->markers[m].Mat]){
      nMarkTex++;
    }
  }

  PetscInt *nMarkTexs;
  PetscMalloc( size*sizeof(PetscInt), &nMarkTexs);CHKERRQ(ierr);

  MPI_Gather( &nMarkTex, 1, MPI_INT, nMarkTexs, 1, MPI_INT, 0, PETSC_COMM_WORLD);

  if(!rank){
    /* first read locally owned textures */
    for(m=0;m<markerset->nMark;m++){
      if(materials->hasTexture[(PetscInt) markerset->markers[m].Mat]){
	fread(&(markerset->markers[m].texture.ctheta[0][0][0]),sizeof(PetscScalar),NT*NT*NT,fh);
	fread(&(markerset->markers[m].texture.cphi[0][0][0]),sizeof(PetscScalar),NT*NT*NT,fh);
      }
    }
    /* loop over other processors */
    for(icpu=1;icpu<size;icpu++){
      PetscScalar *ctheta, *cphi;

      ierr=PetscMalloc( NT3*sizeof(PetscScalar), &ctheta);CHKERRQ(ierr);
      ierr=PetscMalloc( NT3*sizeof(PetscScalar), &cphi);CHKERRQ(ierr);
      printf("restart: [%d] sending %d textures\n",rank,nMarkTexs[icpu]);
      for(m=0;m<nMarkTexs[icpu];m++){
	fread(ctheta,sizeof(PetscScalar),NT3,fh);
	fread(cphi  ,sizeof(PetscScalar),NT3,fh);
	/* NaN checking */
	PetscInt i;
	for(i=0;i<NT3;i++){
	  if(isnan(ctheta[i])){
	    printf("found nan in ctheta on disk\n");
	  }
	}
	/* send ctheta */
	MPI_Send( &ctheta[0], NT3, MPI_DOUBLE, icpu, icpu,PETSC_COMM_WORLD);
	/* send cphi */
	MPI_Send( &cphi[0], NT3, MPI_DOUBLE, icpu, icpu,PETSC_COMM_WORLD);       
      }

      //MPI_Send( ctheta, NT3*nMarkTexs[icpu], MPI_DOUBLE, icpu, 1,PETSC_COMM_WORLD);
      //MPI_Send( cphi, NT3*nMarkTexs[icpu], MPI_DOUBLE, icpu, 2,PETSC_COMM_WORLD);

      PetscFree(ctheta);
      PetscFree(cphi);
    }
  }else{
    /* receive locally owned markers */
    /* first receive token saying that it's time to start */
    MPI_Status status;
    PetscScalar test[NT][NT][NT];
    printf("restart: [%d] receiving %d textures\n",rank,nMarkTex);
    //printf("sizeof scalar array %ld, sizeof ctheta %ld, sizeof NT3 scalars %ld\n",sizeof test,sizeof markerset->markers[0].texture.ctheta,NT3*sizeof(PetscScalar));
    PetscScalar *ctheta, *cphi;
    
    ierr=PetscMalloc( NT3*sizeof(PetscScalar), &ctheta);CHKERRQ(ierr);
    ierr=PetscMalloc( NT3*sizeof(PetscScalar), &cphi);CHKERRQ(ierr);
    
    
    PetscInt mm=0;
    for(m=0;m<markerset->nMark;m++){
      if(materials->hasTexture[(PetscInt) (markerset->markers[m].Mat)]){ 
	MPI_Recv(ctheta, NT3, MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD, &status); 
	MPI_Recv(cphi  , NT3, MPI_DOUBLE, 0, rank, PETSC_COMM_WORLD, &status);
	PetscInt nn=0;
	PetscInt i,j,k;
	for(i=0;i<NT;i++){
	  for(j=0;j<NT;j++){
	    for(k=0;k<NT;k++){
	      markerset->markers[mm].texture.ctheta[i][j][k] = ctheta[nn];
	      markerset->markers[mm].texture.cphi[i][j][k] = cphi[nn];
	      nn++;
	    }
	  }
	}
	mm++;
      }
    }
    /* check for NaN in ctheta */
    for(m=0;m<markerset->nMark;m++){
      if(materials->hasTexture[(PetscInt) (markerset->markers[m].Mat)]){ 
	PetscInt i,j,k;
	for(i=0;i<NT;i++){
	  for(j=0;j<NT;j++){
	    for(k=0;k<NT;k++){
	      if(isnan(markerset->markers[m].texture.ctheta[i][j][k])){
		printf("found nan in texture received from process 0\n");
	      }
	    }
	  }
	}
      }
    }
    PetscFree(ctheta);
    PetscFree(cphi);
    

    if( mm != nMarkTex){
      printf("[%d] restart Error: received %d markers, expected %d\n",rank,mm,nMarkTex);
    }else{
      printf("restart: [%d] done receiving textures\n",rank);
    }
  }


  if( !rank) fclose(fh);


  PetscFunctionReturn(ierr);
}
#endif




PetscErrorCode loadMarkerFieldI( PetscInt *fieldptr, PetscInt *nMarks, FILE *of){
  /* load a scalar-valued marker field*/
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

  PetscInt *tmpField0, *tmpField;
  /* allocate space for temporary field */
  ierr = PetscMalloc( nMarks[rank]*sizeof(PetscInt), &tmpField);CHKERRQ(ierr);
  if(!rank){
    ierr = PetscMalloc( nmax*sizeof(PetscInt), &tmpField0); CHKERRQ(ierr);
  }

  MPI_Status status;
  if(!rank){
    /* write my portion of the field*/
    fread( tmpField, sizeof(PetscInt), nMarks[0], of);/* nMarks[0] is number of markers on processor 0*/
    /* receive marker field from each other process*/
    for(i=1;i<size;i++){/* for better performance, use non-blocking MPI calls here to overlap disk access and communication-NOT IMPLEMENTED YET*/
      fread( tmpField0, sizeof(PetscInt), nMarks[i],of); 
      MPI_Send( tmpField0, nMarks[i], MPI_INT, i, i, PETSC_COMM_WORLD );
    }
  }else{
    /* rank is not zero*/

    MPI_Recv( tmpField, nMarks[rank], MPI_INT, 0, rank,PETSC_COMM_WORLD, &status );
  }

  /* copy tmpField to appropriate location in markers */
  PetscInt m;
  for(m=0;m<nMarks[rank];m++){
    /* cast fieldptr to char * so that we can seek forward by a number of bytes*/
    /* then cast back to PetscScalar * so that we can dereference to a PetscScalar */
    *((PetscInt *)((unsigned char *)fieldptr + stride*m))=tmpField[m];
  }


  if(!rank){  ierr=  PetscFree(tmpField0);CHKERRQ(ierr);}
  ierr=  PetscFree(tmpField);CHKERRQ(ierr);
  PetscFunctionReturn( ierr );
}

PetscErrorCode loadMarkerFieldS( PetscScalar *fieldptr, PetscInt *nMarks, FILE *of){
  /* load a scalar-valued marker field*/
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

  PetscScalar *tmpField0, *tmpField;
  /* allocate space for temporary field */
  ierr = PetscMalloc( nMarks[rank]*sizeof(PetscScalar), &tmpField);CHKERRQ(ierr);
  if(!rank){
    ierr = PetscMalloc( nmax*sizeof(PetscScalar), &tmpField0); CHKERRQ(ierr);
  }

  MPI_Status status;
  if(!rank){
    /* write my portion of the field*/
    fread( tmpField, sizeof(PetscScalar), nMarks[0], of);/* nMarks[0] is number of markers on processor 0*/
    /* receive marker field from each other process*/
    for(i=1;i<size;i++){/* for better performance, use non-blocking MPI calls here to overlap disk access and communication-NOT IMPLEMENTED YET*/
      fread( tmpField0, sizeof(PetscScalar), nMarks[i],of); 
      MPI_Send( tmpField0, nMarks[i], MPI_DOUBLE, i, i, PETSC_COMM_WORLD );
    }
  }else{
    /* rank is not zero*/
    MPI_Recv( tmpField, nMarks[rank], MPI_DOUBLE, 0, rank,PETSC_COMM_WORLD, &status );
  }

  /* copy tmpField to appropriate location in markers */
  PetscInt m;
  for(m=0;m<nMarks[rank];m++){
    /* cast fieldptr to char * so that we can seek forward by a number of bytes*/
    /* then cast back to PetscScalar * so that we can dereference to a PetscScalar */
    *((PetscScalar *)((unsigned char *)fieldptr + stride*m))=tmpField[m];
  }


  if(!rank){  ierr=  PetscFree(tmpField0);CHKERRQ(ierr);}
  ierr=  PetscFree(tmpField);CHKERRQ(ierr);
  PetscFunctionReturn( ierr );
}
