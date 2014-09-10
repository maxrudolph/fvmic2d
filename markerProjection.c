/* This file contains routines to project quantities from the markers onto the nodes. It is designed with modularity in mind so that different projection schemes (i.e. arithmetic, geometric, harmonic averaging) might be tried*/

#include "fdcode.h"
#include "markerProjection.h"
#include "profile.h"

void findCellStag2D( Marker *, GridData *, PetscInt , PetscInt , PetscInt *, PetscInt *);
PetscErrorCode interpolateCellCentersBasicNodes( GridData *, Vec , Vec , PROJECTION_METHOD );

/* project all marker fields to nodes */
PetscErrorCode projectMarkersNodesAll2( MarkerSet *markerset, GridData *grid, NodalFields *nodalFields, Materials *materials, Options *option ){
  /* project all marker fields to nodes */
  /* this should just be a list of commands like */
  /* projectMarkerFieldToNodes( markerField, nodalField, STAGX, STAGY, STAGZ, PROJECTION_METHOD, WEIGHTS ) */
  setLogStage( LOG_PROJECT_M_TO_N );
  PetscErrorCode ierr =0;
  PetscFunctionBegin;
  Marker *markers = markerset->markers;
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].rho, nodalFields->rho, 0, 0, ARITHMETIC);
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].rhodot, nodalFields->rhodot, 0, 0, ARITHMETIC);
  ierr = VecZeroEntries( nodalFields->rhodot );
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].eta, nodalFields->etaN, -1, -1, ARITHMETIC );  
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].eta, nodalFields->etaS,  0,  0, ARITHMETIC );
  //  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].eta, nodalFields->etavx, 0,  1, ARITHMETIC );
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].eta, nodalFields->etavy, 1,  0, ARITHMETIC );
#ifdef TEXTURE
  printf("NOT IMPLEMENTED\n");fflush(stdout);
  abort();
#endif
  // muN     
  //projectMarkerFieldToNodes( grid, markerset, materials->materialMu, PETSC_NULL, nodalFields->muN, -1,  -1, HARMONIC);
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].mu, nodalFields->muN,  -1,  -1, HARMONIC);
  // muS     
  //projectMarkerFieldToNodes( grid, markerset, materials->materialMu, PETSC_NULL, nodalFields->muS,  0,   0, HARMONIC);
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].mu, nodalFields->muS,  0,  0, HARMONIC);
  ierr = VecSet( nodalFields->muS, 1e99 );CHKERRQ(ierr);
  ierr = VecSet( nodalFields->muN, 1e99 );CHKERRQ(ierr);
  
  /* muvx     */
  //projectMarkerFieldToNodes( grid, markerset, materials->materialMu, PETSC_NULL, nodalFields->muvx,  0,  1, HARMONIC);
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].mu, nodalFields->muvx,  0,  1, HARMONIC);
  /* muvy     */
  //projectMarkerFieldToNodes( grid, markerset, materials->materialMu, PETSC_NULL, nodalFields->muvy,  1,  0, HARMONIC);
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].mu, nodalFields->muvy,  1,  0, HARMONIC);
  //  kThermal
  projectMarkerFieldToNodes( grid, markerset, materials->materialkThermal, PETSC_NULL, nodalFields->kThermal,  0,   0, ARITHMETIC);
  //lastT   
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].T, nodalFields->lastT,  0,   0, ARITHMETIC);
  //Cp      
  projectMarkerFieldToNodes( grid, markerset, materials->materialCp, PETSC_NULL, nodalFields->Cp,  0,   0, ARITHMETIC);
  /*   soxx     */
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].s.T11, nodalFields->soxx,  -1,   -1, ARITHMETIC_LOCAL);
  /*   soyy     */
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].s.T22, nodalFields->soyy,  -1,   -1, ARITHMETIC_LOCAL);
  /*   sozz     */
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].s.T33, nodalFields->sozz,  -1,   -1, ARITHMETIC_LOCAL);
  /*   soxy     */
  projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].s.T12, nodalFields->soxy,   0,    0, ARITHMETIC_LOCAL);
  /*   soxz     */
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].s.T13, nodalFields->soxz,   0,    0, ARITHMETIC_LOCAL);
  /*   soyz     */
  //projectMarkerFieldToNodes( grid, markerset, PETSC_NULL, &markers[0].s.T23, nodalFields->soyz,   0,    0, ARITHMETIC_LOCAL);
  /*   ha       */
  //  ierr= interpolateCellCentersBasicNodes( grid, nodalFields->etaN , nodalFields->etaS , HARMONIC );

  PetscLogStagePop();
  PetscFunctionReturn(ierr);
}


PetscErrorCode projectMarkerFieldToNodes(GridData *grid, MarkerSet *markerset, PetscScalar *materialProperty, PetscScalar *markerfield, Vec nodalField, const PetscInt stagx, const PetscInt stagy, const PROJECTION_METHOD pm){
  /* if called with materialProperty non-null, this routine should project a material property instead of a property defined on the marker */
  PetscErrorCode ierr =0;
  PetscFunctionBegin;

  /* make global array for weights */
  Vec wtg, wtl, nfl;
  ierr = DMGetGlobalVector( grid->da, &wtg ); CHKERRQ(ierr);
  ierr = VecZeroEntries( wtg); CHKERRQ(ierr);
  ierr = VecZeroEntries( nodalField ); CHKERRQ(ierr);
  ierr = DMGetLocalVector( grid->da, &wtl ); CHKERRQ(ierr);
  ierr = VecDuplicate( wtl, &nfl); CHKERRQ(ierr);
  /* get local portion of global weights and nodal field */
  ierr=DMGlobalToLocalBegin(grid->da,wtg,INSERT_VALUES,wtl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,wtg,INSERT_VALUES,wtl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalField,INSERT_VALUES,nfl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalField,INSERT_VALUES,nfl);CHKERRQ(ierr);
  
  PetscScalar **wt, **nf;
  ierr=DMDAVecGetArray(grid->da,wtl,&wt);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,nfl,&nf);CHKERRQ(ierr);

  /* extract the markerfield from the markers */
  Marker *markers = markerset->markers;
  PetscScalar *mf;
  ierr = PetscMalloc( markerset->nMark*sizeof(PetscScalar), &mf); CHKERRQ(ierr);
  PetscInt m;
  /* case 1: using a real marker field */
  if( materialProperty == PETSC_NULL){
    PetscInt stride = sizeof( Marker );
    for( m=0; m<markerset->nMark;m++){
      mf[m] = *(PetscScalar *)(((char *) markerfield) + stride*m);
      /* the right way to do this is with memcpy */
      //ierr = PetscMemcpy( &mf[m], markerfield + m*stride, sizeof(PetscScalar) ); CHKERRQ(ierr);

    }
  }
  /* case 2: using something defined based on the marker material properties */
  if( materialProperty != PETSC_NULL){
    for( m=0; m<markerset->nMark;m++){
      mf[m] = materialProperty[ markers[m].Mat ];
    }
  }

  for( m=0; m<markerset->nMark;m++){   /* loop over markers */
    if( markers[m].cellX != -1 ){
      /* do the projection and calculate weights */
      /* Option 1 - arithmetic projection */
      
      /* calculate correct cell */
      /* if stagx = -1, that indicates that this node is offset in the negative direction i.e. as it is for the cell centers throughout the code */
      /* get cellX, cellY */
      PetscInt cellX, cellY;
      findCellStag2D( &markers[m], grid, stagx, stagy, &cellX, &cellY);
      /* calculate fractional distance from upper left corner */
      PetscScalar mdx, mdy;
      switch(stagx){
      case -1:
 	mdx = (markers[m].X - grid->xc[cellX])/(grid->xc[cellX+1]-grid->xc[cellX]); 
	break;
      case 0:
	mdx = (markers[m].X - grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
	break;
      case 1:
	mdx = (markers[m].X - grid->xc[cellX+1])/(grid->xc[ cellX + 2 ]-grid->xc[cellX+1]);
	break;
      default:
	printf("ERROR: Invalid value for stagx\n"); abort();
	break;
      }
      switch(stagy){
      case -1:
	mdy = (markers[m].Y - grid->yc[cellY])/(grid->yc[cellY+1]-grid->yc[cellY]);
	break;
      case 0:
	mdy = (markers[m].Y - grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
	break;
      case 1:
	mdy = (markers[m].Y - grid->yc[cellY+1])/(grid->yc[ cellY + 2 ]-grid->yc[cellY+1]);
	break;
      default:
	printf("ERROR: Invalid value for stagx\n"); abort();
	break;
      }

/*       if( stagx == 1 ){ */
/* 	mdx = (markers[m].X - grid->xc[cellX+1])/(grid->xc[ cellX + 2 ]-grid->xc[cellX+1]); */
/*       }else if( stagx == -1){ */
/* 	mdx = (markers[m].X - grid->xc[cellX])/(grid->xc[cellX+1]-grid->xc[cellX]); */
/*       }else{ */
/* 	mdx = (markers[m].X - grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]); */
/*       } */
/*       if( stagy == 1 ){ */
/* 	mdy = (markers[m].Y - grid->yc[cellY+1])/(grid->yc[ cellY + 2 ]-grid->yc[cellY+1]); */
/*       }else if( stagy == -1){ */
/* 	mdy = (markers[m].Y - grid->yc[cellY])/(grid->yc[cellY+1]-grid->yc[cellY]); */
/*       }else{ */
/* 	mdy = (markers[m].Y - grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]); */
/*       } */
      /* update nodal field and weights */
      switch(pm){

      case ARITHMETIC :
	/*       if( pm == ARITHMETIC ){ */
	{
	const PetscScalar thiswt1 = (1.0-mdx)*(1.0-mdy);
	nf[cellY][cellX] += thiswt1 * mf[m];
	wt[cellY][cellX] += thiswt1;
	const PetscScalar thiswt2 = mdx*(1.0-mdy);
	nf[cellY][cellX+1] += thiswt2 * mf[m];
	wt[cellY][cellX+1] += thiswt2;
	const PetscScalar thiswt3 = (1.0-mdx)*(mdy);
	nf[cellY+1][cellX] += thiswt3 * mf[m];
	wt[cellY+1][cellX] += thiswt3;
	const PetscScalar thiswt4 = mdx*mdy;
	nf[cellY+1][cellX+1] += thiswt4 * mf[m];
	wt[cellY+1][cellX+1] += thiswt4;
	}
	break;
      case ARITHMETIC_LOCAL :
	/*       }else if(pm == ARITHMETIC_LOCAL ){ */
	{
	  PetscScalar thiswt=1.0;
	  if( mdx < 0 ) printf("stagx = %d, stagy = %d, mdx=%e\n",stagx,stagy,mdx);
	  if(mdx > 1.0) printf("mdx=%e\n",mdx);
	  if( mdx >= 0.5 ){ 
	    cellX++;
	    thiswt *= mdx;
	  }else{
	    thiswt *= (1.0-mdx);
	  }
	  if( mdy >= 0.5 ){
	    cellY++;
	    thiswt *= mdy;
	  }else{
	    thiswt *= (1.0-mdy);
	  }
	  nf[cellY][cellX] += thiswt * mf[m];
	  wt[cellY][cellX] += thiswt;
	}
	break;
      case HARMONIC :
	{
	/*       }else if(pm == HARMONIC ){ */
	if( mdx >= 0.5) cellX++;
	if( mdy >= 0.5) cellY++;
	
	nf[cellY][cellX] += 1.0/mf[m];
	wt[cellY][cellX] += 1.0;
	}
	break;
/*       } else { */
      default:
	printf("PROJECTION METHOD NOT IMPLEMENTED\n");
	abort();
      }
      
    }
  }
  ierr=DMDAVecRestoreArray(grid->da,wtl,&wt);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,nfl,&nf);CHKERRQ(ierr);
  /* scatter local to global */
  ierr = DMLocalToGlobalBegin(grid->da,nfl,ADD_VALUES,nodalField);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(  grid->da,nfl,ADD_VALUES,nodalField);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,wtl,ADD_VALUES,wtg);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(  grid->da,wtl,ADD_VALUES,wtg);CHKERRQ(ierr);
  if( pm == ARITHMETIC || pm == ARITHMETIC_LOCAL){
    ierr = VecPointwiseDivide( nodalField, nodalField, wtg);CHKERRQ(ierr);
  } else if( pm == HARMONIC ){
    ierr = VecPointwiseDivide( nodalField, wtg, nodalField); CHKERRQ(ierr);
  } else {
    printf("PROJECTION METHOD NOT IMPLEMENTED\n");
    abort();
  }
  PetscFree(mf);
  ierr = DMRestoreGlobalVector( grid->da, &wtg );CHKERRQ(ierr);
  ierr = DMRestoreLocalVector( grid->da, &wtl );CHKERRQ(ierr);
  ierr = DMRestoreLocalVector( grid->da, &nfl );CHKERRQ(ierr);
  PetscFunctionReturn(ierr); 
}

void findCellStag2D( Marker *marker, GridData *grid, PetscInt stagx, PetscInt stagy, PetscInt *cellX, PetscInt *cellY){
  /* calculate the upper left cell indices */
  cellX[0] = marker->cellX;
  cellY[0] = marker->cellY;
  if( stagx == 0 && stagy == 0) return;
  if( stagx == 0 ){
  } else if( stagx == 1){
    if( marker->X < grid->xc[ marker->cellX+1 ] ) cellX[0] --;
  } else if( stagx == -1){
    if( marker->X >= grid->xc[ marker->cellX+1] ) cellX[0] ++;
  }
  if( stagy == 0 ){
  } else if( stagy == 1){
    if( marker->Y < grid->yc[ marker->cellY+1] ) cellY[0] --;
  } else if( stagy == -1){
    if( marker->Y >= grid->yc[ marker->cellY+1] ) cellY[0] ++;
  }  
}

PetscErrorCode interpolateCellCentersBasicNodes( GridData *grid, Vec from, Vec to, PROJECTION_METHOD pm){
  
  PetscErrorCode ierr;
  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  Vec fl, tl;

  ierr = DMGetLocalVector(grid->da, &fl );
  ierr = DMGetLocalVector(grid->da, &tl );

  ierr=DMGlobalToLocalBegin(grid->da,from,INSERT_VALUES,fl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,  from,INSERT_VALUES,fl);CHKERRQ(ierr);  
  ierr=DMGlobalToLocalBegin(grid->da,to,INSERT_VALUES,tl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,  to,INSERT_VALUES,tl);CHKERRQ(ierr);  
  PetscScalar **f, **t;  
  ierr=DMDAVecGetArray(grid->da,fl,&f);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,tl,&t);CHKERRQ(ierr);

  PetscInt ix, jy;
  for( jy=y;jy>y+n;jy++){
    for( ix=x;ix<x+m;ix++){
      if( grid->xperiodic || ( ix>0 && ix<grid->NX-1)){
	if( jy>0 && jy < grid->NY-1 ){
	  if( pm == HARMONIC ){
	    t[jy][ix] = 4.0/(1.0/f[jy][ix] + 1.0/f[jy][ix+1] + 1.0/f[jy+1][ix] + 1.0/f[jy+1][ix+1]);
	  }else if( pm == ARITHMETIC){
	    t[jy][ix] = (f[jy][ix] + f[jy][ix+1] + f[jy+1][ix] + f[jy+1][ix+1])/4.0;
	  }else{
	    printf("FEATURE NOT IMPLEMENTED\n");
	    abort();
	  }	    
	}
      }
    }
  }
  ierr=DMDAVecRestoreArray(grid->da,fl,&f);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,tl,&t);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,tl,INSERT_VALUES,to);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(  grid->da,tl,INSERT_VALUES,to);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(grid->da, &fl);
  ierr = DMRestoreLocalVector(grid->da, &tl);

  PetscFunctionReturn(ierr);
}
