/* Subroutines pertaining to markers*/
#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include "petscsys.h"
#include "fdcode.h"
#include "markers.h"
#include "viscosity.h"
#include "math.h"
#include "profile.h"
#include "updateDamageViscosity.h"
#include "benchmarkInitialConditions.h"
//#define DEBUG

PetscErrorCode allocateMarkers( PetscInt nMark, MarkerSet *markerset, Options *options){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = PetscMalloc( nMark*sizeof(Marker), &(markerset->markers));CHKERRQ(ierr);
  
#ifdef DEBUG
  printf("Allocated %d markers\n",nMark);
#endif
  /*initialization of fabric would go here*/
  markerset -> maxMark = nMark;
  resetMarkers( markerset, options);
  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyMarkers( MarkerSet *markerset ,Options *options){
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = PetscFree(markerset->markers);CHKERRQ(ierr);
  markerset->markers = PETSC_NULL;
  markerset->nMark = 0;
  markerset->maxMark = 0;
  markerset->nOut = 0;
  PetscFunctionReturn(ierr);
}

void findMarkerCells( MarkerSet *markerset, GridData *grid){
  /* go through list of local markers and figure out which cell I belong to using a bisection method*/
  const PetscInt NX =grid->NX;
  const PetscInt NY =grid->NY;
  const PetscScalar LX =grid->LX;
  const PetscScalar LY =grid->LY;
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt nout=0;
  PetscInt m;

  Marker *markers = markerset->markers;

  for(m=0;m<markerset->nMark;m++){
    if( markers[m].X >= 0.0 && markers[m].Y >= 0.0 && markers[m].X <= LX && markers[m].Y <= LY){
      /* marker is in bounds*/
      PetscInt idxmin=0;
      PetscInt idxmax= grid->xperiodic ? NX : NX-1;
      while( (idxmax - idxmin) >1){
	PetscInt idxc = (idxmax+idxmin)/2;
	if(idxc == idxmin){idxc++;}else if(idxc == idxmax){idxc--;}
	if( markers[m].X >= grid->x[idxc]){ 
	  idxmin=idxc;
	}else{
	  idxmax=idxc;
	}
	//	printf("marker%d X = %e, idxmin %d idxmax %d\n",m,markers->X[m],idxmin,idxmax);
      }
      markers[m].cellX = idxmin; 
      idxmin=0; 
      idxmax=NY-1; 
      while( idxmax - idxmin>1){
	PetscInt idxc = (idxmax+idxmin)/2;
	if(idxc == idxmin){idxc++;}else if(idxc == idxmax){idxc--;}
	if( markers[m].Y >= grid->y[idxc]){ 
	  idxmin=idxc;
	}else{
	  idxmax=idxc;
	}
	//	printf("marker%d X = %e, idxmin %d idxmax %d\n",m,markers->X[m],idxmin,idxmax);
      }
      markers[m].cellY = idxmin;

      idxmax= grid->xperiodic ? NX-1 : NX-2;
      if( markers[m].cellX < 0 || markers[m].cellX > idxmax || markers[m].cellY < 0 || markers[m].cellY > NY-2){
	printf("cell out of bounds marker %d, cellX - %d, cellY = %d\n",m,markers[m].cellX, markers[m].cellY);
	abort();
      }
      


    } else {
      markers[m].cellX = -1;
      markers[m].cellY = -1;
      nout++;
    }
    markers[m].cpu = rank;
  }
  markerset->nOut=nout;
}

PetscErrorCode exchangeMarkers( MarkerSet *markerset, GridData *grid,Options *options){

  PetscErrorCode ierr;
  PetscMPIInt rank,size;
  PetscFunctionBegin;
  setLogStage( LOG_MARK_EXCH );

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /* this routine takes markers that have strayed from the proper cpu and sends them where they belong*/
  /* ignore markers with -1 cellI or cellJ - these are out of bounds */
  /* go through all of my markers. are there any that I no longer should have? */
  /* define local max and min x and y cells*/
  int x,y,z,m,n,p;
  ierr=DMDAGetCorners(grid->da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
  /* if x+m == NX and the grid is periodic, this processor is in the last column of the domain and therefore owns cell NX-1 */
 
  Marker *markers = markerset->markers;

  /* now see if nDiscard > 0*/
  PetscInt nDiscard=0;
  {
    PetscInt m1;
    for(m1=0;m1<markerset->nMark;m1++){
      if(markers[m1].cellX != -1){
	if(markers[m1].cellX<x || markers[m1].cellY <y || markers[m1].cellX>=(x+m) || markers[m1].cellY>=(y+n) ){
	  /* this marker needs to be discarded */
	  //printf("[%d] discards marker in cell x %d y %d bc (x,y,m,n) = (%d,%d,%d,%d)\n",rank,markers[m1].cellX,markers[m1].cellY,x,y,m,n);fflush(stdout);
	  nDiscard++;
	}
      }/* end if in domain*/
    }
  }
#ifdef debug1
  /* calculate the number of markers before exchange */
  PetscInt nGlobal;
  MPI_Allreduce( &markerset->nMark,&nGlobal, 1, MPI_INT , MPI_SUM, PETSC_COMM_WORLD );
  if(!rank) printf("Global marker count before exchange %d\n",nGlobal);
#endif

#ifdef DEBUG
  printf("[%d] found %d markers to move to another cpu\n",rank,nDiscard);fflush(stdout);
#endif
  /* count the ones that I should no longer have and make a list containing their indices*/
  Marker *msend, *mrecv;/* send and receive buffers for exchanging marker info */
  PetscInt nsend,nrecv;
  ierr = PetscMalloc( nDiscard*sizeof(Marker), &msend);
  /* copy all marker fields to a single vector. */
  {
    PetscInt m1;
    PetscInt idx=0;
    for(m1=0;m1<markerset->nMark;m1++){
      if(markers[m1].cellX != -1){
	if(markers[m1].cellX<x || markers[m1].cellY <y || markers[m1].cellX>=(x+m) || markers[m1].cellY>=(y+n) ){
	  /* copy all fields to ougoing markers...*/
	  msend[idx] = markers[m1];

	  /*Flag each outgoing marker -1 and move it outside domain.*/
	  markers[m1].cellX = -1;
	  markers[m1].cellY = -1;
	  markers[m1].X = -1;
	  markers[m1].Y = -1;
	  markerset->nOut++;/* increase count of out-of-domain (i.e. unused) markers*/
	  idx++;
	}
      }
    }
  }
  nsend = nDiscard;
  /* do size-1 times */
  /* tell node to right how many markers I am rejecting, get this number from the node to the left*/
  MPI_Status status;
  PetscInt iex;
  for(iex=1;iex<size;iex++){
    /*if(!rank) {printf("exchange %d of %d\n",iex,size-1);fflush(stdout);}*/
    ierr = MPI_Sendrecv( &nsend, 1, MPI_INT, (rank+1)%size,rank, &nrecv, 1, MPI_INT, (size+rank-1)%size, (size+rank-1)%size, PETSC_COMM_WORLD,&status);CHKERRQ(ierr);
    /*PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] receiving %d markers from cpu %d\n",rank,nrecv,(size+rank-1)%size);*/
    /*PetscSynchronizedFlush(PETSC_COMM_WORLD);*/
    /* allocate an array to hold markers from node to left*/
    ierr = PetscMalloc( nrecv*sizeof(Marker), &mrecv);CHKERRQ(ierr);
    /* pass to right and receive from left, or come up with a more intelligent algorithm to just pass/receive from 8 neighbors*/
    { PetscInt ssend=nsend; PetscInt srecv=nrecv;
      ierr = MPI_Sendrecv( msend, ssend*sizeof(Marker), MPI_BYTE, (rank+1)%size,rank,mrecv, srecv*sizeof(Marker), MPI_BYTE, (size+rank-1)%size,(size+rank-1)%size, PETSC_COMM_WORLD,&status);CHKERRQ(ierr);
    }
    /* look through list of markers that I have received. count number that I will take*/
    PetscInt i;
    PetscInt nTake=0;
    {
      for(i=0;i<nrecv;i++){

	if( mrecv[i].cellX >=x && mrecv[i].cellY >=y && mrecv[i].cellX < (x+m) && mrecv[i].cellY <(y+n)){
	  //printf("[%d] draws marker in cell x %d y %d bc (x,y,m,n) = (%d,%d,%d,%d)\n",rank,mrecv[i].cellX,mrecv[i].cellY,x,y,m,n);fflush(stdout);
	  /* marker is on this cpu. count it*/
	  nTake++;
	}
      }
    }
    /* make sure that there is enough room to take all of the markers available*/
    if( markerset->maxMark-(markerset->nMark-markerset->nOut + nTake) >0){
#ifdef debug1
      if(nTake>0){ printf("%d has room to take requested %d markers\n",rank,nTake);}
#endif
      /* take the markers... first use empty (-1 flagged) spots, then add markers*/
      i=0;
      /*       for(i=0;i<nTake;i++){ */
      while(nTake>0){

	while(  mrecv[i].cellX <x || mrecv[i].cellY <y || mrecv[i].cellX >= (x+m) || mrecv[i].cellY >= (y+n)){
	  i++;
	}
	PetscInt m1=0;
	if(markerset->nOut>0){
	  while(markers[m1].cellX != -1) m1++; /* m1 should now hold first -1 flagged entry*/
	  markerset->nOut--;
	} else{
	  m1=markerset->nMark;
	  markerset->nMark++;
	}
	markers[m1] = mrecv[i];
	nTake--;
	i++;
      }/* end loop over nTake */
    }else{/* not enough room to take the markers */
      printf("%d out of room -cannot take %d markers\n",rank,nTake);
    }


    /* destroy msend*/
    ierr = PetscFree(msend);CHKERRQ(ierr);
    msend = mrecv;/* get ready to pass markers to next cpu*/
    nsend = nrecv;
  }
  ierr = PetscFree(msend);CHKERRQ(ierr);

#ifdef debug1
  /* calculate the number of markers before exchange */
  MPI_Allreduce( &markerset->nMark,&nGlobal, 1, MPI_INT , MPI_SUM, PETSC_COMM_WORLD );
  if(!rank) printf("Global marker count after exchange %d\n",nGlobal);
#endif

  /* end loop */
  PetscLogStagePop();
  PetscFunctionReturn(ierr);

}


void resetMarkers( MarkerSet *markerset,Options *options){
/*   PetscErrorCode ierr=0; */
  PetscInt nMark = markerset->maxMark;
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<nMark;m++){
    markers[m].X = -1;
    markers[m].Y = -1;
    markers[m].Z = -1;
    markers[m].VX = 0;
    markers[m].VY = 0;
    markers[m].VZ = 0;
    markers[m].p = 0.0;
    markers[m].Mat = -1;
    markers[m].T = 0.0;
    markers[m].Tdot = 0.0;
    markers[m].eta =0.0;
    markers[m].s.T11 = 0.0;
    markers[m].s.T22 = 0.0;
    markers[m].s.T33 = 0.0;
    markers[m].s.T13 = 0.0;
    markers[m].s.T23 = 0.0;
    markers[m].s.T12 = 0.0;    
    markers[m].e.T11 = 0.0;
    markers[m].e.T22 = 0.0;
    markers[m].e.T33 = 0.0;
    markers[m].e.T13 = 0.0;
    markers[m].e.T23 = 0.0;
    markers[m].e.T12 = 0.0;

#ifndef TEXTURE
    markers[m].eta = 0.0;
#else
    markers[m].mu = 0.0;
    PetscInt i,j;
    markers[m].texture.N[0][0]=0.0;
    markers[m].texture.N[0][1]=0.0;
    markers[m].texture.N[1][0]=0.0;
    markers[m].texture.N[1][1]=0.0;
    for(i=0;i<6;i++){
      for(j=0;j<6;j++){
	markers[m].texture.M[i][j]=0.0;
      }
    }
#endif
    
    markers[m].p = 0.0;
    /* stress */
    markers[m].s.T11 = 0.0;
    markers[m].s.T22 = 0.0;
    markers[m].s.T33 = 0.0;
    markers[m].s.T23 = 0.0;
    markers[m].s.T13 = 0.0;
    markers[m].s.T12 = 0.0; 
    
    markers[m].rho = 0.0;
    markers[m].rhodot = 0.0;
    /* strain rate */
    markers[m].e.T11 = 0.0;
    markers[m].e.T22 = 0.0;
    markers[m].e.T33 = 0.0;
    markers[m].e.T23 = 0.0;
    markers[m].e.T13 = 0.0;
    markers[m].e.T12 = 0.0;
    /* total strain */
    markers[m].E.T11 = 0.0; 
    markers[m].E.T22 = 0.0; 
    markers[m].E.T33 = 0.0; 
    markers[m].E.T23 = 0.0; 
    markers[m].E.T13 = 0.0; 
    markers[m].E.T12 = 0.0; 
    markers[m].Eii = 0.0;
    /* rotation rate */
    markers[m].wxz = 0.0;
    markers[m].wxy = 0.0;
    markers[m].wyz = 0.0;
    
    markers[m].D = 0.0;
    markers[m].Ddot = 0.0;
    markers[m].isYielding = 0;
  }
}

PetscErrorCode distributeMarkersUniformInDomain( MarkerSet *markerset, Options *options, GridData *grid){
  PetscErrorCode ierr;
  PetscRandom r;
  PetscInt NMX=options->NMX;
  PetscInt NMY=options->NMY;
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  Marker *markers = markerset->markers;

  PetscFunctionBegin;
  ierr=PetscRandomCreate( PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  ierr=PetscRandomSetType(r,PETSCRAND);CHKERRQ(ierr);
  ierr=PetscRandomSetSeed(r, (unsigned long) 0);CHKERRQ(ierr);
  
  /* get local ownership range*/
  int x,y,z,m,n,p;
  ierr=DMDAGetCorners(grid->da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
  /* x,y,z are coordinate indices of upper left node in gerya notation (lower left in petsc notation) */
  
  //PetscInt i,j;
  PetscInt idx=0;/* marker counter*/

  /* x+m should not exceed NX-1 */

  if(x+m >= NX && !grid->xperiodic){ m--;}else{}
  if(x+m > NX+1 && grid->xperiodic){ m--;}else{}
  if(y+n >= NY){ n--;}else{}

  //  PetscScalar mdx = (grid->x[x+m]-grid->x[x])/(NMX*m);
  PetscScalar mdx = grid->LX/(NMX*(NX-1));
  PetscScalar mdy = grid->LY/(NMY*(NY-1));
  //  PetscScalar mdy = (grid->y[y+n]-grid->y[y])/(NMY*n);
  

  /*       PetscInt m1; */
  PetscInt k,l;
  PetscInt NMXl = (grid->x[x+m]-grid->x[x])/mdx;
  PetscInt NMYl = (grid->y[y+n]-grid->y[y])/mdy;


  for(k=0; k<NMXl; k++){
    for(l=0; l<NMYl; l++){
      static PetscScalar tmp;
      ierr=PetscRandomSetInterval(r,-mdx/4.0,mdx/4.0);CHKERRQ(ierr);
      ierr=PetscRandomGetValue(r,&tmp);
      markers[idx].X = grid->x[x] + mdx*(0.5+(PetscScalar) k) + tmp;
      
      ierr=PetscRandomSetInterval(r,-mdy/4.0,mdy/4.0);CHKERRQ(ierr);
      ierr=PetscRandomGetValue(r,&tmp);
      markers[idx].Y = grid->y[y] + mdy*(0.5+(PetscScalar) l) + tmp;
      markers[idx].Z = 0.0;
      idx++;
    }
  }
  
  markerset->nMark = idx;
  ierr=PetscRandomDestroy(&r);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}



PetscErrorCode distributeMarkersUniformInCells( MarkerSet *markerset, Options *options, GridData *grid){
  PetscErrorCode ierr;
  PetscRandom r;
  PetscInt NMX=options->NMX;
  PetscInt NMY=options->NMY;
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  Marker *markers = markerset->markers;

  PetscFunctionBegin;
  ierr=PetscRandomCreate( PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  ierr=PetscRandomSetType(r,PETSCRAND);CHKERRQ(ierr);
  ierr=PetscRandomSetSeed(r, (unsigned long) 0);CHKERRQ(ierr);
  
  /* get local ownership range*/
  int x,y,z,m,n,p;
  ierr=DMDAGetCorners(grid->da,&x,&y,&z,&m,&n,&p); CHKERRQ(ierr);
  /* x,y,z are coordinate indices of upper left node in gerya notation (lower left in petsc notation) */
  
  PetscInt i,j;
  PetscInt idx=0;/* marker counter*/
  
  for(i=x;i<x+m && i<NX-1;i++){
    PetscScalar mdx = (grid->x[i+1]-grid->x[i])/((PetscScalar) NMX);
    for (j=y;j<y+n && j<NY-1;j++){
      PetscScalar mdy = (grid->y[j+1]-grid->y[j])/((PetscScalar) NMY);
/*       PetscInt m1; */
      PetscInt k,l;
      for(k=0;k<NMX;k++){
	for(l=0;l<NMY;l++){
	  PetscScalar tmp=0.0;
	  ierr=PetscRandomSetInterval(r,-mdx/2.0,mdx/2.0);CHKERRQ(ierr);
	  //	  ierr=PetscRandomGetValue(r,&tmp);
	  markers[idx].X = grid->x[i] + mdx*(0.5+(PetscScalar) k) + tmp;
	  
	  ierr=PetscRandomSetInterval(r,-mdy/2.0,mdy/2.0);CHKERRQ(ierr);
	  //ierr=PetscRandomGetValue(r,&tmp);
	  markers[idx].Y = grid->y[j] + mdy*(0.5+(PetscScalar) l) + tmp;
	  markers[idx].Z = 0.0;
	  idx++;
	}
      }
    }
  }
   
  markerset->nMark = idx;
  ierr=PetscRandomDestroy(&r);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}


PetscErrorCode checkMarkerDensity( MarkerSet *markerset, GridData *grid, Options *options, PetscRandom r){
  PetscErrorCode ierr=0;

  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  PetscScalar LX = grid->LX;
  PetscScalar LY = grid->LY;
  const PetscInt minMarkers = options->minMarkers;/* minimum number of markers per quarter-cell*/
  const PetscInt maxMarkers = options->maxMarkers;
  PetscInt nMark = markerset->nMark;
  PetscInt *nMarkers;
  /*   PetscInt *nMarkersp; */
  PetscInt *isout;/* flag for marker out of bounds*/
  Marker *markers = markerset->markers;
  PetscMPIInt rank,size;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscInt x,y,m,n;/* local grid lines owned*/
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  /* adjust size of local domain*/
  /*  1 2      m-2  m-1 */
  /*  x x+1 .. NX-2 NX-1 */
  if(x+m == NX && !grid->xperiodic) m--;
  if(y+n == NY) n--; /* the last cell in any row or column cannot contain any markers because it is outside the domain. do not count or add markers to this cell */


  ierr = PetscMalloc(2*(m)*2*(n)*sizeof(PetscInt), &nMarkers);CHKERRQ(ierr);
  ierr = PetscMalloc(nMark*sizeof(PetscInt), &isout);CHKERRQ(ierr);
  //  ierr = PetscMalloc(2*NX*NY*sizeof(PetscInt), &nMarkersp);CHKERRQ(ierr);

  PetscInt nOut = 0;
  {
    PetscInt i;
    for(i=0;i<2*(n)*2*(m);i++){
      nMarkers[i] = 0;
    }
  }
  PetscInt m1;
  for(m1=0;m1<markerset->nMark;m1++){
    if( markers[m1].cellX != -1){
      /* determine basic quarter-cell that this marker belongs to*/
      PetscInt cellx = (markers[m1].cellX);
      PetscInt celly = (markers[m1].cellY);/* this returns the cell in global indexing*/
      if(markers[m1].X >= grid->xc[cellx+1]){
	cellx = 2*(cellx-x) + 1;     /* this converts global to local indexing*/
      } else {
	cellx = 2*(cellx-x) + 0;
      }
      if(markers[m1].Y >= grid->yc[celly+1]){
	celly=2*(celly-y)+1;
      }else{
	celly=2*(celly-y)+ 0;
      }

      //printf("x %e y %e cellJ %d cellI %d\n",markers[m].X,markers[m].Y,cellJ,cellI);
      if( cellx + celly*2*(m) < 0 || cellx+celly*2*(m) >= (2*(m)*2*(n))) {printf("PANIC %d cellx=%d celly=%d!!!\n", (2*(m-1)*2*(n-1)),cellx,celly);fflush(stdout);}
      nMarkers[cellx+celly*2*(m)] +=1;
      /*   PetscInt pcellJ = floor(markers[m].X/dx-0.5)+1; */
      /*       PetscInt pcellI = floor(markers[m].Y/dy-0.5)+1; */
      /*       if(pcellI < 1){ pcellI = 1;}else if(pcellI > NY-1){ pcellI = NY-1;} */
      /*       if(pcellJ < 1){ pcellJ = 1;}else if(pcellJ > NX-1){ pcellI = NX-1;} */
      
      /* nMarkersp[pcellJ*NY + pcellI] +=1; */
    } else {
      /* mark this marker as being out of bounds */
      isout[nOut] = m1;/* will contain a list of markers that are out of bounds*/
      nOut++; /* total number of markers that have left the domain*/
      //      printf("marker %d out of bounds\n",m);
    }
  }
  {/* print number of markers per cell for debugging*/ 
    /*     PetscInt i,j; */
    /*     for(i=0;i<2*NY;i++){ */
    /*       for(j=0;j<2*NX;j++){ */
    /* 	PetscPrintf(PETSC_COMM_SELF,"%d,",nMarkers[j*2*NY+i]); */
    /*       } */
    /*       printf("\n"); */
    /*     } */
  }
#ifdef DEBUG
  { /* print number of markers per cell for debugging*/ 

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] number of markers:\n",rank);
    PetscInt i,j;
    for(j=0;j<2*(n);j++){
      for(i=0;i<2*(m);i++){
    	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d,",nMarkers[i+2*(m)*j]);
      }
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  }
#endif
  /* count largest number of markers in any quarter-cell */
  PetscInt maxMarkAnyQC = 0;
  {PetscInt i,j;
    for(j=0;j<2*(n);j++){
      for(i=0;i<2*(m);i++){
	if(nMarkers[2*m*j+i]>maxMarkAnyQC) maxMarkAnyQC = nMarkers[2*m*j+i];
      }
    }
  }
  /*check marker density*/
  PetscInt nAdd = 0;
  PetscInt nRemove = 0;
  {
    PetscInt i,j;
    /*   printf("nadd:\n"); */
    for(j=0;j<2*(n);j++){
      for(i=0;i<2*(m);i++){
	if(nMarkers[i+2*(m)*j] < minMarkers){
	  /* add markers to bring number of markers up to original marker density*/
	  /* 	printf("%d,",(PetscInt)(options->NMX*options->NMY/4.0 - nMarkers[2*j*NY+i])); */
	  PetscInt tmp1 = (PetscInt) (ceil((PetscScalar)options->NMX*(PetscScalar)options->NMY/4.0) - nMarkers[2*j*(m)+i]);
	  nAdd += tmp1;
	  nMarkers[2*(m)*j+i] = tmp1;
	} else if(nMarkers[i+2*(m)*j] > maxMarkers){
	  PetscInt tmp1 = (PetscInt) (-maxMarkers + nMarkers[2*j*(m)+i]);
	  nRemove += tmp1;
	  nMarkers[2*(m)*j+i] = -tmp1;
	} else {
	  nMarkers[2*(m)*j+i] = 0;
	}
      }
      /*     printf("\n"); */
    }
  }
#ifdef DEBUG
  { /* print number of markers per cell for debugging*/ 
    
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] adding markers ...\n",rank);
    PetscInt i,j;
    for(j=0;j<2*(n);j++){
      for(i=0;i<2*(m);i++){
    	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d,",nMarkers[i+2*(m)*j]);
      }
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  }
#endif

  if(nRemove > 0){
    printf("[%d] removing %d markers, nOut=%d\n",rank,nRemove,nOut);

    PetscInt *mIdx;
    PetscScalar *mRnd;
    ierr = PetscMalloc( maxMarkAnyQC*sizeof(PetscInt), &mIdx);CHKERRQ(ierr);
    ierr = PetscMalloc( maxMarkAnyQC*sizeof(PetscScalar), &mRnd);CHKERRQ(ierr);

    PetscInt i,j;
    for(j=0;j<2*(n);j++){/* loop over y indices */
      PetscScalar ymin = (!(j%2)) ? grid->y[y+j/2] : grid->yc[y+(j+1)/2];
      PetscScalar ymax = (!(j%2)) ? grid->yc[y+j/2+1] : grid->y[y+(j+1)/2];      
      for(i=0;i<2*(m);i++){/* loop over x indices */
	PetscScalar xmin = (!(i%2)) ? grid->x[x+i/2] : grid->xc[x+(i+1)/2];
	PetscScalar xmax = (!(i%2)) ? grid->xc[x+i/2+1] : grid->x[x+(i+1)/2];
	
	if( nMarkers[2*(m)*j+i] < 0){
	  /* remove some markers from this basic quarter-cell*/
	  /* identify number of markers to remove */
	  PetscInt nToRemove = -nMarkers[2*m*j+i];
	  /* find markers in this quarter cell */
	  PetscInt ii=0;
	  PetscInt im;
	  for( im=0; im<markerset->nMark; im++){
	    /* if marker im is in the current basic quarter cell, set mIdx[ii] == im; ii++; */
	    if( markers[im].X >= xmin && markers[im].X < xmax && markers[im].Y >= ymin && markers[im].Y < ymax){
	      mIdx[ii] = im;
	      ii++;
	    }
	  }
	  PetscInt nCell = ii;
	  //printf("Removing %d markers. quarter cell has %d markers\n",nToRemove,nCell);
	  /* assign each marker in quarter cell a random number between 0 and 1 */
	  for( ii=0;ii<nCell;ii++){
	    mRnd[ii] = drand48();
	  }
	  /* for nRemove lowest random numbers, assign out of bounds location and cell */
	  PetscInt irem;
	  for( irem=0;irem<nToRemove;irem++){
	    PetscScalar minv=1.0;
	    PetscInt minidx=0;
	    for( ii=0;ii<nCell;ii++){
	      if( mRnd[ii] < minv ){
		minv = mRnd[ii];
		minidx = ii;
	      }
	    }
	    mRnd[minidx] = 1.0;
	    /* eliminate marker with mIdx[minidx] */
	    //printf("removing marker %d\n",mIdx[minidx]);
	    markers[ mIdx[minidx] ].X = -1.0;
	    markers[ mIdx[minidx] ].Y = -1.0;
	    markers[ mIdx[minidx] ].cellX = -1;
	    markers[ mIdx[minidx] ].cellY = -1;
	    isout[ nOut ] = mIdx[minidx];
	    nOut++;
	  }
	}	
      }
    }
    ierr = PetscFree( mIdx );CHKERRQ(ierr);
    ierr = PetscFree( mRnd );CHKERRQ(ierr);


  }

  if(nAdd>0) printf("[%d] nAdd=%d\n",rank,nAdd);
  if(nAdd > 0 && (nAdd+markerset->nMark-nOut) < markerset->maxMark){
    printf("adding %d markers, nOut=%d\n",nAdd,nOut); fflush(stdout);
    PetscScalar *newX, *newY;
    ierr = PetscMalloc2(nAdd,&newX,nAdd,&newY);  CHKERRQ(ierr);

    PetscInt iMark =0;
    PetscInt i,j;
    for(j=0;j<2*(n);j++){/* loop over y indices */
      for(i=0;i<2*(m);i++){/* loop over x indices */
	if( nMarkers[2*(m)*j+i] > 0){
	  /* add some markers to this basic quarter-cell*/
	  PetscInt k;
	  for(k=0;k<nMarkers[2*j*(m)+i];k++){
	    /* markers should be distributed about dx*j + dy*i */

	    /* if i is even, this is left half of cell, else right half*/
	    if( !(i%2) ){/*i even*/
	      ierr=PetscRandomSetInterval(r,grid->x[x+i/2],grid->xc[x+i/2+1]);CHKERRQ(ierr);
	    }else{
	      ierr=PetscRandomSetInterval(r,grid->xc[x+(i+1)/2],grid->x[x+(i+1)/2]);CHKERRQ(ierr);
	    }
	    ierr=PetscRandomGetValue(r, &newX[iMark]);CHKERRQ(ierr);
	    if( !(j%2) ){/*j even*/
	      ierr=PetscRandomSetInterval(r,grid->y[y+j/2],grid->yc[y+j/2+1]);CHKERRQ(ierr);
	    }else{
	      ierr=PetscRandomSetInterval(r,grid->yc[y+(j+1)/2],grid->y[y+(j+1)/2]);CHKERRQ(ierr);
	    }
	    ierr=PetscRandomGetValue(r, &newY[iMark]);CHKERRQ(ierr);
	    iMark++;
	  }
	}
      }
    }
    if( iMark != nAdd ) printf("[%d] iMark = %d, nAdd = %d, PROBLEM\n",rank,iMark, nAdd);
    /* several options exist for setting marker properties. */
    /* One would be to weight neighbor values by distance*/
    /* just inherit nearest neighbor's values*/   
    for(i=0;i<nAdd;i++){
      /* find nearest neighbor to the new marker*/
      PetscScalar mindist = LX*LY; /* initialize mindist with something very large*/
      PetscInt minidx = -1;
      for(m1=0;m1<nMark;m1++){/* note that this is the original number of markers*/
	//if( !isout[m]){
	if( markers[m1].cellX != -1){
	  /* require that the marker be within one cell width*/
	  if( fabs(markers[m1].X - newX[i]) <= (grid->x[markers[m1].cellX +1] - grid->x[markers[m1].cellX]) && fabs(markers[m1].Y - newY[i]) <= (grid->y[markers[m1].cellY+1] - grid->y[markers[m1].cellY])){
	    /* this is a candidate. compute euclidean distances*/
	    PetscScalar x1= (markers[m1].X - newX[i]);
	    PetscScalar y1= (markers[m1].Y - newY[i]);
	    PetscScalar r1 = sqrt( x1*x1 + y1*y1 );
	    if(r1 < mindist){ mindist = r1; minidx = m1;}
	  }	
	}
      }
      if(minidx != -1){
	/* take all properties from nearest marker*/
	PetscInt idx;
	if(nOut>0){	/* if there are markers out of the domain, recycle one of them*/
	  idx = isout[nOut-1];
	  nOut--;
	  //printf("recycling marker %d\n",idx);
	}else{
	  idx = markerset->nMark;
	  markerset->nMark++;
	}
        markers[idx] = markers[minidx];
	markers[idx].X = newX[i];
	markers[idx].Y = newY[i];
	if( markers[idx].X <= grid->x[1] ){//enforce slab geotherm
	  markers[idx].T = slab_inflow_temperature( markers[idx].X, markers[idx].Y,options->slabAngle );
	}else if(markers[idx].X >= grid->x[NX-1]){//enforce mantle inflow temperature
	  markers[idx].T = mantle_temperature();
	}

      } else{
	printf("Error: No markers close enough for nearest-neighbor interpolation\n");
      }
      
    } /* end loop over markers to add, finding nearest neighbors*/
    /*     ierr = PetscRandomDestroy(r);CHKERRQ(ierr); */
    ierr = PetscFree2(newX,newY); CHKERRQ(ierr);

  } else if( (nAdd + markerset->nMark-nOut) > markerset->maxMark) {
    printf("[%d] Error: no more markers can be added! nAdd=%d, nMark=%d, nOut=%d, maxMark=%d\n",rank,nAdd,markerset->nMark,nOut,markerset->maxMark); fflush(stdout);
    abort();
  }
  
  ierr = PetscFree(isout);CHKERRQ(ierr);  
  ierr = PetscFree(nMarkers); CHKERRQ(ierr);
  /*   ierr = PetscFree(nMarkersp); CHKERRQ(ierr); */
  
  PetscFunctionReturn(ierr);
}

#ifdef done

PetscErrorCode checkMarkersOutsideDomain( Markers *markers, PetscScalar LX, PetscScalar LY){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  PetscInt i;
  for(i=0;i<markerset->nMark;i++){
    PetscScalar markerEps = 1;
    if( markers->X[i] < 0){ markers->X[i] = markerEps;
    }else if(markers->X[i] > LX){ markers->X[i] = LX-markerEps;
    }else if(markers->Y[i] < 0.0){ markers->Y[i] = markerEps;
    }else if(markers->Y[i] > LY){ markers->Y[i] = LY-markerEps;}
  }
  PetscFunctionReturn(ierr);
}

#endif
/* new routine to project markers onto nodes using DA structures*/
PetscErrorCode projectMarkersNodesAll(MarkerSet *markerset, GridData *grid, NodalFields *nodalFields, Materials *materials, Options *options){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* start by projecting only density*/
  PetscMPIInt rank,size;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /* field variables: rho,etaS,etaN,kThermal,T,Cp,muN,muS,soxx,soxy,soyy,sozz,soxy,soxz,soyz*/
  
  Vec rhol, rhodotl, etaSl, etaNl,etavxl,etavyl, kThermall,Tl,Cpl,muNl,muSl,muvxl,muvyl,soxxl,soyyl,sozzl,soxyl,soxzl,soyzl,hal;
  Vec wtng, wtnl;/* nodal weights for basic nodes*/
  Vec wtvxg, wtvxl, wtvyg, wtvyl; /* nodal weights for velocity nodes */
  Vec wtetasg, wtetasl;/* nodal weights for shear viscosity*/
  Vec wtetang, wtetanl;/* nodal weights for normal viscosity*/
#ifdef TEXTURE
  Vec VPNbl,VPNcl;/* local vectors for viscoplasticity tensor*/
#endif
  Marker *markers = markerset->markers;

  /* create global vectors for weights */
  ierr=DMCreateGlobalVector(grid->da, &wtng);CHKERRQ(ierr);
  ierr=VecDuplicate(wtng, &wtetang);CHKERRQ(ierr);
  ierr=VecDuplicate(wtng, &wtetasg);CHKERRQ(ierr);
  ierr=VecDuplicate(wtng, &wtvxg); CHKERRQ(ierr);
  ierr=VecDuplicate(wtng, &wtvyg); CHKERRQ(ierr);
  
  /* create local vectors for weights */
  ierr=DMCreateLocalVector(grid->da, &wtnl);CHKERRQ(ierr);
  ierr=VecDuplicate(wtnl, &wtetanl);CHKERRQ(ierr);
  ierr=VecDuplicate(wtnl, &wtetasl);CHKERRQ(ierr);
  ierr=VecDuplicate(wtnl, &wtvxl);CHKERRQ(ierr);
  ierr=VecDuplicate(wtnl, &wtvyl);CHKERRQ(ierr);
  
  /* create local vectors for field variables*/
  ierr = VecDuplicate(wtnl, &rhol );CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &rhodotl );CHKERRQ(ierr);
  
  
  ierr = VecDuplicate(wtnl, &etaSl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &etaNl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &etavxl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &etavyl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &muNl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &muSl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &muvxl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &muvyl);CHKERRQ(ierr);
  
  ierr = VecDuplicate(wtnl, &kThermall);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &Tl     );CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &Cpl    );CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &soxxl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &soyyl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &sozzl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &soxyl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &soxzl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &soyzl);CHKERRQ(ierr);
  ierr = VecDuplicate(wtnl, &hal);CHKERRQ(ierr);
#ifdef TEXTURE
  ierr = DMCreateLocalVector(grid->tda, &VPNbl);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(grid->tda, &VPNcl);CHKERRQ(ierr);
#endif
  
  /* zero global fields if necessary*/
  /* zero all weights */
  ierr = VecZeroEntries( wtng ); CHKERRQ(ierr);
  ierr = VecZeroEntries( wtetang ); CHKERRQ(ierr);
  ierr = VecZeroEntries( wtetasg ); CHKERRQ(ierr);
  ierr = VecZeroEntries( wtvxg ); CHKERRQ(ierr);
  ierr = VecZeroEntries( wtvyg ); CHKERRQ(ierr);
  
  ierr = VecZeroEntries(  nodalFields->rho );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->rhodot );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->etaN );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->etaS );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->etavx );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->etavy );CHKERRQ(ierr);

  ierr = VecZeroEntries(  nodalFields->kThermal);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->lastT );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->Cp);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->muN);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->muS);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->muvx );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->muvy );CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->soxx);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->soyy);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->sozz);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->soxy);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->soxz);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->soyz);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->ha);CHKERRQ(ierr);
#ifdef TEXTURE
  ierr = VecZeroEntries(  nodalFields->VPTensorB);CHKERRQ(ierr);
  ierr = VecZeroEntries(  nodalFields->VPTensorC);CHKERRQ(ierr);
#endif


  /* begin scatter - only one of these can be done at a time because the DA only has one scatter context available*/
  /* node weights */
  ierr=DMGlobalToLocalBegin(grid->da,wtng,INSERT_VALUES,wtnl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,wtng,INSERT_VALUES,wtnl);CHKERRQ(ierr);
  /* cell centers */
  ierr=DMGlobalToLocalBegin(grid->da,wtetang,INSERT_VALUES,wtetanl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,wtetang,INSERT_VALUES,wtetanl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,wtetasg,INSERT_VALUES,wtetasl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,wtetasg,INSERT_VALUES,wtetasl);CHKERRQ(ierr);
  /* vx and vy nodes */
  ierr=DMGlobalToLocalBegin(grid->da,wtvxg,INSERT_VALUES,wtvxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,wtvxg,INSERT_VALUES,wtvxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,wtvyg,INSERT_VALUES,wtvyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,wtvyg,INSERT_VALUES,wtvyl);CHKERRQ(ierr);
  
  /* rho */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  /* etaS */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->etaS,INSERT_VALUES,etaSl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->etaS,INSERT_VALUES,etaSl);CHKERRQ(ierr);

  /* scatter global density vector to local density vector*/
  PetscScalar **wtn;   /* weights for basic nodes*/
  PetscScalar **wtetas;/* weights only for markers within dl/2*/
  PetscScalar **wtetan;/* cell center weights*/
  PetscScalar **wtvx;
  PetscScalar **wtvy;
  /* field variables: rho,etaS,etaN,kThermal,T,Cp,muN,muS,soxx,soxy,soyy,sozz,soxy,soxz,soyz*/
  PetscScalar **rho,**rhodot, **etaS,**etavx,**etavy, **etaN, **kThermal,**T,**Cp,**muN,**muS,**muvx,**muvy,**soxx,**soyy,**sozz,**soxy,**soxz,**soyz,**ha;
#ifdef TEXTURE
  /* local viscoplasticity tensors */ 
  Tensor22 **VPNb;
  Tensor22 **VPNc;
#endif

  /* weights*/
  ierr=DMDAVecGetArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,wtvxl,&wtvx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,wtvyl,&wtvy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,wtetasl,&wtetas);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,wtetanl,&wtetan);CHKERRQ(ierr);
  /* field variables */

  ierr=DMDAVecGetArray(grid->da,rhol,&rho);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,rhodotl,&rhodot);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,etaSl,&etaS);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,etaNl,&etaN);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,etavxl,&etavx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,etavyl,&etavy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,kThermall,&kThermal);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,Tl,&T);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,Cpl,&Cp);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,muNl,&muN);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,muSl,&muS);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,muvxl,&muvx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,muvyl,&muvy);CHKERRQ(ierr);

  ierr=DMDAVecGetArray(grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,soyyl,&soyy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,sozzl,&sozz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,soxyl,&soxy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,soxzl,&soxz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,soyzl,&soyz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,hal,&ha);CHKERRQ(ierr);
#ifdef TEXTURE
  ierr=DMDAVecGetArray(grid->tda,VPNbl,&VPNb);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->tda,VPNcl,&VPNc);CHKERRQ(ierr);
#endif


  PetscInt m1;
  /* calculate weight functions first*/
  PetscInt x,y,z,m,n,p;
  DMDAGetCorners(grid->da,&x,&y,&z,&m,&n,&p);
  for(m1=0;m1<markerset->nMark;m1++){
    if(markers[m1].cellX != -1){
      PetscInt cellI = markers[m1].cellY;
      PetscInt cellJ = markers[m1].cellX;
      PetscScalar dx = grid->x[markers[m1].cellX+1] - grid->x[markers[m1].cellX];
      PetscScalar dy = grid->y[markers[m1].cellY+1] - grid->y[markers[m1].cellY];
      PetscScalar mdx = (markers[m1].X - grid->x[markers[m1].cellX])/dx;
      PetscScalar mdy = (markers[m1].Y - grid->y[markers[m1].cellY])/dy;
      if(mdx>1.0 || mdy> 1.0 || mdx<0 || mdy<0) printf("ERROR: mdx=%f,mdy=%f\n",mdx,mdy);
      /*weight*/
      PetscScalar weight = (1.0-mdx)*(1.0-mdy);/* i,j */
      /* node i,j*/
      /*BASIC NODES: rho,rhodot,T,ha,kThermal,Cp*/
      wtn[cellI][cellJ]     += weight;
      rho[cellI][cellJ]     += weight*markers[m1].rho;
      rhodot[cellI][cellJ]  += weight*markers[m1].rhodot;
      T[cellI][cellJ]  += weight*markers[m1].T;
      ha[cellI][cellJ]  += weight*markers[m1].T*materials->materialAlpha[(PetscInt) markers[m1].Mat];
      kThermal[cellI][cellJ]  += weight*materials->materialkThermal[(PetscInt) markers[m1].Mat];
      Cp[cellI][cellJ]  += weight*materials->materialCp[(PetscInt) markers[m1].Mat];
      /* node i+1,j*/
      weight= (1.0-mdx)*(mdy);
      wtn[cellI+1][cellJ]   += weight;
      rho[cellI+1][cellJ]   += weight*markers[m1].rho;
      rhodot[cellI+1][cellJ]  += weight*markers[m1].rhodot;
      T[cellI+1][cellJ]  += weight*markers[m1].T;
      ha[cellI+1][cellJ]  += weight*markers[m1].T*materials->materialAlpha[(PetscInt) markers[m1].Mat];
      kThermal[cellI+1][cellJ]  += weight*materials->materialkThermal[(PetscInt) markers[m1].Mat];
      Cp[cellI+1][cellJ]  += weight*materials->materialCp[(PetscInt) markers[m1].Mat];
      /* node i,j+1*/
      weight= (mdx)*(1.0-mdy);
      wtn[cellI][cellJ+1]   += weight;
      rho[cellI][cellJ+1]   += weight*markers[m1].rho;
      rhodot[cellI][cellJ+1]  += weight*markers[m1].rhodot;
      T[cellI][cellJ+1]  += weight*markers[m1].T;
      ha[cellI][cellJ+1]  += weight*markers[m1].T*materials->materialAlpha[(PetscInt) markers[m1].Mat];
      kThermal[cellI][cellJ+1]  += weight*materials->materialkThermal[(PetscInt) markers[m1].Mat];
      Cp[cellI][cellJ+1]  += weight*materials->materialCp[(PetscInt) markers[m1].Mat];
      /* node i+1,j+1*/
      weight= (mdx)*(mdy);
      wtn[cellI+1][cellJ+1] += weight;
      rho[cellI+1][cellJ+1] += weight*markers[m1].rho;
      rhodot[cellI+1][cellJ+1]  += weight*markers[m1].rhodot;
      T[cellI+1][cellJ+1]  += weight*markers[m1].T;
      ha[cellI+1][cellJ+1]  += weight*markers[m1].T*materials->materialAlpha[(PetscInt) markers[m1].Mat];
      kThermal[cellI+1][cellJ+1]  += weight*materials->materialkThermal[(PetscInt) markers[m1].Mat];
      Cp[cellI+1][cellJ+1]  += weight*materials->materialCp[(PetscInt) markers[m1].Mat];
      /*BASIC NODES with local interpolation*/
      /*etaS,muS,soxy,soxz,soyz*/
      if(mdx<=0.5 && mdy <=0.5){/* contribute to etas(x,y)-upper left*/
	weight=(1.0-mdx)*(1.0-mdy);
	etaS[cellI][cellJ]+= weight*markers[m1].eta;
	muS[cellI][cellJ] += weight/markers[m1].mu;   /* not harmonic averaging for mu as per Gerya's suggestion*/
	/*.mumaterials->materialMu[(PetscInt) markers[m1].Mat];*/ 
#ifdef TEXTURE
	VPNb[cellI][cellJ].T11 += weight*markers[m1].texture.N[0][0];
	VPNb[cellI][cellJ].T12 += weight*markers[m1].texture.N[0][1];
	VPNb[cellI][cellJ].T21 += weight*markers[m1].texture.N[1][0];
	VPNb[cellI][cellJ].T22 += weight*markers[m1].texture.N[1][1];
#endif
	soxy[cellI][cellJ] += weight*markers[m1].s.T12;//sxy;
	soxz[cellI][cellJ] += weight*markers[m1].s.T13;//sxz;
	soyz[cellI][cellJ] += weight*markers[m1].s.T23;//syz;	

	wtetas[cellI][cellJ] += weight;
	/* contribute to xz node (i,j) */
	weight = (1.0-mdx)*mdy;
	/* 	soxz[cellI][cellJ] += weight*markers[m1].s.T13; */ //sxz;
	etavx[cellI][cellJ] += weight*markers[m1].eta;
	muvx[cellI][cellJ] += weight/markers[m1].mu;
	  /* materials->materialMu[(PetscInt) markers[m1].Mat]; */
	wtvx[cellI][cellJ] += weight;
	/* yz node (i,j) */
	weight = mdx*(1.0-mdy);

	muvy[cellI][cellJ] += weight/markers[m1].mu;
	/* materials->materialMu[(PetscInt) markers[m1].Mat]; */
	etavy[cellI][cellJ] += weight*markers[m1].eta;
	wtvy[cellI][cellJ] += weight;
      }else if(mdx>=0.5 && mdy <=0.5){/* contribute to etas(x+1,y) - upper right*/
	weight=(mdx)*(1.0-mdy);
	etaS[cellI][cellJ+1]+= weight*markers[m1].eta;
	muS[cellI][cellJ+1] += weight/markers[m1].mu;
	/* materials->materialMu[(PetscInt) markers[m1].Mat]; */
#ifdef TEXTURE
	VPNb[cellI][cellJ+1].T11 += weight*markers[m1].texture.N[0][0];
	VPNb[cellI][cellJ+1].T12 += weight*markers[m1].texture.N[0][1];
	VPNb[cellI][cellJ+1].T21 += weight*markers[m1].texture.N[1][0];
	VPNb[cellI][cellJ+1].T22 += weight*markers[m1].texture.N[1][1];
#endif
	soxy[cellI][cellJ+1] += weight*markers[m1].s.T12;//sxy;
	soxz[cellI][cellJ+1] += weight*markers[m1].s.T13;//sxz;
	soyz[cellI][cellJ+1] += weight*markers[m1].s.T23;//syz;

	wtetas[cellI][cellJ+1] += weight;
	/* weight for vx node */
	weight = mdx*mdy;

	etavx[cellI][cellJ+1] += weight*markers[m1].eta;
	muvx[cellI][cellJ+1]  += weight/markers[m1].mu;
	/* materials->materialMu[(PetscInt) markers[m1].Mat]; */
	wtvx[cellI][cellJ+1]  += weight;
	/* yz node i,j */
	weight = (1.0-mdx)*(1.0-mdy);

	muvy[cellI][cellJ] += weight/markers[m1].mu;
	/* materials->materialMu[(PetscInt) markers[m1].Mat]; */
	etavy[cellI][cellJ] += weight*markers[m1].eta;
	wtvy[cellI][cellJ]+= weight;
      }else if(mdx<=0.5 && mdy >=0.5){/* contribute to etas(x,y+1) - lower left*/
	weight=(1.0-mdx)*(mdy);
	etaS[cellI+1][cellJ]+= weight*markers[m1].eta;
	muS[cellI+1][cellJ] += weight/markers[m1].mu;
	/* materials->materialMu[(PetscInt) markers[m1].Mat]; */
#ifdef TEXTURE
	VPNb[cellI+1][cellJ].T11 += weight*markers[m1].texture.N[0][0];
	VPNb[cellI+1][cellJ].T12 += weight*markers[m1].texture.N[0][1];
	VPNb[cellI+1][cellJ].T21 += weight*markers[m1].texture.N[1][0];
	VPNb[cellI+1][cellJ].T22 += weight*markers[m1].texture.N[1][1];
#endif
	soxy[cellI+1][cellJ]+= weight*markers[m1].s.T12;//xy;
	soxz[cellI+1][cellJ] += weight*markers[m1].s.T13;//sxz;
	soyz[cellI+1][cellJ] += weight*markers[m1].s.T23;//syz;

	wtetas[cellI+1][cellJ] += weight;

	/* weight for vx node */
	weight = (1.0-mdx)*(1.0-mdy);


	etavx[cellI][cellJ] += weight*markers[m1].eta;
	muvx[cellI][cellJ] += weight/markers[m1].mu;
	  /* materials->materialMu[(PetscInt) markers[m1].Mat]; */
	wtvx[cellI][cellJ] += weight;
	/* yz node (i,j+1) */
	weight = mdx*mdy;

	muvy[cellI+1][cellJ] += weight/markers[m1].mu;
	/* 	  materials->materialMu[(PetscInt) markers[m1].Mat]; */
	etavy[cellI+1][cellJ] += weight*markers[m1].eta;
	wtvy[cellI+1][cellJ] += weight;
      } else if(mdx>=0.5 && mdy >=0.5){/* contribute to etas(x+1,y+1) - lower right*/
	weight=(mdx)*(mdy);
	etaS[cellI+1][cellJ+1]+= weight*markers[m1].eta;
	muS[cellI+1][cellJ+1] += weight/markers[m1].mu;
	/* 	materials->materialMu[(PetscInt) markers[m1].Mat]; */
#ifdef TEXTURE
	VPNb[cellI+1][cellJ+1].T11 += weight*markers[m1].texture.N[0][0];
	VPNb[cellI+1][cellJ+1].T12 += weight*markers[m1].texture.N[0][1];
	VPNb[cellI+1][cellJ+1].T21 += weight*markers[m1].texture.N[1][0];
	VPNb[cellI+1][cellJ+1].T22 += weight*markers[m1].texture.N[1][1];
#endif
	soxy[cellI+1][cellJ+1]  += weight*markers[m1].s.T12;//sxy;
	soxz[cellI+1][cellJ+1]  += weight*markers[m1].s.T13;//sxz;
	soyz[cellI+1][cellJ+1]  += weight*markers[m1].s.T23;//syz;
	wtetas[cellI+1][cellJ+1] += weight;

	/* weight for vx node */
	weight = mdx*(1.0-mdy);

	etavx[cellI][cellJ+1] += weight*markers[m1].eta;
	muvx[cellI][cellJ+1]  += weight/markers[m1].mu;
	/* 	materials->materialMu[(PetscInt) markers[m1].Mat]; */
	wtvx[cellI][cellJ+1]  += weight;
	/* vz node (i,j+1) */
	weight = (1.0-mdx)*mdy;

	muvy[cellI+1][cellJ]  += weight/markers[m1].mu;
	/* materials->materialMu[(PetscInt) markers[m1].Mat]; */
	etavy[cellI+1][cellJ] += weight*markers[m1].eta;
	wtvy[cellI+1][cellJ]  += weight;
      }/* end assignment of etas fields*/
      /* compute marker contributions to 'xz' nodes (colocated with vx nodes in 2.5D */
      
      /* center node*/
      /* etaN,muN,soxx,soyy,sozz */
      {
	PetscScalar weightc = (1.0-fabs(0.5-mdx))*(1.0-fabs(0.5-mdy));
	wtetan[cellI+1][cellJ+1] += weightc;
	etaN[cellI+1][cellJ+1] += weightc*markers[m1].eta;
	muN[cellI+1][cellJ+1] += weightc/markers[m1].mu;
	/* 	materials->materialMu[(PetscInt) markers[m1].Mat]; */
#ifdef TEXTURE
	VPNc[cellI+1][cellJ+1].T11 += weightc*markers[m1].texture.N[0][0];
	VPNc[cellI+1][cellJ+1].T12 += weightc*markers[m1].texture.N[0][1];
	VPNc[cellI+1][cellJ+1].T21 += weightc*markers[m1].texture.N[1][0];
	VPNc[cellI+1][cellJ+1].T22 += weightc*markers[m1].texture.N[1][1];
#endif
	soxx[cellI+1][cellJ+1] += weightc*markers[m1].s.T11;//sxx;
	soyy[cellI+1][cellJ+1] += weightc*markers[m1].s.T22;//syy;
	sozz[cellI+1][cellJ+1] += weightc*markers[m1].s.T33;//szz;
      }
    }
  }
  /* restore arrays for weights*/
  ierr = DMDAVecRestoreArray(grid->da,wtnl,&wtn); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,wtetanl,&wtetan); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,wtetasl,&wtetas); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,wtvxl,&wtvx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,wtvyl,&wtvy); CHKERRQ(ierr);
  /* restore arrays for field variables*/
  ierr=DMDAVecRestoreArray(grid->da,rhol,&rho);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,rhodotl,&rhodot);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etaSl,&etaS);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etaNl,&etaN);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etavxl,&etavx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etavyl,&etavy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,kThermall,&kThermal);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,Tl,&T);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,Cpl,&Cp);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,muNl,&muN);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,muSl,&muS);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,muvxl,&muvx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,muvyl,&muvy);CHKERRQ(ierr);


  ierr=DMDAVecRestoreArray(grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soyyl,&soyy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,sozzl,&sozz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soxyl,&soxy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soxzl,&soxz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soyzl,&soyz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,hal,&ha);CHKERRQ(ierr);
#ifdef TEXTURE
  ierr=DMDAVecRestoreArray(grid->tda,VPNcl,&VPNc);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->tda,VPNbl,&VPNb);CHKERRQ(ierr);
#endif
  
  /* scatter arrays for weights-begin/end syntax is necessary to treat ghosted values properly*/
  ierr = DMLocalToGlobalBegin(grid->da,wtnl,ADD_VALUES,wtng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtnl,ADD_VALUES,wtng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,wtetasl,ADD_VALUES,wtetasg); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtetasl,ADD_VALUES,wtetasg); CHKERRQ(ierr);
  /*   { */
  /*     PetscPrintf(PETSC_COMM_WORLD,"wtetasg:\n");CHKERRQ(ierr); */
  /*     ierr=VecView(wtetasg,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
  /*   } */

  ierr = DMLocalToGlobalBegin(grid->da,wtetanl,ADD_VALUES,wtetang); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtetanl,ADD_VALUES,wtetang); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,wtvxl,ADD_VALUES,wtvxg); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtvxl,ADD_VALUES,wtvxg); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,wtvyl,ADD_VALUES,wtvyg); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtvyl,ADD_VALUES,wtvyg); CHKERRQ(ierr); 

  /* scatter arrays for field variables*/
  ierr = DMLocalToGlobalBegin(grid->da,rhol,ADD_VALUES,nodalFields->rho);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,rhol,ADD_VALUES,nodalFields->rho);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,rhodotl,ADD_VALUES,nodalFields->rhodot);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,rhodotl,ADD_VALUES,nodalFields->rhodot);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,etaSl,ADD_VALUES,nodalFields->etaS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,etaSl,ADD_VALUES,nodalFields->etaS);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,etaNl,ADD_VALUES,nodalFields->etaN);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,etaNl,ADD_VALUES,nodalFields->etaN);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,etavxl,ADD_VALUES,nodalFields->etavx);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,etavxl,ADD_VALUES,nodalFields->etavx);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,etavyl,ADD_VALUES,nodalFields->etavy);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,etavyl,ADD_VALUES,nodalFields->etavy);CHKERRQ(ierr);

#ifdef TEXTURE
  ierr = DMLocalToGlobalBegin(grid->tda,VPNbl,ADD_VALUES,nodalFields->VPTensorB); CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->tda,VPNbl,ADD_VALUES,nodalFields->VPTensorB); CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->tda,VPNcl,ADD_VALUES,nodalFields->VPTensorC); CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->tda,VPNcl,ADD_VALUES,nodalFields->VPTensorC); CHKERRQ(ierr);
#endif
  ierr = DMLocalToGlobalBegin(grid->da,kThermall,ADD_VALUES,nodalFields->kThermal);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,kThermall,ADD_VALUES,nodalFields->kThermal);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,Tl     ,ADD_VALUES,nodalFields->lastT);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,Tl     ,ADD_VALUES,nodalFields->lastT);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,Cpl    ,ADD_VALUES,nodalFields->Cp);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,Cpl    ,ADD_VALUES,nodalFields->Cp);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,muNl,ADD_VALUES,nodalFields->muN);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,muNl,ADD_VALUES,nodalFields->muN);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,muSl,ADD_VALUES,nodalFields->muS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,muSl,ADD_VALUES,nodalFields->muS);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,muvxl,ADD_VALUES,nodalFields->muvx);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,muvxl,ADD_VALUES,nodalFields->muvx);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(grid->da,muvyl,ADD_VALUES,nodalFields->muvy);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,muvyl,ADD_VALUES,nodalFields->muvy);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,soxxl,ADD_VALUES,nodalFields->soxx);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,soxxl,ADD_VALUES,nodalFields->soxx);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,soyyl,ADD_VALUES,nodalFields->soyy);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,soyyl,ADD_VALUES,nodalFields->soyy);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,sozzl,ADD_VALUES,nodalFields->sozz);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,sozzl,ADD_VALUES,nodalFields->sozz);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,soxyl,ADD_VALUES,nodalFields->soxy);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,soxyl,ADD_VALUES,nodalFields->soxy);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,soxzl,ADD_VALUES,nodalFields->soxz);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,soxzl,ADD_VALUES,nodalFields->soxz);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,soyzl,ADD_VALUES,nodalFields->soyz);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,soyzl,ADD_VALUES,nodalFields->soyz);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(grid->da,hal,ADD_VALUES,nodalFields->ha);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,hal,ADD_VALUES,nodalFields->ha);CHKERRQ(ierr);

  /* normalize nodal fields by weightsums*/

  /*BASIC NODES: rho,rhodot,T,ha,kThermal,Cp*/
  ierr = VecPointwiseDivide(nodalFields->rho,nodalFields->rho,wtng);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->rhodot,nodalFields->rhodot,wtng);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->lastT,nodalFields->lastT,wtng);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->ha,nodalFields->ha,wtng);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->kThermal,nodalFields->kThermal,wtng);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->Cp,nodalFields->Cp,wtng);CHKERRQ(ierr);
  /*BASIC NODES with local interpolation: etaS,muS,soxy,soxz,soyz*/
  ierr = VecPointwiseDivide(nodalFields->etaS,nodalFields->etaS,wtetasg);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->etavx,nodalFields->etavx,wtvxg);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->etavy,nodalFields->etavy,wtvyg);CHKERRQ(ierr);

  ierr = VecPointwiseDivide(nodalFields->muS,wtetasg,nodalFields->muS);CHKERRQ(ierr);/* note different form-harmonic average*/
  ierr = VecPointwiseDivide(nodalFields->muvx,wtvxg,nodalFields->muvx);CHKERRQ(ierr);/* note different form-harmonic average*/
  ierr = VecPointwiseDivide(nodalFields->muvy,wtvyg,nodalFields->muvy);CHKERRQ(ierr);/* note different form-harmonic average*/

  ierr = VecPointwiseDivide(nodalFields->soxy,nodalFields->soxy,wtetasg);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->soxz,nodalFields->soxz,wtetasg);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->soyz,nodalFields->soyz,wtetasg);CHKERRQ(ierr);
  /*   ierr = VecPointwiseDivide(nodalFields->soxz,nodalFields->soxz,wtvxg);CHKERRQ(ierr); */  /* use these eventually */
  /*   ierr = VecPointwiseDivide(nodalFields->soyz,nodalFields->soyz,wtvyg);CHKERRQ(ierr); */
#ifdef TEXTURE
  //ierr = VecPointwiseDivide(nodalFields->VPTensorB,nodalFields->VPTensorB,wtetasg);CHKERRQ(ierr);
#endif
  /* cell-center fields*/
  /* etaN,muN,soxx,soyy,sozz */
  ierr = VecPointwiseDivide(nodalFields->etaN,nodalFields->etaN,wtetang);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->muN,wtetang,nodalFields->muN);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->soxx,nodalFields->soxx,wtetang);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->soyy,nodalFields->soyy,wtetang);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(nodalFields->sozz,nodalFields->sozz,wtetang);CHKERRQ(ierr);
#ifdef TEXTURE
  {
    //ierr = VecPointwiseDivide(nodalFields->VPTensorC,nodalFields->VPTensorC,wtetang);CHKERRQ(ierr);
    PetscScalar **wtn, **wts;
    PetscScalar ***Nc, ***Nb;
    ierr = DMDAVecGetArrayDOF(grid->tda, nodalFields->VPTensorC, &Nc);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayDOF(grid->tda, nodalFields->VPTensorB, &Nb);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid->da, wtetang, &wtn); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(grid->da, wtetasg, &wts); CHKERRQ(ierr);
    PetscInt x,y,m,n;
    ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
    PetscInt i,j,k;
    for(j=y;j<y+n;j++){
      for(i=x;i<x+m;i++){
	for(k=0;k<4;k++){
	  Nc[j][i][k] /= wtn[j][i];
	  Nb[j][i][k] /= wts[j][i];
	}
      }
    }
    ierr = DMDAVecRestoreArrayDOF(grid->tda, nodalFields->VPTensorC, &Nc);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayDOF(grid->tda, nodalFields->VPTensorB, &Nb);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid->da, wtetang, &wtn); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(grid->da, wtetasg, &wts); CHKERRQ(ierr);
    
  }
#endif
  /* destroy temporary Vecs for weights*/
  ierr = VecDestroy(&wtng);CHKERRQ(ierr);
  ierr = VecDestroy(&wtnl);CHKERRQ(ierr);
  ierr = VecDestroy(&wtetasl);CHKERRQ(ierr);
  ierr = VecDestroy(&wtetanl);CHKERRQ(ierr);
  ierr = VecDestroy(&wtetasg);CHKERRQ(ierr);
  ierr = VecDestroy(&wtetang);CHKERRQ(ierr);
  ierr = VecDestroy(&wtvxg);CHKERRQ(ierr);
  ierr = VecDestroy(&wtvxl);CHKERRQ(ierr);
  ierr = VecDestroy(&wtvyg);CHKERRQ(ierr);
  ierr = VecDestroy(&wtvyl);CHKERRQ(ierr);
  /* destroy local vecs for field variables*/
  //  Vec rhol, rhodotl, etaSl, etaNl, kThermall,Tl,Cpl,muNl,muSl,soxxl,soyyl,sozzl,soxyl,soxzl,soyzl,hal;
  ierr = VecDestroy(&rhol);CHKERRQ(ierr);
  ierr = VecDestroy(&rhodotl);CHKERRQ(ierr);
  ierr = VecDestroy(&etaSl);CHKERRQ(ierr);
  ierr = VecDestroy(&etaNl);CHKERRQ(ierr);
  ierr = VecDestroy(&etavxl);CHKERRQ(ierr);
  ierr = VecDestroy(&etavyl);CHKERRQ(ierr);
#ifdef TEXTURE
  ierr = VecDestroy(&VPNbl);CHKERRQ(ierr);
  ierr = VecDestroy(&VPNcl);CHKERRQ(ierr);
#endif
  ierr = VecDestroy(&kThermall);CHKERRQ(ierr);
  ierr = VecDestroy(&Tl);CHKERRQ(ierr);
  ierr = VecDestroy(&Cpl);CHKERRQ(ierr);
  ierr = VecDestroy(&muNl);CHKERRQ(ierr);
  ierr = VecDestroy(&muSl);CHKERRQ(ierr);
  ierr = VecDestroy(&muvxl);CHKERRQ(ierr);
  ierr = VecDestroy(&muvyl);CHKERRQ(ierr);
  ierr = VecDestroy(&soxxl);CHKERRQ(ierr);
  ierr = VecDestroy(&soyyl);CHKERRQ(ierr);
  ierr = VecDestroy(&sozzl);CHKERRQ(ierr);
  ierr = VecDestroy(&soxyl);CHKERRQ(ierr);
  ierr = VecDestroy(&soxzl);CHKERRQ(ierr);
  ierr = VecDestroy(&soyzl);CHKERRQ(ierr);
  ierr = VecDestroy(&hal);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}

PetscErrorCode checkPlasticYielding(GridData *grid,MarkerSet *markerset, Materials *materials, PetscScalar dt, PetscInt *isYielding, Options *options, PetscScalar elapsedTime){
  setLogStage( LOG_YIELD );
  const PetscScalar etamin = options->etamin;
  const PetscScalar etamax = options->etamax;
  const PetscScalar fractionalEtamin = options->fractionalEtamin;/* maximum fractional weakening of material relative to materialEta */
  PetscErrorCode ierr=0;
  const PetscInt nMark = markerset->nMark;
  PetscInt m;
  //PetscScalar LX = grid->LX;
  //PetscScalar LY = grid->LY;
  Marker *markers = markerset->markers;

  PetscFunctionBegin;
  for(m=0;m<nMark;m++){
    /*     if(markers[m].X >= 0 && markers[m].Y >= 0 && markers[m].X <= LX && markers[m].Y <= LY){ */
    if(markers[m].cellX != -1){      
      PetscScalar etanew = markers[m].eta;
      PetscScalar X = markers[m].eta/(markers[m].mu*dt + markers[m].eta);
      PetscScalar sxxnew = markers[m].s.T11*X + 2*markers[m].eta*markers[m].e.T11*(1.0-X);
      PetscScalar syynew = markers[m].s.T22*X + 2*markers[m].eta*markers[m].e.T22*(1.0-X);
      PetscScalar szznew = markers[m].s.T33*X + 2*markers[m].eta*markers[m].e.T33*(1.0-X);
      PetscScalar sxynew = markers[m].s.T12*X + 2*markers[m].eta*markers[m].e.T12*(1.0-X);
      PetscScalar sxznew = markers[m].s.T13*X + 2*markers[m].eta*markers[m].e.T13*(1.0-X);
      PetscScalar syznew = markers[m].s.T23*X + 2*markers[m].eta*markers[m].e.T23*(1.0-X);

      /* PetscScalar siinew = sqrt(sxxnew*sxxnew + sxynew*sxynew); */
      PetscScalar siinew = sqrt(0.5*(sxxnew*sxxnew + syynew*syynew+szznew*szznew) + sxynew*sxynew + sxznew*sxznew + syznew*syznew);

      /* this would be the place to re-compute viscosity based on forecast stress and strain-rate */
      updateMarkerViscosity( &markers[m], options, materials, siinew);
      /* correct viscosity for damage */


      if( materials->hasPlasticity[markers[m].Mat] || materials->hasDamage[markers[m].Mat] ){
	PetscInt thismat = (PetscInt) markers[m].Mat;
	/* compute yield stress*/
	//    siiyield = Cohesion + Friction*Pressure
	PetscScalar siiyield;
	/*       switch ( materials->hasPlasticity[(PetscInt) markers[m].Mat] ) { */
	PetscInt plast = materials->hasPlasticity[(PetscInt) markers[m].Mat] ;
	if( plast == 0 ){
	  siiyield = DBL_MAX;
	  /* no plasticity - do nothing */
	  //break;
	}else if( plast == 1){
	  /* simple mohr-coulomb with no strain weakening */
	  siiyield = materials->materialCohesion[(PetscInt) markers[m].Mat]+materials->materialFriction[(PetscInt) markers[m].Mat]*markers[m].p;
	  //break;
	}else if (plast == 2){
	  /* mohr-coulomb with strain weakening */
	  /* look at marker total plastic strain */
	  PetscScalar mEii = markers[m].Eii;
	  PetscScalar c,f;
	  if( mEii < materials->gamma[thismat][0]){
	    /* use initial value */
	    c = materials->C[thismat][0];
	    f = sin(materials->F[thismat][0]/180.0*M_PI);
	  }else if(mEii >= materials->gamma[thismat][0] && mEii <= materials->gamma[thismat][1]){
	    /* use linear interpolation between initial and final */
	    PetscScalar dg = (mEii - materials->gamma[thismat][0])/(materials->gamma[thismat][1]-materials->gamma[thismat][0]);
	    c = (materials->C[thismat][1]-materials->C[thismat][0])*dg + materials->C[thismat][0];
	    f = (sin(materials->F[thismat][1]/180.0*M_PI)-sin(materials->F[thismat][0]/180.0*M_PI))*dg + sin(materials->F[thismat][0]/180.0*M_PI);
	  }else{
	    /* strain greater than gamma1 -use final value */
	    c = materials->C[thismat][1];
	    f = sin(materials->F[thismat][1]/180.0*M_PI);
	  }
	  siiyield = c + f*markers[m].p;	

	}else{
	  siiyield = 0.0;
	}
#ifdef SANDBOX
	/* special frictional boundary condition for numerical sandbox benchmark */
	PetscScalar wallvx = -2.5/100.0/3600.0;
	if(markers[m].X < 0.002 || markers[m].X > grid->LX-0.012 + wallvx*elapsedTime || markers[m].Y > grid->LY-0.002){
	  PetscScalar c=0.0;
	  PetscScalar f=sin(19.0/180.0*M_PI);
	  siiyield = c + f*markers[m].p;
	}
#endif
	if(siiyield < 0.0) siiyield = 0.0;
	if( siinew > siiyield){
	  /* reduce viscosity */
	  PetscScalar exx = markers[m].e.T11;
	  PetscScalar eyy = markers[m].e.T22;
	  PetscScalar ezz = markers[m].e.T33;
	  PetscScalar exy = markers[m].e.T12;
	  PetscScalar exz = markers[m].e.T13;
	  PetscScalar eyz = markers[m].e.T23;
	  PetscScalar eiiold = sqrt(0.5*(exx*exx+eyy*eyy+ezz*ezz) + exy*exy + exz*exz + eyz*eyz);
	  etanew = siiyield/2.0/eiiold;
	  if( etanew/markers[m].eta < fractionalEtamin ){
	    etanew = fractionalEtamin*markers[m].eta ;
	  }
	  if(etanew < etamin){
	    etanew = etamin;
	  }else if(etanew > etamax){
	    etanew = etamax;
	  }
	  PetscScalar siiold = sqrt( 0.5*(markers[m].s.T11*markers[m].s.T11 +markers[m].s.T22*markers[m].s.T22 + markers[m].s.T33*markers[m].s.T33) + markers[m].s.T12*markers[m].s.T12+ markers[m].s.T13*markers[m].s.T13+ markers[m].s.T23*markers[m].s.T23);
	  if(siiold != 0.0){/* if siiold is zero, we expect yielding but cannot adjust old stresses downwards (they are zero already)*/
	    markers[m].s.T11 *= siiyield/siiold;
	    markers[m].s.T22 *= siiyield/siiold;
	    markers[m].s.T33 *= siiyield/siiold;
	    markers[m].s.T12 *= siiyield/siiold;
	    markers[m].s.T13 *= siiyield/siiold;
	    markers[m].s.T23 *= siiyield/siiold;
	  }
	isYielding[0]=1;
	} /* end if yielding */
	
	markers[m].eta = etanew;
      }
    }/* end if in domain*/
  }
  {
    /* do an allreduce on isYielding - if any cpu is still yielding, all should still think that they're yielding*/
    PetscInt tmp1=isYielding[0];
    ierr=MPI_Allreduce( &tmp1,isYielding,1,MPI_INT,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
  }
  
  PetscLogStagePop();
  PetscFunctionReturn(ierr);
}

PetscErrorCode projectMarkersNodesFromScalar(MarkerSet *markerset, GridData *grid, PetscScalar *markerField, Vec nodalField){
  /* projects a markerField onto the nodes (nodalField)*/
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt m;
  /* nodal fields*/
  PetscInt nMark = markerset->nMark; /* number of markers stored locally*/

  Marker *markers = markerset->markers;

  /* initialize weight sum array */
  Vec wtng, wtnl,nfl; /* global, local arrays*/
  ierr=DMCreateGlobalVector(grid->da, &wtng);CHKERRQ(ierr);
  ierr=DMCreateLocalVector(grid->da, &wtnl);CHKERRQ(ierr);
  ierr=VecDuplicate(wtnl, &nfl);CHKERRQ(ierr);
  ierr = VecZeroEntries( wtng);CHKERRQ(ierr);
  ierr = VecZeroEntries( wtnl);CHKERRQ(ierr);
  ierr = VecZeroEntries( nfl );CHKERRQ(ierr);
  ierr = VecZeroEntries( nodalField); CHKERRQ(ierr);
  PetscScalar **wtn;
  PetscScalar **nf;
  ierr = DMDAVecGetArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da,nfl,&nf);CHKERRQ(ierr);

  /*this loop calculates contribution of marker to each node in bounding box.*/
  for(m=0;m<nMark;m++){
    if(markers[m].cellX != -1){/* check to see if marker is in bounds */
      /* find cell i,j that I belong to.*/
      /*       ix          ix+1   */
      /* jy    o           o      */
      /*                          */
      /*            m             */
      /*                          */
      /* jy+1  o           o      */
      
      
      /* calculate contribution to i,j*/
      /*       PetscInt cellI = floor(markers[m].Y/dy); */  /* these are indices from 0 to N-1*/
      /*       PetscInt cellJ = floor(markers[m].X/dx); */
      PetscInt cellX = markers[m].cellX;
      PetscInt cellY = markers[m].cellY;

      /* calculate marker weight*/
      PetscScalar dx = grid->x[cellX+1] - grid->x[cellX];
      PetscScalar dy = grid->y[cellY+1] - grid->y[cellY];
      PetscScalar mdx = (markers[m].X - grid->x[cellX])/dx;
      PetscScalar mdy = (markers[m].Y - grid->y[cellY])/dy;

      /* contribution to cellX, cellY*/
      PetscScalar weight = (1.0 - mdx)*(1.0-mdy);
      wtn[cellY][cellX] += weight;
      nf[cellY][cellX] += weight*markerField[m];
      
      /*cellX+1, cellY*/
      weight = mdx*(1.0-mdy);
      wtn[cellY][cellX+1] += weight;
      nf[cellY][cellX+1]  += weight*markerField[m];
      
      /*cellX  , cellY+1*/
      weight = (1.0-mdx)*mdy;
      wtn[cellY+1][cellX] += weight;
      nf[cellY+1][cellX]  += weight*markerField[m];
      
      /*cellX+1, cellY+1*/
      weight = mdx*mdy;
      wtn[cellY+1][cellX+1] += weight;
      nf[cellY+1][cellX+1] += weight*markerField[m];

    }
  }
  ierr = DMDAVecRestoreArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,nfl,&nf);CHKERRQ(ierr);
  /* scatter arrays for weights-begin/end syntax is necessary to treat ghosted values properly*/
  ierr = DMLocalToGlobalBegin(grid->da,wtnl,ADD_VALUES,wtng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtnl,ADD_VALUES,wtng); CHKERRQ(ierr); 
  /* scatter arrays for field variables*/
  ierr = DMLocalToGlobalBegin(grid->da,nfl,ADD_VALUES,nodalField);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,nfl,ADD_VALUES,nodalField);CHKERRQ(ierr);
  /* normalize by weights */
  ierr = VecPointwiseDivide(nodalField,nodalField,wtng); CHKERRQ(ierr);

  /* destroy weight vectors*/
  ierr = VecDestroy(&wtnl); CHKERRQ(ierr);
  ierr = VecDestroy(&wtng); CHKERRQ(ierr);
  ierr = VecDestroy(&nfl); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode projectMarkersCellCentersFromScalar(MarkerSet *markerset, GridData *grid, PetscScalar *markerField, Vec nodalField){
  /* projects a markerField onto the nodes (nodalField)*/
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt m;
  /* nodal fields*/
  PetscInt nMark = markerset->nMark; /* number of markers stored locally*/

  Marker *markers = markerset->markers;

  /* initialize weight sum array */
  Vec wtng, wtnl,nfl; /* global, local arrays*/
  ierr=DMCreateGlobalVector(grid->da, &wtng);CHKERRQ(ierr);
  ierr=DMCreateLocalVector(grid->da, &wtnl);CHKERRQ(ierr);
  ierr=VecDuplicate(wtnl, &nfl);CHKERRQ(ierr);
  ierr = VecZeroEntries( wtng);CHKERRQ(ierr);
  ierr = VecZeroEntries( wtnl);CHKERRQ(ierr);
  ierr = VecZeroEntries( nfl );CHKERRQ(ierr);
  ierr = VecZeroEntries( nodalField); CHKERRQ(ierr);
  PetscScalar **wtn;
  PetscScalar **nf;
  ierr = DMDAVecGetArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da,nfl,&nf);CHKERRQ(ierr);

  /*this loop calculates contribution of marker to each node in bounding box.*/
  for(m=0;m<nMark;m++){
    if(markers[m].cellX != -1){/* check to see if marker is in bounds */
      /* find cell i,j that I belong to.*/
      /*       ix          ix+1   */
      /* jy    o           o      */
      /*                          */
      /*            m             */
      /*                          */
      /* jy+1  o           o      */
      
      
      /* calculate contribution to i,j*/
      /*       PetscInt cellI = floor(markers[m].Y/dy); */  /* these are indices from 0 to N-1*/
      /*       PetscInt cellJ = floor(markers[m].X/dx); */
      PetscInt cellX = markers[m].cellX;
      PetscInt cellY = markers[m].cellY;

      /* calculate marker weight*/
      PetscScalar dx = grid->x[cellX+1] - grid->x[cellX];
      PetscScalar dy = grid->y[cellY+1] - grid->y[cellY];

      PetscScalar mdx = fabs(markers[m].X - grid->xc[cellX+1])/dx/2.0;
      PetscScalar mdy = fabs(markers[m].Y - grid->yc[cellY+1])/dy/2.0;

      /* contribution to cellX+1, cellY+1*/
      PetscScalar weight = (1.0-mdx)*(1.0-mdy);
      wtn[cellY+1][cellX+1] += weight;
      nf[cellY+1][cellX+1] += weight*markerField[m];      
    }
  }
  ierr = DMDAVecRestoreArray(grid->da,wtnl,&wtn);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,nfl,&nf);CHKERRQ(ierr);
  /* scatter arrays for weights-begin/end syntax is necessary to treat ghosted values properly*/
  ierr = DMLocalToGlobalBegin(grid->da,wtnl,ADD_VALUES,wtng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtnl,ADD_VALUES,wtng); CHKERRQ(ierr); 
  /* scatter arrays for field variables*/
  ierr = DMLocalToGlobalBegin(grid->da,nfl,ADD_VALUES,nodalField);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(grid->da,nfl,ADD_VALUES,nodalField);CHKERRQ(ierr);
  /* normalize by weights */
  ierr = VecPointwiseDivide(nodalField,nodalField,wtng); CHKERRQ(ierr);

  /* destroy weight vectors*/
  ierr = VecDestroy(&wtnl); CHKERRQ(ierr);
  ierr = VecDestroy(&wtng); CHKERRQ(ierr);
  ierr = VecDestroy(&nfl); CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode projectNodalFieldToMarkersS(NodalFields *nodalFields, Vec nodalField, MarkerSet *markerset, PetscScalar *fieldptr, GridData *grid){
  /* this subroutine projects a scalar nodal field (1 dof) to the markers */
  const PetscInt stride = sizeof( Marker );
  /* get the locally owned part of the nodalField, including ghost values */
  Vec fieldl;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr=DMCreateLocalVector(grid->da,&fieldl); CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalField,INSERT_VALUES,fieldl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalField,INSERT_VALUES,fieldl);CHKERRQ(ierr);
  PetscScalar **field;
  ierr=DMDAVecGetArray( grid->da,fieldl,&field);CHKERRQ(ierr);

  Marker *markers = markerset->markers;
  /* begin loop over locally-owned markers*/
  PetscInt mm;
  for(mm=0;mm<markerset->nMark;mm++){
    /* check in bounds? */
    if(markers[mm].cellX != -1){
      PetscInt cellX = markers[mm].cellX;
      PetscInt cellY = markers[mm].cellY;

      /* calculate fracional offset of marker from upper-left corner of cell */
      PetscScalar mdx = (markers[mm].X - grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
      PetscScalar mdy = (markers[mm].Y - grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
      
      PetscScalar val = field[cellY][cellX]*(1-mdx)*(1-mdy) + field[cellY+1][cellX]*(1-mdx)*mdy + field[cellY][cellX+1]*(mdx)*(1-mdy) + field[cellY+1][cellX+1]*mdx*mdy;
      /* pointer to this marker's field */
      *(PetscScalar *)((unsigned char *) fieldptr + stride*mm) = val;
    }/* end if in bounds */  
  }/* end loop over markers */

  ierr=DMDAVecRestoreArray( grid->da,fieldl,&field);CHKERRQ(ierr);
  ierr=VecDestroy(&fieldl);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

PetscErrorCode projectVelocityToMarkers(MarkerSet *markerset, GridData *grid, NodalFields *nodalFields ){
  PetscInt m;

  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  PetscErrorCode ierr=0;
  PetscInt nMark = markerset->nMark;


  /*   PetscScalar *markerVX = markers->VX; */
  /*   PetscScalar *markerVY = markers->VY; */
  /*   PetscScalar *markerVZ = markers->VZ; */
  Marker *markers = markerset->markers;
  Vec vxl, vyl, vzl;
  PetscScalar **vx, **vy, **vz;
  
  PetscFunctionBegin;

  ierr=DMCreateLocalVector(grid->da,&vxl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &vyl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &vzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);

  ierr=DMDAVecGetArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vzl,&vz);CHKERRQ(ierr);


  for(m=0;m<nMark;m++){
    if( markers[m].cellX != -1){
      /*vx node indices*/
      PetscInt vxX, vxY, vyX, vyY, vzX, vzY;
      /*       PetscInt vxI, vxJ, vyI,  vyJ,vzI,vzJ; */
      /*       vxI = floor((markerY[m]-dy/2.0)/dy); */
      /*       vxJ = floor((markerX[m])/dx); */
      /*       vyI = floor((markerY[m])/dy); */
      /*       vyJ = floor((markerX[m]-dx/2.0)/dx); */
      vxX = markers[m].cellX;
      vxY = markers[m].cellY;
      if( vxY > 0 && markers[m].Y < grid->yc[vxY+1] ) vxY--;
      vyX = markers[m].cellX;
      if( vyX > 0 && markers[m].X < grid->xc[vyX+1] ) vyX--;
      vyY = markers[m].cellY;
      vzX = markers[m].cellX;
      /* vz cells are special - offset in both directions */
      if( vzX > 0 && markers[m].X < grid->xc[vzX+1] ) vzX--; 
      vzY = markers[m].cellY;
      if( vzY > 0 && markers[m].Y < grid->yc[vzY+1] ) vzY--;
      
      if( grid->xperiodic ){
	if( vyX < 0 ){ vyX = NX-1;}else if( vyX>NX-1){vyX = 0;}
	if( vxX < 0 ){ vxX = NX-1;}else if( vxX>NX-1){vxX = 0;}
      }else{
	if( vxX < 0){ vxX=0;} else if(vxX>NX-2){ vxX=NX-2;}
	if( vyX < 0){ vyX=0;} else if(vyX>NX-3){ vyX=NX-3;}
      }
      if( vxY < 0){ vxY=0;} else if(vxY>NY-3){ vxY=NY-3;}/* do not let a marker be assigned to one of the ghost cells*/
      if( vyY < 0){ vyY=0;} else if(vyY>NY-2){ vyY=NY-2;}
      
      if( vzY < 1){ vzY=1;} else if(vzY>NY-2){ vzY=NY-2;}
      if( vzX < 1){ vzX=1;} else if(vzX>NX-2){ vzX=NX-2;}
      
      /* compute new x-vel*/
      /* PetscScalar deltaxm=markerX[m] - dx*((PetscScalar) vxJ); */
      PetscScalar deltaxm = markers[m].X - grid->x[vxX];
      /*       PetscScalar deltaym=markerY[m] - dy*(((PetscScalar) vxI)+0.5); */
      PetscScalar deltaym = markers[m].Y - grid->yc[vxY+1];
      
      //      PetscInt vxdof = vxJ*NY + vxI;
      PetscScalar dx = grid->x[vxX+1]-grid->x[vxX];
      PetscScalar dy = grid->yc[vxY+2]-grid->yc[vxY+1];
      markers[m].VX = vx[vxY][vxX]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + vx[vxY][vxX+1]*(deltaxm/dx)*(1-deltaym/dy) + vx[vxY+1][vxX]*(1-deltaxm/dx)*(deltaym/dy) + vx[vxY+1][vxX+1]*(deltaxm*deltaym)/dx/dy; 
      /*compute new y-vel*/
      deltaxm = markers[m].X - grid->xc[vyX+1];
      deltaym = markers[m].Y - grid->y[vyY];
      dx = grid->xc[vyX+2]-grid->xc[vyX+1];
      dy = grid->y[vyY+1]-grid->y[vyY];
      //      PetscInt vydof = vyJ*NY + vyI;
      //      markerVY[m] = vy[vydof]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + vy[vydof+NY]*(deltaxm/dx)*(1-deltaym/dy) + vy[vydof+1]*(1-deltaxm/dx)*(deltaym/dy) + vy[vydof+NY+1]*(deltaxm*deltaym)/dx/dy; 
      markers[m].VY = vy[vyY][vyX]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + vy[vyY][vyX+1]*(deltaxm/dx)*(1-deltaym/dy) + vy[vyY+1][vyX]*(1-deltaxm/dx)*(deltaym/dy) + vy[vyY+1][vyX+1]*(deltaxm*deltaym)/dx/dy; 
      
      /* z interpolation*/
      dx=grid->xc[vzX+1]-grid->xc[vzX];
      dy=grid->yc[vzY+1]-grid->yc[vzY];
      deltaxm = markers[m].X - grid->xc[vzX];
      deltaym = markers[m].Y - grid->yc[vzY];
      markers[m].VZ = vz[vzY][vzX]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + vz[vzY][vzX+1]*(deltaxm/dx)*(1-deltaym/dy) + vz[vzY+1][vzX]*(1-deltaxm/dx)*(deltaym/dy) + vz[vzY+1][vzX+1]*(deltaxm*deltaym)/dx/dy; 
      //      deltaxm = markerX[m] - dx*((PetscScalar) vzJ - 0.5);/* pressure node x is (vzJ+1-0.5)*dx*/
      //      deltaym = markerY[m] - dy*((PetscScalar) vzI - 0.5);
      //      PetscInt vzdof = vzJ*NY+vzI;
      //      markerVZ[m] = vz[vzdof]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + vz[vzdof+NY]*(deltaxm/dx)*(1-deltaym/dy) + vz[vzdof+1]*(1-deltaxm/dx)*(deltaym/dy) + vz[vzdof+NY+1]*(deltaxm*deltaym)/dx/dy; 
    }
  }

  ierr=DMDAVecRestoreArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vzl,&vz);CHKERRQ(ierr);
  ierr=VecDestroy(&vxl);CHKERRQ(ierr);
  ierr=VecDestroy(&vyl);CHKERRQ(ierr);
  ierr=VecDestroy(&vzl);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}

/* marker advection routine - this is first order*/
PetscErrorCode advectMarkers(MarkerSet *markerset, GridData *grid, PetscScalar dt){/* GridData *grid, NodalFields *nodalFields ){ */
      /* 3. advect markers */
  PetscInt m;
  PetscErrorCode ierr=0;
  Marker *markers = markerset->markers;

  PetscFunctionBegin;
  
  for( m=0;m<markerset->nMark;m++){
    if(markers[m].cellX != -1){/* check markers in domain*/
      markers[m].X += markers[m].VX*dt; /* first order advection*/
      if( grid->xperiodic ){
	if( markers[m].X > grid->LX ){
	  markers[m].X -= grid->LX;
	}else if(markers[m].X < 0){
	  markers[m].X += grid->LX;
	}
      }
      markers[m].Y += markers[m].VY*dt; 
      markers[m].Z += markers[m].VZ*dt;
      /*  rotate stresses*/
      PetscScalar wxy = -markers[m].wxy * dt;/* minus sign because my wxy = 0.5*(dvx/dy - dvy/dx) */
      PetscScalar wxz = -markers[m].wxz * dt;/* expressions below are for wxy = 0.5*(dvy/dx-dvx/dy) */
      PetscScalar wyz = -markers[m].wyz * dt;
      
      PetscScalar msxx0 = markers[m].s.T11;//sxx;
      PetscScalar msyy0 = markers[m].s.T22;//yy;
      PetscScalar mszz0 = markers[m].s.T33;//zz;
      PetscScalar msxy0 = markers[m].s.T12;//xy;
      PetscScalar msxz0 = markers[m].s.T13;//xz;
      PetscScalar msyz0 = markers[m].s.T23;//yz;
      markers[m].s.T11 += 2.0*msxy0*-wxy + 2.0*msxz0*-wxz;
      markers[m].s.T22 += 2.0*msxy0*wxy + 2.0*msyz0*-wyz;
      markers[m].s.T33 += 2.0*msxz0*wxz + 2.0*msyz0*wyz;
      markers[m].s.T12 += (msxx0-msyy0)*wxy + msxz0*-wyz - wxz*msyz0;
      markers[m].s.T13 += (msxx0-mszz0)*wxz + msxy0*wyz - wxy*msyz0;
      markers[m].s.T23 += (msyy0-mszz0)*wyz + msxy0*wxz + wxy*msxz0;
      
      /*       markers.sxy[m] =msxx0*sin(2.0*espm) + msxy0*cos(2.0*espm); */
      /*       markers.sxx[m] =msxx0*(cos(espm)*cos(espm) - sin(espm)*sin(espm)) - msxy0*sin(2.0*espm); */
      
      /* 	if(markers.X[m] > grid.LX) markers.X[m] -= grid.LX; */
      /* 	if(markers.Y[m] > grid.LY) markers.Y[m] -= grid.LY; */
    }
  }	
  PetscFunctionReturn(ierr);
}

PetscErrorCode advectMarkersRK(MarkerSet *markerset, NodalFields *nodalFields, GridData *grid, Options *options, BoundaryValues *bv, PetscScalar dt){
  /* advect markers using 4th order runge-kutta method */
  PetscInt NX=grid->NX;
  PetscInt NY=grid->NY;
  PetscInt m;
  PetscErrorCode ierr=0;
  Marker *markers = markerset->markers;
  PetscFunctionBegin;
  setLogStage( LOG_MARK_MOVE );
  Vec vxl,vyl,vzl,wxyl,wxzl,wyzl;
  PetscScalar **vx, **vy, **vz, **wxy, **wxz, **wyz;
  ierr=DMCreateLocalVector(grid->da,&vxl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &vyl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &vzl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &wxyl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &wxzl);CHKERRQ(ierr);
  ierr=VecDuplicate( vxl, &wyzl);CHKERRQ(ierr);
  //ierr=VecDuplicate( vxl, &vzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);


  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->wxy,INSERT_VALUES,wxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->wxy,INSERT_VALUES,wxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->wxz,INSERT_VALUES,wxzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->wxz,INSERT_VALUES,wxzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->wyz,INSERT_VALUES,wyzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->wyz,INSERT_VALUES,wyzl);CHKERRQ(ierr);

  //ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);

  ierr=DMDAVecGetArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vzl,&vz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,wxyl,&wxy); CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,wxzl,&wxz); CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,wyzl,&wyz); CHKERRQ(ierr);

  PetscInt XCMAX, XCMIN;
  if( grid->xperiodic ){
    XCMAX=NX-1;
    XCMIN=-1;
  }else{
    XCMAX=NX-2;
    XCMIN=0;
  }

  for( m=0;m<markerset->nMark;m++){
    if(markers[m].cellX != -1){/* check marker in domain*/
      
      PetscScalar xcur = markers[m].X;
      PetscScalar ycur = markers[m].Y;
      
      PetscInt irk;
      PetscScalar mvx[4], mvy[4], mvz[4], mwxy[4], mwxz[4], mwyz[4];
      for( irk=0; irk<4;irk++){/* 4 sets order of advection scheme */
	
	/* initial guess at current cell */
	PetscInt cellX = markers[m].cellX;
	PetscInt cellY = markers[m].cellY;
	/* check to see if xcur,ycur has moved outside of basic cell (cellX, cellY) */
	
	if( xcur > grid->x[cellX+1]){
	  while(  xcur > grid->x[cellX+1] && cellX < XCMAX){ cellX++;}
	}else if( xcur < grid->x[cellX] ){
	  while( xcur < grid->x[cellX] && cellX > XCMIN){ cellX--; }
	}
	
	
	/* check to see if ycur is outside of basic cell */
	if( ycur > grid->y[cellY+1]){
	  while( ycur > grid->y[cellY+1] && cellY < NY-2){ cellY ++;}
	}else if( ycur < grid->y[cellY]){
	  while( ycur < grid->y[cellY] && cellY > 0){ cellY--;}
	}
	
	PetscInt vxX, vxY, vyX, vyY, vzX, vzY;
	
	vxX = cellX;
	vxY = cellY;
	if( vxY > -1 && ycur < grid->yc[vxY+1] ) vxY--;
	vyX = cellX;
	if( vyX > -1 && xcur < grid->xc[vyX+1] ) vyX--;
	vyY = cellY;
	//vzX = cellX;
	/* vz cells are special - offset in both directions */
	/* 	if( vzX > 0 && markers[m].X < grid->xc[vzX+1] ) vzX--;  */
	/* 	vzY = markers[m].cellY; */
	/* 	if( vzY > 0 && markers[m].Y < grid->yc[vzY+1] ) vzY--; */
	vzX = cellX;
	vzY = cellY;
	
	if( vzX < NX-1 && xcur >= grid->xc[vzX+1]) vzX++;
	if( vzY < NY-1 && ycur >= grid->yc[vzY+1]) vzY++;
	
	if( grid->xperiodic ){
	  // if( vyX < 0 ){ vyX = NX-1;}else if( vyX>NX-1){vyX = 0;}
	  // if( vxX < 0 ){ vxX = NX-1;}else if( vxX>NX-1){vxX = 0;}
	}else{
	  if( vxX < 0){ vxX=0;} else if(vxX>NX-2){ vxX=NX-2;}
	  if( vyX < -1){ vyX=-1;} else if(vyX>NX-2){ vyX=NX-2;}
	}
	
	if( vxY < -1){ vxY=-1;} else if(vxY>NY-2){ vxY=NY-2;}
	if( vyY < 0){ vyY=0;} else if(vyY>NY-2){ vyY=NY-2;}
	
	/* compute new x-vel*/
	/* PetscScalar deltaxm=markerX[m] - dx*((PetscScalar) vxJ); */
	PetscScalar deltaxm = xcur - grid->x[vxX];
	/*       PetscScalar deltaym=markerY[m] - dy*(((PetscScalar) vxI)+0.5); */
	PetscScalar deltaym = ycur - grid->yc[vxY+1];
	
	//      PetscInt vxdof = vxJ*NY + vxI;
	PetscScalar dx = grid->x[vxX+1]-grid->x[vxX];
	PetscScalar dy = grid->yc[vxY+2]-grid->yc[vxY+1];
	
	/* add implicit boundary conditions on top wall */
	PetscScalar vxul;
	PetscScalar vxur;
	
	if( vxY == -1 && bv->mechBCTop.type[0] == 0){/* kinematic bc*/
	  vxul = 2.0*bv->mechBCTop.value[0] - vx[vxY+1][vxX];
	  vxur = 2.0*bv->mechBCTop.value[0] - vx[vxY+1][vxX+1];
	}else if( vxY == -1 && bv->mechBCTop.type[0] == 1){/* free slip */
	  vxul = vx[vxY+1][vxX];
	  vxur = vx[vxY+1][vxX+1];
	}else{
	  vxul = vx[vxY][vxX];
	  vxur = vx[vxY][vxX+1];
	}
	    
	PetscScalar vxll = vx[vxY+1][vxX];
	PetscScalar vxlr = vx[vxY+1][vxX+1];

	mvx[irk]  = vxul*(1.0-deltaxm/dx)*(1.0-deltaym/dy); 
	mvx[irk] += vxur*(deltaxm/dx)*(1.0-deltaym/dy);
	mvx[irk] += vxll*(1.0-deltaxm/dx)*(deltaym/dy); 
	mvx[irk] += vxlr*(deltaxm*deltaym)/dx/dy; 

	/* Use y-ghost BC values */
	PetscScalar vyul;// = vy[vyY][vyX];
	PetscScalar vyur = vy[vyY][vyX+1];
	PetscScalar vyll;// = vy[vyY+1][vyX];
	PetscScalar vylr = vy[vyY+1][vyX+1];
	if( !grid->xperiodic && vyX == -1 && bv->mechBCLeft.type[1] == 0){/* kinematic */
	  vyul = 2.0*bv->mechBCLeft.value[1] - vy[vyY][vyX+1];
	  vyll = 2.0*bv->mechBCLeft.value[1] - vy[vyY+1][vyX+1];
	}else if( !grid->xperiodic && vyX == -1 && bv->mechBCLeft.type[1] == 1 ){/* free slip */
	  vyul = vy[vyY][vyX+1];
	  vyll = vy[vyY+1][vyX+1];
	} else {
	  vyul = vy[vyY][vyX];
	  vyll = vy[vyY+1][vyX];
	}

	/*compute new y-vel*/
	deltaxm = xcur - grid->xc[vyX+1];
	deltaym = ycur - grid->y[vyY];
	dx = grid->xc[vyX+2]-grid->xc[vyX+1];
	dy = grid->y[vyY+1]-grid->y[vyY];

	mvy[irk] =  vyul*(1.0-deltaxm/dx)*(1.0-deltaym/dy);
	mvy[irk] += vyur*(deltaxm/dx)*(1.0-deltaym/dy);
	mvy[irk] += vyll*(1.0-deltaxm/dx)*(deltaym/dy);
	mvy[irk] += vylr*(deltaxm*deltaym)/dx/dy; 

	PetscScalar vzul=vz[vzY][vzX];
	PetscScalar vzur;
	PetscScalar vzll;
	PetscScalar vzlr;
	if( vzX == NX-1 && vzY == NY-1){/* special case for lower right cell */
	  if( bv->mechBCRight.type[2] == 0){/* no-slip */
	    vzur = 2.0*bv->mechBCRight.value[2] - vzul;	    
	  } else if( bv->mechBCRight.type[2] == 1){/* free slip */
	    vzur = vzul;
	  }
	  vzlr = vzur;
	  
	  if( bv->mechBCBottom.type[2] == 0){
	    vzll =  2.0*bv->mechBCBottom.value[2] - vzul;
	  }else if( bv->mechBCBottom.type[2] == 1){
	    vzll = vzul;
	  }
	  
	} else if( vzX == NX-1){/* right boundary but not lower right */
	  vzll = vz[vzY+1][vzX];
	  if( bv->mechBCRight.type[2] == 0){/* no-slip */
	    vzur = 2.0*bv->mechBCRight.value[2] - vzul;	
	    vzlr = 2.0*bv->mechBCRight.value[2] - vzll;	    
	  } else if( bv->mechBCRight.type[2] == 1){/* free slip */
	    vzur = vzul;
	    vzlr = vzll;
	  }
	} else if(vzY == NY-1){/* bottom boundary but not lower right */
	  vzur = vz[vzY][vzX+1];
	  if( bv->mechBCBottom.type[2] == 0 ){
	    vzlr = 2.0*bv->mechBCBottom.value[2] - vzur;
	    vzll = 2.0*bv->mechBCBottom.value[2] - vzul;
	  }else if(bv->mechBCBottom.type[2] == 1){
	    vzlr = vzur;
	    vzll = vzul;
	  }
	} else { /* anywhere else in domain */
	  vzur = vz[vzY][vzX+1];
	  vzll = vz[vzY+1][vzX];
	  vzlr = vz[vzY+1][vzX+1];
	}
	/*compute new z-vel*/
	deltaxm = xcur - grid->xc[vzX];
	deltaym = ycur - grid->yc[vzY];
	dx = grid->xc[vzX+1]-grid->xc[vzX];
	dy = grid->yc[vzY+1]-grid->yc[vzY];

	mvz[irk] =  vzul*(1.0-deltaxm/dx)*(1.0-deltaym/dy);
	mvz[irk] += vzur*(deltaxm/dx)*(1.0-deltaym/dy);
	mvz[irk] += vzll*(1.0-deltaxm/dx)*(deltaym/dy);
	mvz[irk] += vzlr*(deltaxm*deltaym)/dx/dy; 

	/* compute vorticity tensor components by projecting from nodes */
	PetscScalar mdx = (xcur-grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
	PetscScalar mdy = (ycur-grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
	mwxy[irk] =  (1.0-mdx)*(1.0-mdy)*wxy[cellY][cellX];
	mwxy[irk] += (mdx)*(1.0-mdy)*wxy[cellY][cellX+1];
	mwxy[irk] += (1.0-mdx)*(mdy)*wxy[cellY+1][cellX];	
	mwxy[irk] += (mdx)*(mdy)*wxy[cellY+1][cellX+1];

	mwxz[irk] =  (1.0-mdx)*(1.0-mdy)*wxz[cellY][cellX];
	mwxz[irk] += (mdx)*(1.0-mdy)*wxz[cellY][cellX+1];
	mwxz[irk] += (1.0-mdx)*(mdy)*wxz[cellY+1][cellX];	
	mwxz[irk] += (mdx)*(mdy)*wxz[cellY+1][cellX+1];

	mwyz[irk] =  (1.0-mdx)*(1.0-mdy)*wyz[cellY][cellX];
	mwyz[irk] += (mdx)*(1.0-mdy)*wyz[cellY][cellX+1];
	mwyz[irk] += (1.0-mdx)*(mdy)*wyz[cellY+1][cellX];	
	mwyz[irk] += (mdx)*(mdy)*wyz[cellY+1][cellX+1];

	//% Update coordinates for the next Runge-Kutta cycle
	if(irk<3){
	  if(irk<2){
	    xcur=markers[m].X + dt/2*mvx[irk];
	    ycur=markers[m].Y + dt/2*mvy[irk];
	  }else{
	    xcur=markers[m].X + dt*mvx[irk];
	    ycur=markers[m].Y + dt*mvy[irk];
	  }
	}
	if(grid->xperiodic){
	  if( xcur < 0.0){
	    //	    xcur += grid->LX;
	    /* 	    markers[m].cellX = NX-1; */
	  } else if(xcur > grid->LX){
	    //	    xcur -= grid->LX;
	    /* 	    markers[m].cellX = 0; */
	  }
	}
      }/* end loop over rk order */
      //% Recompute velocity using 4-th order Runge_Kutta
      
      markers[m].VX  = (mvx[0]+2*mvx[1]+2*mvx[2]+mvx[3])/6;
      markers[m].VY  = (mvy[0]+2*mvy[1]+2*mvy[2]+mvy[3])/6;
      markers[m].VZ  = (mvz[0]+2*mvz[1]+2*mvz[2]+mvz[3])/6;
      markers[m].wxy = (mwxy[0] + 2*mwxy[1]+2*mwxy[2] + mwxy[3])/6.0;
      markers[m].wxz = (mwxz[0] + 2*mwxz[1]+2*mwxz[2] + mwxz[3])/6.0;
      markers[m].wyz = (mwyz[0] + 2*mwyz[1]+2*mwyz[2] + mwyz[3])/6.0;
      markers[m].X += dt*markers[m].VX;
      markers[m].Y += dt*markers[m].VY;
      markers[m].Z += dt*markers[m].VZ;

      /* rotate stress tensor */

      
      PetscScalar msxx0 = markers[m].s.T11;//sxx;
      PetscScalar msyy0 = markers[m].s.T22;//yy;
      PetscScalar mszz0 = markers[m].s.T33;//zz;
      PetscScalar msxy0 = markers[m].s.T12;//xy;
      PetscScalar msxz0 = markers[m].s.T13;//xz;
      PetscScalar msyz0 = markers[m].s.T23;//yz;
#ifdef JAUMANN
      PetscScalar wxy1 = markers[m].wxy * dt;/* minus sign because my wxy = 0.5*(dvx/dy - dvy/dx) */
      PetscScalar wxz1 = markers[m].wxz * dt;/* expressions below are for wxy = 0.5*(dvy/dx-dvx/dy) */
      PetscScalar wyz1 = markers[m].wyz * dt;
      markers[m].s.T11 += 2.0*msxy0*-wxy1 + 2.0*msxz0*-wxz1;
      markers[m].s.T22 += 2.0*msxy0*wxy1 + 2.0*msyz0*-wyz1;
      markers[m].s.T33 += 2.0*msxz0*wxz1 + 2.0*msyz0*wyz1;
      markers[m].s.T12 += (msxx0-msyy0)*wxy1 + msxz0*-wyz1 - wxz1*msyz0;
      markers[m].s.T13 += (msxx0-mszz0)*wxz1 + msxy0*wyz1 - wxy1*msyz0;
      markers[m].s.T23 += (msyy0-mszz0)*wyz1 + msxy0*wxz1 + wxy1*msxz0;
#else
      PetscScalar thetaz = -markers[m].wxy * dt; /* Angular velocity vector (w_x w_y w_z) <--> (-W_yz, W_xz, -W_xy) */
      PetscScalar thetay =  markers[m].wxz * dt;
      PetscScalar thetax = -markers[m].wyz * dt;
      /* rotate stresses using analytic formulae derived using mathematica */
      markers[m].s.T11 = mszz0*pow(sin(thetay),2.0) + msxz0*cos(thetaz)*sin(2.0*thetay) - 2.0*msyz0*cos(thetay)*sin(thetay)*sin(thetaz) + \
	pow(cos(thetay),2)*(msxx0*pow(cos(thetaz),2) - 2*msxy0*cos(thetaz)*sin(thetaz) + msyy0*pow(sin(thetaz),2));
      markers[m].s.T22 =mszz0*pow(cos(thetay),2)*pow(sin(thetax),2) - 2*msxz0*cos(thetay)*cos(thetaz)*pow(sin(thetax),2)*sin(thetay) + 
	pow(cos(thetaz),2)*sin(thetay)*(msxy0*sin(2*thetax) + msxx0*pow(sin(thetax),2)*sin(thetay)) - \
	2*msxy0*cos(thetaz)*pow(sin(thetax),2)*pow(sin(thetay),2)*sin(thetaz) + \
	pow(sin(thetax),2)*sin(thetaz)*(msyz0*sin(2*thetay) + msyy0*pow(sin(thetay),2)*sin(thetaz)) - \
	2*cos(thetax)*sin(thetax)*(sin(thetay)*sin(thetaz)*((-msxx0 + msyy0)*cos(thetaz) + msxy0*sin(thetaz)) + \
				   cos(thetay)*(msyz0*cos(thetaz) + msxz0*sin(thetaz))) + \
	pow(cos(thetax),2)*(msyy0*pow(cos(thetaz),2) + msxx0*pow(sin(thetaz),2) + msxy0*sin(2*thetaz));
      markers[m].s.T33 = (2*pow(cos(thetax),2)*(mszz0*pow(cos(thetay),2) - 2*msxz0*cos(thetay)*cos(thetaz)*sin(thetay) + \
				     msxx0*pow(cos(thetaz),2)*pow(sin(thetay),2) - 2*msxy0*cos(thetaz)*pow(sin(thetay),2)*sin(thetaz) + \
				     sin(thetaz)*(msyz0*sin(2*thetay) + msyy0*pow(sin(thetay),2)*sin(thetaz))) + \
	     pow(sin(thetax),2)*(msxx0 + msyy0 + (-msxx0 + msyy0)*cos(2*thetaz) + 2*msxy0*sin(2*thetaz)) + \
	     sin(2*thetax)*(2*cos(thetay)*(msyz0*cos(thetaz) + msxz0*sin(thetaz)) + \
			    sin(thetay)*(-2*msxy0*cos(2*thetaz) + (-msxx0 + msyy0)*sin(2*thetaz))))/2.;
      markers[m].s.T12 =(sin(thetax)*(4*cos(2*thetay)*(-(msxz0*cos(thetaz)) + msyz0*sin(thetaz)) + \
			 sin(2*thetay)*(msxx0 + msyy0 - 2*mszz0 + (msxx0 - msyy0)*cos(2*thetaz) - 2*msxy0*sin(2*thetaz))) + \
	    2*cos(thetax)*(2*sin(thetay)*(msyz0*cos(thetaz) + msxz0*sin(thetaz)) + \
			   cos(thetay)*(2*msxy0*cos(2*thetaz) + (msxx0 - msyy0)*sin(2*thetaz))))/4.;
      markers[m].s.T13 = sin(thetax)*(sin(thetay)*(msyz0*cos(thetaz) + msxz0*sin(thetaz)) + \
			 cos(thetay)*(msxy0*cos(2*thetaz) + (msxx0 - msyy0)*cos(thetaz)*sin(thetaz))) + \
	(cos(thetax)*(4*cos(2*thetay)*(msxz0*cos(thetaz) - msyz0*sin(thetaz)) - \
		      sin(2*thetay)*(msxx0 + msyy0 - 2*mszz0 + (msxx0 - msyy0)*cos(2*thetaz) - 2*msxy0*sin(2*thetaz))))/4.;
      markers[m].s.T23 = (2*cos(thetax)*sin(thetax)*(-(mszz0*pow(cos(thetay),2)) + pow(cos(thetaz),2)*(msyy0 - msxx0*pow(sin(thetay),2)) + \
					msxz0*cos(thetaz)*sin(2*thetay) - 2*msyz0*cos(thetay)*sin(thetay)*sin(thetaz) + \
					(msxx0 - msyy0*pow(sin(thetay),2))*pow(sin(thetaz),2) + msxy0*(1 + pow(sin(thetay),2))*sin(2*thetaz)) + \
	     pow(sin(thetax),2)*(-2*cos(thetay)*(msyz0*cos(thetaz) + msxz0*sin(thetaz)) + \
				   sin(thetay)*(2*msxy0*cos(2*thetaz) + (msxx0 - msyy0)*sin(2*thetaz))) + \
	     pow(cos(thetax),2)*(2*cos(thetay)*(msyz0*cos(thetaz) + msxz0*sin(thetaz)) + \
				   sin(thetay)*(-2*msxy0*cos(2*thetaz) + (-msxx0 + msyy0)*sin(2*thetaz))))/2.;
#endif
#ifdef TEXTURE
      /* calculate marker rotation rate in terms of crystal theta, phi */
      /* goal: compute thetadot and phidot for macroscopic rotation in terms of theta, phi (crystal) and marker wij */
      /* this code ONLY works if JAUMANN is not set */
      PetscInt c1,c2,c3;
      for(c1=0;c1<NT;c1++){
	for(c2=0;c2<NT;c2++){
	  for(c3=0;c3<NT;c3++){
	    /* calculate vector c */
	    PetscScalar ct =  markers[m].texture.ctheta[c1][c2][c3] ;
	    PetscScalar cp =  markers[m].texture.cphi[c1][c2][c3] ;
	    /* cpx, cpy, cpz are new orientations after rotation through theta_i */
	    PetscScalar cpx = cos(thetay)*cos(ct + thetaz)*sin(cp) + cos(cp)*sin(thetay);
	    PetscScalar cpy = -(cos(cp)*cos(thetay)*sin(thetax)) + sin(cp)*(cos(ct + thetaz)*sin(thetax)*sin(thetay) + cos(thetax)*sin(ct + thetaz));
	    PetscScalar cpz = cos(cp)*cos(thetax)*cos(thetay) + sin(cp)*(-(cos(thetax)*cos(ct + thetaz)*sin(thetay)) + sin(thetax)*sin(ct + thetaz));
	    /* recover new theta and phi from cp */
	    cp = acos( cpz );
	    ct = atan2( cpy, cpx );
	    /* ctheta and cphi are guaranteed to be within [-pi,pi] and [0, pi] */
	    markers[m].texture.ctheta[c1][c2][c3] = ct;
	    markers[m].texture.cphi[c1][c2][c3]  = cp;
	  }
	}
      }
#endif
      /* update total plastic strain */
      markers[m].Eii += dt*sqrt( 0.5*(markers[m].e.T11*markers[m].e.T11 + markers[m].e.T22*markers[m].e.T22 + markers[m].e.T33*markers[m].e.T33) + markers[m].e.T12*markers[m].e.T12 + markers[m].e.T13*markers[m].e.T13 + markers[m].e.T23*markers[m].e.T23 );

      if(grid-> xperiodic){
	if( markers[m].X > grid->LX ){ 
	  markers[m].X -= grid->LX;
	}else if(markers[m].X < 0){
	  markers[m].X += grid->LX;
	}
      }
    }
  }

  ierr=DMDAVecRestoreArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vzl,&vz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,wxyl,&wxy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,wxzl,&wxz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,wyzl,&wyz);CHKERRQ(ierr);
  ierr=VecDestroy(&vxl);
  ierr=VecDestroy(&vyl);
  ierr=VecDestroy(&vzl);
  ierr=VecDestroy(&wxyl);
  ierr=VecDestroy(&wxzl);
  ierr=VecDestroy(&wyzl);
  PetscLogStagePop();
  PetscFunctionReturn(ierr);

}



#ifdef done
void printMarkersAllASCIIMatlab( Markers *markers, PetscInt iTime){
  /* OUTPUT OF MARKERS*/
  FILE *mf;
  char markName[80];
  PetscInt i;
  sprintf(markName,"output/loadMarkers%d.m",iTime);
  mf=fopen(markName,"w");
  fprintf(mf,"%%X\tY\tVX\tVY\tT\tEta\texx\texy\tsxx\tsxy\tp\twxy\tMaterialId\n");
  fprintf(mf,"markerInfo=[");
  for(i=0; i<markers->nMark;i++){
    fprintf(mf,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d;\n",markers->X[i],markers->Y[i],markers->VX[i],markers->VY[i],markers->T[i],markers->eta[i],markers->exx[i],markers->exy[i],markers->sxx[i],markers->sxy[i],markers->p[i],markers->wxy[i],markers->Mat[i]);
  }
  fprintf(mf,"];");
  fclose(mf);
}
#endif
