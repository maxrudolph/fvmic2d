#include "fdcode.h"

void gridMarkerFieldS( FILE *, PetscInt , PetscScalar *, PetscInt, PetscInt, PetscInt *, PetscInt *, PetscScalar *);
void gridMarkerFieldI( FILE *, PetscInt , PetscScalar *, PetscInt, PetscInt, PetscInt *, PetscInt *, PetscScalar *);
void getWeights( PetscInt, PetscScalar *, PetscScalar *, PetscInt *, PetscInt *, PetscScalar *, PetscScalar , PetscScalar );


/* calling syntax should be input_file, LX, LY, NX, NY */
int main(int argc, char **args){
  PetscErrorCode ierr;
  if(argc != 6){/* Check for correct calling syntax */
    printf("Wrong number of arguments\n");
    return(-2);
  }
  char *ifn = args[1];
  printf("reading from %s\n",ifn);

  char ofn[160];
  sprintf( &ofn[0], "%s.gridded", ifn);

  PetscScalar LX, LY, elapsedTime;
  PetscInt NX, NY, nMark;
  sscanf(args[2],"%le",&LX);/* Convert text arguments to numbers */
  sscanf(args[3],"%le",&LY);
  sscanf(args[4],"%d",&NX);
  sscanf(args[5],"%d",&NY);
  printf("using LX = %e, LY = %e, NX = %d, NY = %d\n",LX,LY,NX,NY);

  FILE *ifile, *ofile;
  ifile = fopen( ifn, "rb"); /* Open input and output files */
  ofile = fopen( ofn, "wb");
  printf("writing to %s\n",ofn);
  /* read the file */
  
  fread( &nMark, sizeof(PetscInt),1,ifile);/* number of markers*/
  fread( &elapsedTime, sizeof(PetscScalar),1,ifile); /*elapsed time */
  printf("reading %d markers, elapsedTime = %e\n",nMark,elapsedTime);
  /* write NX, NY, elapsedTime */
  fwrite( &NX, sizeof(PetscInt), 1, ofile);
  fwrite( &NY, sizeof(PetscInt), 1, ofile);
  fwrite( &elapsedTime, sizeof(PetscScalar), 1, ofile);

  /* allocate scalar field to store marker information */
  PetscScalar *sfield;
  ierr = PetscMalloc( NX*NY*sizeof(PetscScalar), &sfield);

  fseek( ifile, 3*sizeof(PetscInt)*nMark, SEEK_CUR);/* cpu, cellx, celly*/

  /* allocate arrays for X and Y */
  PetscScalar *x, *y;
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &x);
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &y);
  fread( x, sizeof(PetscScalar), nMark, ifile);
  fread( y, sizeof(PetscScalar), nMark, ifile);
  /* compute node nearest to each marker */
  PetscInt *nodeX, *nodeY;
  ierr = PetscMalloc( nMark*sizeof(PetscInt), &nodeX); /* stores number of closest node */
  ierr = PetscMalloc( nMark*sizeof(PetscInt), &nodeY);
  PetscInt m;
  PetscScalar dx = LX/(NX-1);
  PetscScalar dy = LY/(NY-1);
  printf("calculating marker nodes, dx, dy = %e, %e\n",dx,dy);
  for(m=0;m<nMark;m++){
    if( x[m] < 0.0 || y[m] < 0.0 || x[m] > LX || y[m] > LY){
      nodeX[m] = -1;
      nodeY[m] = -1;
    }else{      
      PetscInt cellx = (PetscInt)(x[m]/dx);
      PetscInt celly = (PetscInt)(y[m]/dy);/* round down */  
      if( x[m] > (cellx+0.5)*dx ) cellx++;/* push to right */
      if( y[m] > (celly+0.5)*dy ) celly++;/* push down */
      if( cellx < 0 ){ cellx = 0;}else if(cellx > NX-1) cellx = NX-1;
      if( celly < 0 ){ celly = 0;}else if(celly > NY-1) celly = NY-1;
      nodeX[m] = cellx;
      nodeY[m] = celly;
      //    printf("m %d, x = %e, y = %e, cellx %d, celly, %d\n",m,x[m], y[m],cellx, celly);
    }
  }

  PetscScalar *wt;
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &wt);
  getWeights( nMark, x, y, nodeX, nodeY, wt, dx, dy);

  /* Z */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* vx */
  printf("reading VX\n");
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  printf("writing VX\n");
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* vY */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* vZ */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* material */
  gridMarkerFieldI( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);/* T */
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);/* Tdot */
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* eta */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* mu */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* TEXTURE */
#ifdef TEXTURE
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
#endif
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);/* D */
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);/* Ddot */
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* exx */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* exy */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* Eii */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* sxx */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* syy */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* szz */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* sxz */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* syz */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* sxy */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* p */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* rho */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);
  /* rhodot */
  gridMarkerFieldS( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);

  /* end reading fields in sequence */

  /* CPU # */
  /* go back to beginning of file and seek to cpu # */
  fseek( ifile, 3*sizeof(PetscInt) + 2*sizeof(PetscScalar), SEEK_SET);/* cellx, celly*/

  gridMarkerFieldI( ifile, nMark, sfield, NX, NY, nodeX, nodeY, wt);
  fwrite( sfield, sizeof(PetscScalar), NX*NY, ofile);

  fclose(ifile);
  fclose(ofile);
  ierr = PetscFree(x);
  ierr = PetscFree(y);
  ierr = PetscFree(sfield);
  ierr = PetscFree(wt);

}

void getWeights( PetscInt nMark, PetscScalar *x, PetscScalar *y, PetscInt *cellx, PetscInt *celly, PetscScalar *wt, PetscScalar dx, PetscScalar dy){
  PetscInt m;
  for(m=0;m<nMark;m++){
    PetscScalar xn = cellx[m]*dx;
    PetscScalar yn = celly[m]*dy;
    PetscScalar mdx = fabs(x[m]-xn)/(dx/2.0);
    PetscScalar mdy = fabs(y[m]-yn)/(dy/2.0); 
    wt[m] = sqrt( mdx*mdx + mdy*mdy );    
  }
}


/* routine to read marker field and grid it */
void gridMarkerFieldS( FILE *ifile, PetscInt nMark, PetscScalar *field, PetscInt NX, PetscInt NY, PetscInt *cellx, PetscInt *celly, PetscScalar *wt){
  /* read the field*/
  PetscScalar *mfield;
  PetscMalloc( nMark*sizeof(PetscScalar), &mfield);

  /* read marker field from input file*/
  fread( mfield, sizeof(PetscScalar), nMark, ifile);

  PetscInt m;
  /* zero out nodal field array */
  for(m=0;m<NX*NY;m++) field[m] = 0.0; 
  for(m=0;m<NX*NY;m++) wt[m] = 0.0; 
  /* loop through markers, project marker value onto nearest node */
  for(m=0;m<nMark;m++){
    if( cellx[m] != -1 ){
      
      field[ celly[m]*NX + cellx[m] ] += mfield[m];
      wt[ celly[m]*NX + cellx[m]] += 1.0;
    }
  }
  for(m=0;m<NX*NY;m++){
    field[m] /= wt[m];
  }
  /* check for NaN */
  {
    int ix,jy;
    for(ix=0;ix<NX;ix++){
      for(jy=0;jy<NY;jy++){
	int idx = ix+jy*NX;
	if(isnan(field[idx])){
	     int sx=ix;
	     while(ix<NX-1){
	       if( !isnan( field[sx+jy*NX] ) ){
		 field[idx] = field[sx+jy*NX];
		 break;
	       }
	       sx++;
	     }   
	}
      }
    }
  }
  
  PetscFree(mfield);
  
}

/* routine to read marker field and grid it */
void gridMarkerFieldI( FILE *ifile, PetscInt nMark, PetscScalar *field, PetscInt NX, PetscInt NY, PetscInt *cellx, PetscInt *celly, PetscScalar *wt){
  /* read the field*/
  PetscInt *mfield;
  PetscMalloc( nMark*sizeof(PetscInt), &mfield);

  /* read marker field from input file*/
  fread( mfield, sizeof(PetscInt), nMark, ifile);

  PetscInt m;
  /* zero out nodal field array */
  for(m=0;m<NX*NY;m++) field[m] = 0.0; 
  for(m=0;m<NX*NY;m++) wt[m] = 0.0; 
  /* loop through markers, project marker value onto nearest node */
  for(m=0;m<nMark;m++){
    if(cellx[m] != -1){
      field[ celly[m]*NX + cellx[m] ] += (PetscScalar) mfield[m];
      wt[ celly[m]*NX + cellx[m]] += 1.0;
    }
  }
  for(m=0;m<NX*NY;m++){
    field[m] /= wt[m];
  }

  PetscFree(mfield);

}
