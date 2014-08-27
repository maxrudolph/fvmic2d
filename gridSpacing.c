
#include "fdcode.h"
#include "gridSpacing.h"
#include "math.h"

/* This file contains subroutines that generate uniform and non-uniform coordinates to use in constructing grids */

void gridSpacingUniform( PetscScalar *x, PetscScalar LX, PetscInt NX){
  PetscInt ix;
  PetscScalar dx = LX/((PetscScalar) NX - 1.0);
  x[0] = 0.0;
  for(ix=1;ix<NX-1;ix++){
    x[ix] = 0.0 + dx*((PetscScalar) ix);
  }
  x[NX-1] = LX;
}

void smoothGrid( PetscScalar *x, PetscScalar NX, PetscInt nSmooth){
  PetscInt ix;
  for(ix=0;ix<NX;ix++){
    PetscInt jx;
    /* average x value over [ix-nSmooth ix+nSmooth] */
    PetscInt nxm = nSmooth;
    PetscInt nxp = nSmooth;
    if( ix < nxm ) nxm = ix;
    if( ix + nxp >= NX ) nxp = NX-1 - ix;
    PetscScalar val = 0.0;
    for(jx=-nxm;jx<=nxp;jx++){
      val += x[ix+jx]/(nxm+nxp+1.0);
    }
    x[ix] = val;
  }
}

void circShift( PetscScalar *x, PetscScalar LX, PetscInt NX, PetscInt shift){
  /* shift vector of gridline locations circularly*/
  /* cache last value in array */
  PetscScalar *dx; 

  PetscMalloc( (NX-1)*sizeof(PetscScalar), &dx);
  PetscInt ix;
  for(ix = 0; ix < NX-1; ix ++){
    dx[ix] = x[ix+1]-x[ix];
  }
  if( shift < 0){
    PetscInt ishift;
    for(ishift =1;ishift <= -shift;ishift++){
      PetscInt ix;
      PetscScalar tmp = dx[0];
      for(ix=0;ix<NX-2;ix++){
	dx[ix] = dx[ix+1];
      }
      dx[NX-2] = tmp;
    }
  }else{
    PetscInt ishift;
    for(ishift =1;ishift <= shift;ishift++){
      PetscInt ix;
      PetscScalar tmp = dx[NX-2];
      for( ix = NX-2; ix>0; ix--){
	dx[ix] = dx[ix-1];
      }
      dx[0] = tmp;
    }
  }
  x[0] = 0.0;
  for(ix=1;ix<NX-1;ix++){
    x[ix] = x[ix-1] + dx[ix];
  }
  x[NX-1] = LX;
  PetscFree( dx );
}

void gridSpacingRefinedCenter( PetscScalar *x, PetscScalar LX, PetscInt NX, PetscScalar hcontrast){
  if(!( NX % 2 )){
    printf("N must be odd for this mesh\n");
    abort();
  }
  PetscScalar L = LX/2.0;
  PetscInt NE = (NX-1)/2;
  PetscInt Nr=NE/hcontrast;/* number of elements in refined region */
  PetscInt Nc=NE-Nr-5;
  PetscScalar a = L/(10.0+Nr+Nc*hcontrast);
  PetscInt ix;
  PetscScalar f=1.2407234;

  x[0] = 0.0;
  for(ix=1;ix<=Nc;ix++){
    x[ix] = x[ix-1] + a*hcontrast;
  }
  for(ix=Nc+1;ix<=Nc+5;ix++){
    x[ix] = x[ix-1] + a*pow(f,6.0-(ix-Nc));
  }
  for(ix=Nc+6;ix<=Nc+Nr+6;ix++){
    x[ix] = x[ix-1] + a;
  }
  
  /* mirror other half of grid */
  x[NE] = LX/2;
  for(ix=(NX+1)/2;ix<NX-1;ix++){
    x[ix] = LX - x[NX-1 - ix];
  }
  x[NX-1] = LX;
}

void gridSpacingConstantInnerOuter( PetscScalar *x, PetscScalar LX, PetscInt NX, PetscScalar hcontrast){
  /* Make a grid spaced so that the outer and inner portions have constant spacing set by hcontrast with a smooth transition between the two regions */
  PetscScalar f=1.2407234;
  PetscScalar L = LX/2.0;
  PetscInt NE = (NX-1)/2;
  PetscInt Nr=NE/hcontrast;/* number of elements in refined region */
  PetscInt Nc=NE-Nr-5;
  PetscScalar a = L/(10.0+Nr+Nc*hcontrast);

  PetscInt ix;
  x[0] = 0.0;
  for(ix=1;ix<=Nr;ix++){
    x[ix] = x[ix-1] + a;
  }
  for(ix=Nr+1;ix<=Nr+5;ix++){
    x[ix] = x[ix-1] + a*pow(f,ix-Nr);
  }
  for(ix=Nr+6;ix<=Nc+Nr+6;ix++){
    x[ix] = x[ix-1] + hcontrast*a;
  }
  /* mirror other half of grid */
  x[NE] = LX/2;
  for(ix=(NX+1)/2;ix<NX-1;ix++){
    x[ix] = LX - x[NX-1 - ix];
  }
  x[NX-1] = LX; 
}

void getCellCenters( PetscScalar *x, PetscScalar *xc, PetscScalar LX, PetscInt NX){
  PetscInt ix;
  for(ix=1;ix<NX;ix++){
    xc[ix] = (x[ix-1] + x[ix])/2.0;
  }
  xc[0] = -xc[1];
  xc[NX] = 2.0*LX - xc[NX-1];  
}

void makePeriodic( PetscScalar *x, PetscScalar *xc, PetscScalar LX, PetscInt NX){
  x[NX] = LX;
  x[NX+1] = LX+x[1];
  x[-1] = x[NX-1]-LX;
  x[-2] = x[NX-2]-LX;
  xc[NX] = (LX + x[NX-1])/2.0;
  xc[NX+1] = xc[1]+LX;
  xc[0] = 0.0-(LX-xc[NX]);
  xc[-1] = 0.0-(LX-xc[NX-1]);
}
