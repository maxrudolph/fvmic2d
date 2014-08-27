#include "fdcode.h"
#include "texture.h"
#include "math.h"


PetscScalar psmax( PetscScalar a, PetscScalar b ){
  PetscScalar val= ((a) > (b)) ? (a) : (b) ;
  return(val);
}

PetscErrorCode limitAnisotropicViscosity(MarkerSet *markerset,Materials *materials, Options *options){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  /* make a mask array to identify markers with viscoplasticity tensors that are too large */
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar etamax = options->etamax;
  for(m=0;m<markerset->nMark;m++){
    if( markers[m].cellX != -1 ){
      if(materials->hasTexture[(PetscInt) markers[m].Mat]){
	/* look at all viscoplasticity tensor components*/
	PetscScalar maxN = psmax(fabs(markers[m].texture.N[0][0]),fabs(markers[m].texture.N[0][1]));
	maxN = psmax(maxN,fabs(markers[m].texture.N[1][0]));
	maxN = psmax(maxN,fabs(markers[m].texture.N[1][1]));
	if(isnan(maxN) || isinf(maxN) || maxN == 0.0){
	  printf("problem with maxN = %e\n",maxN); abort();
	}
	if(maxN > etamax){
#ifdef debug1
	  printf("Max N before viscosity adjustment %e\n",maxN);
#endif
	  /* scale M, the viscoplasticity tensor for this marker */
	  maxN = etamax/maxN;
	  PetscInt i,j;
	  for(i=0;i<6;i++){
	    for(j=0;j<6;j++){
	      /* modify the viscoplasticity tensor by scaling all components*/
	      markers[m].texture.M[i][j] /= maxN;
	      /* modify the single crystal rotation rates by scaling */
	      PetscInt ix,iy,iz;
	      for(ix=0;ix<NT;ix++){
		for(iy=0;iy<NT;iy++){
		  for(iz=0;iz<NT;iz++){
		    markers[m].texture.cthetadot[ix][iy][iz] *= maxN;
		    markers[m].texture.cphidot[ix][iy][iz] *= maxN;
		  }
		}
	      }
	    }
	  }
	  /* invert the scaled M to N */
	  invertMtoN(&markers[m].texture);
#ifdef debug1
	  maxN = psmax(fabs(markers[m].texture.N[0][0]),fabs(markers[m].texture.N[0][1]));
	  maxN = psmax(maxN,fabs(markers[m].texture.N[1][0]));
	  maxN = psmax(maxN,fabs(markers[m].texture.N[1][1]));
	  printf("Max N after viscosity adjustment %e\n",maxN);
#endif
	  
	  
	  /*       markers[m].texture.N[0][0] /= maxN; */
	  /*       markers[m].texture.N[0][1] /= maxN; */
	  /*       markers[m].texture.N[1][0] /= maxN; */
	  /*       markers[m].texture.N[1][1] /= maxN; */		    
	}
      }
    }
  }
  PetscFunctionReturn(ierr);
}
