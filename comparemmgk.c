#include "fdcode.h"
#include "texture.h"
#include "petscksp.h"
#include "petsctime.h"
#include "time.h"
#include "tensorOps.h"
#include "viscosity.h"

#define MM 10

/* compare micro macro rheology with goldsby and kohlstedt */

int main(int argc, char **args){
  PetscInitialize(&argc,&args,(char *)0,"help message");
  srand48( (unsigned) time( NULL));
  TextureInfo tex[MM];
  PetscErrorCode ierr;
  const PetscScalar pi = M_PI;
  printf("made it here\n");
  MarkerSet markerset;
  Marker *markers;
  Options options;
  Materials materials;
  PetscInt nTime = 5000;
  /* set grain size */
  options.grainSize = 1e-3;/* grain size in meters*/
  options.thermalBCBottom.value[0] = 273;

  materials.hasTexture[0] = 1;
  materials.hasTexture[1] = 0;

  materials.hasEtaT[1] = 4;
  options.textureDevelopment = 1;
  /* initialize markers*/
  PetscFunctionBegin;
  ierr=PetscMalloc( MM*sizeof(Marker), &markers );CHKERRQ(ierr);

  markerset.markers=&markers[0];
  markerset.nMark = 2;

  PetscScalar Tmin = 100; PetscScalar Tmax = 260;
  PetscInt nT = 50;

  PetscInt istress;
  PetscInt nstress = 50;
  PetscInt iT;
  printf("grainsize: %e m\n",options.grainSize);
  printf("T\t sii\t\t eii_mm\t\t eii_gk\t\t eta_mm\t\t eta_gk\n");

  /* loop over temperature */
  for(iT = 0; iT< nT; iT++){
    PetscScalar T = (Tmax-Tmin)*((PetscScalar) iT)/((PetscScalar) nT-1) + Tmin;
    
    for(istress=0;istress<nstress;istress++){  /* loop over stresses */
      PetscScalar s = pow(10.0,((PetscScalar) istress)/((PetscScalar)nstress-1.0)*5.0);
      /* set stress state*/
      PetscInt m;
      for(m=0;m<markerset.nMark;m++){
	markers[m].s.T11 = s;
	markers[m].s.T12 = 0.0;
	markers[m].s.T22 = s;
	markers[m].s.T23 = 0.0;
	markers[m].s.T13 = 0.0;
	markers[m].s.T33 = 0.0;
	markers[m].T = T;
	markers[m].p = 1e6;
	markers[m].wxy = 0.0;
	markers[m].wxz = 0.0;
	markers[m].wyz = 0.0;
	markers[m].Mat = 0;
      }
      markers[1].Mat = 1;
      PetscScalar sii = IIdev(&( markers[0].s ));
      /* randomize marker texture distribution*/
      initializeIsotropicFabric( &markerset );
      
      PetscInt iTime=0;
      PetscScalar elapsedTime = 0.0;
      /* form gk viscosity */ 
      updateMarkerViscosity( &markers[1], &options, &materials, sii );
      Tensor33s D;
      /* calculate mm viscosity */
      formViscoplasticMMarker( &markers[0], &options, &D, 1 );      
      PetscScalar eii = sqrt(0.5*( D.T11*D.T11 + D.T22*D.T22 + D.T33*D.T33 ) + D.T12*D.T12 + D.T23*D.T23 + D.T13*D.T13);
      printf("%e\t%e\t %e\t %e\t %e\t %e\n",T,sii,eii,sii/2.0/markers[1].eta,sii/2.0/eii,markers[1].eta);
    }
  }
  
  

  
  /* compute sii*/
  
  
  
  return(0);
}

PetscScalar drand(){
  return ((PetscScalar) rand()) / ((PetscScalar) RAND_MAX);
}
