#include "fdcode.h"
#include "texture.h"
#include "petscksp.h"
#include "petsctime.h"
#include "time.h"
#include "tensorOps.h"

#define MM 1

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
  PetscInt nTime = 10001;
  /* set grain size */
  options.grainSize = 10.0;/* grain size in meters*/
  options.thermalBCBottom.value[0] = 260;

  materials.hasTexture[0] = 1;
  options.textureDevelopment = 1;

  /* initialize markers*/
  PetscFunctionBegin;
  ierr=PetscMalloc( MM*sizeof(Marker), &markers );CHKERRQ(ierr);

  markerset.markers=&markers[0];
  markerset.nMark = 1;

  /* set stress state*/
  PetscInt m;
  for(m=0;m<markerset.nMark;m++){
    markers[m].s.T11 = 1e5;
    markers[m].s.T12 = 0.0;
    markers[m].s.T22 = -0.5e5;
    markers[m].s.T23 = 0.0;
    markers[m].s.T13 = 0.0;
    markers[m].s.T33 = -0.5e5;
    markers[m].T = 200.0;
    markers[m].wxy = 0.0;
    markers[m].wxz = 0.0;
    markers[m].wyz = 0.0;
    markers[m].Mat = 0;
    markers[m].cellX = 1;
    markers[m].cellY = 1;
    markers[m].texture.d = options.grainSize;
  }
  markers[0].Mat = 0;
  PetscScalar sii = IIdev(&( markers[0].s ));
  /* randomize marker texture distribution*/
  initializeIsotropicFabric( &markerset );

  PetscInt iTime;
  PetscScalar elapsedTime = 0.0;
  PetscScalar Eii = 0.0; /* total plastic strain */
  
  {
    /* dump the fabric to a file */
    char filename[80];
    sprintf(filename,"output/texture.%d.dat",-1);
    FILE *fh;
    fh = fopen(filename,"w");
    Marker M = markerset.markers[0];
    /* print stress second invariant, strain rate second invariant */
    fprintf(fh,"NT=%d\n",NT);
    fprintf(fh,"eii=%e\n",0.0);
    fprintf(fh,"Eii=%e\n",0.0);
    fprintf(fh,"eta=%e\n",0.0);
    PetscInt i,j,k;
    for(i=0;i<NT;i++){
      for(j=0;j<NT;j++){
	  for(k=0;k<NT;k++){
	    fprintf(fh,"%e\t%e\n",M.texture.ctheta[i][j][k],M.texture.cphi[i][j][k]);
	  }
      }
    }
    fclose(fh);
  }    
  
  for(iTime=0;iTime < nTime;iTime++){
    printf("beginning timestep %d\n",iTime);
    //formViscoplasticMNewton( &markerset, &materials, &options);
    /* compute crystal rotation rates */
    Tensor33s D;/* deviatoric strain rate tensor */    
    for(m=0;m<markerset.nMark;m++){
      //PetscErrorCode formViscoplasticM(Marker *,Options *,Tensor33s *, PetscInt );
      Marker *M = &(markerset.markers[m]);
      formViscoplasticMMarker( M, &options, &D, 1 );
      
      
#ifdef verbose
      PetscInt i,j,k;
      for(i=0;i<NT;i++){
	for(j=0;j<NT;j++){
	  for(k=0;k<NT;k++){
	    printf("rotation: [%d][%d][%d] = (%e,%e)\n",i,j,k,M->texture.cthetadot[i][j][k],M->texture.cphidot[i][j][k]);
	  }
	}
      }
#endif
    }
    
    /* update the fabric */
    PetscScalar dt = 3.15e11;
    /* determine maximum timestep for fabric change */
    PetscScalar dtf = getFabricDT( &markerset );
    printf("maximum dt for fabric = %e\n",dtf);
    //if( dtf < dt) dt = dtf;
    updateFabric( &markerset, &materials, dt);
    elapsedTime += dt;
    PetscScalar eii = sqrt(0.5*( D.T11*D.T11 + D.T22*D.T22 + D.T33*D.T33 ) + D.T12*D.T12 + D.T23*D.T23 + D.T13*D.T13);

    Eii += eii*dt;
    PetscScalar Emax = 10.0;
    //#ifdef verbose
    if( !(iTime %100) || Eii>Emax ){
      for(m=0;m<markerset.nMark;m++){
	invertMtoN( &( markers[m].texture ) );
	printf("Marker %d:\n",m);
	printf("M tensor:\n");
	PetscInt i,j;
	for(i=0;i<6;i++){
	  for(j=0;j<6;j++){
	    printf("%e\t",markers[m].texture.M[i][j]);
	  }
	  printf("\n");
	}
	printf("N tensor (2d):\n");
	for(i=0;i<2;i++){
	  for(j=0;j<2;j++){
	    printf("%e\t",markers[m].texture.N[i][j]);
	  }
	  printf("\n");
	}
      }
    }
    //#endif


    if( !(iTime % 20) || Eii > Emax ){
      /* dump the fabric to a file */
      char filename[80];
      sprintf(filename,"output/texture.%d.dat",iTime);
      FILE *fh;
      fh = fopen(filename,"w");
      Marker M = markerset.markers[0];
      /* print stress second invariant, strain rate second invariant */
      fprintf(fh,"NT=%d\n",NT);
      fprintf(fh,"eii=%e\n",eii);
      fprintf(fh,"Eii=%e\n",Eii);
      fprintf(fh,"eta=%e\n",sii/2.0/eii);
      PetscInt i,j,k;
      for(i=0;i<NT;i++){
	for(j=0;j<NT;j++){
	  for(k=0;k<NT;k++){
	    fprintf(fh,"%e\t%e\n",M.texture.ctheta[i][j][k],M.texture.cphi[i][j][k]);
	  }
	}
      }
      fclose(fh);

    }
    if( Eii > Emax )break;
  }/* end timestep loop */

  /* calculate effective deviatoric strain-rate*/
  //  Tensor33s *M = &(markers[m].texture.M);
  PetscScalar sxx1 = markers[m].s.T11;
  PetscScalar syy1 = markers[m].s.T22;
  PetscScalar sxy1 = markers[m].s.T12;

  /* compute sii*/

  printf("sii = %e Pa, effective viscosity = sii/dii = \n",sii);


  /* change orientation of stresses and deform again */
  for(m=0;m<markerset.nMark;m++){
    markers[m].s.T11 = 0.0;
    markers[m].s.T12 = 1e5;
    markers[m].s.T22 = 0.0;
    markers[m].s.T23 = 0.0;
    markers[m].s.T13 = 0.0;
    markers[m].s.T33 = 0.0;
    markers[m].T = 200.0;
    markers[m].wxy = 0.0;
    markers[m].wxz = 0.0;
    markers[m].wyz = 0.0;
    markers[m].Mat = 0;
    markers[m].cellX = 1;
    markers[m].cellY = 1;
    markers[m].texture.d = options.grainSize;
  }
  Tensor33s D;/* deviatoric strain rate tensor */    
  for(m=0;m<markerset.nMark;m++){
    //PetscErrorCode formViscoplasticM(Marker *,Options *,Tensor33s *, PetscInt );
    Marker *M = &(markerset.markers[m]);
    formViscoplasticMMarker( M, &options, &D, 1 );
    invertMtoN( &( M->texture ) );
  }
  PetscScalar eii = sqrt(0.5*( D.T11*D.T11 + D.T22*D.T22 + D.T33*D.T33 ) + D.T12*D.T12 + D.T23*D.T23 + D.T13*D.T13);
  printf("Simple shear step\n");
  for(m=0;m<markerset.nMark;m++){
    printf("Marker %d:\n",m);
    printf("M tensor:\n");
    PetscInt i,j;
    for(i=0;i<6;i++){
      for(j=0;j<6;j++){
	printf("%e\t",markers[m].texture.M[i][j]);
      }
      printf("\n");
    }
    printf("N tensor (2d):\n");
    for(i=0;i<2;i++){
      for(j=0;j<2;j++){
	printf("%e\t",markers[m].texture.N[i][j]);
      }
      printf("\n");
    }
  }

  ierr = PetscFinalize();CHKERRQ(ierr);
  
}

PetscScalar drand(){
  return ((PetscScalar) rand()) / ((PetscScalar) RAND_MAX);
}
