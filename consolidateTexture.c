#include "fdcode.h"

#ifdef TEXTURE
/* read the c-axis orientations (theta,phi) */
int main(int argc, char ** argv){
  PetscInt rank=0;
  PetscInt size;
  PetscErrorCode ierr=0;
  PetscInt m;

  if( argc != 2){
    printf("error: wrong number of arguments\n");
    abort();
  }

  char *ifn = argv[1];

  PetscInt iMonte;
  PetscInt iTime;

  char fnt[160];

  FILE *ofh;

  PetscInt irank=0;
  size = 1; /* initial guess. */

  PetscInt nt = NT;

  ofh = fopen(ifn, "wb");
  PetscInt nTotal = 0;
  while(irank < size){
    FILE *fh;
    sprintf(fnt,"%s.%d",ifn,irank);
    printf("Reading Texture from %s\n",fnt);
  
    fh = fopen(fnt,"r");
      
    fread(&size,sizeof(PetscInt),1,fh);
    printf("size=%d\n",size);
    PetscInt nMarkTex = 0; /* number of markers in this file */
    fread(&nMarkTex,sizeof(PetscInt),1,fh);    

    fread(&nt,sizeof(PetscInt),1,fh);
    
    if( nt != NT ){ printf("Error: File NT=%d is not compatible with compiled-in setting %d\n",nt, NT); abort(); }
    
    PetscInt NT3 = NT*NT*NT;
    
    PetscScalar *ctheta;
    PetscScalar *cphi;
    ierr = PetscMalloc(nMarkTex*NT3*sizeof(PetscScalar), &ctheta );CHKERRQ(ierr);
    ierr = PetscMalloc(nMarkTex*NT3*sizeof(PetscScalar), &cphi   );CHKERRQ(ierr);
    
    for(m=0;m<nMarkTex;m++){
      fread(&(ctheta[NT3*m]),sizeof(PetscScalar),NT3,fh);
      fread(&(cphi[NT3*m]),sizeof(PetscScalar),NT3,fh);
    }
    
    /* write to the output file */
    if( irank == 0 ){
      fwrite(&size,sizeof(PetscInt),1,ofh);
      fwrite(&nMarkTex,sizeof(PetscInt),1,ofh);
      fwrite(&nt,sizeof(PetscInt),1,ofh);      
    }
    for(m=0;m<nMarkTex;m++){
      fwrite(&(ctheta[NT3*m]),sizeof(PetscScalar),NT3,ofh);
      fwrite(&(cphi[NT3*m]),sizeof(PetscScalar),NT3,ofh);
    }
    ierr = PetscFree(ctheta); CHKERRQ(ierr);
    ierr = PetscFree(cphi); CHKERRQ(ierr);
    fclose(fh);
    nTotal += nMarkTex;
    irank++;
  }
  
  fseek( ofh, sizeof(PetscInt), SEEK_SET );
  fwrite(&nTotal,sizeof(PetscInt),1,ofh);
  printf("nTotal=%d\n",nTotal);
  
  fclose(ofh);


  
  return 0;
}
#endif
