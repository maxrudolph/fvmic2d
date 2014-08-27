/* This subroutine is designed to read in Markers from the isotropic code and create a texture file for each CPU */
#include "fdcode.h"
#include "math.h"

static char help[] = "To run, type progname infile outfile\n";

int main(int argc, char **args){
  PetscInitialize(&argc,&args,PETSC_NULL,help);
  PetscErrorCode ierr;
  /* first argument should be Marker file */
  /* load markers - these should be markers as defined in isotropic code */
    char *ifn = args[1];
    printf("reading from %s\n",ifn);
    char *ofn = args[2];
    
    FILE *ifile, *ofile;
    ifile = fopen( ifn, "rb"); /* Open input and output files */

    PetscInt nMark;
    PetscScalar elapsedTime;
    fread( &nMark, sizeof(PetscInt),1,ifile);/* number of markers*/
    fread( &elapsedTime, sizeof(PetscScalar),1,ifile); /*elapsed time */

    printf("reading %d markers, elapsedTime = %e\n",nMark,elapsedTime);
    
    //fseek( ifile, 3*sizeof(PetscInt)*nMark, SEEK_CUR);/* cellx, celly*/
    
    PetscInt *cpu;
    ierr = PetscMalloc( nMark*sizeof(PetscInt), &cpu);
    fread( cpu, sizeof(PetscInt), nMark, ifile);
    /* find maximum cpu number */
    PetscInt m;
    PetscInt cpumax = 0;
    for(m=0;m<nMark;m++){
      if( cpu[m] > cpumax) cpumax = cpu[m];
    }
    PetscInt ncpu = cpumax+1;
    printf("found %d cpus\n",ncpu);
    
    /* now count number of markers on each cpu */
    PetscInt *nMarks;
    ierr = PetscMalloc( ncpu*sizeof(PetscInt), &nMarks);
    for( m=0;m<ncpu;m++) nMarks[m] = 0;
    for( m=0 ; m<nMark ; m++){
      nMarks[ cpu[m] ]++;
    }
    for( m = 0; m < ncpu; m++){
      printf("cpu %d : %d markers\n", m ,nMarks[m]); 
    }
    /* Look up the all-important NT parameter and print it out to be safe */
    printf("!!! NT = %d !!!\n",NT);
    PetscInt NT3 = NT*NT*NT;

    /* initialize texture arrays */
    PetscInt icpu;
    for( icpu=0;icpu<ncpu;icpu++){
      /* allocate memory for ctheta, cphi arrays */
      PetscScalar *ctheta,*cphi;
      ierr=PetscMalloc( NT3*nMarks[icpu]*sizeof(PetscScalar), &ctheta);CHKERRQ(ierr);
      ierr=PetscMalloc( NT3*nMarks[icpu]*sizeof(PetscScalar), &cphi); CHKERRQ(ierr);
      for( m=0; m<nMarks[icpu]; m++){/* loop over markers for this cpu */
	PetscInt ic;
	for( ic=0 ; ic< NT3; ic++){
	  /* make isotropic fabric */
	  ctheta[NT3*m+ic] =  2.0*M_PI*drand48();
	  PetscScalar z = -1.0 + 2.0*drand48();
	  cphi[NT3*m+ic] = acos(z);
	}
      }     
      /* write texture file */
      char ofn1[80];
      sprintf(ofn1,"%s.%d",ofn,icpu);
      ofile = fopen( ofn1, "wb");
      printf("writing to %s\n",ofn1);
      fwrite( &ncpu, sizeof(PetscInt), 1, ofile);
      PetscInt nt = NT;

      fwrite( &nMarks[icpu], sizeof(PetscInt), 1, ofile);
      fwrite( &nt, sizeof(PetscInt),1,ofile);
      for( m=0;m<nMarks[icpu];m++){
	fwrite( &ctheta[m*NT3],sizeof(PetscScalar),NT3,ofile);
	fwrite( &cphi[m*NT3],sizeof(PetscScalar),NT3,ofile);
      }
      fclose(ofile);
      printf("made it here\n");fflush(stdout);
      /* free texture arrays */
      ierr=PetscFree( ctheta );CHKERRQ(ierr);
      ierr=PetscFree( cphi );CHKERRQ(ierr);
    }

    fclose( ifile );
    
    ierr = PetscFree( cpu);
    ierr = PetscFree( nMarks );
        
    /* create texture for these markers */
    ierr = PetscFinalize();CHKERRQ(ierr);

}
