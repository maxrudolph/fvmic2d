#include<stdio.h>
#include<stdlib.h>
#include "petscksp.h"

static char help[] = "Help statement goes here\n";

int main(int argc, char **args){
  PetscMPIInt rank,size;

  PetscInitialize(&argc,&args,(char *)0,help);  //INITIALIZE PETSC

  MPI_Comm_size(PETSC_COMM_WORLD,&size);  //Get MPI rank
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscRandom r;
  PetscErrorCode ierr;

  ierr=  PetscRandomCreate(PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  ierr = PetscRandomSetType(r,PETSCRAND48);CHKERRQ(ierr);

  {
    /* seed random number generator with cpu time */
    PetscLogDouble time;
    unsigned long seed;
    if(!rank) ierr = PetscGetTime(&time); CHKERRQ(ierr);
    ierr = MPI_Bcast( &time, sizeof(PetscLogDouble),MPI_BYTE , 0, PETSC_COMM_WORLD);
    seed = (unsigned long) time;
    //    seed = 1;
    printf("seeding random number generator with %ld\n",seed);fflush(stdout);
    ierr = PetscRandomSetSeed(r,seed);CHKERRQ(ierr);
    ierr = PetscRandomSeed(r);CHKERRQ(ierr);/* seed the generator*/
  }
  
  
  ierr = PetscRandomDestroy( r );CHKERRQ(ierr);

  ierr = PetscFinalize();
}
