#include "fdcode.h"
#include "memuse.h"

void getMemoryUse(){
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  /* get memory use */
  PetscLogDouble memuse;
  PetscScalar memuses;
  PetscScalar memuset;
  PetscMemoryGetCurrentUsage( &memuse );
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"[%d] memory use %e \n",rank,memuse);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  memuses = (PetscScalar) memuse;
  MPI_Allreduce( &memuses, &memuset, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  if(!rank) PetscPrintf(PETSC_COMM_SELF,"Global memory use %e\n",memuset); 
}  
