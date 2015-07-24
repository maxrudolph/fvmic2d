
#include<stdio.h>
#include<stdlib.h>
#include "petscksp.h"
#include "petsctime.h"
#include "fdcode.h"
#include "initialize_problem.h"
#include "version.h"
#include "pressureNullSpace.h"

PetscErrorCode initialize_problem(Problem *problem){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;

  problem->mech_system.pc = PETSC_NULL; /* nullify preconditioner */
  /* random number generator context*/
  ierr = PetscRandomCreate(PETSC_COMM_WORLD, &(problem->random) );CHKERRQ(ierr);
  ierr = PetscRandomSetType(problem->random,PETSCRAND48);CHKERRQ(ierr);

  MPI_Comm_size(PETSC_COMM_WORLD,&problem->parallel.size);  //Get MPI rank, size
  MPI_Comm_rank(PETSC_COMM_WORLD,&problem->parallel.rank);

  {
    /* seed random number generator with cpu time */
    PetscLogDouble time;
    unsigned long seed;
    if(!problem->parallel.rank) ierr = PetscTime(&time); CHKERRQ(ierr);
    ierr = MPI_Bcast( &time, sizeof(PetscLogDouble),MPI_BYTE , 0, PETSC_COMM_WORLD);
    seed = (unsigned long) time;
    seed = 1;
    printf("seeding random number generator with %ld\n",seed);fflush(stdout);
    ierr = PetscRandomSetSeed(problem->random,seed);CHKERRQ(ierr);
    ierr = PetscRandomSeed(problem->random);CHKERRQ(ierr);/* seed the generator*/
  }
  if(!problem->parallel.rank) printf("Version Information : %s\n",MARKERCODE_HG_VERSION);

  PetscFunctionReturn(ierr);
}

PetscErrorCode initialize_matrices(Problem *problem){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  /* set up local to global mapping*/
  ISLocalToGlobalMapping ltogm;
  ierr=DMGetLocalToGlobalMapping(problem->grid.da,&ltogm);CHKERRQ(ierr);
  
  /* initialize LHS, matrix and RHS, matrix for mechanical problem */
  ierr = DMCreateMatrix(problem->grid.vda,&problem->mech_system.lhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem->mech_system.lhs,"CmechLHS");
  ierr = createPressureNullSpace( &problem->grid, &problem->mech_system.ns );CHKERRQ(ierr); 

  problem->Kbond = 0.0;
  problem->Kcont = 0.0;
  {
    PetscInt m,n; 
    ierr=MatGetSize(problem->mech_system.lhs,&m,&n);CHKERRQ(ierr);
    PetscInt ndof=problem->grid.NX * problem->grid.NY * 3;
    if(!problem->parallel.rank) printf("LHS is %d by %d, ndof=%d\n",m,n,ndof); 
    ierr=MatGetOwnershipRange(problem->mech_system.lhs,&m,&n);CHKERRQ(ierr);
    printf("[%d] has rows %d-%d\n",problem->parallel.rank,m,n);CHKERRQ(ierr);
  }

  ierr=DMGetGlobalVector(problem->grid.vda,&problem->mech_system.rhs);CHKERRQ(ierr);
  ierr = PetscObjectSetName( (PetscObject) problem->mech_system.rhs,"CmechRHS"); 

  /* initialize LHS matrix, RHS Vec and nodalHeating Vec for thermal problem*/
  ierr=DMGetGlobalVector(problem->grid.da,&problem->thermal_system.rhs);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem->thermal_system.rhs,"CthermalRHS");
  ierr=DMCreateMatrix(problem->grid.da,&problem->thermal_system.lhs);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem->thermal_system.lhs,"CthermalLHS");
  ierr=VecDuplicate(problem->thermal_system.rhs,&problem->nodal_heating);CHKERRQ(ierr);
  ierr=PetscObjectSetName( (PetscObject) problem->nodal_heating,"nodalHeating");
  {
    PetscBool dop = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-stokes_ksp_build_pc_diag",&dop,0);
    if( dop ){
      printf("[%d] Allocating Preconditioner\n",problem->parallel.rank);
      /*       ierr=MatDuplicate( LHS, MAT_DO_NOT_COPY_VALUES ,&P ); CHKERRQ(ierr); */
      ierr = DMCreateMatrix(problem->grid.vda, &problem->mech_system.pc); CHKERRQ(ierr);
      ierr = PetscObjectSetName( (PetscObject) problem->mech_system.pc, "MechanicalPreconditioner");
    }
  }

  ierr = VecDuplicate(problem->mech_system.rhs, &problem->mech_system.solution);  
  ierr = VecZeroEntries(problem->mech_system.solution);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);

}
