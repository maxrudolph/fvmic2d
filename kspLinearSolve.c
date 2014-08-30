#include<stdlib.h>
#include<stdio.h>

#include "petscksp.h"
#include "mpi.h"
#include "kspLinearSolve.h"
#include "profile.h"
#include "fdcode.h"

//#define DEBUG
PetscErrorCode KSPMonitorStokesBlocks(KSP ,PetscInt ,PetscReal,void *);

PetscErrorCode kspLinearSolve( Mat LHS, Mat P, Vec RHS, Vec X,const char *label, MatNullSpace ns ){
  PetscErrorCode ierr;
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  setLogStage( LOG_SOLVE );
  /* This code creates a linear solver context and solves the linear system. In order to allow thermal and mechanical problems to share this solve routine while specifying different solvers on the command line, the char *label needs to be passed in and also prefixed to any ksp options on the command line. */
      KSP ksp;
      PC prec;
      Vec zerovec;
      KSPConvergedReason reason;
#ifdef DEBUG
      PetscViewer viewer;
#endif

      PetscFunctionBegin;

      ierr = VecDuplicate(RHS, &zerovec);CHKERRQ(ierr);
      ierr = VecZeroEntries(zerovec);CHKERRQ(ierr);
      ierr = MatDiagonalSet(LHS,zerovec,ADD_VALUES);CHKERRQ(ierr);
      
      ierr = VecDestroy(&zerovec);CHKERRQ(ierr);
      ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr); 
      ierr = KSPSetOptionsPrefix(ksp,label); CHKERRQ(ierr);
      if( P == PETSC_NULL ){
	ierr = KSPSetOperators(ksp,LHS,LHS);CHKERRQ(ierr); 
      }else{
	ierr = KSPSetOperators(ksp,LHS,P);CHKERRQ(ierr); 
      }
      ierr = KSPGetPC(ksp,&prec);CHKERRQ(ierr); 
      ierr = PCSetFromOptions(prec); CHKERRQ(ierr); 
      ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr); 
      const char prefix[] = "stokes_";
      if(!strncmp( prefix, label, 6)){
	
	PetscBool fs = PETSC_FALSE;
	ierr = PetscObjectTypeCompare((PetscObject) prec, PCFIELDSPLIT, &fs);CHKERRQ(ierr);
	if( fs ){
	  const PetscInt ufields[] = {DOF_U,DOF_V}, pfields[]={DOF_P};
	  ierr = PCFieldSplitSetBlockSize(prec,3); CHKERRQ(ierr);
	  ierr = PCFieldSplitSetFields(prec,"u",2,ufields,ufields);CHKERRQ(ierr);
	  ierr = PCFieldSplitSetFields(prec,"p",1,pfields,pfields);CHKERRQ(ierr);
	}

	PetscBool stokes_monitor = PETSC_FALSE;
	PetscOptionsGetBool(PETSC_NULL,"-stokes_ksp_monitor_blocks",&stokes_monitor,0);
	if( stokes_monitor){
	  KSPMonitorSet(ksp,KSPMonitorStokesBlocks,PETSC_NULL,PETSC_NULL);
	}
	if( !FIX_PRESSURE && (ns != PETSC_NULL)){
	  if(!rank) printf("Setting null space \n");
	  ierr = KSPSetNullSpace( ksp, ns ); CHKERRQ(ierr);
	}
      }
      
      ierr = KSPSolve( ksp,  RHS,  X);//CHKERRQ(ierr);

#ifdef DEBUG
      char fn[80];
      sprintf(fn,"./output/%slinearsystem.petscbin",label);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fn,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);  
      ierr = MatView( LHS, viewer);CHKERRQ(ierr);
      if( P != PETSC_NULL){
	ierr = MatView( P, viewer);CHKERRQ(ierr);	
      }
      ierr = VecView( X  ,viewer);CHKERRQ(ierr);    
      ierr = VecView( RHS  ,viewer);CHKERRQ(ierr);    

      Vec v,w,Br;
      MatGetVecs( LHS, &w,&v );
      KSPBuildResidual( ksp, v, w, &Br );
      ierr = VecView( Br  ,viewer);CHKERRQ(ierr);    
      VecDestroy(&v);
      VecDestroy(&w);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
#endif
      ierr = KSPGetConvergedReason( ksp, &reason);
      if(reason < 1){
	printf("converged reason %d\n",reason); 
	printf("failure to converge!\n");
	abort();
      }
      ierr = KSPDestroy( &ksp );CHKERRQ(ierr);
      PetscLogStagePop();
      PetscFunctionReturn(ierr);
}

/* KSPMonitorStokesBlocks is lifted from petsc example ex42.c by Dave May */
PetscErrorCode KSPMonitorStokesBlocks(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy)
{
  const PetscInt bs = 3;/* block size */
  PetscReal norms[bs];
  Vec Br, v,w;
  Mat A;

  KSPGetOperators( ksp, &A, PETSC_NULL);
  MatGetVecs( A, &w,&v );
  VecSetBlockSize(v,bs);
  
  KSPBuildResidual( ksp, v, w, &Br );

  VecStrideNorm(Br,0,NORM_2,&norms[0]);
  VecStrideNorm(Br,1,NORM_2,&norms[1]);
  VecStrideNorm(Br,2,NORM_2,&norms[2]);
   
  VecDestroy(&v);
  VecDestroy(&w);
  
  PetscPrintf(PETSC_COMM_WORLD,"  %d KSP Component U,V,P residual norm [ %1.12e, %1.12e, %1.12e ]\n",n,norms[0],norms[1],norms[2]);
  
  return(0);
  
}
