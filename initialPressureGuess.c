/* This file contains subroutines that will generate an initial guess at the pressure field. The way that they work is by imposing zero pressure at the top left corner of the domain and then using the solver specified for the diffusion equation to obtain an initial pressure field, which is then copied into the correct spot in the mechanical solution vector */
#include "fdcode.h"
#include "initialPressureGuess.h"
#include "kspLinearSolve.h"

PetscErrorCode initialPressureGuess( GridData *grid , NodalFields *nodalFields, Options *options, Vec S ){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  const PetscScalar gy = options->gy;/* acceleration due to gravity in x, y directions */
  const PetscScalar gx = 0.0;
  const PetscInt NX = grid->NX;
  const PetscInt NY = grid->NY;

  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  /* create LHS matrix for pressure poisson equation, solve it for the correct pressure in this cell */
  Mat LHS;
  Vec RHS;
  ierr= DMCreateMatrix(grid->da,&LHS);CHKERRQ(ierr); /* for petsc-dev */
  ierr= MatZeroEntries(LHS);CHKERRQ(ierr);
  ierr = DMGetGlobalVector(grid->da,&RHS);CHKERRQ(ierr);
  ierr=VecZeroEntries(RHS); CHKERRQ(ierr);
  /* get density arrays */
  Vec rhol;
  PetscScalar **rho;
  ierr = DMCreateLocalVector(grid->da,&rhol);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  
  ierr = DMDAVecGetArray(grid->da,rhol, &rho);CHKERRQ(ierr);

  PetscInt ix,jy;
  PetscInt ixl, jyl;

  /* Set up DOF indexing */
  PetscInt nlocal;
  ISLocalToGlobalMapping ltogm;
  PetscInt *globalIdx;
  ierr=DMGetLocalToGlobalMapping(grid->da,&ltogm); CHKERRQ(ierr);
  ierr=ISLocalToGlobalMappingGetIndices(ltogm,&globalIdx);CHKERRQ(ierr);

  PetscInt xg,yg,mg,ng;
  ierr=DMDAGetGhostCorners(grid->da,&xg,&yg,PETSC_NULL,&mg,&ng,PETSC_NULL);CHKERRQ(ierr);
  PetscInt **pdof;
  PetscMalloc( ng*sizeof(PetscInt *), &pdof);
  for(jy=0;jy<ng;jy++){
    ierr=PetscMalloc( mg*sizeof(PetscInt), &pdof[jy]);CHKERRQ(ierr);
  }
  for(jy=0;jy<ng;jy++){
    for(ix=0;ix<mg;ix++){
      pdof[jy][ix] = globalIdx[ix+mg*jy];/* vz goes into a separate matrix*/
    }
  }
  ierr = ISLocalToGlobalMappingRestoreIndices(ltogm,&globalIdx);CHKERRQ(ierr);
  PetscScalar dx = grid->x[1]-grid->x[0];
  PetscScalar dy = grid->y[1]-grid->y[0];
  PetscScalar Kb = 4.0/(dx*dx+dy*dy); /* scaling to condition LHS matrix properly */

  for(ix=x;ix<x+m;ix++){
    ixl = ix-xg;
    for(jy=y;jy<y+n;jy++){
      jyl = jy-yg;
      /* div(grad p) = grad(rho).g = d(rho)/dy * gy + d(rho)/dx * gx */
      /* i+1,j term */
      if( (!grid->xperiodic && ix == 0) || jy == 0){
	PetscScalar dx = grid->x[1]-grid->x[0];
	PetscScalar dy = grid->y[1]-grid->y[0];
	PetscScalar lval = 4.0/(dx*dx+dy*dy); /* approximate scaling to condition LHS matrix properly */
	
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],lval,ADD_VALUES); CHKERRQ(ierr);

      }else if( ix == 3 && jy == 3){/* specify zero pressure in one cell */
	PetscScalar dx = grid->x[ix]-grid->x[ix-1];
	PetscScalar dy = grid->y[jy]-grid->y[jy-1];
	PetscScalar lval = 4.0/(dx*dx+dy*dy); /* approximate scaling to condition LHS matrix properly */
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],lval,ADD_VALUES); CHKERRQ(ierr);

      }else{
	PetscInt colidx[5];
	PetscScalar vals[5];
	PetscInt rowidx = pdof[jyl][ixl];
	/* right */
	colidx[0] = (!grid->xperiodic && ix==NX-1) ? -1 : pdof[jyl][ixl+1];
	vals[0] = 1.0/(grid->x[ix]-grid->x[ix-1])/(grid->xc[ix+1]-grid->xc[ix]);
	
	/* center */
	colidx[1] = pdof[jyl][ixl];
	vals[1] =  1.0/(grid->x[ix]-grid->x[ix-1])*( -1.0/(grid->xc[ix+1]-grid->xc[ix]) -1.0/(grid->xc[ix]-grid->xc[ix-1])) \
	  +  1.0/(grid->y[jy]-grid->y[jy-1])*( -1.0/(grid->yc[jy+1]-grid->yc[jy]) -1.0/(grid->yc[jy]-grid->yc[jy-1]));

	/* left */
	colidx[2] = (!grid->xperiodic && ix==1) ? -1 : pdof[jyl][ixl-1];
	vals[2] =  1.0/(grid->x[ix]-grid->x[ix-1])/(grid->xc[ix]-grid->xc[ix-1]);

	/* down */
	colidx[3] = jy == NY-1 ? -1 : pdof[jyl+1][ixl];
	vals[3] = 1.0/(grid->y[jy]-grid->y[jy-1])/(grid->yc[jy+1]-grid->yc[jy]);
	
	/* up */
	colidx[4] = jy == 1 ? -1 : pdof[jyl-1][ixl];
     	vals[4] = 1.0/(grid->y[jy]-grid->y[jy-1])/(grid->yc[jy]-grid->yc[jy-1]);

	ierr = MatSetValues( LHS, 1,&rowidx,5,colidx,vals,ADD_VALUES);CHKERRQ(ierr);	
	PetscScalar rval = gy * ( (rho[jy][ix] + rho[jy][ix-1])/2.0-(rho[jy-1][ix]+rho[jy-1][ix-1])/2.0 )/(grid->y[jy]-grid->y[jy-1]) \
	  + gx * ((rho[jy][ix] + rho[jy-1][ix])/2.0-(rho[jy][ix-1]+rho[jy-1][ix-1])/2.0 )/(grid->x[ix]-grid->x[ix-1]);
	ierr = VecSetValue( RHS, pdof[jyl][ixl], rval, ADD_VALUES); CHKERRQ(ierr);
      }
      
      /* deal with special boundary condition cases */
      if(!grid->xperiodic ){
	if( ix == 1 ){
	  PetscScalar p_im_j_coeff =  1.0/(grid->x[ix]-grid->x[ix-1])/(grid->xc[ix]-grid->xc[ix-1]);
	  ierr = MatSetValue( LHS,pdof[jyl][ixl], pdof[jyl][ixl], p_im_j_coeff, ADD_VALUES);
	}else if(ix == NX-1){
	  PetscScalar p_ip_j_coeff =  1.0/(grid->x[ix]-grid->x[ix-1])/(grid->xc[ix+1]-grid->xc[ix]);
	  ierr = MatSetValue( LHS,pdof[jyl][ixl], pdof[jyl][ixl], p_ip_j_coeff, ADD_VALUES);
	}
      }
      if(jy == 1 ){/* force dpdy = rho*g */
	PetscScalar p_i_jm_coeff = 1.0/(grid->y[jy]-grid->y[jy-1])/(grid->yc[jy]-grid->yc[jy-1]);
	PetscScalar rhoc = (rho[jy][ix] + rho[jy][ix-1] + rho[jy-1][ix] + rho[jy-1][ix-1])/4.0;
	PetscScalar dy = grid->yc[jy]-grid->yc[jy-1];
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],p_i_jm_coeff, ADD_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl], p_i_jm_coeff*rhoc*gy*dy ,ADD_VALUES);CHKERRQ(ierr);
      }else if(jy == NY-1){
	PetscScalar p_i_jp_coeff =  1.0/(grid->y[jy]-grid->y[jy-1])/(grid->yc[jy+1]-grid->yc[jy]);
	PetscScalar rhoc = (rho[jy][ix] + rho[jy][ix-1] + rho[jy-1][ix] + rho[jy-1][ix-1])/4.0;
	PetscScalar dy =  grid->yc[jy]-grid->yc[jy-1];
	ierr = MatSetValue( LHS, pdof[jyl][ixl],pdof[jyl][ixl],p_i_jp_coeff,ADD_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue( RHS, pdof[jyl][ixl], -p_i_jp_coeff*rhoc*gy*dy, ADD_VALUES);CHKERRQ(ierr);
	
      }

    }
  }






  /* assemble matrix and vector */
  ierr = VecAssemblyBegin(RHS);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(RHS);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  Vec P;/* solution vector */
  ierr = VecDuplicate( RHS, &P); CHKERRQ(ierr);
  ierr = VecZeroEntries(P); CHKERRQ(ierr);
  
  /* solve the linear system */
  ierr = kspLinearSolve(LHS, PETSC_NULL , RHS, P,"energy_",PETSC_NULL);CHKERRQ(ierr);/* note that I use the energy label here because we are solving a poisson equation for pressure, the same form as the energy equation so whatever solver was chosen for energy will be fine for this as well */
  PetscViewer viewer;
  ierr=  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"output/pressureGuess.petscbin",FILE_MODE_WRITE,&viewer);
  ierr= MatView( LHS, viewer);CHKERRQ(ierr);
  ierr=  VecView( RHS,viewer );CHKERRQ(ierr);
  ierr=  VecView( P,viewer );CHKERRQ(ierr);
  ierr=  PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  
  /* get the dof corresponding to the pressure in the mechanical solution vector */
  PetscScalar ***Sa;/* solution array */
  PetscScalar **p;/* pressure solution array */

  ierr = DMDAVecGetArray(grid->da,P,&p); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(grid->vda,S,&Sa);CHKERRQ(ierr);

  for( jy = y; jy<y+n;jy++){
    for(ix = x; ix<x+m;ix++){
      Sa[jy][ix][DOF_P] = p[jy][ix];
    }
  }

  ierr = DMDAVecRestoreArrayDOF(grid->vda,S,&Sa);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,P,&p); CHKERRQ(ierr);
  ierr = VecDestroy( &P);CHKERRQ(ierr);
  
  ierr=DMDAVecRestoreArray(grid->da,rhol, &rho);CHKERRQ(ierr);
  ierr = VecDestroy(&rhol); CHKERRQ(ierr);
  ierr = MatDestroy(&LHS);CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(grid->da,&RHS);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
