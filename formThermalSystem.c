#include "fdcode.h"
#include "thermalSystem.h"
#include "initialConditions.h"
#include "profile.h"

PetscErrorCode formThermalSystem(Problem *problem, GridData *grid, NodalFields *nodalFields, Mat thermalLHS, Vec thermalRHS, PetscScalar dt, Options *options, PetscInt doSteady){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  setLogStage( LOG_THERM_ASM );
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscInt NX = grid-> NX;
  PetscInt NY = grid-> NY;
  /* Shorthand form of grid->x, y*/
  PetscScalar *gx = grid->x;
  PetscScalar *gy = grid->y;


  /* Find local ownership range */
  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  PetscInt xg,yg,mg,ng;/* get the ghost corners too */
  ierr=DMDAGetGhostCorners(grid->da,&xg,&yg,PETSC_NULL,&mg,&ng,PETSC_NULL);CHKERRQ(ierr);
  /* get the global indices of the local nodes*/
  ISLocalToGlobalMapping ltogm;
  const PetscInt *globalIdx;
  ierr=DMGetLocalToGlobalMapping(grid->da,&ltogm); CHKERRQ(ierr);
  ierr=ISLocalToGlobalMappingGetIndices(ltogm,&globalIdx);CHKERRQ(ierr);

  /* Get vectors with ghost information for kThermal, rho, Cp */
  /* create local vectors containing ghost information */
  Vec kThermall, rhol, Cpl, lastTl, vxl, vyl;
  ierr = DMCreateLocalVector( grid->da, &kThermall); CHKERRQ(ierr);
  ierr = VecDuplicate( kThermall, &rhol); CHKERRQ(ierr);
  ierr = VecDuplicate( kThermall, &Cpl); CHKERRQ(ierr);
  ierr = VecDuplicate( kThermall, &lastTl); CHKERRQ(ierr);
  ierr = VecDuplicate( kThermall, &vxl); CHKERRQ(ierr);
  ierr = VecDuplicate( kThermall, &vyl); CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->kThermal, INSERT_VALUES, kThermall);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->kThermal, INSERT_VALUES, kThermall);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->rho, INSERT_VALUES, rhol);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->rho, INSERT_VALUES, rhol);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->Cp, INSERT_VALUES, Cpl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->Cp, INSERT_VALUES, Cpl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->lastT, INSERT_VALUES, lastTl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->lastT, INSERT_VALUES, lastTl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->vx, INSERT_VALUES, vxl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->vx, INSERT_VALUES, vxl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin( grid->da, nodalFields->vy, INSERT_VALUES, vyl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd( grid->da, nodalFields->vy, INSERT_VALUES, vyl);CHKERRQ(ierr);
  /* get local arrays for matrix assembly */
  PetscScalar **kThermal, **rho, **Cp, **lastT, **vx, **vy;
  ierr = DMDAVecGetArray(grid->da, kThermall, &kThermal); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da, rhol, &rho); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da, Cpl, &Cp); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da, lastTl, &lastT); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da, vxl, &vx); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da, vyl, &vy); CHKERRQ(ierr);
  
  /* make a DOF mapping */
  PetscInt **Tdof;
  PetscInt ix, jy; /* loop indices */
  PetscMalloc( ng*sizeof(PetscInt *), &Tdof);
  for(jy=0;jy<ng;jy++){
    ierr=PetscMalloc( mg*sizeof(PetscInt), &Tdof[jy]);CHKERRQ(ierr);
  }

  for(jy=0;jy<ng;jy++){
    for(ix=0;ix<mg;ix++){
      Tdof[jy][ix] = globalIdx[ix+mg*jy];
    }
  }
  ierr=ISLocalToGlobalMappingRestoreIndices(ltogm,&globalIdx);CHKERRQ(ierr);
  /* assemble the matrix */

  /*   PetscInt y1 = y; */ /* these are the DA corners from above */
  /*   PetscInt x1 = x; */ /* while currently unused (hence commented), these would be useful for building an array containing global subscripts (DOFs) in local indexing*/
  
  for( jy=y; jy<y+n;jy++){
    PetscInt jyl = jy-yg; 
    for(ix=x;ix<x+m;ix++){
      PetscInt ixl = ix-xg; 
      PetscInt idxnode = Tdof[jyl][ixl];
      /* boundary conditions*/
      if( !grid->xperiodic && ix == 0 ){/* left wall */
	/* fixed temperature */
	ierr = MatSetValue( thermalLHS, idxnode, idxnode, 1.0, INSERT_VALUES);CHKERRQ(ierr);
       
      } else if(ix == NX-1 && (grid->y[jy] <= plate_depth(grid->LX,options) || (jy>0 && (vx[jy][ix]+vx[jy-1][ix])/2.0 < 0.0  )) ){/* RIGHT WALL */
	// enforce overriding plate geotherm
	ierr = MatSetValue( thermalLHS, idxnode, idxnode, 1.0, INSERT_VALUES);CHKERRQ(ierr);
      }else if(ix == NX-1){
	ierr = MatSetValue( thermalLHS, idxnode, idxnode, 1.0, INSERT_VALUES);CHKERRQ(ierr);

      } else if( jy ==0){
	/* Top boundary */
	if( options->thermalBCTop.type[0] == 0){
	  /* fixed temperature */
	  ierr = MatSetValue( thermalLHS, idxnode, idxnode, 1.0, INSERT_VALUES);CHKERRQ(ierr);	  
	}else if(options->thermalBCTop.type[0] == 1){
	  /* fixed temperature gradient */
	  PetscScalar dyt1 = grid->y[jy+1]-grid->y[jy];
	  ierr = MatSetValue( thermalLHS, idxnode, idxnode, -1.0/dyt1, INSERT_VALUES);CHKERRQ(ierr);
	  ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl+1][ixl], 1.0/dyt1, INSERT_VALUES);CHKERRQ(ierr);
	}
      } else if( jy == NY-1){/* Bottom */
	/* Bottom boundary */
	ierr = MatSetValue( thermalLHS, idxnode, idxnode, 1.0, INSERT_VALUES);CHKERRQ(ierr);
      } else{ /* interior */
	if(ix == NX-1) printf("ix, jy %d %d\n",ix,jy);
	PetscScalar ka,kb,kc,kd;
	
	//       ix-1 dx1 ix dx2 ix+1
	// jy-1  x        x      x
	// dy1            c 
	// jy    x   a    x   b  x
	// dy2            d 
	// jy+1  x        x      x
	//

	ka = (kThermal[jy][ix] + kThermal[jy][ix-1])/2.0;
	kb = (kThermal[jy][ix] + kThermal[jy][ix+1])/2.0;
	kc = (kThermal[jy][ix] + kThermal[jy-1][ix])/2.0;
	kd = (kThermal[jy][ix] + kThermal[jy+1][ix])/2.0;

	PetscScalar dx2 = gx[ix+1]-gx[ix];
	PetscScalar dx1 = gx[ix] - gx[ix-1];
	PetscScalar dy1 = gy[jy] - gy[jy-1];
	PetscScalar dy2 = gy[jy+1] - gy[jy];
	/* right hand side value */
	if( !doSteady ){/* Normal RHS and i,j term for transient problems */
	  PetscScalar rho0 = problem->materials.materialRho[0];
	  ierr = VecSetValue(thermalRHS,idxnode, (rho0+rho[jy][ix])*Cp[jy][ix]/dt*lastT[jy][ix],ADD_VALUES); CHKERRQ(ierr); /* this is add values because it is assumed that the RHS already contains internal heating terms */
	  ierr = MatSetValue( thermalLHS, idxnode, idxnode, (rho0+rho[jy][ix])*Cp[jy][ix]/dt + 2.0/(dx2+dx1)*(kb/dx2 + ka/dx1) + 2.0/(dy1+dy2)*(kd/dy2+kc/dy1) , INSERT_VALUES); CHKERRQ(ierr); 
	}else{ /* special RHS for steady problem */
	  // do nothing to RHS - it contains only internal heating
	  ierr = MatSetValue( thermalLHS, idxnode, idxnode, 2.0/(dx2+dx1)*(kb/dx2 + ka/dx1) + 2.0/(dy1+dy2)*(kd/dy2+kc/dy1) , INSERT_VALUES); CHKERRQ(ierr);
	}
	/* left hand side values */
	/* i,j*/
	  // 	ierr = MatSetValue( thermalLHS, idxnode, idxnode, rho[jy][ix]*Cp[jy][ix]/dt + 2.0/(dx2+dx1)*(kb/dx2 + ka/dx1) + 2.0/(dy1+dy2)*(kd/dy2+kc/dy1) , INSERT_VALUES); CHKERRQ(ierr); 
	/* i,j+1 */
	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl][ixl+1], -2.0/(dx1+dx2)*kb/dx2, INSERT_VALUES); CHKERRQ(ierr);
	/* i,j-1 */
	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl][ixl-1], -2.0/(dx1+dx2)*ka/dx1, INSERT_VALUES); CHKERRQ(ierr);
	/* i-1,j */
	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl-1][ixl], -2.0/(dy1+dy2)*kc/dy1, INSERT_VALUES); CHKERRQ(ierr);
	/* i+1,j */
	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl+1][ixl], -2.0/(dy1+dy2)*kd/dy2, INSERT_VALUES); CHKERRQ(ierr);	
      }
      
      
    }/* end loop over x*/
  }/* end loop over y */
  /* assemble the system */
  ierr = MatAssemblyBegin( thermalLHS, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( thermalLHS, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin( thermalRHS);
  ierr = VecAssemblyEnd(thermalRHS);
  /* go back through RHS and enforce BCs */
  for( jy=y; jy<y+n;jy++ ){
    PetscInt jyl = jy-yg; 
    for(ix=x;ix<x+m;ix++){
      PetscInt ixl = ix-xg; 
      PetscInt idxnode = Tdof[jyl][ixl];
      /* boundary conditions*/
      if( ix == 0 ){/* left wall */
	
	
	//slab temperature
	PetscScalar thisT = slab_inflow_temperature( grid->x[ix], grid->y[jy], options->slabAngle);
	ierr = VecSetValue( thermalRHS, idxnode, thisT , INSERT_VALUES); CHKERRQ(ierr);


      } else if(ix == NX-1 && (grid->y[jy] <= plate_depth(grid->LX,options) || (jy>0 && (vx[jy][ix]+vx[jy-1][ix])/2.0 < 0.0  )) ){/* RIGHT WALL */
	PetscScalar thisT;
	if( grid->y[jy] <= plate_depth( grid->x[ix],options ) ){
	  thisT = plate_geotherm( grid->x[ix], grid->y[jy],options );
	}else{
	  thisT = 1573.0;
	}
	ierr = VecSetValue( thermalRHS, idxnode, thisT , INSERT_VALUES); CHKERRQ(ierr);
	
      } else if(ix == NX-1){
	/* enforce zero curvature on outflow section */
	
	/* 	PetscScalar dx2 = gx[ix]-gx[ix-1]; */
	/* 	PetscScalar dx1 = gx[ix] - gx[ix-1]; */
	/* 	PetscScalar kb = (kThermal[jy][ix]); */
	/* 	PetscScalar coef = -2.0/(dx1+dx2)*kb/dx2; */
	/* 	printf("ix, jy, dx1 dx2 kb coef = %d %d %le, %le, %le, %le\n",ix,jy,dx1,dx2,kb,coef); */
	/* 	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl][ixl],   2.0*coef, ADD_VALUES);CHKERRQ(ierr); */
	/* 	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl][ixl-1],-1.0*coef, ADD_VALUES);CHKERRQ(ierr); */
	ierr = VecSetValue(thermalRHS, idxnode, lastT[jy][ix], INSERT_VALUES); CHKERRQ(ierr);

      } else if( jy ==0){/* TOP - constant T*/
	if(options->thermalBCTop.type[0] == 0){
	  ierr = VecSetValue( thermalRHS, idxnode, options->thermalBCTop.value[0], INSERT_VALUES); CHKERRQ(ierr);
	}else if(options->thermalBCTop.type[0] == 1){
	  ierr = VecSetValue( thermalRHS, idxnode, 0.0 , INSERT_VALUES); CHKERRQ(ierr);
	}
      } else if( jy == NY-1){/* Bottom */
	// this boundary is 100% outflow.
	/* 	PetscScalar dy1 = gy[jy] - gy[jy-1]; */
	/* 	PetscScalar dy2 = gy[jy] - gy[jy-1]; */
	/* 	PetscScalar kd = (kThermal[jy][ix]); */
	/* 	PetscScalar coef =  -2.0/(dy1+dy2)*kd/dy2; */
 	//T[i,j+1] = T[i,j] + T[i,j] - T[i,j-1]
	  /* 	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl][ixl], 2.0*coef, ADD_VALUES);CHKERRQ(ierr); */
	  /* 	ierr = MatSetValue( thermalLHS, idxnode, Tdof[jyl-1][ixl], -1.0*coef, ADD_VALUES);CHKERRQ(ierr); */
	ierr = VecSetValue( thermalRHS, idxnode, lastT[jy][ix], INSERT_VALUES); CHKERRQ(ierr);
      } else{
	/* do nothing */
      }
    }
  }
  ierr = MatAssemblyBegin( thermalLHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( thermalLHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecAssemblyBegin( thermalRHS);
  ierr = VecAssemblyEnd( thermalRHS);
  /* free the DOF Map*/
  for(jy=0;jy<ng;jy++){
    ierr=PetscFree(Tdof[jy]);CHKERRQ(ierr);
  }
  ierr = PetscFree(Tdof);CHKERRQ(ierr);

  /* restore borrowed arrays*/
  ierr = DMDAVecRestoreArray( grid->da, kThermall, &kThermal); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray( grid->da, rhol, &rho); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray( grid->da, Cpl, &Cp); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray( grid->da, lastTl, &lastT); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray( grid->da, vxl, &vx); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray( grid->da, vyl, &vy); CHKERRQ(ierr);
  
  /* destroy local vectors*/
  ierr = VecDestroy(& kThermall );CHKERRQ(ierr);
  ierr = VecDestroy(& rhol ); CHKERRQ(ierr);
  ierr = VecDestroy(& Cpl ); CHKERRQ(ierr);
  ierr = VecDestroy(& lastTl ); CHKERRQ(ierr);
  ierr = VecDestroy(& vxl ); CHKERRQ(ierr);
  ierr = VecDestroy(& vyl ); CHKERRQ(ierr);


  PetscLogStagePop();
  PetscFunctionReturn(ierr);
}
  
