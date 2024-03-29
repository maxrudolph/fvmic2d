/* subroutine to form LHS and RHS for stokes flow problem */
/* forms LHS and RHS for visco-elasto-plastic problem with INcompressible continuity*/
/* assumes REGULAR GRID!!! */

#include "fdcode.h"
#include "vep_system.h"
#include "profile.h"
#include "benchmarkInitialConditions.h"
#define NADD 27 /* number of entries to add to matrix at a time*/



PetscErrorCode formVEPSystem(NodalFields *nodalFields, GridData *grid, Mat LHS,Mat P,  Vec RHS, PetscScalar *Kbond, PetscScalar *Kcont, PetscScalar gy, PetscScalar dt, Options *options, BoundaryValues *bv){
  /* Form system of equations for viscoelastoplastic problem */
  /* P is a preconditioner matrix, which can be optionally populated during the LHS assembly process */


  PetscErrorCode ierr;
  PetscMPIInt rank;//,size;
  PetscInt NX = grid->NX;
  PetscInt NY = grid-> NY;

  PetscScalar LX = grid-> LX; 
  PetscScalar LY = grid->LY; 

  PetscScalar gx = options->gx; /* this can be changed to add horizontal gravity*/
  PetscScalar gz = 0.0;

  PetscScalar svx = options->slabVelocity * cos( options->slabAngle );
  PetscScalar svy = options->slabVelocity * sin( options->slabAngle );
  

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  setLogStage( LOG_MECH_ASM );
  ierr = VecZeroEntries( RHS ); CHKERRQ(ierr);
  ierr = MatZeroEntries( LHS ); CHKERRQ(ierr);



  /* soxxZ etc should be global vectors complete with ghost information*/
  Vec soxxZg,soyyZg,sozzZg,soxyZg,etaSZg,etaNZg,X;
  PetscFunctionBegin;
  ierr=VecDuplicate( nodalFields->soxx, &soxxZg);CHKERRQ(ierr);
  ierr=VecDuplicate( nodalFields->soxx, &soyyZg);CHKERRQ(ierr);
  //  ierr=VecDuplicate( nodalFields->soxx, &sozzZg);CHKERRQ(ierr);
  ierr=VecDuplicate( nodalFields->soxx, &soxyZg);CHKERRQ(ierr);
  //ierr=VecDuplicate( nodalFields->soxx, &soxzZg);CHKERRQ(ierr);
  //ierr=VecDuplicate( nodalFields->soxx, &soyzZg);CHKERRQ(ierr);

  ierr=VecDuplicate( nodalFields->soxx, &etaSZg);CHKERRQ(ierr);
  ierr=VecDuplicate( nodalFields->soxx, &etaNZg);CHKERRQ(ierr);
  //ierr=VecDuplicate( nodalFields->soxx, &etavxZg);CHKERRQ(ierr);
  //ierr=VecDuplicate( nodalFields->soxx, &etavyZg);CHKERRQ(ierr);

  ierr=VecDuplicate( nodalFields->soxx, &X);CHKERRQ(ierr);/*viscoelasticity factor*/

/*   PetscScalar *soxxZ, *soxyZ, *soxzZ, *soyzZ, *soyyZ; */   /*last stresses modified by viscoelasticity factor*/
  /* the visco-elasto-plastic implementation looks just like a viscous implementation but with additional terms on the RHS and a modified effective viscosity on the LHS*/
  /* I group together the viscosity (eta_vp) with Z, the 'viscoelasticity factor'*/
  /* Z = mu*dt/(mu*dt+etavp) */

  /* copy etaN to X*/
  ierr=VecCopy(nodalFields->etaN,X);
  ierr=VecAXPY(X,dt,nodalFields->muN);/* X now holds etaN+dt*muN */
  ierr=VecPointwiseDivide(X,nodalFields->etaN,X);/* X now holds etaN/(etaN+dt*muN)*/
  /* scale soxx and soyy*/
  ierr=VecCopy(nodalFields->soxx,soxxZg);CHKERRQ(ierr);
  ierr=VecPointwiseMult(soxxZg,X,soxxZg);CHKERRQ(ierr);
  ierr=VecCopy(nodalFields->soyy,soyyZg);CHKERRQ(ierr);
  ierr=VecPointwiseMult(soyyZg,X,soyyZg);CHKERRQ(ierr);
  /* now negate X and add 1.0*/
  ierr=VecScale(X,-1.0);CHKERRQ(ierr);
  ierr=VecShift(X,1.0);CHKERRQ(ierr);
  ierr=VecCopy(nodalFields->etaN,etaNZg);CHKERRQ(ierr);
  ierr=VecPointwiseMult(etaNZg,X,etaNZg);CHKERRQ(ierr);

  ierr=VecCopy(nodalFields->etaS,X);CHKERRQ(ierr);
  ierr=VecAXPY(X,dt,nodalFields->muS);CHKERRQ(ierr);
  ierr=VecPointwiseDivide(X,nodalFields->etaS,X);CHKERRQ(ierr);
/*   ierr=VecView(X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
  ierr=VecCopy(nodalFields->soxy,soxyZg);CHKERRQ(ierr);
  ierr=VecPointwiseMult(soxyZg,X,soxyZg);CHKERRQ(ierr);
  /*   ierr=VecCopy(nodalFields->soxz,soxzZg);CHKERRQ(ierr); */
  /*   ierr=VecPointwiseMult(soxzZg,X,soxzZg);CHKERRQ(ierr); */

  /*   ierr=VecCopy(nodalFields->soyz,soyzZg);CHKERRQ(ierr); */
  /*   ierr=VecPointwiseMult(soyzZg,X,soyzZg);CHKERRQ(ierr); */
  /* now negate X and add 1.0*/
  ierr=VecScale(X,-1.0);CHKERRQ(ierr);
  ierr=VecShift(X,1.0);CHKERRQ(ierr);

  ierr=VecCopy(nodalFields->etaS,etaSZg); CHKERRQ(ierr);
  ierr=VecPointwiseMult(etaSZg,X,etaSZg); CHKERRQ(ierr);

  /* out-of-plane flow uses etavx for viscosity when computing sxz */
  //ierr=VecCopy(nodalFields->etavx,X);
  //ierr=VecAXPY(X,dt,nodalFields->muvx);/* X now holds etaN+dt*muN */
  //ierr=VecPointwiseDivide(X,nodalFields->etavx,X);/* X now holds etaN/(etaN+dt*muN)*/
  /* soxz adjustment would go here if soxz was evaluated at vx nodes */
  //ierr=VecScale(X,-1.0);CHKERRQ(ierr);
  //ierr=VecShift(X,1.0);CHKERRQ(ierr);
  //ierr=VecCopy(nodalFields->etavx,etavxZg);CHKERRQ(ierr);
  //ierr=VecPointwiseMult(etavxZg,X,etavxZg);CHKERRQ(ierr);
  
  //ierr=VecCopy(nodalFields->etavy,X);
  //ierr=VecAXPY(X,dt,nodalFields->muvy);/* X now holds etaN+dt*muN */
  //  ierr=VecPointwiseDivide(X,nodalFields->etavy,X);/* X now holds etaN/(etaN+dt*muN)*/

  //ierr=VecScale(X,-1.0);CHKERRQ(ierr);
  //ierr=VecShift(X,1.0);CHKERRQ(ierr);
  //ierr=VecCopy(nodalFields->etavy,etavyZg);CHKERRQ(ierr);
  //ierr=VecPointwiseMult(etavyZg,X,etavyZg);CHKERRQ(ierr);


  /* compute Kbond and Kcont*/
  PetscScalar etamin;
  {
    PetscScalar dx=grid->x[1]-grid->x[0];
    PetscScalar dy=grid->y[1]-grid->y[0];
    ierr=VecMin(etaSZg,PETSC_NULL,&etamin);CHKERRQ(ierr);    
    Kbond[0] = 4*etamin/((dx+dy)*(dx+dy));
    Kcont[0] = 2.0*etamin/(dx+dy);
  }

  if(!rank) printf("found etamin = %e\n, Kcont = %e, Kbond = %e\n",etamin,Kcont[0],Kbond[0]);

  /* now calculate mass loss from system */
  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  Vec rhol,rhodotl,etaNZl,etaSZl,soxxZl,soyyZl,soxyZl;/* declare other local vecs here*/
  PetscScalar **rho,**rhodot, **etaNZ, **etaSZ, **soxxZ,**soyyZ, **soxyZ;
  ierr=DMCreateLocalVector(grid->da, &rhol);CHKERRQ(ierr);
  ierr=VecDuplicate(rhol,&rhodotl);CHKERRQ(ierr);
  ierr=VecDuplicate(rhol,&etaNZl);CHKERRQ(ierr);
  ierr=VecDuplicate(rhol,&etaSZl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(rhol,&etavxZl);CHKERRQ(ierr); */
/*   ierr=VecDuplicate(rhol,&etavyZl);CHKERRQ(ierr); */
  ierr=VecDuplicate(rhol,&soxxZl);CHKERRQ(ierr);
  ierr=VecDuplicate(rhol,&soxyZl);CHKERRQ(ierr);
  ierr=VecDuplicate(rhol,&soyyZl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(rhol,&soxzZl);CHKERRQ(ierr); */
/*   ierr=VecDuplicate(rhol,&soyzZl);CHKERRQ(ierr); */

  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->rhodot,INSERT_VALUES,rhodotl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->rhodot,INSERT_VALUES,rhodotl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,etaNZg,INSERT_VALUES,etaNZl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,etaNZg,INSERT_VALUES,etaNZl);CHKERRQ(ierr);

  ierr=DMGlobalToLocalBegin(grid->da,etaSZg,INSERT_VALUES,etaSZl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,etaSZg,INSERT_VALUES,etaSZl);CHKERRQ(ierr);

/*   ierr=DMGlobalToLocalBegin(grid->da,etavxZg,INSERT_VALUES,etavxZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,etavxZg,INSERT_VALUES,etavxZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalBegin(grid->da,etavyZg,INSERT_VALUES,etavyZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,etavyZg,INSERT_VALUES,etavyZl);CHKERRQ(ierr); */

  ierr=DMGlobalToLocalBegin(grid->da,soxxZg,INSERT_VALUES,soxxZl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,soxxZg,INSERT_VALUES,soxxZl);CHKERRQ(ierr);

/*   ierr=DMGlobalToLocalBegin(grid->da,soxyZg,INSERT_VALUES,soxyZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,soxyZg,INSERT_VALUES,soxyZl);CHKERRQ(ierr); */

/*   ierr=DMGlobalToLocalBegin(grid->da,soxzZg,INSERT_VALUES,soxzZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,soxzZg,INSERT_VALUES,soxzZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalBegin(grid->da,soyzZg,INSERT_VALUES,soyzZl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,soyzZg,INSERT_VALUES,soyzZl);CHKERRQ(ierr); */


  /* get arrays of values*/
  ierr=DMDAVecGetArray(grid->da,rhol,&rho);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,rhodotl,&rhodot);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,etaSZl,&etaSZ);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,etaNZl,&etaNZ);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray(grid->da,etavxZl,&etavxZ);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray(grid->da,etavyZl,&etavyZ);CHKERRQ(ierr); */

  ierr=DMDAVecGetArray(grid->da,soxxZl,&soxxZ);CHKERRQ(ierr);
  ierr=DMDAVecGetArray(grid->da,soxyZl,&soxyZ);CHKERRQ(ierr);
  /*   ierr=DMDAVecGetArray(grid->da,soxzZl,&soxzZ);CHKERRQ(ierr); */
  /*   ierr=DMDAVecGetArray(grid->da,soyzZl,&soyzZ);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray(grid->da,soyyZl,&soyyZ);CHKERRQ(ierr);


  PetscInt y1=y;
  PetscInt x1=x;
  if(y==0) y1=1;
  if(x==0) x1=1;/* do not loop over ghost cells*/


  /* now we want to assemble the global LHS, which requires knowing the mapping from local to global indices*/
  ISLocalToGlobalMapping ltogm;
  const  PetscInt *globalIdx;
  ierr=DMGetLocalToGlobalMapping(grid->da,&ltogm); CHKERRQ(ierr);
  ierr=ISLocalToGlobalMappingGetIndices(ltogm,&globalIdx);CHKERRQ(ierr);

  ierr=DMDAGetGhostCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);

  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  PetscInt xg,yg,mg,ng;
  ierr=DMDAGetGhostCorners(grid->da,&xg,&yg,PETSC_NULL,&mg,&ng,PETSC_NULL);CHKERRQ(ierr);


  /*make 2d arrays to easily reference local vx,vy,p dofs*/
  PetscInt **vxdof;
  PetscInt **vydof;
  PetscInt **pdof;
  //  PetscInt **vzdof;
  PetscInt ix,jy;
  PetscMalloc( (ng+1)*sizeof(PetscInt *), &vxdof);
  PetscMalloc( (ng+1)*sizeof(PetscInt *), &vydof);
  PetscMalloc( ng*sizeof(PetscInt *), &pdof);
  //  PetscMalloc( ng*sizeof(PetscInt *), &vzdof);
  
  vxdof++; vydof++;
  ierr=PetscMalloc( (mg+1)*sizeof(PetscInt), &vxdof[-1]);CHKERRQ(ierr); vxdof[-1]++;
  ierr=PetscMalloc( (mg+1)*sizeof(PetscInt), &vydof[-1]);CHKERRQ(ierr); vydof[-1]++;
  for(jy=0;jy<ng;jy++){
    ierr=PetscMalloc( (mg+1)*sizeof(PetscInt), &vxdof[jy]);CHKERRQ(ierr); vxdof[jy]++;
    ierr=PetscMalloc( (mg+1)*sizeof(PetscInt), &vydof[jy]);CHKERRQ(ierr); vydof[jy]++;
    ierr=PetscMalloc( mg*sizeof(PetscInt), &pdof[jy]);CHKERRQ(ierr);
    //ierr=PetscMalloc( mg*sizeof(PetscInt), &vzdof[jy]);CHKERRQ(ierr);
  }

  for(jy=-1;jy<ng;jy++){
    for(ix=-1;ix<mg;ix++){
      if( ix == -1 || jy == -1){
	vxdof[jy][ix] = -1;
	vydof[jy][ix] = -1;
      }else{
	pdof[jy][ix] = 3*globalIdx[ix+mg*jy]+DOF_P;
	vxdof[jy][ix] = 3*globalIdx[ix+mg*jy]+DOF_U;
	vydof[jy][ix] = 3*globalIdx[ix+mg*jy]+DOF_V;
	//vzdof[jy][ix] = globalIdx[ix+mg*jy];/* vz goes into a separate matrix*/
      }
    }
  }
  ierr = ISLocalToGlobalMappingRestoreIndices(ltogm,&globalIdx); CHKERRQ(ierr);

  /* print out the local vxdofs*/

  for(jy=y;jy<y+n;jy++){
    PetscInt jyl=jy-yg;
    for(ix=x;ix<x+m;ix++){
      PetscInt ixl=ix-xg;
      /* arrays to hold values to add.*/
      PetscScalar vals[NADD];
      PetscInt rowidx;
      PetscInt colidx[NADD];

      /* Continuity equation */
      if( ix>0 && jy > 0 && 
	  in_either( grid->xc[ix], grid->y[jy-1], options) && 
	  in_either( grid->xc[ix], grid->y[jy], options) && 
	  in_either( grid->x[ix] , grid->yc[jy], options) && 
	  in_either( grid->x[ix-1], grid->yc[jy], options) ){
	
	ierr = MatSetValue(LHS,pdof[jyl][ixl], pdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl], 0.0,INSERT_VALUES); CHKERRQ(ierr);
      }else if( jy == 0 ){/* make ghost nodes have same pressure as first row so that null space will be constant */
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],   1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl+1][ixl],-1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES); CHKERRQ(ierr);
      }else if( !grid->xperiodic && ix == 0){
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],   1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl+1],-1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES); CHKERRQ(ierr);	       
	
      } else if( FIX_PRESSURE && ix==NX-1 && jy > 0 && grid->yc[jy] > plate_depth(grid->LX,options) && grid->yc[jy-1] <= plate_depth(grid->LX,options) ){/* pressure specified in one cell */	
	printf("Fixing pressure\n");
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES);CHKERRQ(ierr);
      } else {
	PetscScalar dx=grid->x[ix]-grid->x[ix-1];
	PetscScalar dy=grid->y[jy]-grid->y[jy-1];

	ierr = MatSetValue(LHS,pdof[jyl][ixl],vxdof[jyl-1][ixl],Kcont[0]/dx,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],vxdof[jyl-1][ixl-1],-Kcont[0]/dx,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],vydof[jyl][ixl-1],Kcont[0]/dy,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],vydof[jyl-1][ixl-1],-Kcont[0]/dy,INSERT_VALUES);CHKERRQ(ierr);


	/* this term adds volumetric changes due to the equation of state (phase changes)*/
	/* this is not the same thing as compressibility - that would require additional terms in the momentum equations too */
	PetscScalar rhobar = (rho[jy-1][ix-1] + rho[jy-1][ix] + rho[jy][ix-1] + rho[jy][ix])*0.25;
	PetscScalar rhodotbar =  (rhodot[jy-1][ix-1] + rhodot[jy-1][ix] + rhodot[jy][ix-1] + rhodot[jy][ix])*0.25; 
	//PetscScalar Rval = -1.0/rhobar*rhodotbar*Kcont[0];
	PetscScalar Rval = 0.0;

	ierr = VecSetValue(RHS,pdof[jyl][ixl],Rval,INSERT_VALUES);CHKERRQ(ierr);
      }/* end continuity equation*/
      
      
      /*         %x-momentum */
      if( in_slab( grid->x[ix], grid->yc[jy+1], options ) ){   //slab
	ierr = MatSetValue(LHS, vxdof[jyl][ixl],vxdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);	
	ierr = VecSetValue(RHS, vxdof[jyl][ixl],Kbond[0]*svx,INSERT_VALUES);CHKERRQ(ierr);       
      }else if( in_plate( grid->x[ix], grid->yc[jy+1], options) ){ //rigid over-riding plate
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vxdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);	
	ierr = VecSetValue(RHS, vxdof[jyl][ixl], 0.0 ,INSERT_VALUES);CHKERRQ(ierr);       		  
      }else if( ix == NX-1 && jy < NY-1){
	PetscScalar eta = etaNZ[jy+1][ix];
	PetscScalar dx = grid->x[ix] - grid->x[ix-1];
	ierr = MatSetValue(LHS,vxdof[jyl][ixl], vxdof[jyl][ixl], 2.0*eta/dx, INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,vxdof[jyl][ixl], vxdof[jyl][ixl-1], -2.0*eta/dx, INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,vxdof[jyl][ixl], pdof[jyl+1][ixl], -Kcont[0], INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS,vxdof[jyl][ixl], 0.0, INSERT_VALUES);CHKERRQ(ierr);

      }else if( jy == NY-1 && ix < NX-1){ /* BOTTOM BOUNDARY */
	  PetscScalar dxc = grid->xc[ix+1] - grid->xc[ix];
	  PetscScalar dyc = grid->yc[jy] - grid->yc[jy-1];
	  ierr=MatSetValue(LHS,vxdof[jyl][ixl], vxdof[jyl][ixl], 1.0/dyc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vxdof[jyl][ixl], vxdof[jyl-1][ixl],  -1.0/dyc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vxdof[jyl][ixl], vydof[jyl  ][ixl], 1.0/dxc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vxdof[jyl][ixl], vydof[jyl  ][ixl-1],  -1.0/dxc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=VecSetValue(RHS,vxdof[jyl][ixl], 0.0, INSERT_VALUES);CHKERRQ(ierr);	
      }else if (jy == NY-1 && ix == NX-1){
	PetscScalar dyc = grid->yc[NY]-grid->yc[NY-1];	
	printf("dyc = %le\n",dyc);
	ierr=MatSetValue(LHS,vxdof[jyl][ixl], vxdof[jyl][ixl],    1.0/dyc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr=MatSetValue(LHS,vxdof[jyl][ixl], vxdof[jyl-1][ixl], -1.0/dyc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr=VecSetValue(RHS,vxdof[jyl][ixl], 0.0, INSERT_VALUES);CHKERRQ(ierr);	
      }else{
	if( jy == NY-1 ){ fprintf(stderr,"jy==NY-1 !!!\n"); }
	/* normal x-stokes stencil */
	rowidx = vxdof[jyl][ixl];
	PetscScalar v0000 = grid->xc[ix];
	PetscScalar v0001 = Kcont[0];
	PetscScalar v0002 = grid->x[ix];
	PetscScalar v0003 = grid->y[jy];
	PetscInt v0004 = 1 + ix;
	PetscInt v0005 = 1 + jy;
	PetscScalar v0006 = etaSZ[jy][ix];
	PetscScalar v0007 = grid->xc[v0004];
	PetscScalar v0008 = grid->y[v0005];
	PetscScalar v0009 = grid->yc[v0005];
	PetscScalar v0010 = etaSZ[v0005][ix];
	PetscScalar v0011 = etaNZ[v0005][ix];
	PetscScalar v0012 = etaNZ[v0005][v0004];
	if( ix == NX-1 ) v0012 = v0011;
	PetscScalar v0013 = 1.0/(v0000 - v0007);
	PetscScalar v0014 = 1.0/(v0002 - grid->x[v0004]);
	if( ix == NX-1){
	  v0014 = 1.0/(grid->x[NX-2] - grid->x[NX-1]);
	}

	PetscScalar v0015 = 1/(v0003 - v0008);
	PetscScalar v0016 = 1/(-v0009 + grid->yc[jy]);
	PetscScalar v0017 = 1/(-v0002 + grid->x[-1 + ix]);
	PetscScalar v0018 = 1/(v0009 - grid->yc[2 + jy]);
	colidx[0] = -1;
	colidx[1] = -1;
	colidx[2] = -1;
	vals[0]=0;
	vals[1]=0;
	vals[2]=0;
	colidx[3] = vxdof[jyl+0][ixl-1];
	colidx[4] = vydof[jyl+0][ixl-1];
	colidx[5] = -1;
	vals[3]=(4*v0011*v0013*v0017)/3.;
	vals[4]=((3*v0006 - 2*v0011)*v0013*v0015)/3.;
	vals[5]=0;
	colidx[6] = -1;
	colidx[7] = vydof[jyl+1][ixl-1];
	colidx[8] = -1;
	vals[6]=0;
	vals[7]=((-3*v0010 + 2*v0011)*v0013*v0015)/3.;
	vals[8]=0;
	colidx[9] = vxdof[jyl-1][ixl+0];
	colidx[10] = -1;
	colidx[11] = -1;
	vals[9]=v0006*v0015*v0016;
	vals[10]=0;
	vals[11]=0;
	colidx[12] = vxdof[jyl+0][ixl+0];
	colidx[13] = vydof[jyl+0][ixl+0];
	colidx[14] = -1;
	vals[12]=(-4*v0012*v0013*v0014)/3. - (4*v0011*v0013*v0017)/3. + \
	  (v0006*v0016 + v0010*v0018)/(-v0003 + v0008);
	vals[13]=((-3*v0006 + 2*v0012)*v0013*v0015)/3.;
	vals[14]=0;
	colidx[15] = vxdof[jyl+1][ixl+0];
	colidx[16] = vydof[jyl+1][ixl+0];
	colidx[17] = pdof[jyl+1][ixl+0];
	vals[15]=v0010*v0015*v0018;
	vals[16]=((3*v0010 - 2*v0012)*v0013*v0015)/3.;
	vals[17]=v0001/(-v0000 + v0007);
	colidx[18] = -1;
	colidx[19] = -1;
	colidx[20] = -1;
	vals[18]=0;
	vals[19]=0;
	vals[20]=0;
	colidx[21] = vxdof[jyl+0][ixl+1];
	colidx[22] = -1;
	colidx[23] = -1;
	vals[21]=(4*v0012*v0013*v0014)/3.;
	vals[22]=0;
	vals[23]=0;
	colidx[24] = -1;
	colidx[25] = -1;
	colidx[26] = pdof[jyl+1][ixl+1];
	vals[24]=0;
	vals[25]=0;
	vals[26]=v0001*v0013;

	PetscScalar Rval = -gx*(rho[jy][ix]+rho[jy+1][ix])/2.0;//-(soxxZ[jy+1][ix+1]-soxxZ[jy+1][ix])/(grid->xc[ix+1]-grid->xc[ix])-(soxyZ[jy+1][ix]-soxyZ[jy][ix])/(grid->y[jy+1]-grid->y[jy]);
	
	ierr = VecSetValue( RHS, vxdof[jyl][ixl],Rval,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues( LHS, 1,&rowidx,27,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);	
      } 

      /* %y-momentum */
      /* y-equation has implicit bcs on left, explicit on top, bottom, and right */
      if( in_slab( grid->xc[ix+1], grid->y[jy], options ) ){/* slab internal BC */
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl][ixl], Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vydof[jyl][ixl], Kbond[0] * svy,INSERT_VALUES);CHKERRQ(ierr);	

      }else if( in_plate( grid->xc[ix+1], grid->y[jy], options ) ){/* plate internal BC */
	ierr = MatSetValue(LHS, vydof[jyl][ixl],vydof[jyl][ixl], Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vydof[jyl][ixl], 0.0,INSERT_VALUES);CHKERRQ(ierr);     	
      }else if( jy == NY-1 && ix < NX-1){
	PetscScalar eta = etaNZ[jy][ixl+1];
	PetscScalar dy = grid->y[jy]-grid->y[jy-1];
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl][ixl], 2.0*eta/dy, INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl-1][ixl],-2.0*eta/dy, INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], pdof[jyl][ixl+1],-Kcont[0], INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS,vydof[jyl][ixl], 0.0, INSERT_VALUES);CHKERRQ(ierr);
	
      }else if( ix == NX-1 && !in_plate(grid->x[ix],grid->y[jy],options) /*&& jy < NY-1 */){// not in plate - zero shear stress on RIGHT BOUNDARY
	  PetscScalar dxc = grid->xc[ix+1] - grid->xc[ix];
	  PetscScalar dyc = grid->yc[jy+1] - grid->yc[jy];
	  ierr=MatSetValue(LHS,vydof[jyl][ixl], vxdof[jyl][ixl],   1.0/dyc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vydof[jyl][ixl], vxdof[jyl-1][ixl],-1.0/dyc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vydof[jyl][ixl], vydof[jyl][ixl],   1.0/dxc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vydof[jyl][ixl], vydof[jyl][ixl-1],-1.0/dxc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=VecSetValue(RHS,vydof[jyl][ixl], 0.0, INSERT_VALUES);CHKERRQ(ierr);
      }else if( ix == NX-1 && jy == NY-1 ){//fix bottom boundary
	PetscScalar dxc = grid->xc[NX]-grid->xc[NX-1];
	printf("dxc=%le\n",dxc);
	ierr = MatSetValue(LHS, vydof[jyl][ixl],vydof[jyl][ixl],    1.0/dxc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl],vydof[jyl][ixl-1], -1.0/dxc*Kcont[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vydof[jyl][ixl], 0.0,INSERT_VALUES);CHKERRQ(ierr);     	

      }else {
	rowidx = vydof[jyl][ixl];
	/* normal y-stokes stencil*/
	PetscScalar v0000 = grid->yc[jy];
	PetscScalar v0001 = Kcont[0];
	PetscScalar v0002 = grid->y[jy];
	PetscInt v0003 = 1 + jy;
	PetscInt v0004 = 1 + ix;
	PetscScalar v0005 = etaSZ[jy][ix];
	PetscScalar v0006 = grid->yc[v0003];
	PetscScalar v0007 = grid->xc[v0004];
	PetscScalar v0008 = etaSZ[jy][v0004];
	PetscScalar v0009 = etaNZ[jy][v0004];
	PetscScalar v0010 = etaNZ[v0003][v0004];
	if( jy == NY-1 ){
	  v0010 = v0009;
	}
	PetscScalar v0011 = 1/(v0000 - v0006);
	PetscScalar v0012 = 1/(-v0000 + v0006);
	PetscScalar v0013 = 1/(grid->x[ix] - grid->x[v0004]);
	PetscScalar v0014 = 1/(v0002 - grid->y[v0003]);
	if( jy == NY-1 ){
	  v0014 = v0013;
	}
	PetscScalar v0015 = 1/(-v0002 + grid->y[-1 + jy]);
	PetscScalar v0016 = 1/(-v0007 + grid->xc[ix]);
	PetscScalar v0017 = 1/(v0007 - grid->xc[2 + ix]);
	//	if( jy == NY-1) printf("grid->xc[2+ix] = %le\n",grid->xc[2+ix]);
	colidx[0] = -1;
	colidx[1] = -1;
	colidx[2] = -1;
	vals[0]=0;
	vals[1]=0;
	vals[2]=0;
	colidx[3] = -1;
	colidx[4] = vydof[jyl+0][ixl-1];
	colidx[5] = -1;
	vals[3]=0;
	vals[4]=v0005*v0013*v0016;
	vals[5]=0;
	colidx[6] = -1;
	colidx[7] = -1;
	colidx[8] = -1;
	vals[6]=0;
	vals[7]=0;
	vals[8]=0;
	colidx[9] = vxdof[jyl-1][ixl+0];
	colidx[10] = vydof[jyl-1][ixl+0];
	colidx[11] = -1;
	vals[9]=((3*v0005 - 2*v0009)*v0011*v0013)/3.;
	vals[10]=(4*v0009*v0011*v0015)/3.;
	vals[11]=0;
	colidx[12] = vxdof[jyl+0][ixl+0];
	colidx[13] = vydof[jyl+0][ixl+0];
	colidx[14] = -1;
	vals[12]=((-3*v0005 + 2*v0010)*v0011*v0013)/3.;
	vals[13]=(4*v0012*(v0010*v0014 + v0009*v0015))/3. - v0005*v0013*v0016 \
	  - v0008*v0013*v0017;
	vals[14]=0;
	colidx[15] = -1;
	colidx[16] = vydof[jyl+1][ixl+0];
	colidx[17] = -1;
	vals[15]=0;
	vals[16]=(4*v0010*v0011*v0014)/3.;
	vals[17]=0;
	colidx[18] = vxdof[jyl-1][ixl+1];
	colidx[19] = -1;
	colidx[20] = -1;
	vals[18]=((-3*v0008 + 2*v0009)*v0011*v0013)/3.;
	vals[19]=0;
	vals[20]=0;
	colidx[21] = vxdof[jyl+0][ixl+1];
	colidx[22] = vydof[jyl+0][ixl+1];
	colidx[23] = pdof[jyl+0][ixl+1];
	vals[21]=((3*v0008 - 2*v0010)*v0011*v0013)/3.;
	vals[22]=v0008*v0013*v0017;
	vals[23]=v0001*v0012;
	colidx[24] = -1;
	colidx[25] = -1;
	colidx[26] = pdof[jyl+1][ixl+1];
	vals[24]=0;
	vals[25]=0;
	vals[26]=v0001*v0011;

	PetscScalar Rval = -gy*(rho[jy][ix]+rho[jy][ix+1])/2.0;//-(soyyZ[jy+1][ix+1]-soyyZ[jy][ix+1])/(grid->yc[jy+1]-grid->yc[jy])-(soxyZ[jy][ix+1]-soxyZ[jy][ix])/(grid->x[ix+1]-grid->x[ix]); 
	
	ierr = VecSetValue( RHS, vydof[jyl][ixl],Rval,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues( LHS, 1,&rowidx,27,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
      }

    }/* end loop over x grid-line*/

  }/* end loop over y grid-line*/
  ierr = VecAssemblyBegin(RHS);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(RHS);CHKERRQ(ierr);

  ierr = MatAssemblyBegin( LHS, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( LHS, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);


  /* go back through RHS and modify to account for boundary conditions using ADD_VALUES */
  for(jy=y;jy<y+n;jy++){
    PetscInt jyl=jy-yg;
    for(ix=x;ix<x+m;ix++){
      PetscInt ixl=ix-xg;
      
      /* x-stokes equation */     
      /* top boundary - adjust x-stokes equation to account for ghosted values*/
      
      if(0 && ix == NX-1 && jy < NY-1 && !in_plate( grid->x[ix], grid->yc[jy+1], options) ){// modify for outflow BC
	/* get coefficient for vx[ix-1][jy] */
	PetscScalar dx = grid->x[ix]-grid->x[ix-1];
	PetscScalar dy = grid->y[jy+1]-grid->y[jy];

	PetscScalar Pcoef = Kcont[0]/(grid->xc[ix-1] - grid->xc[ix]);// what would normally multiply Pi+1,j+1
	PetscScalar eta = etaNZ[jy+1][ix];
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vxdof[jyl][ixl],    Pcoef * 2.0*eta/dx/Kcont[0] , ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vydof[jyl+1][ixl], -Pcoef * 2.0*eta/dy/Kcont[0] , ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vydof[jyl][ixl],    Pcoef * 2.0*eta/dy/Kcont[0] , ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vxdof[jyl][ixl-1], -Pcoef * 2.0*eta/dx/Kcont[0] , ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], pdof[jyl+1][ixl],  -Pcoef , ADD_VALUES); CHKERRQ(ierr);	

	
      }
      
      /* y-stokes */
      // bottom boundary
      if( 0 &&  jy == NY-1 && ix < NX-1 && in_slab(grid->xc[ix+1], 2.0*grid->y[jy] - grid->y[jy-1], options) && !in_slab(grid->xc[ix+1],grid->y[jy],options) ){
	// special case where ghost node is in slab
	PetscScalar eta = etaNZ[jy][ix+1];
	PetscScalar Pcoef = Kcont[0]/(grid->yc[jy] - grid->yc[1 + jy]);
	PetscScalar dx = grid->x[ix+1]-grid->x[ix];
	PetscScalar dy = grid->y[jy]-grid->y[jy-1];
	PetscScalar vyg = svy;

	//ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl][ixl],    Pcoef*2.0*eta/dy/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	//ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl][ixl+1], -Pcoef*2.0*eta/dx/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	//ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl][ixl],    Pcoef*2.0*eta/dx/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl-1][ixl], -Pcoef*2.0*eta/dy/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], pdof[jyl][ixl+1],  -Pcoef, ADD_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vydof[jyl][ixl], -Pcoef*2.0*eta/dy/Kcont[0] * vyg, ADD_VALUES);CHKERRQ(ierr);
      } else if( 0 && jy == NY-1 && ix < NX-1 && !in_slab( grid->xc[ix+1], grid->y[jy], options) ){
	/* get coefficient for vx[ix-1][jy] */
	PetscScalar eta = etaNZ[jy][ix+1];
	PetscScalar Pcoef = Kcont[0]/(grid->yc[jy] - grid->yc[1 + jy]);
	PetscScalar dx = grid->x[ix+1]-grid->x[ix];
	PetscScalar dy = grid->y[jy]-grid->y[jy-1];
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl][ixl],    Pcoef*2.0*eta/dy/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl][ixl+1], -Pcoef*2.0*eta/dx/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl][ixl],    Pcoef*2.0*eta/dx/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl-1][ixl], -Pcoef*2.0*eta/dy/Kcont[0] , ADD_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS, vydof[jyl][ixl], pdof[jyl][ixl+1],  -Pcoef, ADD_VALUES);CHKERRQ(ierr);
      }
      

      
    }/* end x-loop */
  }/* end y-loop */

  ierr = VecAssemblyBegin(RHS);

  ierr = MatAssemblyBegin( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatAssemblyEnd( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecAssemblyEnd(RHS);CHKERRQ(ierr);





  /* assemble preconditioner */
  if( P != PETSC_NULL ){
    ierr = MatZeroEntries( P );CHKERRQ(ierr);
    ierr = MatCopy( LHS, P, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
    PetscBool pc_diag = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,"-stokes_ksp_build_pc_diag",&pc_diag,0);
    /*if( pc_diag )*/{
      printf("assembling preconditioner\n");
      for(jy=y;jy<y+n;jy++){
	PetscInt jyl=jy-yg;
	for(ix=x;ix<x+m;ix++){
	  PetscInt ixl=ix-xg;
	  PetscScalar bcval = 0.0*1.0*Kbond[0];
	  /**              %continuity*/		       
	  if( jy == 0 ){/* make ghost nodes have same pressure as first row so that null space will be constant */
	    ierr = MatSetValue(P,pdof[jyl][ixl],pdof[jyl][ixl],bcval,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(P,pdof[jyl][ixl],pdof[jyl+1][ixl],-bcval,INSERT_VALUES); CHKERRQ(ierr);
	  }else if( !grid->xperiodic && ix == 0){
	    ierr = MatSetValue(P,pdof[jyl][ixl],pdof[jyl][ixl],bcval,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(P,pdof[jyl][ixl],pdof[jyl][ixl+1],-bcval,INSERT_VALUES); CHKERRQ(ierr);	    
	  } else if( FIX_PRESSURE && jy==1 && ix==2 ){/* pressure specified in one cell */       
	    ierr = MatSetValue(P,pdof[jyl][ixl],pdof[jyl][ixl],bcval,INSERT_VALUES);CHKERRQ(ierr);

	  }else{
	    /* 	    PetscScalar dx = grid->x[ix+1]-grid->x[ix]; */
	    /* 	    PetscScalar dy = grid->y[jy+1]-grid->y[jy]; */
	    
	    ierr = MatSetValue( P, pdof[jyl][ixl], pdof[jyl][ixl], 0.0*etaNZ[jy][ix]/Kcont[0], INSERT_VALUES );CHKERRQ(ierr);
	  }
	}
      }
    }
    ierr = MatAssemblyBegin( P, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  ( P, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }

  ierr=DMDAVecRestoreArray(grid->da,rhol,&rho);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,rhodotl,&rhodot);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etaSZl,&etaSZ);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etaNZl,&etaNZ);CHKERRQ(ierr);
  //  ierr=DMDAVecRestoreArray(grid->da,etavxZl,&etavxZ);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray(grid->da,etavyZl,&etavyZ);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soxxZl,&soxxZ);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soyyZl,&soyyZ);CHKERRQ(ierr);

  ierr=DMDAVecRestoreArray(grid->da,soxyZl,&soxyZ);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray(grid->da,soxzZl,&soxzZ);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray(grid->da,soyzZl,&soyzZ);CHKERRQ(ierr);

  ierr=VecDestroy(& rhol);CHKERRQ(ierr);
  ierr=VecDestroy(&  rhodotl);CHKERRQ(ierr);
  ierr=VecDestroy(&  etaNZl);CHKERRQ(ierr);
  ierr=VecDestroy(&  etaSZl);CHKERRQ(ierr);
  //ierr=VecDestroy(&  etavxZl);CHKERRQ(ierr);
  //ierr=VecDestroy(&  etavyZl);CHKERRQ(ierr);
  ierr=VecDestroy(&  soxxZl);CHKERRQ(ierr);
  ierr=VecDestroy(&  soxyZl);CHKERRQ(ierr);
  //ierr=VecDestroy(&  soxzZl);CHKERRQ(ierr);
  //ierr=VecDestroy(&  soyzZl);CHKERRQ(ierr);
  ierr=VecDestroy(&  soyyZl);CHKERRQ(ierr);
  ierr=VecDestroy(&  soxxZg);CHKERRQ(ierr);
  ierr=VecDestroy(&  soyyZg);CHKERRQ(ierr);
  //ierr=VecDestroy(&  sozzZg);CHKERRQ(ierr);
  ierr=VecDestroy(&  soxyZg);CHKERRQ(ierr);
  //ierr=VecDestroy(&  soxzZg);CHKERRQ(ierr);
  //ierr=VecDestroy(&  soyzZg);CHKERRQ(ierr);

  ierr=VecDestroy(&  etaSZg);CHKERRQ(ierr);
  ierr=VecDestroy(&  etaNZg);CHKERRQ(ierr);
  //ierr=VecDestroy(&  etavxZg);CHKERRQ(ierr);
  //ierr=VecDestroy(&  etavyZg);CHKERRQ(ierr);

  ierr=VecDestroy(&  X);CHKERRQ(ierr);


  /* free dof maps*/
  //ierr=PetscFree(vxdof[-1]);CHKERRQ(ierr);
  //ierr=PetscFree(vydof[-1]);CHKERRQ(ierr);
  for(jy=0;jy<ng;jy++){
    vxdof[jy] --;
    vydof[jy] --;
    ierr=PetscFree(vxdof[jy]);CHKERRQ(ierr);
    ierr=PetscFree(vydof[jy]);CHKERRQ(ierr);
    ierr=PetscFree(pdof[jy]); CHKERRQ(ierr);
    //ierr=PetscFree(vzdof[jy]);CHKERRQ(ierr);
  }
  vxdof[-1] --;
  vydof[-1] --;
  ierr=PetscFree(vxdof[-1]);CHKERRQ(ierr);
  ierr=PetscFree(vydof[-1]);CHKERRQ(ierr);

  vxdof--; vydof--;
  ierr=PetscFree(vxdof);CHKERRQ(ierr);
  ierr=PetscFree(vydof);CHKERRQ(ierr);
  ierr=PetscFree(pdof);CHKERRQ(ierr);
  //  ierr=PetscFree(vzdof);CHKERRQ(ierr);

/*   PetscScalar rhsnorm; */
/*   ierr = VecNorm( RHS, NORM_2, &rhsnorm); CHKERRQ(ierr); */
/*   rhsnorm = 1.0/rhsnorm; */
/*   ierr = VecScale( RHS, rhsnorm ); */
/*   ierr = MatScale( LHS, rhsnorm ); */
/*   ierr = MatScale( P, rhsnorm); */



  PetscLogStagePop();
  PetscFunctionReturn(ierr);

}
