/* subroutine to form LHS and RHS for stokes flow problem */
/* forms LHS and RHS for visco-elasto-plastic problem with INcompressible continuity*/
/* assumes REGULAR GRID!!! */

#include "fdcode.h"
#include "anisotropic_system.h"
#define NADD 75 /* number of entries to add to matrix at a time*/

PetscErrorCode formAnisotropicSystem(NodalFields *nodalFields, GridData *grid, Mat LHS, Vec RHS, Mat LHSz, Vec RHSz, PetscScalar *Kbond, PetscScalar *Kcont, PetscScalar gy, PetscScalar dt, BoundaryValues *bv, Options *options){
  PetscErrorCode ierr;
  PetscMPIInt rank;

  PetscInt NX = grid->NX;
  PetscInt NY = grid-> NY;
  PetscScalar dx=grid->x[1]-grid->x[0];/* characteristic grid dimensions*/
  PetscScalar dy=grid->y[1]-grid->y[0];
  PetscScalar *xc = grid->xc;
  PetscScalar *yc = grid->yc;

  PetscScalar LX = grid->LX; 
  PetscScalar LY = grid->LY; 

  PetscScalar gx = 0.0; /* this can be changed to add horizontal gravity*/
  PetscInt NDOF = 3*NX*NY;

  //  PetscScalar gz = 0.0;
  Tensor22 **Nc, **Nb;/* viscoplasticity tensors*/

  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  Vec rhol;
  ierr=DMCreateLocalVector(grid->da, &rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->rho,INSERT_VALUES,rhol);CHKERRQ(ierr);

  PetscScalar **rho;
  ierr=DMDAVecGetArray(grid->da,rhol,&rho);CHKERRQ(ierr);

  /* get local viscoplasticity tensors */
  Vec Nbl, Ncl;
  ierr = DMCreateLocalVector(grid->tda, &Nbl);CHKERRQ(ierr);
  ierr = VecDuplicate(Nbl, &Ncl); CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(grid->tda,nodalFields->VPTensorB,INSERT_VALUES,Nbl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(grid->tda,nodalFields->VPTensorB,INSERT_VALUES,Nbl);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(grid->tda,nodalFields->VPTensorC,INSERT_VALUES,Ncl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(grid->tda,nodalFields->VPTensorC,INSERT_VALUES,Ncl);CHKERRQ(ierr);
  /* get viscoplasticity tensors */

  ierr = DMDAVecGetArray(grid->tda, Nbl, &Nb);
  ierr = DMDAVecGetArray(grid->tda, Ncl, &Nc);

  /* now calculate mass loss from system */
  PetscInt x,y,m,n;
  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);


  PetscScalar outflow = 0.0;
  if( bv->mechBCLeft.type[0] == 0){
    outflow -= LY*bv->mechBCLeft.value[0];
  }
  if( bv->mechBCRight.type[0] == 0){
    outflow += LY*bv->mechBCRight.value[0];
  }

  PetscScalar inflow = outflow/LX; 
 



  /* compute Kbond and Kcont*/
  PetscScalar etamin=1e99;  
  /* calculate global minimum effective viscosity */
  {
    PetscInt i,j;
    if(x<1) x=1;
    if(y<1) y=1;
    if(x+m>NX-2)m=NX-2-x;
    if(y+n>NY-2)n=NY-2-y;
    for(j=y;j<y+n;j++){
      for(i=x;i<x+m;i++){
	PetscScalar a1=Nb[j][i].T11;
	PetscScalar a2=Nb[j][i].T12;
	PetscScalar a3=Nb[j][i].T21;
	PetscScalar a4=Nb[j][i].T22;
	PetscScalar eta = sqrt(a1*a1+a2*a2+a3*a3+a4*a4);
	if(eta<etamin) etamin=eta;
      }
    }
  }
  {
    PetscScalar sbuf=etamin;
    ierr=MPI_Allreduce(&sbuf,&etamin,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
  }

  
  Kbond[0] = 4*etamin/((dx+dy)*(dx+dy));
  Kcont[0] = 2.0*etamin/(dx+dy);
  printf("found etamin = %e\n, Kcont = %e\n",etamin,Kcont[0]);

  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);/* restore correct corners*/

  /* now we want to assemble the global LHS, which requires knowing the mapping from local to global indices*/
  PetscInt nlocal;
  PetscInt *globalIdx;
  ierr=DMDAGetGlobalIndices(grid->da,&nlocal,&globalIdx);CHKERRQ(ierr);
  //  ierr=DAGetGhostCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);

  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  PetscInt xg,yg,mg,ng;
  ierr=DMDAGetGhostCorners(grid->da,&xg,&yg,PETSC_NULL,&mg,&ng,PETSC_NULL);CHKERRQ(ierr);

  /*make 2d arrays to easily reference local vx,vy,p dofs*/
  PetscInt **vxdof;
  PetscInt **vydof;
  PetscInt **pdof;

  PetscInt ix,jy;
  PetscMalloc( (ng+1)*sizeof(PetscInt *), &vxdof);
  PetscMalloc( (ng+1)*sizeof(PetscInt *), &vydof);
  PetscMalloc( ng*sizeof(PetscInt *), &pdof);
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
	pdof[jy][ix]  = 3*globalIdx[ix+mg*jy]+2;
	vxdof[jy][ix] = 3*globalIdx[ix+mg*jy]+0;
	vydof[jy][ix] = 3*globalIdx[ix+mg*jy]+1;
      //vzdof[jy][ix] = globalIdx[ix+mg*jy];/* vz goes into a separate matrix*/
      }
    }
  }
  /* print out the local vxdofs*/

  for(jy=y;jy<y+n;jy++){
    PetscInt jyl=jy-yg;
    for(ix=x;ix<x+m;ix++){
      PetscInt ixl=ix-xg;
      /* arrays to hold values to add.*/
      PetscScalar vals[NADD];
      PetscInt rowidx;
      PetscInt colidx[NADD];
      
      /**              %continuity*/		       
      if( (grid->xperiodic && jy == 0) || (!grid->xperiodic && jy == 0 && ix == 0) ){
	/* %pressure ghost nodes. Set equal to zero. */
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES); CHKERRQ(ierr);
	
      } else if( 0 && ((jy==1 && ix == 1) || (jy==NY-1 && ix==1))){
	//upper and lower left corner - dp/dx = 0;	//LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl], 1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl+1],-1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES); CHKERRQ(ierr);
	
      } else if( 0 && ((jy==1 && ix==NX-1) || (jy==NY-1 && ix== NX-1))) {
	//upper and lower right corners - dp/dx = 0;

	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);

	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl-1],-1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);

	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES); CHKERRQ(ierr);
	
      } else if( 0 && jy==1 && ix==2 ){
	/* specify pressure = 0.0 in one cell */
	ierr = MatSetValue(LHS,pdof[jyl][ixl],pdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof[jyl][ixl],0.0,INSERT_VALUES);CHKERRQ(ierr);
      } else {
	dx=grid->x[ix]-grid->x[ix-1];
	dy=grid->y[jy]-grid->y[jy-1];

	ierr = MatSetValue(LHS,pdof[jyl][ixl],vxdof[jyl-1][ixl],Kcont[0]/dx,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],vxdof[jyl-1][ixl-1],-Kcont[0]/dx,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],vydof[jyl][ixl-1],Kcont[0]/dy,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValue(LHS,pdof[jyl][ixl],vydof[jyl-1][ixl-1],-Kcont[0]/dy,INSERT_VALUES);CHKERRQ(ierr);
	//      end

	/* this term adds volumetric changes due to the equation of state (phase changes)*/
	/* this is not the same thing as compressibility - that would require additional terms in the momentum equations too */
	/* 	PetscScalar rhobar = (rho[jy-1][ix-1] + rho[jy-1][ix] + rho[jy][ix-1] + rho[jy][ix])*0.25; */
	/* 	PetscScalar rhodotbar =  (rhodot[jy-1][ix-1] + rhodot[jy-1][ix] + rhodot[jy][ix-1] + rhodot[jy][ix])*0.25;  */
	/* 	PetscScalar Rval = -1.0/rhobar*rhodotbar*Kcont[0]; */ 
	PetscScalar Rval = 0.0; /* incompressible */
	ierr = VecSetValue(RHS,pdof[jyl][ixl],Rval,INSERT_VALUES);CHKERRQ(ierr);
      }/* end continuity equation*/
      
      
      /*         %x-momentum */
      /* TOP BOUNDARY IS TREATED IMPLICITLY, ALL OTHERS ARE EXPLICIT */
      if(jy == NY-1){ /* BOTTOM BOUNDARY */
	/* Last row of nodes are used explicitly for boundary conditions. Note that the boundary conditions are treated implicitly for the top and left but explicitly for the right and bottom*/
	if( bv->mechBCBottom.type[0] == 0 ){
	  /* prescribed velocity vx[j,i] + vx[j-1,i] = 2.0*vb */
	  ierr =  MatSetValue(LHS,vxdof[jyl][ixl],vxdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr =  MatSetValue(LHS,vxdof[jyl][ixl],vxdof[jyl-1][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecSetValue(RHS, vxdof[jyl][ixl],2.0*Kbond[0]*bv->mechBCBottom.value[0],INSERT_VALUES);CHKERRQ(ierr);
	}else if( bv->mechBCBottom.type[0] == 1){
	  /* prescribed velocity gradient */
	  /* (vx[j,i] - vx[j-1,i])/dy = bcv */
	  PetscScalar dy1 = 2.0*(grid->LY - grid->yc[NY-1]);
	  ierr =  MatSetValue(LHS,vxdof[jyl][ixl],vxdof[jyl][ixl],1.0*Kbond[0]/dy1,INSERT_VALUES);CHKERRQ(ierr);
	  ierr =  MatSetValue(LHS,vxdof[jyl][ixl],vxdof[jyl-1][ixl],-1.0*Kbond[0]/dy1,INSERT_VALUES);CHKERRQ(ierr);
	  ierr =  VecSetValue(RHS, vxdof[jyl][ixl], Kbond[0]*bv->mechBCBottom.value[0],INSERT_VALUES);CHKERRQ(ierr);
	}
	
      } else if( !grid->xperiodic && (ix==0 && jy < NY-1)){ /* LEFT BOUNDARY */
	ierr = MatSetValue(LHS,vxdof[jyl][ixl],  vxdof[jyl][ixl],-1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vxdof[jyl][ixl], Kbond[0]*bv->mechBCLeft.value[0],INSERT_VALUES);CHKERRQ(ierr);
	
      } else if( !grid->xperiodic && (ix == NX-1 && jy < NY-1)){/* RIGHT BOUNDARY */
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vxdof[jyl][ixl],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vxdof[jyl][ixl], Kbond[0]*bv->mechBCRight.value[0],INSERT_VALUES);CHKERRQ(ierr);
         
      } else if( grid->xperiodic && bv->mechBCTop.type[0] == 1 && bv->mechBCBottom.type[0] == 1 &&  ix == 0 && jy == 0){
	/* left top boundary - fix x-velocity at one point */
	/* this is only needed when the domain is periodic and free slip boundary conditions */
	ierr = MatSetValue(LHS, vxdof[jyl][ixl], vxdof[jyl][ixl], 1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vxdof[jyl][ixl], 0.0 ,INSERT_VALUES); CHKERRQ(ierr);	

      } else {
	/* Normal anisotropic stencil - generated in Mathematica*/

	rowidx = vxdof[jyl][ixl];
	/* BEGIN MACHINE-GENERATED CODE */
	PetscScalar v0000 = grid->xc[ix];
	PetscScalar v0001 = Kcont[0];
	PetscScalar v0002 = grid->x[ix];
	PetscScalar v0003 = grid->y[jy];
	PetscScalar v0004 = grid->yc[jy];
	PetscInt v0005 = 1 + ix;
	PetscInt v0006 = 1 + jy;
	PetscInt v0007 = 2 + jy;
	PetscInt v0008 = -1 + ix;
	PetscScalar v0009 = Nb[jy][ix].T22;
	PetscScalar v0010 = jy == 0 ? Nc[1][ix].T21 : Nc[jy][ix].T21;
	PetscScalar v0011 = grid->xc[v0005];
	PetscScalar v0012 = grid->yc[v0007];
	PetscScalar v0013 = grid->yc[v0006];
	PetscScalar v0014 = grid->y[v0006];
	PetscScalar v0015 = -v0000;
	PetscScalar v0016 = grid->x[v0005];
	PetscScalar v0017 = -v0004;
	PetscScalar v0018 = -v0002;
	PetscScalar v0019 = grid->x[v0008];
	PetscScalar v0020 = Nb[jy][v0005].T12;
	PetscScalar v0021 = jy == 0 ? Nc[1][v0005].T21 : Nc[jy][v0005].T21;
	PetscScalar v0022 = Nb[v0006][ix].T22;
	PetscScalar v0023 = v0007 > NY-1 ? Nc[jy+1][ix].T21 : Nc[v0007][ix].T21;
	PetscScalar v0024 = Nc[v0006][ix].T21;
	PetscScalar v0025 = Nc[v0006][ix].T11;
	PetscScalar v0026 = Nb[jy][v0008].T12;
	PetscScalar v0027 = -v0011;
	PetscScalar v0028 = -v0013;
	PetscScalar v0029 = Nb[v0006][v0005].T12;
	PetscScalar v0030 = v0002 + v0015;
	PetscScalar v0031 = v0007 > NY-1 ? Nc[jy+1][ix].T21 : Nc[v0007][v0005].T21;
	PetscScalar v0032 = Nc[v0006][v0005].T21;
	PetscScalar v0033 = Nc[v0006][v0005].T11;
	PetscScalar v0034 = Nb[v0006][v0008].T12;
	PetscScalar v0035 = v0000 + v0027;
	PetscScalar v0036 = v0003 - v0014;
	PetscScalar v0037 = v0002 - v0016;
	PetscScalar v0038 = v0004 + v0028;
	PetscScalar v0039 = v0003 + v0028;
	PetscScalar v0040 = v0002 + v0027;
	PetscScalar v0041 = v0018 + v0019;
	PetscScalar v0042 = v0013 + v0017;
	PetscScalar v0043 = 1/v0035;
	PetscScalar v0044 = -v0012 + v0013;
	PetscScalar v0045 = v0014 + v0028;
	PetscScalar v0046 = 1/v0036;
	PetscScalar v0047 = 1/v0037;
	PetscScalar v0048 = 1/v0038;
	PetscScalar v0049 = 1/(v0011 + v0015);
	PetscScalar v0050 = 1/(-v0003 + v0014);
	PetscScalar v0051 = 1/v0041;
	PetscScalar v0052 = 1/(v0002 - v0019);
	PetscScalar v0053 = 1/(v0015 + grid->xc[v0008]);
	PetscScalar v0054 = 1/(v0011 - grid->xc[2 + ix]);
	PetscScalar v0055 = 1/v0044;
	PetscScalar v0056 = v0013*(v0012 + v0017) + v0014*v0038 + \
	  v0003*v0044;
	colidx[0] = -1;
	colidx[1] = -1;
	colidx[2] = -1;
	vals[0]=0;
	vals[1]=0;
	vals[2]=0;
	colidx[3] = -1;
	colidx[4] = -1;
	colidx[5] = -1;
	vals[3]=0;
	vals[4]=0;
	vals[5]=0;
	colidx[6] = -1;
	colidx[7] = vydof[jyl+0][ixl+-2];
	colidx[8] = -1;
	vals[6]=0;
	vals[7]=(v0026*v0043*v0053)/8.;
	vals[8]=0;
	colidx[9] = -1;
	colidx[10] = vydof[jyl+1][ixl+-2];
	colidx[11] = -1;
	vals[9]=0;
	vals[10]=(v0034*v0043*v0053)/8.;
	vals[11]=0;
	colidx[12] = -1;
	colidx[13] = -1;
	colidx[14] = -1;
	vals[12]=0;
	vals[13]=0;
	vals[14]=0;
	colidx[15] = -1;
	colidx[16] = -1;
	colidx[17] = -1;
	vals[15]=0;
	vals[16]=0;
	vals[17]=0;
	colidx[18] = (ix > 0 && jy > 0) ? vxdof[jyl+-1][ixl+-1] : -1;
	colidx[19] = -1;
	colidx[20] = -1;
	vals[18]=((8*v0010*v0039*v0040 + \
		   v0026*v0036*v0041)*v0043*v0046*v0048*v0051)/8.;
	vals[19]=0;
	vals[20]=0;
	colidx[21] = vxdof[jyl+0][ixl+-1];
	colidx[22] = vydof[jyl+0][ixl+-1];
	colidx[23] = -1;
	vals[21]=(v0043*v0046*v0048*v0051*v0055*(v0034*v0036*v0038*v0041 - \
						 v0036*(v0026*v0041 + 8*v0025*v0042)*v0044 - 8*v0024*v0040*v0056))/8.;
	vals[22]=(v0049*(4*v0009*v0050 + v0026*v0053))/8.;
	vals[23]=0;
	colidx[24] = vxdof[jyl+1][ixl+-1];
	colidx[25] = vydof[jyl+1][ixl+-1];
	colidx[26] = -1;
	vals[24]=(v0043*(-(v0034*v0036*v0041) + \
			 8*v0023*v0040*v0045)*v0046*v0051*v0055)/8.;
	vals[25]=(v0049*(4*v0022*v0046 + v0034*v0053))/8.;
	vals[26]=0;
	colidx[27] = -1;
	colidx[28] = -1;
	colidx[29] = -1;
	vals[27]=0;
	vals[28]=0;
	vals[29]=0;
	colidx[30] = -1;
	colidx[31] = -1;
	colidx[32] = -1;
	vals[30]=0;
	vals[31]=0;
	vals[32]=0;
	colidx[33] = jy > 0 ? vxdof[jyl+-1][ixl+0] : -1;
	colidx[34] = -1;
	colidx[35] = -1;
	vals[33]=((v0009*v0035*v0037*v0041 - 2*v0039*(v0010*v0037*v0040 + \
						      v0021*v0030*v0041))*v0043*v0046*v0047*v0048*v0051)/2.;
	vals[34]=0;
	vals[35]=0;
	colidx[36] = vxdof[jyl+0][ixl+0];
	colidx[37] = vydof[jyl+0][ixl+0];
	colidx[38] = -1;
	vals[36]=-(v0033*v0043*v0047) + v0002*v0032*v0043*v0046*v0047 - \
	  (v0009*v0046*v0048)/2. + v0002*v0003*v0032*v0043*v0046*v0047*v0048 + \
	  v0000*v0004*v0032*v0043*v0046*v0047*v0048 + \
	  v0003*v0015*v0032*v0043*v0046*v0047*v0048 + \
	  v0002*v0017*v0032*v0043*v0046*v0047*v0048 + \
	  (v0000*v0032*v0049*v0050)/(v0016 + v0018) - v0025*v0043*v0051 + \
	  v0024*v0046*v0051 + v0002*v0024*v0043*v0046*v0051 + \
	  v0003*v0024*v0046*v0048*v0051 + \
	  v0002*v0003*v0024*v0043*v0046*v0048*v0051 + \
	  v0000*v0004*v0024*v0043*v0046*v0048*v0051 + \
	  v0003*v0015*v0024*v0043*v0046*v0048*v0051 + \
	  v0002*v0017*v0024*v0043*v0046*v0048*v0051 + \
	  (v0013*v0024*v0050*v0052)/(v0012 + v0028) + \
	  (v0004*v0024*v0050*v0052)/v0042 + v0000*v0024*v0049*v0050*v0052 - \
	  (v0022*v0046*v0055)/2. + \
	  (v0013*v0018*v0024*v0043*v0046*v0055)/(-v0002 + v0019) + \
	  v0000*v0013*v0032*v0043*v0046*v0047*v0055 + \
	  v0002*v0014*v0032*v0043*v0046*v0047*v0055 + \
	  v0014*v0015*v0032*v0043*v0046*v0047*v0055 + \
	  v0013*v0018*v0032*v0043*v0046*v0047*v0055 + \
	  v0014*v0024*v0046*v0051*v0055 + \
	  v0000*v0013*v0024*v0043*v0046*v0051*v0055 + \
	  v0002*v0014*v0024*v0043*v0046*v0051*v0055 + \
	  v0014*v0015*v0024*v0043*v0046*v0051*v0055;
	vals[37]=(v0049*(4*v0009*v0046 + v0020*v0054))/8.;
	vals[38]=0;
	colidx[39] = vxdof[jyl+1][ixl+0];
	colidx[40] = vydof[jyl+1][ixl+0];
	colidx[41] = pdof[jyl+1][ixl+0];
	vals[39]=(v0043*(v0022*v0035*v0037*v0041 - 2*(v0023*v0037*v0040 + \
						      v0030*v0031*v0041)*v0045)*v0046*v0047*v0051*v0055)/2.;
	vals[40]=(v0049*(4*v0022*v0050 + v0029*v0054))/8.;
	vals[41]=v0001*v0049;
	colidx[42] = -1;
	colidx[43] = -1;
	colidx[44] = -1;
	vals[42]=0;
	vals[43]=0;
	vals[44]=0;
	colidx[45] = -1;
	colidx[46] = -1;
	colidx[47] = -1;
	vals[45]=0;
	vals[46]=0;
	vals[47]=0;
	colidx[48] = (jy > 0 && ix > 0) ? vxdof[jyl+-1][ixl+1] : -1;
	colidx[49] = -1;
	colidx[50] = -1;
	vals[48]=((-(v0020*v0036*v0037) + \
		   8*v0021*v0030*v0039)*v0043*v0046*v0047*v0048)/8.;
	vals[49]=0;
	vals[50]=0;
	colidx[51] = vxdof[jyl+0][ixl+1];
	colidx[52] = vydof[jyl+0][ixl+1];
	colidx[53] = -1;
	vals[51]=(v0043*v0046*v0047*v0048*v0055*(-(v0029*v0036*v0037*v0038) + \
						 v0036*(v0020*v0037 + 8*v0033*v0038)*v0044 - 8*v0030*v0032*v0056))/8.;
	vals[52]=(v0020*v0043*v0054)/8.;
	vals[53]=0;
	colidx[54] = vxdof[jyl+1][ixl+1];
	colidx[55] = vydof[jyl+1][ixl+1];
	colidx[56] = pdof[jyl+1][ixl+1];
	vals[54]=(v0043*(v0029*v0036*v0037 + \
			 8*v0030*v0031*v0045)*v0046*v0047*v0055)/8.;
	vals[55]=(v0029*v0043*v0054)/8.;
	vals[56]=v0001*v0043;
	colidx[57] = -1;
	colidx[58] = -1;
	colidx[59] = -1;
	vals[57]=0;
	vals[58]=0;
	vals[59]=0;
	colidx[60] = -1;
	colidx[61] = -1;
	colidx[62] = -1;
	vals[60]=0;
	vals[61]=0;
	vals[62]=0;
	colidx[63] = -1;
	colidx[64] = -1;
	colidx[65] = -1;
	vals[63]=0;
	vals[64]=0;
	vals[65]=0;
	colidx[66] = -1;
	colidx[67] = -1;
	colidx[68] = -1;
	vals[66]=0;
	vals[67]=0;
	vals[68]=0;
	colidx[69] = -1;
	colidx[70] = -1;
	colidx[71] = -1;
	vals[69]=0;
	vals[70]=0;
	vals[71]=0;
	colidx[72] = -1;
	colidx[73] = -1;
	colidx[74] = -1;
	vals[72]=0;
	vals[73]=0;
	vals[74]=0;
	/* END MACHINE GENERATED CODE */

	/* Some coefficients may correspond to terms that are out of domain. Zero these values */
	{
	  PetscInt ii;
	  for(ii = 0; ii<75; ii++){
	    if(colidx[ii] >= NDOF ) colidx[ii] = -1;/* do not add this value to the matrix */
	  }
	}
	PetscScalar Rval = -gx*(rho[jy][ix]+rho[jy+1][ix])/2.0;/* -(soxxZ[jy+1][ix+1]-soxxZ[jy+1][ix])/dxc-(soxyZ[jy+1][ix]-soxyZ[jy][ix])/dy;*/
	ierr = VecSetValue( RHS, vxdof[jyl][ixl],Rval,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues( LHS, 1,&rowidx,75,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
      }


      /*         %y-momentum */
      /* TOP, BOTTOM, RIGHT ARE TREATED EXPLICITLY, LEFT IS IMPLICIT */
      if(  jy == 0 && (grid->xperiodic || ix < NX-1) ) {/* TOP */
	/* fixed velocity */
	ierr = MatSetValue(LHS, vydof[jyl][ixl],vydof[jyl][ixl], Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vydof[jyl][ixl], Kbond[0] * inflow,INSERT_VALUES);CHKERRQ(ierr);
      }else if(  jy == NY-1 && (grid->xperiodic ||  ix < NX-1) ){/* BOTTOM */
	/* vy = 0; */
	ierr=MatSetValue(LHS,vydof[jyl][ixl],vydof[jyl][ixl],Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecSetValue(RHS, vydof[jyl][ixl], 0.0 ,INSERT_VALUES);CHKERRQ(ierr);
      }else if( !grid->xperiodic && ix == NX-1 ){/* RIGHT */
	if( bv->mechBCRight.type[1] == 0){/* prescribed velocity case */
	  ierr=MatSetValue(LHS,vydof[jyl][ixl], vydof[jyl][ixl],  1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vydof[jyl][ixl], vydof[jyl][ixl-1],1.0*Kbond[0],INSERT_VALUES);CHKERRQ(ierr);
	  ierr=VecSetValue(RHS,vydof[jyl][ixl], 2.0*Kbond[0]*bv->mechBCRight.value[1], INSERT_VALUES); CHKERRQ(ierr);	  
	}else if(bv->mechBCRight.type[1] ==1){ /* free slip case */
	  ierr=MatSetValue(LHS,vydof[jyl][ixl],vydof[jyl][ixl],   1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	  ierr=MatSetValue(LHS,vydof[jyl][ixl],vydof[jyl][ixl-1],-1.0*Kbond[0],INSERT_VALUES); CHKERRQ(ierr);
	  ierr=VecSetValue(RHS,vydof[jyl][ixl], 0.0, INSERT_VALUES); CHKERRQ(ierr);
	}

      } else {
	/* anisotropic stokes stencil */
	rowidx = vydof[jyl][ixl];
	/* BEGIN MACHINE-GENERATED CODE */
	PetscScalar v0000 = grid->yc[jy];
	PetscScalar v0001 = Kcont[0];
	PetscScalar v0002 = grid->y[jy];
	PetscScalar v0003 = grid->x[ix];
	PetscScalar v0004 = grid->xc[ix];
	PetscInt v0005 = 1 + jy;
	PetscInt v0006 = 1 + ix;
	PetscInt v0007 = 2 + ix;
	PetscInt v0008 = -1 + jy;
	PetscScalar v0009 = Nb[jy][ix].T22;
	PetscScalar v0010 = ix == 0 ? Nc[jy][1].T21 : Nc[jy][ix].T21;
	PetscScalar v0011 = grid->yc[v0005];
	PetscScalar v0012 = -v0000;
	PetscScalar v0013 = grid->xc[v0006];
	PetscScalar v0014 = grid->x[v0006];
	PetscScalar v0015 = grid->y[v0005];
	PetscScalar v0016 = grid->xc[v0007];
	PetscScalar v0017 = grid->yc[2 + jy];
	PetscScalar v0018 = -v0002;
	PetscScalar v0019 = grid->y[v0008];
	PetscScalar v0020 = -v0003;
	PetscScalar v0021 = grid->yc[v0008];
	PetscScalar v0022 = Nb[jy][v0006].T22;
	PetscScalar v0023 = ix == NX-2 ? Nc[jy][ix+1].T21 : Nc[jy][v0007].T21;
	PetscScalar v0024 = Nb[v0005][ix].T12;
	PetscScalar v0025 = Nc[jy][v0006].T21;
	PetscScalar v0026 = Nc[jy][v0006].T11;
	PetscScalar v0027 = Nb[v0008][ix].T12;
	PetscScalar v0028 = ix == 0 ? Nc[v0005][ix+1].T21 : Nc[v0005][ix].T21;
	PetscScalar v0029 = -v0011;
	PetscScalar v0030 = -v0013;
	PetscScalar v0031 = v0002 + v0012;
	PetscScalar v0032 = ix == NX-2 ? Nc[v0005][ix+1].T21 : Nc[v0005][v0007].T21;
	PetscScalar v0033 = Nb[v0005][v0006].T12;
	PetscScalar v0034 = Nb[v0008][v0006].T12;
	PetscScalar v0035 = Nc[v0005][v0006].T21;
	PetscScalar v0036 = Nc[v0005][v0006].T11;
	PetscScalar v0037 = v0000 + v0029;
	PetscScalar v0038 = v0002 - v0015;
	PetscScalar v0039 = v0003 - v0014;
	PetscScalar v0040 = v0018 + v0019;
	PetscScalar v0041 = v0002 + v0029;
	PetscScalar v0042 = v0004 + v0030;
	PetscScalar v0043 = v0003 + v0030;
	PetscScalar v0044 = 1/v0037;
	PetscScalar v0045 = v0014 + v0030;
	PetscScalar v0046 = 1/v0038;
	PetscScalar v0047 = v0013 - v0016;
	PetscScalar v0048 = 1/v0039;
	PetscScalar v0049 = 1/(v0011 + v0012);
	PetscScalar v0050 = 1/v0040;
	PetscScalar v0051 = 1/(v0014 + v0020);
	PetscScalar v0052 = 1/(v0000 - v0021);
	PetscScalar v0053 = 1/(v0012 + v0021);
	PetscScalar v0054 = 1/v0042;
	PetscScalar v0055 = 1/(v0002 - v0019);
	PetscScalar v0056 = 1/v0047;
	PetscScalar v0057 = 1/(v0011 - v0017);
	PetscScalar v0058 = 1/(v0017 + v0029);
	colidx[0] = -1;
	colidx[1] = -1;
	colidx[2] = -1;
	vals[0]=0;
	vals[1]=0;
	vals[2]=0;
	colidx[3] = -1;
	colidx[4] = -1;
	colidx[5] = -1;
	vals[3]=0;
	vals[4]=0;
	vals[5]=0;
	colidx[6] = -1;
	colidx[7] = -1;
	colidx[8] = -1;
	vals[6]=0;
	vals[7]=0;
	vals[8]=0;
	colidx[9] = -1;
	colidx[10] = -1;
	colidx[11] = -1;
	vals[9]=0;
	vals[10]=0;
	vals[11]=0;
	colidx[12] = -1;
	colidx[13] = -1;
	colidx[14] = -1;
	vals[12]=0;
	vals[13]=0;
	vals[14]=0;
	colidx[15] = -1;
	colidx[16] = -1;
	colidx[17] = -1;
	vals[15]=0;
	vals[16]=0;
	vals[17]=0;
	colidx[18] = -1;
	colidx[19] = vydof[jyl+-1][ixl+-1];
	colidx[20] = -1;
	vals[18]=0;
	vals[19]=-((v0027*v0039*v0040 + \
		    8*v0010*v0041*v0043)*v0044*v0048*v0050*v0054)/8.;
	vals[20]=0;
	colidx[21] = -1;
	colidx[22] = vydof[jyl+0][ixl+-1];
	colidx[23] = -1;
	vals[21]=0;
	vals[22]=((2*v0028*v0031*v0040*v0043 + v0038*(v0009*v0037*v0040 + \
						      2*v0010*v0041*v0043))*v0044*v0046*v0048*v0050*v0054)/2.;
	vals[23]=0;
	colidx[24] = -1;
	colidx[25] = vydof[jyl+1][ixl+-1];
	colidx[26] = -1;
	vals[24]=0;
	vals[25]=((v0024*v0038*v0039 - \
		   8*v0028*v0031*v0043)*v0044*v0046*v0048*v0054)/8.;
	vals[26]=0;
	colidx[27] = -1;
	colidx[28] = -1;
	colidx[29] = -1;
	vals[27]=0;
	vals[28]=0;
	vals[29]=0;
	colidx[30] = (jy > 1) ? vxdof[jyl+-2][ixl+0] : -1;
	colidx[31] = -1;
	colidx[32] = -1;
	vals[30]=-(v0027*v0044*v0053)/8.;
	vals[31]=0;
	vals[32]=0;
	colidx[33] = vxdof[jyl+-1][ixl+0];
	colidx[34] = vydof[jyl+-1][ixl+0];
	colidx[35] = -1;
	vals[33]=(v0049*(4*v0009*v0051 + v0027*v0052))/8.;
	vals[34]=(v0044*(8*v0002*v0003*v0013*v0025 - \
			 8*v0002*v0004*v0013*v0025 - 8*v0003*v0011*v0013*v0025 + \
			 8*v0004*v0011*v0013*v0025 + 8*v0002*v0004*v0014*v0025 - \
			 8*v0004*v0011*v0014*v0025 - 8*v0002*v0013*v0014*v0025 + \
			 8*v0011*v0013*v0014*v0025 - 8*v0002*v0003*v0016*v0025 + \
			 8*v0003*v0011*v0016*v0025 + 8*v0002*v0013*v0016*v0025 - \
			 8*v0011*v0013*v0016*v0025 + v0002*v0013*v0014*v0027 + \
			 v0002*v0003*v0016*v0027 + v0003*v0013*v0018*v0027 + \
			 v0014*v0016*v0018*v0027 + v0003*v0013*v0019*v0027 + \
			 v0014*v0016*v0019*v0027 + v0016*v0019*v0020*v0027 + \
			 v0014*v0019*v0027*v0030 - v0034*v0039*v0040*v0042 + \
			 8*v0026*v0039*v0042*v0047)*v0048*v0050*v0054*v0056)/8.;
	vals[35]=0;
	colidx[36] = vxdof[jyl+0][ixl+0];
	colidx[37] = vydof[jyl+0][ixl+0];
	colidx[38] = -1;
	vals[36]=(v0049*(4*v0009*v0048 + v0024*v0058))/8.;
	vals[37]=-(v0036*v0044*v0046) + v0000*v0035*v0044*v0046*v0048 - \
	  v0026*v0044*v0050 - v0025*v0048*v0050 + v0000*v0025*v0044*v0048*v0050 \
	  + (v0002*v0035*v0049*v0051)/(v0015 + v0018) - (v0009*v0048*v0054)/2. \
	  + (v0003*v0018*v0025*v0044*v0048*v0054)/(-v0002 + v0019) + \
	  v0000*v0003*v0035*v0044*v0046*v0048*v0054 + \
	  v0002*v0004*v0035*v0044*v0046*v0048*v0054 + \
	  v0004*v0012*v0035*v0044*v0046*v0048*v0054 + \
	  v0003*v0018*v0035*v0044*v0046*v0048*v0054 + \
	  v0004*v0025*v0048*v0050*v0054 + \
	  v0000*v0003*v0025*v0044*v0048*v0050*v0054 + \
	  v0002*v0004*v0025*v0044*v0048*v0050*v0054 + \
	  v0004*v0012*v0025*v0044*v0048*v0050*v0054 + \
	  (v0003*v0025*v0051*v0055)/(-v0004 + v0013) + \
	  (v0014*v0025*v0051*v0055)/(v0016 + v0030) + \
	  v0002*v0025*v0049*v0051*v0055 - (v0022*v0048*v0056)/2. + \
	  (v0014*v0018*v0025*v0044*v0048*v0056)/(-v0002 + v0019) + \
	  v0002*v0013*v0035*v0044*v0046*v0048*v0056 + \
	  v0012*v0013*v0035*v0044*v0046*v0048*v0056 + \
	  v0000*v0014*v0035*v0044*v0046*v0048*v0056 + \
	  v0014*v0018*v0035*v0044*v0046*v0048*v0056 + \
	  v0013*v0025*v0048*v0050*v0056 + \
	  v0002*v0013*v0025*v0044*v0048*v0050*v0056 + \
	  v0012*v0013*v0025*v0044*v0048*v0050*v0056 + \
	  v0000*v0014*v0025*v0044*v0048*v0050*v0056;
	vals[38]=0;
	colidx[39] = vxdof[jyl+1][ixl+0];
	colidx[40] = vydof[jyl+1][ixl+0];
	colidx[41] = -1;
	vals[39]=-(v0024*v0044*v0057)/8.;
	vals[40]=(v0044*v0046*(v0002*v0013*v0014*v0024 + \
			       v0003*v0013*v0015*v0024 + v0002*v0003*v0016*v0024 + \
			       v0014*v0015*v0016*v0024 + v0003*v0013*v0018*v0024 + \
			       v0014*v0016*v0018*v0024 + v0015*v0016*v0020*v0024 + \
			       v0014*v0015*v0024*v0030 - 8*v0000*v0003*v0013*v0035 + \
			       8*v0002*v0003*v0013*v0035 + 8*v0000*v0004*v0013*v0035 - \
			       8*v0002*v0004*v0013*v0035 - 8*v0000*v0004*v0014*v0035 + \
			       8*v0002*v0004*v0014*v0035 + 8*v0000*v0013*v0014*v0035 - \
			       8*v0002*v0013*v0014*v0035 + 8*v0000*v0003*v0016*v0035 - \
			       8*v0002*v0003*v0016*v0035 - 8*v0000*v0013*v0016*v0035 + \
			       8*v0002*v0013*v0016*v0035 + v0033*v0038*v0039*v0042 + \
			       8*v0036*v0039*v0042*v0047)*v0048*v0054*v0056)/8.;
	vals[41]=0;
	colidx[42] = -1;
	colidx[43] = -1;
	colidx[44] = -1;
	vals[42]=0;
	vals[43]=0;
	vals[44]=0;
	colidx[45] = jy > 1 ? vxdof[jyl+-2][ixl+1] : -1;
	colidx[46] = -1;
	colidx[47] = -1;
	vals[45]=-(v0034*v0044*v0053)/8.;
	vals[46]=0;
	vals[47]=0;
	colidx[48] = vxdof[jyl+-1][ixl+1];
	colidx[49] = vydof[jyl+-1][ixl+1];
	colidx[50] = -1;
	vals[48]=(v0049*(4*v0022*v0048 + v0034*v0052))/8.;
	vals[49]=(v0044*(v0034*v0039*v0040 - \
			 8*v0023*v0041*v0045)*v0048*v0050*v0056)/8.;
	vals[50]=0;
	colidx[51] = vxdof[jyl+0][ixl+1];
	colidx[52] = vydof[jyl+0][ixl+1];
	colidx[53] = pdof[jyl+0][ixl+1];
	vals[51]=(v0049*(4*v0022*v0051 + v0033*v0058))/8.;
	vals[52]=(v0044*(2*v0031*v0032*v0040*v0045 + v0038*(v0022*v0037*v0040 \
							    + 2*v0023*v0041*v0045))*v0046*v0048*v0050*v0056)/2.;
	vals[53]=v0001*v0049;
	colidx[54] = vxdof[jyl+1][ixl+1];
	colidx[55] = vydof[jyl+1][ixl+1];
	colidx[56] = pdof[jyl+1][ixl+1];
	vals[54]=-(v0033*v0044*v0057)/8.;
	vals[55]=-(v0044*(v0033*v0038*v0039 + \
			  8*v0031*v0032*v0045)*v0046*v0048*v0056)/8.;
	vals[56]=v0001*v0044;
	colidx[57] = -1;
	colidx[58] = -1;
	colidx[59] = -1;
	vals[57]=0;
	vals[58]=0;
	vals[59]=0;
	colidx[60] = -1;
	colidx[61] = -1;
	colidx[62] = -1;
	vals[60]=0;
	vals[61]=0;
	vals[62]=0;
	colidx[63] = -1;
	colidx[64] = -1;
	colidx[65] = -1;
	vals[63]=0;
	vals[64]=0;
	vals[65]=0;
	colidx[66] = -1;
	colidx[67] = -1;
	colidx[68] = -1;
	vals[66]=0;
	vals[67]=0;
	vals[68]=0;
	colidx[69] = -1;
	colidx[70] = -1;
	colidx[71] = -1;
	vals[69]=0;
	vals[70]=0;
	vals[71]=0;
	colidx[72] = -1;
	colidx[73] = -1;
	colidx[74] = -1;
	vals[72]=0;
	vals[73]=0;
	vals[74]=0;

	/* END MACHINE-GENERATED CODE */

	/* Some coefficients may correspond to terms that are out of domain. Zero these values */
	{
	  PetscInt ii;
	  for(ii = 0; ii<75; ii++){
	    if(colidx[ii] >= NDOF ) colidx[ii] = -1;/* do not add this value to the matrix */
	  }
	}

	/* code generated by Mathematica */
	PetscScalar Rval = -gy*(rho[jy][ix]+rho[jy][ix+1])/2.0;
	  
	ierr = VecSetValue( RHS, vydof[jyl][ixl],Rval,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues( LHS, 1,&rowidx,75,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
      }
	  
    }/* end loop over x grid-line*/
	
  }/* end loop over y grid-line*/
    
  /* flush LHS and RHS before adding implicit boundary condition values */    
  ierr = VecAssemblyBegin(RHS);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(RHS);CHKERRQ(ierr);
  ierr = MatAssemblyBegin( LHS, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( LHS, MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);


  /* Implicit boundary condition values */
  for(jy=y;jy<y+n;jy++){
    PetscInt jyl=jy-yg;
    for(ix=x;ix<x+m;ix++){
      PetscInt ixl=ix-xg;
      /* arrays to hold values to add.*/
      PetscScalar vals[NADD];

      /* X Stokes */
      /* Left Boundary */
      if( !grid->xperiodic && ix == 1 && jy < NY-1 ){
	/* put in ghost values for vy[i-2][j] and vy[i-2][j+1] */
	/* free slip - vy[i-2][j] = vy[i-1][j] */
	PetscScalar vyim2jcoef = Nb[jy][-1 + ix].T12/(8.*(-grid->xc[-1 + ix] + grid->xc[ix])*(-grid->xc[ix] + grid->xc[1 + ix]));
	PetscScalar vyim2jp1coef = Nb[1 + jy][-1 + ix].T12/(8.*(-grid->xc[-1 + ix] + grid->xc[ix])*(-grid->xc[ix] + grid->xc[1 + ix]));
	if( bv->mechBCLeft.type[1] == 1){
	  ierr = MatSetValue( LHS, vxdof[jyl][ixl], vydof[jyl][ixl-1], vyim2jcoef, ADD_VALUES );CHKERRQ(ierr);
	  ierr = MatSetValue( LHS, vxdof[jyl][ixl], vydof[jyl+1][ixl-1], vyim2jp1coef, ADD_VALUES );CHKERRQ(ierr);
	}
	/* if j == 0, also put in vx[i-1][j-1], vx[i,j-1], vx[i+1,j-1] */
	if( jy ==0 ){
	  PetscScalar vxim1jm1coef = (Nb[jy][-1 + ix].T12*(grid->x[-1 + ix] - grid->x[ix])*(grid->y[jy] - \
				     grid->y[1 + jy]) + 8*Nc[jy+1][ix].T21*(grid->x[ix] - grid->xc[1 + \
					ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[-1 + ix] - \
						grid->x[ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - \
						grid->y[1 + jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  PetscScalar vxijm1coef = (Nb[jy][ix].T22*(grid->x[-1 + ix] - grid->x[ix])*(grid->x[ix] - \
grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 + ix]) - 2*(Nc[jy+1][1 + \
ix].T21*(grid->x[-1 + ix] - grid->x[ix])*(grid->x[ix] - grid->xc[ix]) \
+ Nc[jy+1][ix].T21*(grid->x[ix] - grid->x[1 + ix])*(grid->x[ix] - \
grid->xc[1 + ix]))*(grid->y[jy] - grid->yc[1 + jy]))/(2.*(grid->x[-1 \
+ ix] - grid->x[ix])*(grid->x[ix] - grid->x[1 + ix])*(grid->xc[ix] - \
grid->xc[1 + ix])*(grid->y[jy] - grid->y[1 + jy])*(grid->yc[jy] - \
						   grid->yc[1 + jy]));
	  PetscScalar vxip1jm1coef = (-(Nb[jy][1 + ix].T12*(grid->x[ix] - grid->x[1 + ix])*(grid->y[jy] - \
				     grid->y[1 + jy])) + 8*Nc[jy+1][1 + ix].T21*(grid->x[ix] - \
				    grid->xc[ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[ix] - \
				    grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - \
				    grid->y[1 + jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  if( bv->mechBCTop.type[0] == 1){
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl-1], vxim1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl], vxijm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl+1], vxip1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	  }	  
	}/* end special case for top left */
      }

      /* Right Boundary*/
      if( (!grid->xperiodic && ix == NX-2 && jy == 0)){
	/* put in ghost values only if j ==0 : vx[i-1][j-1], vx[i][j-1], vx[i+1][j-1] */
	PetscScalar vxim1jm1coef = (Nb[jy][-1 + ix].T12*(grid->x[-1 + ix] - grid->x[ix])*(grid->y[jy] - \
				     grid->y[1 + jy]) + 8*Nc[jy+1][ix].T21*(grid->x[ix] - grid->xc[1 + \
					ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[-1 + ix] - \
					grid->x[ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - \
					grid->y[1 + jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  PetscScalar vxijm1coef = (Nb[jy][ix].T22*(grid->x[-1 + ix] - grid->x[ix])*(grid->x[ix] - \
grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 + ix]) - 2*(Nc[jy+1][1 + \
ix].T21*(grid->x[-1 + ix] - grid->x[ix])*(grid->x[ix] - grid->xc[ix]) \
+ Nc[jy+1][ix].T21*(grid->x[ix] - grid->x[1 + ix])*(grid->x[ix] - \
grid->xc[1 + ix]))*(grid->y[jy] - grid->yc[1 + jy]))/(2.*(grid->x[-1 \
+ ix] - grid->x[ix])*(grid->x[ix] - grid->x[1 + ix])*(grid->xc[ix] - \
grid->xc[1 + ix])*(grid->y[jy] - grid->y[1 + jy])*(grid->yc[jy] - \
						   grid->yc[1 + jy]));
	  PetscScalar vxip1jm1coef = (-(Nb[jy][1 + ix].T12*(grid->x[ix] - grid->x[1 + ix])*(grid->y[jy] - \
				     grid->y[1 + jy])) + 8*Nc[jy+1][1 + ix].T21*(grid->x[ix] - \
				    grid->xc[ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[ix] - \
				    grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - \
				    grid->y[1 + jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  if( bv->mechBCTop.type[0] == 1){
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl-1], vxim1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl], vxijm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl+1], vxip1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	  }
      }
      /* Top Boundary BUT NOT SIDES UNLESS XPERIODIC*/
      if( jy ==0 && (grid->xperiodic || (ix > 1 && ix < NX-2)) ){
	  PetscScalar vxim1jm1coef = (Nb[jy][-1 + ix].T12*(grid->x[-1 + ix] - grid->x[ix])*(grid->y[jy] - \
				     grid->y[1 + jy]) + 8*Nc[jy+1][ix].T21*(grid->x[ix] - grid->xc[1 + \
					ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[-1 + ix] - \
					grid->x[ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - \
					grid->y[1 + jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  PetscScalar vxijm1coef = (Nb[jy][ix].T22*(grid->x[-1 + ix] - grid->x[ix])*(grid->x[ix] - \
grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 + ix]) - 2*(Nc[jy+1][1 + \
ix].T21*(grid->x[-1 + ix] - grid->x[ix])*(grid->x[ix] - grid->xc[ix]) \
+ Nc[jy+1][ix].T21*(grid->x[ix] - grid->x[1 + ix])*(grid->x[ix] - \
grid->xc[1 + ix]))*(grid->y[jy] - grid->yc[1 + jy]))/(2.*(grid->x[-1 \
+ ix] - grid->x[ix])*(grid->x[ix] - grid->x[1 + ix])*(grid->xc[ix] - \
grid->xc[1 + ix])*(grid->y[jy] - grid->y[1 + jy])*(grid->yc[jy] - \
						   grid->yc[1 + jy]));
	  PetscScalar vxip1jm1coef = (-(Nb[jy][1 + ix].T12*(grid->x[ix] - grid->x[1 + ix])*(grid->y[jy] - \
				     grid->y[1 + jy])) + 8*Nc[jy+1][1 + ix].T21*(grid->x[ix] - \
				    grid->xc[ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[ix] - \
				    grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - \
				    grid->y[1 + jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  if( bv->mechBCTop.type[0] == 1){
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl-1], vxim1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl], vxijm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue( LHS, vxdof[jyl][ixl], vxdof[jyl][ixl+1], vxip1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	  }
      }
      /* Bottom Boundary - no changes needed! */
      

      /* Y Stokes */
      /* Left Boundary */
      if( !grid->xperiodic && (ix == 0 && jy > 0 && jy < NY-1)){
	/* add ghost values for vy[i-1][j-1], vy[i-1][j], vy[i-1][j+1] */
	PetscScalar vyim1jm1coef = -(Nb[-1 + jy][ix].T12*(grid->x[ix] - grid->x[1 + ix])*(grid->y[-1 + \
jy] - grid->y[jy]) + 8*Nc[jy][ix+1].T21*(grid->x[ix] - grid->xc[1 + \
ix])*(grid->y[jy] - grid->yc[1 + jy]))/(8.*(grid->x[ix] - grid->x[1 + \
ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[-1 + jy] - \
					grid->y[jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  PetscScalar vyim1jcoef = (2*Nc[1 + jy][ix+1].T21*(grid->x[ix] - grid->xc[1 + ix])*(grid->y[-1 + \
jy] - grid->y[jy])*(grid->y[jy] - grid->yc[jy]) + (grid->y[jy] - \
grid->y[1 + jy])*(2*Nc[jy][ix+1].T21*(grid->x[ix] - grid->xc[1 + \
ix])*(grid->y[jy] - grid->yc[1 + jy]) + Nb[jy][ix].T22*(grid->y[-1 + \
jy] - grid->y[jy])*(grid->yc[jy] - grid->yc[1 + \
jy])))/(2.*(grid->x[ix] - grid->x[1 + ix])*(grid->xc[ix] - grid->xc[1 \
+ ix])*(grid->y[-1 + jy] - grid->y[jy])*(grid->y[jy] - grid->y[1 + \
							       jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  PetscScalar vyim1jp1coef = (Nb[1 + jy][ix].T12*(grid->x[ix] - grid->x[1 + ix])*(grid->y[jy] - \
grid->y[1 + jy]) - 8*Nc[1 + jy][ix+1].T21*(grid->x[ix] - grid->xc[1 + \
ix])*(grid->y[jy] - grid->yc[jy]))/(8.*(grid->x[ix] - grid->x[1 + \
ix])*(grid->xc[ix] - grid->xc[1 + ix])*(grid->y[jy] - grid->y[1 + \
							      jy])*(grid->yc[jy] - grid->yc[1 + jy]));
	  if( bv->mechBCLeft.type[1] == 1 ){/* free slip */
	    ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl-1][ixl], vyim1jm1coef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl][ixl], vyim1jcoef, ADD_VALUES);CHKERRQ(ierr);
	    ierr = MatSetValue(LHS, vydof[jyl][ixl], vydof[jyl+1][ixl], vyim1jp1coef, ADD_VALUES);CHKERRQ(ierr);
	  }
	  
	  /* if also at top boundary, add ghost values for vx[i][j-2] and vx[i+1][j-2] */
	  if( jy == 1){
	    PetscScalar vxijm2coef = Nb[-1 + jy][ix].T12/(8.*(grid->yc[-1 + jy] - grid->yc[jy])*(-grid->yc[jy] + grid->yc[1 + jy]));
	    PetscScalar vxip1jm2coef = Nb[-1 + jy][1 + ix].T12/(8.*(grid->yc[-1 + jy] - grid->yc[jy])*(-grid->yc[jy] + grid->yc[1 + jy]));
	    if( bv->mechBCTop.type[0] == 1){
	      ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl-1][ixl], vxijm2coef, ADD_VALUES);CHKERRQ(ierr);
	      ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl-1][ixl+1], vxip1jm2coef, ADD_VALUES);CHKERRQ(ierr);
	    }	   
	  }
      }
      /* Right Boundary */
      /* only need to make changes for ix==NX-2 and jy == 1 */

      /* Top Boundary */
      if( jy == 1 && (!grid->xperiodic || (ix > 0 && ix < NX-1) )){
	/* need to set implicit ghost values for vx[i][j-2] and vx[i+1][j-2] */
	PetscScalar vxijm2coef = Nb[-1 + jy][ix].T12/(8.*(grid->yc[-1 + jy] - grid->yc[jy])*(-grid->yc[jy] + grid->yc[1 + jy]));
	PetscScalar vxip1jm2coef = Nb[-1 + jy][1 + ix].T12/(8.*(grid->yc[-1 + jy] - grid->yc[jy])*(-grid->yc[jy] + grid->yc[1 + jy]));
	if( bv->mechBCTop.type[0] == 1){
	  ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl-1][ixl], vxijm2coef, ADD_VALUES);CHKERRQ(ierr);
	  ierr = MatSetValue(LHS, vydof[jyl][ixl], vxdof[jyl-1][ixl+1], vxip1jm2coef, ADD_VALUES);CHKERRQ(ierr);
	}      

      }
      /* Bottom Boundary */
      /* Bottom boundary only special at jy == NY-1 and ix == 0, treated above */
      
    }
  }


  ierr = DMDAVecRestoreArray(grid->tda,Nbl,&Nb);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->tda,Ncl,&Nc);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(grid->da,rhol,&rho);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(RHS);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(RHS);CHKERRQ(ierr);
  ierr = MatAssemblyBegin( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  



  ierr = VecDestroy(& Nbl ); CHKERRQ(ierr);
  ierr = VecDestroy(& Ncl ); CHKERRQ(ierr);

  ierr=VecDestroy(& rhol);CHKERRQ(ierr);


  /* free dof maps*/
  //ierr=PetscFree(vxdof[-1]);CHKERRQ(ierr);
  //ierr=PetscFree(vydof[-1]);CHKERRQ(ierr);
  for(jy=0;jy<ng;jy++){
    vxdof[jy] --;
    vydof[jy] --;
    ierr=PetscFree(vxdof[jy]);CHKERRQ(ierr);
    ierr=PetscFree(vydof[jy]);CHKERRQ(ierr);
    ierr=PetscFree(pdof[jy]);CHKERRQ(ierr); 
    /*     ierr=PetscFree(vzdof[jy]);CHKERRQ(ierr); */
  }
  vxdof--; vydof--;
  ierr=PetscFree(vxdof);CHKERRQ(ierr);
  ierr=PetscFree(vydof);CHKERRQ(ierr);
  ierr=PetscFree(pdof);CHKERRQ(ierr);
  /*   ierr=PetscFree(vzdof);CHKERRQ(ierr); */

  /*   ierr = VecAssemblyEnd(RHSz); CHKERRQ(ierr); */

  /*   ierr = MatAssemblyEnd( LHSz, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); */

  PetscFunctionReturn(ierr);

}
