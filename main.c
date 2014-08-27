//#include "../include/fem.h"
//
/*
 *Main routine of FEM code. This copy is specifically for ridge damage problem
 *Do not put any compiler macros here. They should be in fem.h.
 */
static char help[] = "Help statement goes here\n\
 petscmpiexec -n 2 ./Fem ./input/fault.1 output/test -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps >log1\n\
It can be one line or many \n";

#include<stdio.h>
#include<stdlib.h>
//#include<malloc.h>
#include "petscksp.h"
#include "mpi.h"
#include "structures.h"
#include "kspLinearSolve.h"

//Include file for PETSc linear solve commands

//extern struct ELEMENT linearHexahedralElementCreate();
//extern double evaluateNodalBasisFunction();

//extern struct ELEMENT quadTriElementCreate();
//extern double evaluateNodalBasisFunction();

#define NADD 12 
//maximum number of entries to add at once.

#undef __FUNCT__
#define __FUNCT__ "main"

int main(int argc, char **args)
{
  PetscMPIInt rank,size;
  PetscErrorCode ierr;
  KSP	ksp; /* linear solver context */
  PetscViewer viewer;
  PetscLogStage stages[4];  
  /*end debugging stuff*/
  
  PetscInt NX,NY;/* number of nodes in x,y,directions*/
  //  PetscScalar eta;

  PetscScalar *etaS,*etaN;/* center and nodal viscosity */
  
  PetscScalar *rho; /*nodal density dist.*/
  PetscScalar LX,LY;/* model dimensions (x,y)*/
  PetscScalar dx,dy;/*x,y grid spacing*/
  PetscScalar dx2,dy2,dxdy;/*dx*dx, dy*dy*/
  Mat LHS;
  Vec RHS;

  /* markers*/
  PetscScalar *markerX, *markerY;/* marker nodal positions*/
  PetscScalar *markerVX, *markerVY;/* marker velocities*/
  char *markerMat;/*marker material*/

  //INITIALIZE PETSC
  PetscInitialize(&argc,&args,(char *)0,help); //copied from petsc/src/ksp/ksp/examles/tutorials
  //Set MPI rank
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  /* grid parameters */
  NX=21;
  NY=31;
  PetscInt NMX=10;/* number of markers per cell, x*/
  PetscInt NMY=10;/* markers per cell y*/
  PetscInt ndof = NX*NY*3;
  PetscInt nTime = 100;


  /* Model Parameters */
  LX=1000000.0;
  LY=1500000.0;
  //  eta=1.0e21;/* viscosity, pa-s*/
  /*   PetscScalar rholeft = 3200; */
  /*   PetscScalar rhoright = 3300; */
  /*   PetscScalar etaleft=1.0e20; */
  /*   PetscScalar etaright=1.0e22; */

  PetscScalar gy = 10;
  PetscScalar materialEta[2] = {1.0e20, 1.0e22};
  PetscScalar materialRho[2] = {3200.0, 3300.0};
  PetscScalar etamin = 1.0e20;
  /* grid setup*/
  PetscInt nMark = NX*NMX*NY*NMY;
  PetscScalar dtMax = 3.16e13;/* timestep in seconds*/
  PetscScalar elapsedTime=0.0;
  dx=LX/((PetscScalar) NX - 1);
  dy=LY/((PetscScalar) NY - 1);
  dx2 = dx*dx;
  dy2 = dy*dy;
  dxdy=dx*dy;

  /* scaling parameters*/
  PetscScalar Kbond = 4*etamin/((dx+dy)*(dx+dy));
  PetscScalar Kcont = 2*etamin/(dx+dy);

  /* initialize markers*/
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &markerX);CHKERRQ(ierr);
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &markerY);CHKERRQ(ierr);
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &markerVX);CHKERRQ(ierr);
  ierr = PetscMalloc( nMark*sizeof(PetscScalar), &markerVY);CHKERRQ(ierr);
  ierr = PetscMalloc( nMark*sizeof(char), &markerMat); CHKERRQ(ierr);
  PetscRandom r;
  PetscRandomCreate( PETSC_COMM_WORLD, &r);CHKERRQ(ierr);
  PetscRandomSetType(r,PETSCRAND);CHKERRQ(ierr);
  ierr=PetscRandomSetSeed(r, (unsigned long) 0);CHKERRQ(ierr);

  PetscScalar tmp;
  PetscInt i,j;
  PetscScalar mdx = LX/(((PetscScalar) NMX)*((PetscScalar) NX)-1.0);
  PetscScalar mdy = LY/(((PetscScalar) NMY)*((PetscScalar) NY)-1.0);
  for(i=1;i<=NMX*NX;i++){
    for (j=1;j<=NMY*NY;j++){
      PetscRandomGetValue(r, &tmp);/* random number in [0,1]*/
      markerX[ (i-1)*NMY*NY + j-1]= mdx*((PetscScalar)i - 1.0)+ (tmp-0.5)*mdx;
      PetscRandomGetValue(r, &tmp);
      markerY[ (i-1)*NMY*NY + j-1]= mdy*((PetscScalar)j - 1.0)+ (tmp-0.5)*mdy;
    }
  }
  /* make sure that markers are inside domain*/
  PetscScalar markerEps = 1;
  for(i=0;i<nMark;i++){
    if( markerX[i] < 0){ markerX[i] = markerEps;
    }else if(markerX[i] > LX){ markerX[i] = LX-markerEps;
    }else if(markerY[i] < 0.0){ markerY[i] = markerEps;
    }else if(markerY[i] > LY){ markerY[i] = LY-markerEps;}
  }

  for(i=1;i<=nMark;i++){
    if(markerX[ i-1] < LX/2){
      markerMat[i-1] = 0;
    } else{
      markerMat[i-1] = 1;
    }
  }


  /*initial density distribution*/
  ierr = PetscMalloc( NX*NY*sizeof(PetscScalar), &rho);CHKERRQ(ierr);
  /* arrays to hold nodal and mid-cell viscosity */
  ierr = PetscMalloc( NX*NY*sizeof(PetscScalar), &etaS);CHKERRQ(ierr);
  ierr = PetscMalloc( NX*NY*sizeof(PetscScalar), &etaN);CHKERRQ(ierr);

  PetscInt iTime;

  for(iTime=0;iTime<nTime;iTime++){
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Beginning Timestep %d\n",iTime);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);    
    
/*     for(i=0;i<nMark;i++){ */
/*       if( markerX[i] < 0.0){ markerX[i] = markerEps; */
/*       }else if(markerX[i] > LX){ markerX[i] = LX-markerEps; */
/*       }else if(markerY[i] < 0.0){ markerY[i] = markerEps; */
/*       }else if(markerY[i] > LY){ markerY[i] = LY-markerEps;} */
/*     } */

    /* 1. project marker quantities onto nodes */
    PetscInt m;
    /* zero out nodal viscosity and density*/
    for(i=0;i<NX*NY;i++){
      rho[i] = 0.0;
      etaS[i] = 0.0;
      etaN[i] = 0.0;
    }
    
    PetscScalar *weightSum;/* sum of nodal weights*/
    ierr = PetscMalloc( NX*NY*sizeof(PetscScalar), &weightSum);CHKERRQ(ierr);
    for(i=0;i<NX*NY;i++){
      weightSum[i]=0.0;
    }
    
    /*this loop calculates contribution of marker to each node in bounding box.*/
    for(m=0;m<nMark;m++){
      /* find cell i,j that I belong to.*/
      /*      j           j+1   */
      /* i    o           o     */
      /*                        */
      /*            m           */
      /*                        */
      /* i+1  o           o     */

      /* calculate contribution to i,j*/
      PetscInt cellI = floor(markerY[m]/dy);/* these are indices from 0 to N-1*/
      PetscInt cellJ = floor(markerX[m]/dx);
      //printf("cellI = %d, cellJ = %d\n",cellI,cellJ);
      /* calculate marker weight*/
      PetscScalar deltaxm = markerX[m] - ((PetscScalar) cellJ)*dx;
      PetscScalar deltaym = markerY[m] - ((PetscScalar) cellI)*dy;
      /*weight*/
      PetscScalar weight = (1 - deltaxm/dx)*(1-deltaym/dy);
      rho[ cellJ*NY + cellI] += weight*materialRho[(PetscInt) markerMat[m]];
      etaS[cellJ*NY + cellI] += weight*materialEta[(PetscInt) markerMat[m]];
      weightSum[  cellJ*NY + cellI] += weight;
      
      /* calculate contribution to node i+1,j+1*/
      if(cellI<NY-1 && cellJ<NX-1){
	deltaxm = ((PetscScalar) cellJ+1)*dx - markerX[m];
	deltaym = ((PetscScalar) cellI+1)*dy - markerY[m];
	weight = (1 - deltaxm/dx)*(1-deltaym/dy);
	rho[ (cellJ+1)*NY + cellI+1] += weight*materialRho[(PetscInt) markerMat[m]];
	etaS[(cellJ+1)*NY + cellI+1] += weight*materialEta[(PetscInt) markerMat[m]];
	weightSum[  (cellJ+1)*NY + cellI+1] += weight;
      }
      /* contribution to node i,j+1*/
      if(cellJ < NX-1){
	deltaxm = ((PetscScalar) cellJ+1)*dx - markerX[m];
	deltaym = markerY[m]-((PetscScalar) cellI)*dy;
	weight = (1 - deltaxm/dx)*(1-deltaym/dy);
	rho[ (cellJ+1)*NY + cellI] += weight*materialRho[(PetscInt) markerMat[m]];
	etaS[(cellJ+1)*NY + cellI] += weight*materialEta[(PetscInt) markerMat[m]];
	weightSum[  (cellJ+1)*NY + cellI] += weight;
      }
      /* contribution to node i+1,j*/
      if(cellI<NY-1){
	deltaxm = markerX[m]-((PetscScalar) cellJ)*dx;
	deltaym = ((PetscScalar) cellI+1)*dy - markerY[m];
	weight = (1 - deltaxm/dx)*(1-deltaym/dy);
	rho[ cellJ*NY + cellI+1 ] += weight*materialRho[(PetscInt) markerMat[m]];
	etaS[cellJ*NY + cellI+1 ] += weight*materialEta[(PetscInt) markerMat[m]];
	weightSum[  cellJ*NY + cellI+1] += weight;
      }
    }

    
    /* normalize nodal quantities*/
    for(i=0;i<NX*NY;i++){
      rho[i] = rho[i]/weightSum[i];
      etaS[i] = etaS[i]/weightSum[i];
    }
    ierr = PetscFree(weightSum);
    
    /* calculate etaN by nodal averaging */
    for(j=1;j<NX-1;j++){/* the limits of iteration seem wrong to me but are consistent with TGs code*/
      for(i=1;i<NY-1;i++){
	PetscInt idxnode = j*(NY) + i;
	/*       etaN(i,j) = (etaS(i-1,j-1) + etaS(i,j-1) + etaS(i,j) + etaS(i-1,j))/4; */
	etaN[idxnode] = (etaS[idxnode-NY-1] + etaS[idxnode-NY] + etaS[idxnode] + etaS[idxnode-1])/4.0;
      }
    }
    
    /* initialize matrix for left hand side and vector for right hand side*/
    ierr = MatCreate( PETSC_COMM_WORLD, &LHS);CHKERRQ(ierr);
    ierr = MatSetSizes( LHS, PETSC_DECIDE,PETSC_DECIDE,ndof,ndof); CHKERRQ(ierr);
    ierr = MatSetType( LHS, MATAIJ); CHKERRQ(ierr);
    ierr = PetscObjectSetName( (PetscObject) LHS,"CLHS");
    ierr = VecCreate( PETSC_COMM_WORLD, &RHS); CHKERRQ(ierr);
    ierr = VecSetSizes( RHS, PETSC_DECIDE, ndof);CHKERRQ(ierr);
    ierr = VecSetType( RHS, VECMPI);CHKERRQ(ierr);
    ierr = PetscObjectSetName( (PetscObject) RHS,"CRHS");
    
    
    
    for (i=1;i<=NY;i++){
      for(j=1;j<=NX;j++){
	PetscInt idxnode = (j-1)*(NY)+(i-1);
	PetscInt pdof=3*idxnode+0;
	PetscInt vxdof=3*idxnode+1;
	PetscInt vydof=3*idxnode+2;
	
	PetscScalar vals[NADD];
	PetscInt rowidx[NADD];
	PetscInt colidx[NADD];
	
	
	/**              %continuity*/		       
	if( (i == 1 || j == 1) ){
	  /*       %pressure ghost nodes. Set equal to zero. */
	  ierr = MatSetValue(LHS,pdof,pdof,1.0*Kbond,INSERT_VALUES); CHKERRQ(ierr);
	  ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES); CHKERRQ(ierr);
	} else if( (i==2 && j == 2) || (i==NY && j==2)){
	  //upper and lower left corner - dp/dx = 0;
	  //LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	  ierr = MatSetValue(LHS,pdof,pdof, 1.0*Kbond,INSERT_VALUES); CHKERRQ(ierr);
	  //      LHS(p(i,j),p(i,j+1)) = -1.0*Kbond;
	  ierr = MatSetValue(LHS,pdof,pdof + 3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  //      RHS(p(i,j)) = 0.0;
	  ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES); CHKERRQ(ierr);
	} else if( (i==2 && j==NX) || (i==NY && j== NX)) {
	  //upper and lower right corners - dp/dx = 0;
	  //LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	  ierr = MatSetValue(LHS,pdof,pdof,1.0*Kbond,INSERT_VALUES); CHKERRQ(ierr);
	  //LHS(p(i,j),p(i,j-1)) = -1.0*Kbond;
	  ierr = MatSetValue(LHS,pdof,pdof - 3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  //	RHS(p(i,j)) = 0.0;
	  ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES); CHKERRQ(ierr);
	  
	} else if( i==2 && j==3){
	  //	LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	  ierr = MatSetValue(LHS,pdof,pdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  //RHS(p(i,j)) = 0.0;
	  ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES);CHKERRQ(ierr);
	} else {
	  //LHS(p(i,j), vx(i-1,j)) =  Kcont/dx;
	  ierr = MatSetValue(LHS,pdof,vxdof - 3,Kcont/dx,INSERT_VALUES);CHKERRQ(ierr);
	  //LHS(p(i,j), vx(i-1,j-1) ) = -Kcont/dx;
	  ierr = MatSetValue(LHS,pdof,vxdof-3-3*NY,-Kcont/dx,INSERT_VALUES);CHKERRQ(ierr);
	  //LHS(p(i,j), vy(i,j-1) ) =  Kcont/dy;
	  ierr = MatSetValue(LHS,pdof,vydof-3*NY,Kcont/dy,INSERT_VALUES);CHKERRQ(ierr);
	  //LHS(p(i,j), vy(i-1,j-1) ) = -Kcont/dy;
	  ierr = MatSetValue(LHS,pdof,vydof-3-3*NY,-Kcont/dy,INSERT_VALUES);CHKERRQ(ierr);
	  //      end
	}
	
	/*         %x-momentum */
	if(i == NY){
	  /*             %ghost vx nodes. set to zero. */
	  /*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; RHS(vx(i,j)) = 0.0; */
	  ierr =  MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /* RHS is automatic*/
	  
	} else if( (j==1 || j == NX) && i < NY){
	  /*             %left boundary. vx=0; */
	  /*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; RHS(vx(i,j)) = 0.0; */
	  ierr =  MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	}else if( (i ==1) && j>1 && j<NX ){
	  /*             %top boundary - free slip, dvx/y = 0 */
	  /*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; */
	  ierr=MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             LHS(vx(i,j),vx(i+1,j)) = -1.0*Kbond; */
	  ierr=MatSetValue(LHS,vxdof,vxdof+3,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             RHS(vx(i,j)) = 0.0; */
	}else if( i==NY-1 && j>1 && j<NX ) {
	  /*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; */
	  ierr=MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             LHS(vx(i,j),vx(i-1,j)) = -1.0*Kbond; */
	  ierr=MatSetValue(LHS,vxdof,vxdof-3,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             RHS(vx(i,j)) = 0.0; */
	}else{
	  
	  /*   LHS(vx(i,j),vx(i,j+1)) =  2*etaN(i+1,j+1)/dx2; */
	  rowidx[0] = vxdof; colidx[0] = vxdof+3*NY; vals[0] = 2*etaN[idxnode+NY+1]/dx2;
	  /*             LHS(vx(i,j),vx(i,j)) =   -2*etaN(i+1,j+1)/dx2- 2*etaN(i+1,j)/dx2 -etaS(i+1,j)/dy2 - etaS(i,j)/dy2; */
	  rowidx[1] = vxdof; colidx[1] = vxdof; vals[1] = -2.0*etaN[idxnode+NY+1]/dx2 -2.0*etaN[idxnode+1]/dx2 - etaS[idxnode+1]/dy2 - etaS[idxnode]/dy2;
	  /*             LHS(vx(i,j),vx(i,j-1)) = 2*etaN(i+1,j)/dx2; */
	  rowidx[2] = vxdof;
	  colidx[2] = vxdof-3*NY;
	  vals[2] = 2*etaN[idxnode+1]/dx2;
	  /*             LHS(vx(i,j),vx(i+1,j)) = etaS(i+1,j)/dy2;  */
	  rowidx[3] = vxdof;
	  colidx[3] = vxdof+3;
	  vals[3] = etaS[idxnode+1]/dy2;
	  /*             LHS(vx(i,j),vy(i+1,j)) = etaS(i+1,j)/dxdy; */
	  rowidx[4] = vxdof;
	  colidx[4] = vydof+3;
	  vals[4] = etaS[idxnode+1]/dxdy;
	  /*             LHS(vx(i,j),vy(i+1,j-1)) = -etaS(i+1,j)/dxdy; */
	  rowidx[5] = vxdof;
	  colidx[5] = vydof+3-3*NY;
	  vals[5] = -1.0*etaS[idxnode+1]/dxdy;
	  /*             LHS(vx(i,j),vx(i-1,j)) = +etaS(i,j)/dy2; */
	  rowidx[6] = vxdof;
	  colidx[6] = vxdof-3;
	  vals[6] = 1.0*etaS[idxnode]/dy2;
	  /*             LHS(vx(i,j),vy(i,j)) = -etaS(i,j)/dxdy; */
	  rowidx[7] = vxdof;
	  colidx[7] = vydof;
	  vals[7] = -1.0*etaS[idxnode]/dxdy;
	  /*             LHS(vx(i,j),vy(i,j-1)) = etaS(i,j)/dxdy; */
	  rowidx[8] = vxdof;
	  colidx[8] = vydof-3*NY;
	  vals[8] = etaS[idxnode]/dxdy;
	  /*             LHS(vx(i,j),p(i+1,j+1)) = - Kcont/dx; */
	  rowidx[9] = vxdof;
	  colidx[9] = pdof+3*NY+3;
	  vals[9] = -Kcont/dx;
	  /*             LHS(vx(i,j),p(i+1,j)) =  Kcont/dx; */
	  rowidx[10] = vxdof;
	  colidx[10] = pdof+3;
	  vals[10] = Kcont/dx;
	  /*             RHS(vx(i,j)) =  0; */
	  ierr = MatSetValues( LHS, 1,&rowidx[0],11,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
	  
	  
	} 
	
	/*         %y-momentum */
	if( j==NX ) {
	  /*             %ghost vy nodes. set to zero. */
	  /*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; RHS(vy(i,j)) =0.0; */
	  ierr=MatSetValue(LHS,vydof,vydof,Kbond,INSERT_VALUES);CHKERRQ(ierr);
	}else if( i==1 || i == NY) {
	  /*             %top and bottom boundary. vy=0 */
	  /*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; RHS(vy(i,j)) =0.0; */
	  ierr=MatSetValue(LHS,vydof,vydof,Kbond,INSERT_VALUES);CHKERRQ(ierr);
	}else if( j==1 && i>1 && i < NY){
	  /*             %left wall - free slip - dvy/dx = 0; */
	  /*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; */
	  ierr=MatSetValue(LHS,vydof,vydof, Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             LHS(vy(i,j),vy(i,j+1)) = -1.0*Kbond; */
	  ierr=MatSetValue(LHS,vydof,vydof+3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             RHS(vy(i,j)) = 0.0; */
	}else if( j == NX-1 && i >1 && i < NY) {
	  /*             %right wall - free slip - dvy/dx =0; */
	  /*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; */
	  ierr=MatSetValue(LHS,vydof,vydof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             LHS(vy(i,j),vy(i,j-1)) = -1.0*Kbond; */
	  ierr=MatSetValue(LHS,vydof,vydof-3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	  /*             RHS(vy(i,j)) = 0.0; */
	}else {
	  
	  /*             LHS(vy(i,j),vy(i+1,j)) = 2*etaN(i+1,j+1)/dy2; */
	  rowidx[0] = vydof; colidx[0] = vydof+3; vals[0] = 2*etaN[idxnode+NY+1]/dy2;
	  
	  /* LHS(vy(i,j),vy(i,j)) = -2*etaN(i+1,j+1)/dy2 - 2*etaN(i,j+1)/dy2 -etaS(i,j+1)/dx2 - etaS(i,j)/dx2 ; */
	  rowidx[1] = vydof;
	  colidx[1] = vydof;
	  vals[1] = -2.0*etaN[idxnode+NY+1]/dy2 -2.0*etaN[idxnode+NY]/dy2 - etaS[idxnode+NY]/dx2 - etaS[idxnode]/dx2;
	  /* LHS(vy(i,j),vy(i-1,j)) = 2*etaN(i,j+1)/dy2; */
	  rowidx[2] = vydof;
	  colidx[2] = vydof-3;
	  vals[2] = 2*etaN[idxnode+NY]/dy2;
	  /* LHS(vy(i,j),vy(i,j+1)) = etaS(i,j+1)/dx2; */
	  rowidx[3] = vydof;
	  colidx[3] = vydof+3*NY;
	  vals[3] = etaS[idxnode+NY]/dx2;
	  /*             LHS(vy(i,j),vx(i,j+1)) = etaS(i,j+1)/dxdy; */
	  rowidx[4] = vydof;
	  colidx[4] = vxdof+3*NY;
	  vals[4] = etaS[idxnode+NY]/dxdy;
	  /*             LHS(vy(i,j),vx(i-1,j+1)) = -etaS(i,j+1)/dxdy; */
	  rowidx[5] = vydof;
	  colidx[5] = vxdof-3+3*NY;
	  vals[5] = -1.0*etaS[idxnode+NY]/dxdy;
	  /*  LHS(vy(i,j),vy(i,j-1)) = etaS(i,j)/dx2; */
	  rowidx[6] = vydof;
	  colidx[6] = vydof-3*NY;
	  vals[6] = etaS[idxnode]/dx2;
	  /*  LHS(vy(i,j),vx(i,j))   = - etaS(i,j)/dxdy; */
	  rowidx[7] = vydof;
	  colidx[7] = vxdof;
	  vals[7] = -1.0*etaS[idxnode]/dxdy;
	  /* LHS(vy(i,j),vx(i-1,j)) = etaS(i,j)/dxdy; */
	  rowidx[8] = vydof;
	  colidx[8] = vxdof-3;
	  vals[8] = etaS[idxnode]/dxdy;
	  /* LHS(vy(i,j),p(i+1,j+1)) = -Kcont/dy; */
	  rowidx[9] = vydof;
	  colidx[9] = pdof+3*NY+3;
	  vals[9] = -Kcont/dy;
	  /* LHS(vy(i,j),p(i,j+1)) = Kcont/dy; */
	  rowidx[10] = vydof;
	  colidx[10] = pdof+3*NY;
	  vals[10] = Kcont/dy;
	  /*   RHS(vy(i,j)) = -gy*(RHO(i,j)+RHO(i,j+1))/2; */
	  ierr = VecSetValue( RHS, vydof, -gy*(rho[idxnode]+rho[idxnode+NY])/2.0,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = MatSetValues( LHS, 1,&rowidx[0],11,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
	}
      }/* end loop over i (Y)*/
    }/* end loop over j (x) */
    
    /* finalize assembly*/
    ierr = VecAssemblyBegin(RHS);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(RHS);CHKERRQ(ierr);
    
    ierr = MatAssemblyBegin( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    
    
    /* solve linear system */
    //  ierr 
    Vec S;
    ierr = VecDuplicate( RHS, &S);
    ierr = PetscObjectSetName( (PetscObject) S,"CS");
    ierr = kspLinearSolve( LHS, RHS, S);CHKERRQ(ierr);
    
    /*get solution as an array */
    PetscScalar *Sa;
    ierr = VecGetArray(S,&Sa); CHKERRQ(ierr);
    
    PetscScalar vxmax=0.0;
    PetscScalar vymax=0.0;
    
    /* Project nodal velocity solution onto markers*/
    for(m=0;m<nMark;m++){
      /*vx node indices*/
      PetscInt vxI, vxJ, vyI,  vyJ;
      vxI = floor((markerY[m]-dy/2.0)/dy);
      vxJ = floor((markerX[m])/dx);
      vyI = floor((markerY[m])/dy);
      vyJ = floor((markerX[m]-dx/2.0)/dx);
      
      if( vxI < 0){ vxI=0;} else if(vxI>NY-3){ vxI=NY-3;}/* do not let a marker be assigned to one of the ghost cells*/
      if( vxJ < 0){ vxJ=0;} else if(vxJ>NX-2){ vxJ=NX-2;}
      if( vyI < 0){ vyI=0;} else if(vyI>NY-2){ vyI=NY-2;}
      if( vyJ < 0){ vyJ=0;} else if(vyJ>NX-3){ vyJ=NX-3;}

      /* compute new x-vel*/
      PetscScalar deltaxm=markerX[m] - dx*((PetscScalar) vxJ);
      PetscScalar deltaym=markerY[m] - dy*(((PetscScalar) vxI)+0.5);
      PetscInt vxdof = 3*vxJ*NY + 3*vxI + 1;
      /* markerVX(myMarkers) = S(vx(i,j))*(1-deltaXM/dx).*(1-deltaYM/dy) + S(vx(i,j+1))*(deltaXM/dx).*(1-deltaYM/dy) + S(vx(i+1,j))*(1-deltaXM/dx).*(deltaYM/dy) + S(vx(i+1,j+1))*(deltaXM.*deltaYM)/dx/dy; */
      markerVX[m] = Sa[vxdof]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + Sa[vxdof+3*NY]*(deltaxm/dx)*(1-deltaym/dy) + Sa[vxdof+3]*(1-deltaxm/dx)*(deltaym/dy) + Sa[vxdof+3*NY+3]*(deltaxm*deltaym)/dx/dy; 
      /*compute new y-vel*/
      deltaxm=markerX[m] - dx*(((PetscScalar) vyJ)+0.5);
      deltaym=markerY[m] - dy*((PetscScalar) vyI);
      PetscInt vydof = 3*vyJ*NY+3*vyI + 2;
      markerVY[m] = Sa[vydof]*(1.0-deltaxm/dx)*(1.0-deltaym/dy) + Sa[vydof+3*NY]*(deltaxm/dx)*(1-deltaym/dy) + Sa[vydof+3]*(1-deltaxm/dx)*(deltaym/dy) + Sa[vydof+3*NY+3]*(deltaxm*deltaym)/dx/dy; 
      
    }
    
    ierr = VecRestoreArray(S,&Sa);CHKERRQ(ierr); 
    /* determine correct time step*/
    for(m=0;m<nMark;m++){
      if(fabs(markerVY[m]) > vymax) vymax=fabs(markerVY[m]);
      if(fabs(markerVX[m]) > vxmax) vxmax=fabs(markerVX[m]);
    }
    PetscScalar dt;
    dt = dx/2.0/vxmax;
    if( dt > dy/2.0/vymax) dt=dy/2.0/vymax;
    if(dt > dtMax) dt=dtMax;
    printf("Using Timestep %e years\n",dt/3.1557e7);
    elapsedTime += dt;
    printf("elapsedTime %e years\n",elapsedTime/3.1557e7);
    /* move markers*/
    for(m=0;m<nMark;m++){
      markerX[m] += markerVX[m]*dt;
      markerY[m] += markerVY[m]*dt;
    }
    
    /* save marker information */
    FILE *mf;
    char markName[80];
    sprintf(markName,"output/loadMarkers%d.m",iTime);

    mf=fopen(markName,"w");
    fprintf(mf,"markerInfo=[");
    for(i=0; i<nMark;i++){
      fprintf(mf,"%e,%e,%e,%e,%d;\n",markerX[i],markerY[i],markerVX[i],markerVY[i],markerMat[i]);
    }
    fprintf(mf,"];");
    fclose(mf);
    /* dump nodal density*/
    FILE *of,*of2,*of3;
    char name1[80];
    char name2[80];
    char name3[80];
    sprintf(name1,"output/loadrho%d.m",iTime);
    sprintf(name2,"output/loadetaS%d.m",iTime);
    sprintf(name3,"output/loadetaN%d.m",iTime);
    

    /*save nodal values of rho, etaS, etaN*/
    of=fopen(name1,"w");
    of2=fopen(name2,"w");
    of3=fopen(name3,"w");
    
    fprintf(of,"Crho=[");
    fprintf(of2,"CetaS=[");
    fprintf(of3,"CetaN=[");
    for(i=0;i<NY;i++){
      for(j=0;j<NX;j++){
	fprintf(of,"%e,",rho[j*NY+i]);
	fprintf(of2,"%e,",etaS[j*NY+i]);
	fprintf(of3,"%e,",etaN[j*NY+i]);
      }
      fprintf(of,";\n");
      fprintf(of2,";\n");
      fprintf(of3,";\n");
    }
    fprintf(of,"];\n");
    fprintf(of2,"];\n");
    fprintf(of3,"];\n");
    

    fclose(of);
    fclose(of2);
    fclose(of3);
    
  }/* end time step*/


  /* clean up*/
  ierr = MatDestroy( LHS);CHKERRQ(ierr);
  ierr = VecDestroy(RHS);CHKERRQ(ierr);
  ierr = PetscFree(rho);CHKERRQ(ierr);
  ierr = PetscFree(etaS);CHKERRQ(ierr);
  ierr = PetscFree(etaN);CHKERRQ(ierr);
  /* clean up markers*/
  ierr = PetscFree(markerX);CHKERRQ(ierr);
  ierr = PetscFree(markerY);CHKERRQ(ierr);
  ierr = PetscFree(markerVX);CHKERRQ(ierr);
  ierr = PetscFree(markerVY);CHKERRQ(ierr);
  ierr = PetscFree(markerMat);CHKERRQ(ierr);



  PetscLogStageRegister("Main",&stages[0]);
  PetscLogStageRegister("MechNewt",&stages[1]);
  PetscLogStageRegister("Thermal Problem",&stages[2]);
  PetscLogStageRegister("Thermal Problem Debug",&stages[3]);
  PetscLogStagePush(stages[0]);
  


  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}




