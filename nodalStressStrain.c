#include "fdcode.h"
#include "nodalStressStrain.h"
#include "profile.h"
/* calculate stress and strain at nodes*/

PetscErrorCode nodalStressStrain( GridData *grid, NodalFields *nodalFields, Options *options, BoundaryValues *bv, PetscScalar dt, Vec nodalHeating, PetscScalar gy){
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  //PetscScalar dx = grid->dx;
  //PetscScalar dy = grid->dy;
  PetscInt x,y,m,n;
  PetscErrorCode ierr;
  PetscInt ix,jy;  
  
  PetscFunctionBegin;
  setLogStage( LOG_NODAL_STR );

  /* retrieve local arrays of nodal fields, including ghost node values*/
  Vec vxl,vyl;
  Vec exxl,eyyl,exyl;

  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  
  ierr=DMCreateLocalVector(grid->da,&vxl);
  ierr=VecDuplicate(vxl,&vyl);
  //  ierr=VecDuplicate(vxl,&vzl);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);

/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr); */

  ierr=VecDuplicate(vxl,&exxl);CHKERRQ(ierr);
  ierr=VecDuplicate(vxl,&eyyl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(vxl,&ezzl);CHKERRQ(ierr); */
  ierr=VecDuplicate(vxl,&exyl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(vxl,&exzl);CHKERRQ(ierr); */
/*   ierr=VecDuplicate(vxl,&eyzl);CHKERRQ(ierr); */

  PetscScalar **vx, **vy, **exx, **eyy, **exy;  
  ierr=DMDAVecGetArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  //  ierr=DMDAVecGetArray( grid->da,vzl,&vz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,exyl,&exy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,exzl,&exz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,eyzl,&eyz);CHKERRQ(ierr); */


  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){
      if( (!grid->xperiodic && ix==0) || jy == 0){
	/* account for ghost cells */
	exx[jy][ix] = 0.0;
	eyy[jy][ix] = 0.0;
/* 	ezz[jy][ix] = 0.0; */
      }else{
	/* this is a cell-centered quantity*/
	//      PetscScalar divv = (vx[i-1+j*NY]-vx[i-1+(j-1)*NY])/dx +(vy[i+NY*(j-1)]-vy[i-1+NY*(j-1)])/dy;
	PetscScalar dx=grid->x[ix]-grid->x[ix-1];
	PetscScalar dy=grid->y[jy]-grid->y[jy-1];
	PetscScalar divv = (vx[jy-1][ix]-vx[jy-1][ix-1])/dx + (vy[jy][ix-1]-vy[jy-1][ix-1])/dy;/* velocity divergence == trace of dev. strain rate*/
	
	/* exx[j*NY+i]=(vx[i-1+j*NY]-vx[i-1+(j-1)*NY])/dx - divv/3.0; */ /* should be deviatoric strain rate*/
	exx[jy][ix]=(vx[jy-1][ix]-vx[jy-1][ix-1])/dx - divv/3.0; /* should be deviatoric strain rate*/
	eyy[jy][ix]=(vy[jy][ix-1]-vy[jy-1][ix-1])/dy - divv/3.0; /* should be deviatoric strain rate*/
	//ezz[jy][ix]=-divv/3.0; /* deviatoric strain rate in z direction. ezz = 0 but e'zz is nonzero*/
      }
    }
  }
  /* scatter deviatoric strain rates back to global vectors*/
  ierr=DMDAVecRestoreArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
/*   ierr=DMDAVecRestoreArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); */
  ierr=DMLocalToGlobalBegin(grid->da,exxl,INSERT_VALUES,nodalFields->edotxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,exxl,INSERT_VALUES,nodalFields->edotxx);CHKERRQ(ierr);

  ierr=DMLocalToGlobalBegin(grid->da,eyyl,INSERT_VALUES,nodalFields->edotyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,eyyl,INSERT_VALUES,nodalFields->edotyy);CHKERRQ(ierr);

/*   ierr=DMLocalToGlobalBegin(grid->da,ezzl,INSERT_VALUES,nodalFields->edotzz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalEnd(grid->da,ezzl,INSERT_VALUES,nodalFields->edotzz);CHKERRQ(ierr); */


  /* Compute cell-centered quantities: ds_{ii}, s_{ii} */
  Vec X,X1;
  ierr=VecDuplicate(nodalFields->soxx,&X);CHKERRQ(ierr);
  /*   PetscScalar X = etaN[i+NY*j]/(etaN[i+NY*j]+dt*muN[i+NY*j]); */
  ierr = VecCopy(nodalFields->etaN,X);CHKERRQ(ierr); /* X = etaN */
  ierr = VecAXPY(X,dt,nodalFields->muN);CHKERRQ(ierr); /* X = X+dt*muN */
  ierr = VecPointwiseDivide(X,nodalFields->etaN,X);CHKERRQ(ierr);/* X = etaN./X  */
  ierr = VecDuplicate(X, &X1);CHKERRQ(ierr);
  ierr = VecCopy(X,X1);CHKERRQ(ierr); /* X1 = X */
  ierr = VecScale(X1,-2.0);CHKERRQ(ierr); /* X1 = X1*-2.0 */
  ierr = VecShift(X1,2.0);CHKERRQ(ierr); /* X1 = X1 + 2.0 */
  ierr = VecPointwiseMult(X1,nodalFields->etaN,X1);CHKERRQ(ierr); /* X1 = etaN.*X1 -> X1 = etaN.*(-2.0*X+2.0) */

  Vec stemp,etemp;
  ierr=VecDuplicate(X, &stemp);CHKERRQ(ierr);
  ierr=VecDuplicate(X, &etemp);CHKERRQ(ierr);

  /* calculate new sxx*/
  /*       PetscScalar sxxnew = (1.0-X)*2.0*etaN[i+NY*j]*exx[i+NY*j] + X*soxx[i+NY*j]; */
  ierr=VecPointwiseMult(etemp,X1,nodalFields->edotxx);CHKERRQ(ierr); /* etemp = X1.*exx */
  ierr=VecPointwiseMult(stemp,X,nodalFields->soxx);  /* stemp = X.*soxx */
  ierr=VecAXPY(stemp,1.0,etemp);CHKERRQ(ierr);/* stemp = stemp+1.0*etmp */
  /*dsxx = sxxnew-soxx*/
  ierr=VecCopy(nodalFields->soxx,nodalFields->dsxx);CHKERRQ(ierr); /* dsxx = soxx */
  ierr=VecAYPX(nodalFields->dsxx,-1.0,stemp);CHKERRQ(ierr);/*dsxx = sxxnew - 1.0*soxx*/
  ierr=VecCopy(stemp,nodalFields->soxx);CHKERRQ(ierr); /* soxx[i+NY*j]= sxxnew; */
  /* new syy */
  /* PetscScalar syynew = (1.0-X)*2.0*etaN[i+NY*j]*eyy[i+NY*j] + X*soyy[i+NY*j]; */
  ierr = VecPointwiseMult( etemp, X1, nodalFields->edotyy);CHKERRQ(ierr);
  ierr = VecPointwiseMult( stemp, X, nodalFields->soyy);CHKERRQ(ierr);
  ierr = VecAXPY( stemp,1.0,etemp);CHKERRQ(ierr);
  ierr = VecCopy(nodalFields->soyy,nodalFields->dsyy);CHKERRQ(ierr);
  ierr = VecAYPX(nodalFields->dsyy,-1.0,stemp);CHKERRQ(ierr);
  ierr = VecCopy(stemp,nodalFields->soyy);CHKERRQ(ierr);
  /* new szz */
  /*PetscScalar szznew = (1.0-X)*2.0*etaN[i+NY*j]*ezz[i+NY*j] + X*sozz[i+NY*j]; */
  /*   ierr = VecPointwiseMult( etemp, X1, nodalFields->edotzz);CHKERRQ(ierr); */
  /*   ierr = VecPointwiseMult( stemp, X, nodalFields->sozz);CHKERRQ(ierr); */
  /*   ierr = VecAXPY( stemp,1.0,etemp);CHKERRQ(ierr); */
  /*   ierr = VecCopy(nodalFields->sozz,nodalFields->dszz);CHKERRQ(ierr); */
  /*   ierr = VecAYPX(nodalFields->dszz,-1.0,stemp);CHKERRQ(ierr); */
  /*   ierr = VecCopy(stemp,nodalFields->sozz);CHKERRQ(ierr); */

  /* end cell-centered quantities */
  Vec wxyl;
  PetscScalar **wxy;
  ierr=DMCreateLocalVector(grid->da,&wxyl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(wxyl,&wxzl);CHKERRQ(ierr); */
/*   ierr=VecDuplicate(wxyl,&wyzl);CHKERRQ(ierr); */

  ierr=DMDAVecGetArray( grid->da,wxyl,&wxy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,wxzl,&wxz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,wyzl,&wyz);CHKERRQ(ierr); */

  /* compute internal basic node quantities */

  /* vorticity */
  PetscInt x1=x;
  PetscInt y1=y;
  PetscInt m1=m;
  PetscInt n1=n;
  //if(x1 == 0 && !grid->xperiodic){ x1++; m1--;}
  //if(y1 == 0){ y1++; n1--;}
  //if( !grid->xperiodic && x1+m1 >= NX) m1=NX-x1-1;
  /* restrict loop indices to interior nodes if this cpu owns cells on boundary of domain*/
  //if( y1+n1 >= NY) n1=NY-y1-1;
  
  for(jy=y1;jy<y1+n1;jy++){
    for(ix=x1;ix<x1+m1;ix++){
      PetscScalar vxm, vym, vxp, vyp;/* ghost values - these depend on type of boundary condition */

      //      PetscScalar vzul, vzur, vzll, vzlr; /* z values depend on whether we are at a boundary */

      if(!grid->xperiodic && ix == 0){/* if grid is periodic in x-direction, do not worry about left boundary */
	if( bv->mechBCLeft.type[1] == 0 ){/* prescribed velocity */
	  vym = 2.0*bv->mechBCLeft.value[1] - vy[jy][ix];
	}else if(bv->mechBCLeft.type[1] == 1){/* free slip */
	  vym = vy[jy][ix];
	}	
      } else {
	vym = vy[jy][ix-1];
      }
      vyp = vy[jy][ix];

      if(jy == 0){
	if( bv->mechBCTop.type[0] == 0){/* prescribed velocity */
	  vxm = 2.0*bv->mechBCTop.value[0] - vx[jy][ix];
	}else if(bv->mechBCTop.type[0] == 1){
	  vxm = vx[jy][ix];
	}        
      }else{
	vxm = vx[jy-1][ix];
      }
      vxp = vx[jy][ix];

      if( ix == NX-1 || jy == NY-1){
	if( ix == NX-1 && jy < NY-1 ){
	  /* implicit ghost z-velocity nodes */
	  if( bv->mechBCRight.type[2] == 0){	  
	    //vzur = 2.0*bv->mechBCRight.value[2]-vz[jy][ix];
	    //vzlr = 2.0*bv->mechBCRight.value[2]-vz[jy+1][ix];
	  } else if( bv->mechBCRight.type[2] == 1){
	    //vzur = vz[jy][ix];
	    //vzlr = vz[jy+1][ix];
	  }
	  //	  vzul = vz[jy][ix];
	  //vzll = vz[jy+1][ix];
	}
	if( jy == NY-1 ){
	  if( ix < NX-1 ){
	    if(bv->mechBCBottom.type[2] == 1){
	      //vzll = vz[jy][ix];
	      //vzlr = vz[jy][ix+1];
	    }
	    //vzul = vz[jy][ix];
	    //vzur = vz[jy][ix+1];
	  }else{
	    /* special case for lower right corner */
	    if(bv->mechBCRight.type[2] == 0){/* prescribed velocity */
	      //vzur = 2.0*bv->mechBCRight.value[2]-vz[jy][ix];
	    }else if(bv->mechBCRight.type[2] == 1){
	      //vzur = vz[jy][ix];
	    }
	    if( bv->mechBCBottom.type[2] == 1){
	      //vzlr = vzur; /* this imposes free slip on lower boundary */
	      //vzll = vz[jy][ix];
	    }
	    //vzul = vz[jy][ix];
	  }
	}
      }else{
/* 	vzul = vz[jy][ix]; */
/* 	vzur = vz[jy][ix+1]; */
/* 	vzll = vz[jy+1][ix]; */
/* 	vzlr = vz[jy+1][ix+1]; */
      }
      

      PetscScalar dx=grid->xc[ix+1]-grid->xc[ix];
      PetscScalar dy=grid->yc[jy+1]-grid->yc[jy];
      
      wxy[jy][ix] = 0.5*((vxp-vxm)/dy-(vyp-vym)/dx);
      exy[jy][ix] = 0.5*((vxp-vxm)/dy+(vyp-vym)/dx);


    }
  }
  ierr=DMDAVecRestoreArray( grid->da,wxyl,&wxy);CHKERRQ(ierr);
/*   ierr=DMDAVecRestoreArray( grid->da,wxzl,&wxz);CHKERRQ(ierr); */
/*   ierr=DMDAVecRestoreArray( grid->da,wyzl,&wyz);CHKERRQ(ierr); */
  ierr=DMDAVecRestoreArray( grid->da,exyl,&exy);CHKERRQ(ierr);
/*   ierr=DMDAVecRestoreArray( grid->da,exzl,&exz);CHKERRQ(ierr); */
/*   ierr=DMDAVecRestoreArray( grid->da,eyzl,&eyz);CHKERRQ(ierr); */
  ierr=DMLocalToGlobalBegin(grid->da,wxyl,INSERT_VALUES,nodalFields->wxy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,wxyl,INSERT_VALUES,nodalFields->wxy);CHKERRQ(ierr);
/*   ierr=DMLocalToGlobalBegin(grid->da,wxzl,INSERT_VALUES,nodalFields->wxz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalEnd(grid->da,wxzl,INSERT_VALUES,nodalFields->wxz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalBegin(grid->da,wyzl,INSERT_VALUES,nodalFields->wyz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalEnd(grid->da,wyzl,INSERT_VALUES,nodalFields->wyz);CHKERRQ(ierr); */
  ierr=DMLocalToGlobalBegin(grid->da,exyl,INSERT_VALUES,nodalFields->edotxy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,exyl,INSERT_VALUES,nodalFields->edotxy);CHKERRQ(ierr);
/*   ierr=DMLocalToGlobalBegin(grid->da,exzl,INSERT_VALUES,nodalFields->edotxz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalEnd(grid->da,exzl,INSERT_VALUES,nodalFields->edotxz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalBegin(grid->da,eyzl,INSERT_VALUES,nodalFields->edotyz);CHKERRQ(ierr); */
/*   ierr=DMLocalToGlobalEnd(grid->da,eyzl,INSERT_VALUES,nodalFields->edotyz);CHKERRQ(ierr); */
  
  /* PetscScalar X = etaS[idxnode]/(etaS[idxnode]+dt*muS[idxnode]); */
  ierr=VecCopy(nodalFields->etaS,X); CHKERRQ(ierr);
  ierr=VecAXPY(X,dt,nodalFields->muS); CHKERRQ(ierr);
  ierr=VecPointwiseDivide(X,nodalFields->etaS,X);CHKERRQ(ierr);
  ierr = VecCopy(X,X1);CHKERRQ(ierr); /* X1 = X */
  ierr = VecScale(X1,-2.0);CHKERRQ(ierr); /* X1 = X1*-2.0 */
  ierr = VecShift(X1,2.0);CHKERRQ(ierr); /* X1 = X1 + 2.0 */
  ierr = VecPointwiseMult(X1,nodalFields->etaS,X1);CHKERRQ(ierr); /* X1 = etaS.*X1 -> X1 = etaS.*(-2.0*X+2.0) */

  /* compute new sxy*/
  ierr=VecPointwiseMult(etemp,X1,nodalFields->edotxy);CHKERRQ(ierr); /* etemp = X1.*exx */
  ierr=VecPointwiseMult(stemp,X,nodalFields->soxy);  /* stemp = X.*soxx */
  ierr=VecAXPY(stemp,1.0,etemp);CHKERRQ(ierr);/* stemp = stemp+1.0*etmp */  /*dsxx = sxxnew-soxx*/
  ierr=VecCopy(nodalFields->soxy,nodalFields->dsxy);CHKERRQ(ierr); /* dsxx = soxx */
  ierr=VecAYPX(nodalFields->dsxy,-1.0,stemp);CHKERRQ(ierr);/*dsxx = sxxnew - 1.0*soxx*/
  ierr=VecCopy(stemp,nodalFields->soxy);CHKERRQ(ierr); /* soxx[i+NY*j]= sxxnew; */
#ifdef DEBUG
  {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"./loadsxy.m",&viewer);CHKERRQ(ierr);   
    ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr); 
    ierr = VecView( nodalFields->soxy, viewer);CHKERRQ(ierr);
    ierr = VecView( nodalFields->edotxy, viewer);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) X, "X");CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) X1, "X1");CHKERRQ(ierr);
    ierr = VecView(X, viewer);CHKERRQ(ierr);
    ierr = VecView(X1, viewer);CHKERRQ(ierr);

    ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
  }
#endif


/*   PetscScalar sxynew = (1-X)*2.0*etaS[idxnode]*exy[idxnode]+X*soxy[idxnode]; */
  
/*   PetscScalar sxznew = (1-X)*2.0*etaS[idxnode]*exz[idxnode]+X*soxz[idxnode]; */
/*   PetscScalar syznew = (1-X)*2.0*etaS[idxnode]*eyz[idxnode]+X*soyz[idxnode]; */
  
/*   dsxy[idxnode] = sxynew-soxy[idxnode]; */
/*   soxy[idxnode] = sxynew; */
/*   dsxz[idxnode] = sxznew-soxz[idxnode]; */
/*   soxz[idxnode] = sxznew; */
/*   dsyz[idxnode] = syznew-soyz[idxnode]; */
/*   soyz[idxnode] = syznew; */


  /* shear heating*/

  /* get all of the stress tensor components and viscosity as locally-accessible arrays */
  Vec soxxl,soyyl,soxyl,etaNl,etaSl,eiil,siil;
  
  ierr=DMCreateLocalVector(grid->da,&soxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&soyyl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(soxxl,&sozzl);CHKERRQ(ierr); */
  ierr=VecDuplicate(soxxl,&soxyl);CHKERRQ(ierr);
/*   ierr=VecDuplicate(soxxl,&soxzl);CHKERRQ(ierr); */
/*   ierr=VecDuplicate(soxxl,&soyzl);CHKERRQ(ierr); */
  ierr=VecDuplicate(soxxl,&etaNl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&etaSl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&eiil);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&siil);CHKERRQ(ierr);
  /* get local portions of all vectors with ghost data*/
  /*stress*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->sozz,INSERT_VALUES,sozzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->sozz,INSERT_VALUES,sozzl);CHKERRQ(ierr); */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxz,INSERT_VALUES,soxzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxz,INSERT_VALUES,soxzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyz,INSERT_VALUES,soyzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyz,INSERT_VALUES,soyzl);CHKERRQ(ierr); */
  /* strain rate*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxx,INSERT_VALUES,exxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxx,INSERT_VALUES,exxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotyy,INSERT_VALUES,eyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotyy,INSERT_VALUES,eyyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotzz,INSERT_VALUES,ezzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotzz,INSERT_VALUES,ezzl);CHKERRQ(ierr); */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxy,INSERT_VALUES,exyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxy,INSERT_VALUES,exyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxz,INSERT_VALUES,exzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxz,INSERT_VALUES,exzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotyz,INSERT_VALUES,eyzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotyz,INSERT_VALUES,eyzl);CHKERRQ(ierr); */
  /* viscosity*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->etaN,INSERT_VALUES,etaNl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->etaN,INSERT_VALUES,etaNl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->etaS,INSERT_VALUES,etaSl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->etaS,INSERT_VALUES,etaSl);CHKERRQ(ierr);
  /* get local arrays*/
  
  PetscScalar **soxx, **soyy, **soxy, **etaN, **etaS, **sii, **eii;  
  ierr=DMDAVecGetArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);

  ierr=DMDAVecGetArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,soxzl,&soxz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,soyzl,&soyz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,etaNl,&etaN);CHKERRQ(ierr);  
  ierr=DMDAVecGetArray( grid->da,etaSl,&etaS);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,exyl,&exy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,exzl,&exz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,eyzl,&eyz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,eiil,&eii);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,siil,&sii);CHKERRQ(ierr);

  /*  ha[idxnode] = ha[idxnode]*gy*rho[idxnode]*vy[idxnode]; */
/*   ierr=VecScale(nodalFields->ha,gy);CHKERRQ(ierr); */
/*   ierr=VecPointwiseMult(nodalFields->ha,nodalFields->ha,nodalFields->rho);CHKERRQ(ierr); */
/*   ierr=VecPointwiseMult(nodalFields->ha,nodalFields->ha,nodalFields->vy);CHKERRQ(ierr); */
/*   ierr=VecCopy(nodalFields->ha,nodalHeating);CHKERRQ(ierr); */
  ierr=VecZeroEntries(nodalHeating);CHKERRQ(ierr);

  PetscScalar **nodalHeatingA;
  Vec nodalHeatingl;
  ierr=VecDuplicate(soxxl,&nodalHeatingl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalHeating,INSERT_VALUES,nodalHeatingl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalHeating,INSERT_VALUES,nodalHeatingl);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,nodalHeatingl,&nodalHeatingA);CHKERRQ(ierr);

  /* nodalHeating is what is ultimately important in the energy equation */

  x1=x;
  y1=y;
  m1=m;
  n1=n;
  if(!grid->xperiodic && x1 == 0){ x1++; m1--;}
  if(y1 == 0){ y1++; n1--;}
  if(!grid->xperiodic && x1+m1 >= NX) m1=NX-x1-1;
  if( y1+n1 >= NY) n1=NY-y1-1;/* restrict loop indices to interior basic nodes if this cpu owns cells on boundary of domain*/


  for(jy=y1;jy<y1+n1;jy++){
    for(ix=x1;ix<x1+m1;ix++){
      PetscScalar hs;
      hs=(soxx[jy][ix]*soxx[jy][ix]/etaN[jy][ix]+soxx[jy+1][ix]*soxx[jy+1][ix]/etaN[jy+1][ix]+soxx[jy][ix+1]*soxx[jy][ix+1]/etaN[jy][ix+1]+soxx[jy+1][ix+1]*soxx[jy+1][ix+1]/etaN[jy+1][ix+1])/8.0;
      hs+=(soyy[jy][ix]*soyy[jy][ix]/etaN[jy][ix]+soyy[jy+1][ix]*soyy[jy+1][ix]/etaN[jy+1][ix]+soyy[jy][ix+1]*soyy[jy][ix+1]/etaN[jy][ix+1]+soyy[jy+1][ix+1]*soyy[jy+1][ix+1]/etaN[jy+1][ix+1])/8.0 ;
      //hs+=(sozz[jy][ix]*sozz[jy][ix]/etaN[jy][ix]+sozz[jy+1][ix]*sozz[jy+1][ix]/etaN[jy+1][ix]+sozz[jy][ix+1]*sozz[jy][ix+1]/etaN[jy][ix+1]+sozz[jy+1][ix+1]*sozz[jy+1][ix+1]/etaN[jy+1][ix+1])/8.0;
      hs+= soxy[jy][ix]*soxy[jy][ix]/etaS[jy][ix] ;
      //hs+= soxz[jy][ix]*soxz[jy][ix]/etaS[jy][ix] ;
      //hs+= soyz[jy][ix]*soyz[jy][ix]/etaS[jy][ix] ; /* factors of 8 are to do nodal averaging (1/4) and then to assume soxx = exx/(2*eta) */
      
      nodalHeatingA[jy][ix]+=hs;/* adiabatic heating is already in nodalHeating. add shear heating*/
      
      
      eii[jy][ix] =     sqrt((exx[jy][ix]*exx[jy][ix]+exx[jy+1][ix]*exx[jy+1][ix]+exx[jy][ix+1]*exx[jy][ix+1]+exx[jy+1][ix+1]*exx[jy+1][ix+1])/8.0 + \
			     (eyy[jy][ix]*eyy[jy][ix]+eyy[jy+1][ix]*eyy[jy+1][ix]+eyy[jy][ix+1]*eyy[jy][ix+1]+eyy[jy+1][ix+1]*eyy[jy+1][ix+1])/8.0 + \
			     exy[jy][ix]*exy[jy][ix]) ;
      
      sii[jy][ix] =     sqrt((soxx[jy][ix]*soxx[jy][ix]+soxx[jy+1][ix]*soxx[jy+1][ix]+soxx[jy][ix+1]*soxx[jy][ix+1]+soxx[jy+1][ix+1]*soxx[jy+1][ix+1])/8.0 + \
			     (soyy[jy][ix]*soyy[jy][ix]+soyy[jy+1][ix]*soyy[jy+1][ix]+soyy[jy][ix+1]*soyy[jy][ix+1]+soyy[jy+1][ix+1]*soyy[jy+1][ix+1])/8.0 + \
			     soxy[jy][ix]*soxy[jy][ix]);		
      
      
    }
  }
  
  /* restore local arrays */
  ierr=DMDAVecRestoreArray(grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soyyl,&soyy);CHKERRQ(ierr);
  /*   ierr=DMDAVecRestoreArray(grid->da,sozzl,&sozz);CHKERRQ(ierr); */
  ierr=DMDAVecRestoreArray(grid->da,soxyl,&soxy);CHKERRQ(ierr);
  /*   ierr=DMDAVecRestoreArray(grid->da,soxzl,&soxz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray(grid->da,soyzl,&soyz);CHKERRQ(ierr); */
  ierr=DMDAVecRestoreArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
  /*   ierr=DMDAVecRestoreArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray( grid->da,wyzl,&exy);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray( grid->da,wyzl,&exz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray( grid->da,wyzl,&eyz);CHKERRQ(ierr); */
  
  ierr=DMDAVecRestoreArray(grid->da,nodalHeatingl,&nodalHeatingA);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,siil,&sii);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,eiil,&eii);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etaNl,&etaN);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,etaSl,&etaS);CHKERRQ(ierr);
  /* scatter local arrays to global arrays*/
  ierr=DMLocalToGlobalBegin(grid->da,eiil,INSERT_VALUES,nodalFields->eii);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,eiil,INSERT_VALUES,nodalFields->eii);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,siil,INSERT_VALUES,nodalFields->sii);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,siil,INSERT_VALUES,nodalFields->sii);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,nodalHeatingl,INSERT_VALUES,nodalHeating);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,nodalHeatingl,INSERT_VALUES,nodalHeating);CHKERRQ(ierr);
  
  
  ierr=VecDestroy(&nodalHeatingl);CHKERRQ(ierr);
  ierr=VecDestroy(&wxyl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&wxzl);CHKERRQ(ierr); */
  /*   ierr=VecDestroy(&wyzl);CHKERRQ(ierr); */
  ierr=VecDestroy(&exxl);CHKERRQ(ierr);
  ierr=VecDestroy(&eyyl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&ezzl);CHKERRQ(ierr); */
  ierr=VecDestroy(&exyl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&exzl);CHKERRQ(ierr); */
  /*   ierr=VecDestroy(&eyzl);CHKERRQ(ierr); */
  ierr=VecDestroy(&eiil);CHKERRQ(ierr);
  /* stress*/
  ierr=VecDestroy(&soxxl);CHKERRQ(ierr);
  ierr=VecDestroy(&soyyl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&sozzl);CHKERRQ(ierr); */
  ierr=VecDestroy(&soxyl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&soxzl);CHKERRQ(ierr); */
  /*   ierr=VecDestroy(&soyzl);CHKERRQ(ierr); */
  ierr=VecDestroy(&siil);CHKERRQ(ierr);
  /* viscosity */
  ierr=VecDestroy(&etaNl);CHKERRQ(ierr);
  ierr=VecDestroy(&etaSl);CHKERRQ(ierr);
  
  ierr = VecDestroy(&stemp);CHKERRQ(ierr);
  ierr = VecDestroy(&etemp);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&X1);CHKERRQ(ierr);
  
  
  ierr=DMDAVecRestoreArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  /*   ierr=DMDAVecRestoreArray( grid->da,vzl,&vz);CHKERRQ(ierr); */
  ierr = VecDestroy(&vxl);CHKERRQ(ierr);
  ierr = VecDestroy(&vyl);CHKERRQ(ierr);
  //  ierr = VecDestroy(&vzl);CHKERRQ(ierr);
  PetscLogStagePop();
  PetscFunctionReturn(ierr);
}

