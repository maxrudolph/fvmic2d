#include "fdcode.h"
#include "nodalStressStrain.h"

/* calculate stress and strain at nodes*/
/* this version uses the anisotropic constitutive equation to calculate predicted nodal stresses */

PetscErrorCode nodalStressStrain( GridData *grid, NodalFields *nodalFields,Options *options, BoundaryValues *bv, PetscScalar dt, Vec nodalHeating, PetscScalar gy){
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  //PetscScalar dx = grid->dx;
  //PetscScalar dy = grid->dy;
  PetscInt x,y,m,n;
  PetscErrorCode ierr;
  PetscInt ix,jy;  
  
  PetscFunctionBegin;

  /* retrieve local arrays of nodal fields, including ghost node values*/
  Vec vxl,vyl,vzl,pl;
  Vec exxl,eyyl,ezzl,exyl,exzl,eyzl,VPNbl,VPNcl;

  ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  
  ierr=DMCreateLocalVector(grid->da,&vxl);
  ierr=VecDuplicate(vxl,&vyl);
  ierr=VecDuplicate(vxl,&vzl);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vx,INSERT_VALUES,vxl);CHKERRQ(ierr);
  /*y*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vy,INSERT_VALUES,vyl);CHKERRQ(ierr);
  /*z*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->vz,INSERT_VALUES,vzl);CHKERRQ(ierr);
  /*get viscoplasticity tensors*/
  ierr = DMCreateLocalVector(grid->tda,&VPNbl);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(grid->tda,&VPNcl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(grid->tda,nodalFields->VPTensorB,INSERT_VALUES,VPNbl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(grid->tda,nodalFields->VPTensorB,INSERT_VALUES,VPNbl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(grid->tda,nodalFields->VPTensorC,INSERT_VALUES,VPNcl);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(grid->tda,nodalFields->VPTensorC,INSERT_VALUES,VPNcl);CHKERRQ(ierr);

  ierr=VecDuplicate(vxl,&exxl); CHKERRQ(ierr);
  ierr=VecDuplicate(vxl,&eyyl); CHKERRQ(ierr);
  ierr=VecDuplicate(vxl,&ezzl); CHKERRQ(ierr);
  ierr=VecDuplicate(vxl,&exyl); CHKERRQ(ierr);
  ierr=VecDuplicate(vxl,&exzl); CHKERRQ(ierr);
  ierr=VecDuplicate(vxl,&eyzl); CHKERRQ(ierr);

  PetscScalar **vx, **vy, **vz, **exx, **eyy, **ezz, **exy, **exz, **eyz;
  Tensor22 **VPNb, **VPNc;  
  ierr=DMDAVecGetArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,vzl,&vz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->tda,VPNbl,&VPNb);CHKERRQ(ierr);/* get viscoplasticity N */
  ierr=DMDAVecGetArray( grid->tda,VPNcl,&VPNc);CHKERRQ(ierr);

  ierr=DMDAVecGetArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,eyyl,&eyy);CHKERRQ(ierr); 
  ierr=DMDAVecGetArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); 


  /* get all stress tensor components */

  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){/* note that these loops start at 1 to account for ghost cells*/
      if( (!grid->xperiodic && ix==0) || jy == 0){
	/* account for ghost cells */
	exx[jy][ix] = 0.0;
	eyy[jy][ix] = 0.0;
	ezz[jy][ix] = 0.0;
      }else{
	/* this is a cell-centered quantity*/
	//      PetscScalar divv = (vx[i-1+j*NY]-vx[i-1+(j-1)*NY])/dx +(vy[i+NY*(j-1)]-vy[i-1+NY*(j-1)])/dy;
	PetscScalar dx=grid->x[ix]-grid->x[ix-1];
	PetscScalar dy=grid->y[jy]-grid->y[jy-1];
	PetscScalar divv = (vx[jy-1][ix]-vx[jy-1][ix-1])/dx + (vy[jy][ix-1]-vy[jy-1][ix-1])/dy;/* velocity divergence == trace of dev. strain rate*/
	
	exx[jy][ix]=(vx[jy-1][ix]-vx[jy-1][ix-1])/dx - divv/3.0; /* should be deviatoric strain rate*/
	eyy[jy][ix]=(vy[jy][ix-1]-vy[jy-1][ix-1])/dy - divv/3.0; /* should be deviatoric strain rate*/
	ezz[jy][ix]=-divv/3.0; /* deviatoric strain rate in z direction. ezz = 0 but e'zz is nonzero if there is divergence in the plane*/
      }
    }
  }
  /* scatter deviatoric strain rates back to global vectors*/
  ierr=DMDAVecRestoreArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,eyyl,&eyy);CHKERRQ(ierr); 
  ierr=DMDAVecRestoreArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); 
  ierr=DMLocalToGlobalBegin(grid->da,exxl,INSERT_VALUES,nodalFields->edotxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,exxl,INSERT_VALUES,nodalFields->edotxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,eyyl,INSERT_VALUES,nodalFields->edotyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,eyyl,INSERT_VALUES,nodalFields->edotyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,ezzl,INSERT_VALUES,nodalFields->edotzz);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,ezzl,INSERT_VALUES,nodalFields->edotzz);CHKERRQ(ierr);

  /* end cell-centered quantities */

  /* compute internal basic node quantities */
  Vec wxyl;/* wxzl, wyzl;*/
  PetscScalar **wxy; /*, **wxz, **wyz;*/
  ierr=DMCreateLocalVector(grid->da,&wxyl);CHKERRQ(ierr);
  /*   ierr=VecDuplicate(wxyl,&wxzl);CHKERRQ(ierr); */
  /*   ierr=VecDuplicate(wxyl,&wyzl);CHKERRQ(ierr); */
  
  ierr=DMDAVecGetArray( grid->da,wxyl,&wxy);CHKERRQ(ierr);
  /*   ierr=DMDAVecGetArray( grid->da,wxzl,&wxz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecGetArray( grid->da,wyzl,&wyz);CHKERRQ(ierr); */
  
  ierr=DMDAVecGetArray( grid->da,exyl,&exy);CHKERRQ(ierr); 
  /*   ierr=DMDAVecGetArray( grid->da,exzl,&exz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecGetArray( grid->da,eyzl,&eyz);CHKERRQ(ierr); */
  
  PetscInt x1=x;
  PetscInt y1=y;
  PetscInt m1=m;
  PetscInt n1=n;

  for(jy=y1;jy<y1+n1;jy++){
    for(ix=x1;ix<x1+m1;ix++){
      
      PetscScalar vxm, vym, vxp, vyp;/* ghost values - these depend on type of boundary condition */
      if(!grid->xperiodic && ix == 0){
	if( bv->mechBCLeft.type[1] == 0 ){/* prescribed velocity */
	  vym = 2.0*bv->mechBCLeft.value[1] - vy[jy][ix];
	}else if(bv->mechBCLeft.type[1] == 1){
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

      PetscScalar dx=grid->xc[ix+1]-grid->xc[ix];
      PetscScalar dy=grid->yc[jy+1]-grid->yc[jy];

      wxy[jy][ix] = -0.5*((vyp-vym)/dx-(vxp-vxm)/dy);
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
  /*   ierr=DALocalToGlobal(grid->da,wxzl,INSERT_VALUES,nodalFields->wxz);CHKERRQ(ierr); */
  /*   ierr=DALocalToGlobal(grid->da,wyzl,INSERT_VALUES,nodalFields->wyz);CHKERRQ(ierr); */
  ierr=DMLocalToGlobalBegin(grid->da,exyl,INSERT_VALUES,nodalFields->edotxy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd(grid->da,exyl,INSERT_VALUES,nodalFields->edotxy);CHKERRQ(ierr);
  /*   ierr=DALocalToGlobal(grid->da,exzl,INSERT_VALUES,nodalFields->edotxz);CHKERRQ(ierr); */
  /*   ierr=DALocalToGlobal(grid->da,eyzl,INSERT_VALUES,nodalFields->edotyz);CHKERRQ(ierr); */

  /* NODAL STRESS */

  /* get local arrays for strain-rate with ghost values */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxx,INSERT_VALUES,exxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxx,INSERT_VALUES,exxl);CHKERRQ(ierr);

  /*   ierr=DAGlobalToLocalBegin(grid->da,nodalFields->edotyy,INSERT_VALUES,eyyl);CHKERRQ(ierr); */
  /*   ierr=DAGlobalToLocalEnd(grid->da,nodalFields->edotyy,INSERT_VALUES,eyyl);CHKERRQ(ierr); */
  
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxy,INSERT_VALUES,exyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxy,INSERT_VALUES,exyl);CHKERRQ(ierr);

  ierr=DMDAVecGetArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  /*   ierr=DMDAVecGetArray( grid->da,eyyl,&eyy);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,exyl,&exy);CHKERRQ(ierr); 
  
  /* get all of the stress tensor components and viscosity as locally-accessible arrays */
  Vec soxxl,soyyl,sozzl,soxyl,soxzl,soyzl,dsxxl,dsxyl,dsyyl;
  ierr=DMCreateLocalVector(grid->da,&soxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&soyyl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&sozzl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&soxyl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&soxzl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&soyzl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&dsxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&dsyyl);CHKERRQ(ierr);
  ierr=VecDuplicate(soxxl,&dsxyl);CHKERRQ(ierr);
  /* get local portions of all vectors with ghost data*/
  /*stress*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr); 
  /*   ierr=DAGlobalToLocalBegin(grid->da,nodalFields->sozz,INSERT_VALUES,sozzl);CHKERRQ(ierr); */
  /*   ierr=DAGlobalToLocalEnd(grid->da,nodalFields->sozz,INSERT_VALUES,sozzl);CHKERRQ(ierr); */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
  /*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxz,INSERT_VALUES,soxzl);CHKERRQ(ierr); */
  /*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxz,INSERT_VALUES,soxzl);CHKERRQ(ierr); */
  /*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyz,INSERT_VALUES,soyzl);CHKERRQ(ierr); */
  /*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyz,INSERT_VALUES,soyzl);CHKERRQ(ierr); */
  /* get local arrays for stress tensor components*/  
  PetscScalar **soxx, **soyy, **sozz, **soxy, **soxz, **soyz, **sii, **eii, **dsxx, **dsyy, **dsxy;  
  ierr=DMDAVecGetArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);
  /*   ierr=DMDAVecGetArray( grid->da,soyyl,&soyy);CHKERRQ(ierr); */
  /*   ierr=DMDAVecGetArray( grid->da,sozzl,&sozz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
  /*   ierr=DMDAVecGetArray( grid->da,soxzl,&soxz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecGetArray( grid->da,soyzl,&soyz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,dsxxl,&dsxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsyyl,&dsyy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsxyl,&dsxy);CHKERRQ(ierr);

  /* calculate cell-centered stress components using assumed constitutive relationship */

  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){/* note that these loops start at 1 to account for ghost cells*/
      if( (!grid->xperiodic && ix == 0) || jy == 0){
	dsxx[jy][ix] = 0.0;
	dsyy[jy][ix] = 0.0;
	soxx[jy][ix] = 0.0;
	soyy[jy][ix] = 0.0;
      }else{
	/* sxx = N11*exx + N12*exy */
	       
	PetscScalar exyc = (exy[jy-1][ix-1] + exy[jy-1][ix] + exy[jy][ix-1] + exy[jy][ix])/4.0;
	dsxx[jy][ix] = VPNc[jy][ix].T11*exx[jy][ix] + VPNc[jy][ix].T12*exyc - soxx[jy][ix];
	
	soxx[jy][ix] += dsxx[jy][ix];
	dsyy[jy][ix] = -dsxx[jy][ix];
	soyy[jy][ix] -= dsxx[jy][ix];
	
	if(isnan(soxx[jy][ix])) printf("NAN in soxx\n");
      }
    }
  }


  for(jy=y;jy<y+n;jy++){
    for(ix=x;ix<x+m;ix++){
      /* calculate shear stresses sxy, dsxy */
      /* account for ghost boundary nodes along left and top */       
	PetscScalar exxul, exxur, exxll, exxlr;
	PetscInt ixul, ixur, ixll, ixlr;
	PetscInt jyul, jyur, jyll, jylr;

	ixul = ix;
	ixur = ix+1;
	ixll = ix;
	ixlr = ix+1;
	jyul = jy;
	jyur = jy;
	jyll = jy+1;
	jylr = jy+1;
	if( ix == 0 && !grid->xperiodic){/* doing this is exactly correct for free slip boundaries */
	  if( bv->mechBCLeft.type[1] == 1){
	    ixul++;
	    ixll++;
	  }
	} else if(ix == NX-1 && !grid->xperiodic ){
	  if( bv->mechBCRight.type[1] == 1){
	    ixur--;
	    ixlr--;
	  }
	}
	if( jy == 0){
	  if( bv->mechBCTop.type[0] == 1){
	    jyul++;
	    jyur++;
	  }
	}else if(jy == NY-1){
	  if( bv->mechBCBottom.type[0] == 1){
	    jyll--;
	    jylr--;
	  }
	}
	exxul = exx[jyul][ixul];
	exxur = exx[jyur][ixur];
	exxll = exx[jyll][ixll];
	exxlr = exx[jylr][ixlr];
	
	PetscScalar exxb = (exxul+exxur+exxll+exxlr)/4.0;
	dsxy[jy][ix] = VPNb[jy][ix].T21*exxb + VPNb[jy][ix].T22*exy[jy][ix]-soxy[jy][ix];
	soxy[jy][ix] += dsxy[jy][ix];
	
    }
  }
  
  ierr=DMDAVecRestoreArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  /*   ierr=DMDAVecRestoreArray( grid->da,eyyl,&eyy);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); */
  ierr=DMDAVecRestoreArray( grid->da,exyl,&exy);CHKERRQ(ierr);
  /*   ierr=DMDAVecRestoreArray( grid->da,exzl,&exz);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray( grid->da,eyzl,&eyz);CHKERRQ(ierr); */

  ierr=DMDAVecRestoreArray(grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soyyl,&soyy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,soxyl,&soxy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,dsxxl,&dsxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,dsyyl,&dsyy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray(grid->da,dsxyl,&dsxy);CHKERRQ(ierr);

  ierr=DMLocalToGlobalBegin(grid->da,soxxl,INSERT_VALUES,nodalFields->soxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,soxxl,INSERT_VALUES,nodalFields->soxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,soyyl,INSERT_VALUES,nodalFields->soyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,soyyl,INSERT_VALUES,nodalFields->soyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,soxyl,INSERT_VALUES,nodalFields->soxy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,soxyl,INSERT_VALUES,nodalFields->soxy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,dsxxl,INSERT_VALUES,nodalFields->dsxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,dsxxl,INSERT_VALUES,nodalFields->dsxx);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,dsyyl,INSERT_VALUES,nodalFields->dsyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,dsyyl,INSERT_VALUES,nodalFields->dsyy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalBegin(grid->da,dsxyl,INSERT_VALUES,nodalFields->dsxy);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,dsxyl,INSERT_VALUES,nodalFields->dsxy);CHKERRQ(ierr);

  /* scatter global arrays back to local arrays */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);/* stress */
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr); 
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr); 
  /* strain rate */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxx,INSERT_VALUES,exxl);CHKERRQ(ierr);/* stress */
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxx,INSERT_VALUES,exxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotyy,INSERT_VALUES,eyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotyy,INSERT_VALUES,eyyl);CHKERRQ(ierr); 
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->edotxy,INSERT_VALUES,exyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->edotxy,INSERT_VALUES,exyl);CHKERRQ(ierr); 
  /* get the arrays */
  ierr=DMDAVecGetArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,exyl,&exy);CHKERRQ(ierr);



  ierr = VecZeroEntries(nodalHeating);CHKERRQ(ierr);
  PetscScalar **nodalHeatingA;
  Vec nodalHeatingl;
  ierr=DMCreateLocalVector(grid->da,&nodalHeatingl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalHeating,INSERT_VALUES,nodalHeatingl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,  nodalHeating,INSERT_VALUES,nodalHeatingl);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,nodalHeatingl,&nodalHeatingA);CHKERRQ(ierr);

  /* nodalHeating is what is ultimately important in the energy equation */
/*   PetscInt x1,y1,m1,n1; */
  x1=x;
  y1=y;
  m1=m;
  n1=n;
  if(x1 == 0){ x1++; m1--;}
  if(y1 == 0){ y1++; n1--;}
  if( x1+m1 >= NX) m1=NX-x1-1;
  if( y1+n1 >= NY) n1=NY-y1-1; /* restrict loop indices to interior basic nodes if this cpu owns cells on boundary of domain*/
  
  for(jy=y1;jy<y1+n1;jy++){
    for(ix=x1;ix<x1+m1;ix++){
      PetscScalar hs;
      /* hs=(soxx[jy][ix]*soxx[jy][ix]/etaN[jy][ix]+soxx[jy+1][ix]*soxx[jy+1][ix]/etaN[jy+1][ix]+soxx[jy][ix+1]*soxx[jy][ix+1]/etaN[jy][ix+1]+soxx[jy+1][ix+1]*soxx[jy+1][ix+1]/etaN[jy+1][ix+1])/8.0; */
/*       hs+=(soyy[jy][ix]*soyy[jy][ix]/etaN[jy][ix]+soyy[jy+1][ix]*soyy[jy+1][ix]/etaN[jy+1][ix]+soyy[jy][ix+1]*soyy[jy][ix+1]/etaN[jy][ix+1]+soyy[jy+1][ix+1]*soyy[jy+1][ix+1]/etaN[jy+1][ix+1])/8.0 ; */
/*       hs+=(sozz[jy][ix]*sozz[jy][ix]/etaN[jy][ix]+sozz[jy+1][ix]*sozz[jy+1][ix]/etaN[jy+1][ix]+sozz[jy][ix+1]*sozz[jy][ix+1]/etaN[jy][ix+1]+sozz[jy+1][ix+1]*sozz[jy+1][ix+1]/etaN[jy+1][ix+1])/8.0; */
      hs =  (soxx[jy][ix]*exx[jy][ix]+soxx[jy+1][ix]*exx[jy+1][ix]+soxx[jy][ix+1]*exx[jy][ix+1]+soxx[jy+1][ix+1]*exx[jy+1][ix+1])/4;
      hs += (soyy[jy][ix]*eyy[jy][ix]+soyy[jy+1][ix]*eyy[jy+1][ix]+soyy[jy][ix+1]*eyy[jy][ix+1]+soyy[jy+1][ix+1]*eyy[jy+1][ix+1])/4;
      
      /*       hs+= soxy[jy][ix]*soxy[jy][ix]/etaS[jy][ix] ; */
      hs += 2.0*soxy[jy][ix]*exy[jy][ix]; /* factor of 2.0 to account for soxy and soyx */

      /* in 2D code, exz = eyz = 0 */

      
      //nodalHeatingA[jy][ix]+=hs;/* adiabatic heating is already in nodalHeating. add shear heating*/
      
      nodalHeatingA[jy][ix]+=hs;
      
      /*eii[jy][ix] =     sqrt((exx[jy][ix]*exx[jy][ix]+exx[jy+1][ix]*exx[jy+1][ix]+exx[jy][ix+1]*exx[jy][ix+1]+exx[jy+1][ix+1]*exx[jy+1][ix+1])/8.0 + \
	(eyy[jy][ix]*eyy[jy][ix]+eyy[jy+1][ix]*eyy[jy+1][ix]+eyy[jy][ix+1]*eyy[jy][ix+1]+eyy[jy+1][ix+1]*eyy[jy+1][ix+1])/8.0 + \
	(ezz[jy][ix]*ezz[jy][ix]+ezz[jy+1][ix]*ezz[jy+1][ix]+ezz[jy][ix+1]*ezz[jy][ix+1]+ezz[jy+1][ix+1]*ezz[jy+1][ix+1])/8.0 + \
	exy[jy][ix]*exy[jy][ix] +					\
	exz[jy][ix]*exz[jy][ix] +					\
	eyz[jy][ix]*eyz[jy][ix]);
	sii[jy][ix] =     sqrt((soxx[jy][ix]*soxx[jy][ix]+soxx[jy+1][ix]*soxx[jy+1][ix]+soxx[jy][ix+1]*soxx[jy][ix+1]+soxx[jy+1][ix+1]*soxx[jy+1][ix+1])/8.0 + \
	(soyy[jy][ix]*soyy[jy][ix]+soyy[jy+1][ix]*soyy[jy+1][ix]+soyy[jy][ix+1]*soyy[jy][ix+1]+soyy[jy+1][ix+1]*soyy[jy+1][ix+1])/8.0 + \
	(sozz[jy][ix]*sozz[jy][ix]+sozz[jy+1][ix]*sozz[jy+1][ix]+sozz[jy][ix+1]*sozz[jy][ix+1]+sozz[jy+1][ix+1]*sozz[jy+1][ix+1])/8.0 + \
	soxy[jy][ix]*soxy[jy][ix] +					\
	soxz[jy][ix]*soxz[jy][ix] +					\
	soyz[jy][ix]*soyz[jy][ix]);
      */
    }
  }
  /* restore local arrays */
  ierr=DMDAVecRestoreArray( grid->da,exxl,&exx);CHKERRQ(ierr); 
  ierr=DMDAVecRestoreArray( grid->da,eyyl,&eyy);CHKERRQ(ierr); 
  ierr=DMDAVecRestoreArray( grid->da,exyl,&exy);CHKERRQ(ierr);   
  ierr=DMDAVecRestoreArray( grid->da,soxxl,&soxx);CHKERRQ(ierr); 
  ierr=DMDAVecRestoreArray( grid->da,soyyl,&soyy);CHKERRQ(ierr); 
  ierr=DMDAVecRestoreArray( grid->da,soxyl,&soxy);CHKERRQ(ierr); 
  ierr=DMDAVecRestoreArray(grid->da,nodalHeatingl,&nodalHeatingA);CHKERRQ(ierr);

  ierr=DMDAVecRestoreArray( grid->tda,VPNbl,&VPNb);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->tda,VPNcl,&VPNc);CHKERRQ(ierr);
  ierr=VecDestroy(&VPNbl);CHKERRQ(ierr);
  ierr=VecDestroy(&VPNcl);CHKERRQ(ierr);


  /*   ierr=DMDAVecRestoreArray(grid->da,siil,&sii);CHKERRQ(ierr); */
  /*   ierr=DMDAVecRestoreArray(grid->da,eiil,&eii);CHKERRQ(ierr); */
  
  /* scatter local arrays to global arrays*/
  /*   ierr=DALocalToGlobal(grid->da,eiil,INSERT_VALUES,nodalFields->eii);CHKERRQ(ierr); */
  /*   ierr=DALocalToGlobal(grid->da,siil,INSERT_VALUES,nodalFields->sii);CHKERRQ(ierr); */
  ierr=DMLocalToGlobalBegin(grid->da,nodalHeatingl,INSERT_VALUES,nodalHeating);CHKERRQ(ierr);
  ierr=DMLocalToGlobalEnd  (grid->da,nodalHeatingl,INSERT_VALUES,nodalHeating);CHKERRQ(ierr);
  ierr=VecDestroy(&nodalHeatingl); CHKERRQ(ierr);
  ierr=VecDestroy(&wxyl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&wxzl);CHKERRQ(ierr);  */
  /*   ierr=VecDestroy(&wyzl);CHKERRQ(ierr);  */
  ierr=VecDestroy(&exxl);CHKERRQ(ierr);
  ierr=VecDestroy(&eyyl);CHKERRQ(ierr);
  ierr=VecDestroy(&ezzl);CHKERRQ(ierr);
  ierr=VecDestroy(&exyl);CHKERRQ(ierr);
  ierr=VecDestroy(&exzl);CHKERRQ(ierr);
  ierr=VecDestroy(&eyzl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&eiil);CHKERRQ(ierr); */
  /* stress*/
  ierr=VecDestroy(&soxxl);CHKERRQ(ierr);
  ierr=VecDestroy(&soyyl);CHKERRQ(ierr);
  ierr=VecDestroy(&sozzl);CHKERRQ(ierr);
  ierr=VecDestroy(&soxyl);CHKERRQ(ierr);
  ierr=VecDestroy(&soxzl);CHKERRQ(ierr);
  ierr=VecDestroy(&soyzl);CHKERRQ(ierr);
  /*   ierr=VecDestroy(&siil);CHKERRQ(ierr); */

  ierr=VecDestroy(&dsxxl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsyyl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsxyl);CHKERRQ(ierr);

  ierr=DMDAVecRestoreArray( grid->da,vxl,&vx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vyl,&vy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,vzl,&vz);CHKERRQ(ierr);
  ierr = VecDestroy(&vxl);CHKERRQ(ierr);
  ierr = VecDestroy(&vyl);CHKERRQ(ierr);
  ierr = VecDestroy(&vzl);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}

