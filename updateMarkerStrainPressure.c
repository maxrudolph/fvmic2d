#include "fdcode.h"
#include "updateMarkerStrainPressure.h"
//#include "texture.h"
#include "updateDamageViscosity.h"
#include "viscosity.h"
#include "profile.h"

/* subroutine to update the marker strain rate and pressure fields*/
PetscErrorCode updateMarkerStrainPressure( GridData *grid, NodalFields *nodalFields,MarkerSet *markerset, Materials *materials, Options *options,PetscScalar dt){
  PetscErrorCode ierr = 0;
/*   PetscScalar LX = grid->LX; */
/*   PetscScalar LY=grid->LY; */
/*  PetscScalar dx = grid->dx;*/
/*  PetscScalar dy = grid->dy;*/
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  Vec dsxxl, dsyyl, dsxyl, soxxl,soyyl,soxyl,exxl,eyyl,exyl,wxyl,pl;
  PetscScalar **dsxx, **dsyy, **dsxy, **soxx,**soyy,**soxy,**exx,**eyy,**exy,**wxy,**p;
  Marker *markers = markerset->markers;

  PetscFunctionBegin;
  setLogStage( LOG_MARK_STR_PRESS );
  /* allocate all of the local vectors*/
  /* stress changes*/
  ierr=DMCreateLocalVector(grid->da,&dsxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&dsyyl);CHKERRQ(ierr);
  //  ierr=VecDuplicate(dsxxl,&dszzl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&dsxyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&dsxzl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&dsyzl);CHKERRQ(ierr);
  /*stress*/
  ierr=VecDuplicate(dsxxl,&soxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&soyyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&sozzl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&soxyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&soxzl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&soyzl);CHKERRQ(ierr);
  /*strain rate*/
  ierr=VecDuplicate(dsxxl,&exxl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&eyyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&ezzl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&exyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&exzl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&eyzl);CHKERRQ(ierr);
  /*vorticity*/
  ierr=VecDuplicate(dsxxl,&wxyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&wxzl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&wyzl);CHKERRQ(ierr);
  /*pressure*/
  ierr=VecDuplicate(dsxxl,&pl);CHKERRQ(ierr);
  
  /* populate the local arrays with global values*/
  /* stress changes */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsxx,INSERT_VALUES,dsxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsxx,INSERT_VALUES,dsxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsyy,INSERT_VALUES,dsyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsyy,INSERT_VALUES,dsyyl);CHKERRQ(ierr);
  /*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dszz,INSERT_VALUES,dszzl);CHKERRQ(ierr); */
  /*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dszz,INSERT_VALUES,dszzl);CHKERRQ(ierr); */
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsxy,INSERT_VALUES,dsxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsxy,INSERT_VALUES,dsxyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsxz,INSERT_VALUES,dsxzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsxz,INSERT_VALUES,dsxzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsyz,INSERT_VALUES,dsyzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsyz,INSERT_VALUES,dsyzl);CHKERRQ(ierr); */
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
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->p,INSERT_VALUES,pl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->p,INSERT_VALUES,pl);CHKERRQ(ierr);
  /* vorticity*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->wxy,INSERT_VALUES,wxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->wxy,INSERT_VALUES,wxyl);CHKERRQ(ierr);
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->wxz,INSERT_VALUES,wxzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->wxz,INSERT_VALUES,wxzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalBegin(grid->da,nodalFields->wyz,INSERT_VALUES,wyzl);CHKERRQ(ierr); */
/*   ierr=DMGlobalToLocalEnd(grid->da,nodalFields->wyz,INSERT_VALUES,wyzl);CHKERRQ(ierr); */
  /* now get local work vectors for all of these fields*/
  ierr=DMDAVecGetArray( grid->da,dsxxl,&dsxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsyyl,&dsyy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,dszzl,&dszz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,dsxyl,&dsxy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,dsxzl,&dsxz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,dsyzl,&dsyz);CHKERRQ(ierr); */
  /* stress*/
  ierr=DMDAVecGetArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,sozzl,&sozz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,soxzl,&soxz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,soyzl,&soyz);CHKERRQ(ierr); */
  /*strain rate*/
  ierr=DMDAVecGetArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,ezzl,&ezz);CHKERRQ(ierr); */
  ierr=DMDAVecGetArray( grid->da,exyl,&exy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,exzl,&exz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,eyzl,&eyz);CHKERRQ(ierr); */
  /*vorticity*/
  ierr=DMDAVecGetArray( grid->da,wxyl,&wxy);CHKERRQ(ierr);
/*   ierr=DMDAVecGetArray( grid->da,wxzl,&wxz);CHKERRQ(ierr); */
/*   ierr=DMDAVecGetArray( grid->da,wyzl,&wyz);CHKERRQ(ierr); */
  /* pressure*/
  ierr=DMDAVecGetArray( grid->da,pl,&p);CHKERRQ(ierr);

  /* project all fields onto markers*/
  PetscInt m;
  for(m=0;m<markerset->nMark;m++){
    if( markers[m].cellX != -1){
      /* marker is in bounds. */
      PetscInt cellX = markers[m].cellX;
      PetscInt cellY = markers[m].cellY;
      /*       cellI = floor((markerY[m]-dy/2.0)/dy)+1; */
      PetscInt cellXMax = grid->xperiodic ? NX : NX-2; /* periodic case is NX, NOT NX-1 */
      if(cellX < cellXMax && markers[m].X >= grid->xc[cellX+1]) cellX++;
      if(cellY < NY-1 && markers[m].Y >= grid->yc[cellY+1]) cellY++;
      
      /*       cellJ = floor((markerX[m]-dx/2.0)/dx)+1; */
      /* first do fields interpolated on cell centers*/
      
      /* make sure that marker is assigned to a sensible cell. cellI,J=0 -> ghost*/
      if( cellY < 1){ cellY=1;} else if(cellY>NY-2){ cellY=NY-2;}/* do not let a marker be assigned to one of the ghost cells */
      cellXMax = grid->xperiodic ? NX-1 : NX - 2;
      if( cellX < 1){ cellX=1;} else if(cellX>cellXMax){ cellX=cellXMax;}
      /*pressure*/
      PetscScalar mdx = (markers[m].X - grid->xc[cellX])/(grid->xc[cellX+1]-grid->xc[cellX]);
      PetscScalar mdy = (markers[m].Y - grid->yc[cellY])/(grid->yc[cellY+1]-grid->yc[cellY]);
      
      PetscScalar pm=0.0;
      //PetscInt idxc = cellI+cellJ*NY;
      pm+=(1.0-mdx)*(1.0-mdy)*p[cellY][cellX];
      pm+=(1.0-mdx)*mdy*p[cellY+1][cellX];
      pm+=mdx*(1.0-mdy)*p[cellY][cellX+1];
      pm+=mdx*mdy*p[cellY+1][cellX+1];
      markers[m].p = pm;
      /* exx */
      PetscScalar mexx=0.0;
      mexx+=(1.0-mdx)*(1.0-mdy)*exx[cellY][cellX];
      mexx+=(1.0-mdx)*mdy*exx[cellY+1][cellX];
      mexx+=mdx*(1.0-mdy)*exx[cellY][cellX+1];
      mexx+=mdx*mdy*exx[cellY+1][cellX+1];
      /* eyy */
      PetscScalar meyy=0.0;
      meyy+=(1.0-mdx)*(1.0-mdy)*eyy[cellY][cellX];
      meyy+=(1.0-mdx)*mdy*eyy[cellY+1][cellX];
      meyy+=mdx*(1.0-mdy)*eyy[cellY][cellX+1];
      meyy+=mdx*mdy*eyy[cellY+1][cellX+1];
      /*       ezz */
/*       PetscScalar mezz=0.0; */
/*       mezz+=(1.0-mdx)*(1.0-mdy)*ezz[cellY][cellX]; */
/*       mezz+=(1.0-mdx)*mdy*ezz[cellY+1][cellX]; */
/*       mezz+=mdx*(1.0-mdy)*ezz[cellY][cellX+1]; */
/*       mezz+=mdx*mdy*ezz[cellY+1][cellX+1]; */
      PetscScalar msxx=0.0;
      msxx+=(1.0-mdx)*(1.0-mdy)*soxx[cellY][cellX];
      msxx+=(1.0-mdx)*mdy*soxx[cellY+1][cellX];
      msxx+=mdx*(1.0-mdy)*soxx[cellY][cellX+1];
      msxx+=mdx*mdy*soxx[cellY+1][cellX+1];
      PetscScalar msyy=0.0;
      msyy+=(1.0-mdx)*(1.0-mdy)*soyy[cellY][cellX];
      msyy+=(1.0-mdx)*mdy*soyy[cellY+1][cellX];
      msyy+=mdx*(1.0-mdy)*soyy[cellY][cellX+1];
      msyy+=mdx*mdy*soyy[cellY+1][cellX+1];
/*       PetscScalar mszz=0.0; */
/*       mszz+=(1.0-mdx)*(1.0-mdy)*sozz[cellY][cellX]; */
/*       mszz+=(1.0-mdx)*mdy*sozz[cellY+1][cellX]; */
/*       mszz+=mdx*(1.0-mdy)*sozz[cellY][cellX+1]; */
/*       mszz+=mdx*mdy*sozz[cellY+1][cellX+1]; */
      PetscScalar mdsxx=0.0;
      mdsxx+=(1.0-mdx)*(1.0-mdy)*dsxx[cellY][cellX];
      mdsxx+=(1.0-mdx)*mdy*dsxx[cellY+1][cellX];
      mdsxx+=mdx*(1.0-mdy)*dsxx[cellY][cellX+1];
      mdsxx+=mdx*mdy*dsxx[cellY+1][cellX+1];
      PetscScalar mdsyy=0.0;
      mdsyy+=(1.0-mdx)*(1.0-mdy)*dsyy[cellY][cellX];
      mdsyy+=(1.0-mdx)*mdy*dsyy[cellY+1][cellX];
      mdsyy+=mdx*(1.0-mdy)*dsyy[cellY][cellX+1];
      mdsyy+=mdx*mdy*dsyy[cellY+1][cellX+1];
/*       PetscScalar mdszz=0.0; */
/*       mdszz+=(1.0-mdx)*(1.0-mdy)*dszz[cellY][cellX]; */
/*       mdszz+=(1.0-mdx)*mdy*dszz[cellY+1][cellX]; */
/*       mdszz+=mdx*(1.0-mdy)*dszz[cellY][cellX+1]; */
/*       mdszz+=mdx*mdy*dszz[cellY+1][cellX+1]; */
      /* now do quantities defined on internal basic nodes*/
      cellX=markers[m].cellX;
      cellY=markers[m].cellY;
      
      cellXMax = grid->xperiodic ? NX-1 : NX-2;
      if( cellY < 0){ cellY=0;} else if(cellY>NY-2){ cellY=NY-2;}/* do not let a marker be assigned to one of the boundary cells*/
      if( cellX < 0){ cellX=0;} else if(cellX>cellXMax){ cellX=cellXMax;}
      mdx= (markers[m].X-grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
      mdy =(markers[m].Y-grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);

      /* vorticity*/
      PetscScalar mwxy=0.0;
      mwxy+=(1.0-mdx)*(1.0-mdy)*wxy[cellY][cellX];
      mwxy+=(1.0-mdx)*mdy*wxy[cellY+1][cellX];
      mwxy+=mdx*(1.0-mdy)*wxy[cellY][cellX+1];
      mwxy+=mdx*mdy*wxy[cellY+1][cellX+1];
      markers[m].wxy = mwxy;/* spin*/
      //PetscScalar mwxz=0.0;
      //mwxz+=(1.0-mdx)*(1.0-mdy)*wxz[cellY][cellX];
      //mwxz+=(1.0-mdx)*mdy*wxz[cellY+1][cellX];
      //mwxz+=mdx*(1.0-mdy)*wxz[cellY][cellX+1];
      //mwxz+=mdx*mdy*wxz[cellY+1][cellX+1];
      //      markers[m].wxy = mwxz;/* spin*/

      //PetscScalar mwyz=0.0;
      //mwyz+=(1.0-mdx)*(1.0-mdy)*wyz[cellY][cellX];
      //mwyz+=(1.0-mdx)*mdy*wyz[cellY+1][cellX];
      //mwyz+=mdx*(1.0-mdy)*wyz[cellY][cellX+1];
      //mwyz+=mdx*mdy*wyz[cellY+1][cellX+1];
      //markers[m].wyz = mwyz;/* spin*/
      
      PetscScalar mexy=0.0;
      mexy+=(1.0-mdx)*(1.0-mdy)*exy[cellY][cellX];
      mexy+=(1.0-mdx)*mdy*exy[cellY+1][cellX];
      mexy+=mdx*(1.0-mdy)*exy[cellY][cellX+1];
      mexy+=mdx*mdy*exy[cellY+1][cellX+1];
      
      //PetscScalar mexz=0.0;
      //mexz+=(1.0-mdx)*(1.0-mdy)*exz[cellY][cellX];
      //mexz+=(1.0-mdx)*mdy*exz[cellY+1][cellX];
      //mexz+=mdx*(1.0-mdy)*exz[cellY][cellX+1];
      //mexz+=mdx*mdy*exz[cellY+1][cellX+1];
      
      //PetscScalar meyz=0.0;
      //meyz+=(1.0-mdx)*(1.0-mdy)*eyz[cellY][cellX];
      //meyz+=(1.0-mdx)*mdy*eyz[cellY+1][cellX];
      //meyz+=mdx*(1.0-mdy)*eyz[cellY][cellX+1];
      //meyz+=mdx*mdy*eyz[cellY+1][cellX+1];
      
      /* marker stress */
      PetscScalar msxy=0.0;
      msxy+=(1.0-mdx)*(1.0-mdy)*soxy[cellY][cellX];
      msxy+=(1.0-mdx)*mdy*soxy[cellY+1][cellX];
      msxy+=mdx*(1.0-mdy)*soxy[cellY][cellX+1];
      msxy+=mdx*mdy*soxy[cellY+1][cellX+1];
      //PetscScalar msxz=0.0;
      //msxz+=(1.0-mdx)*(1.0-mdy)*soxz[cellY][cellX];
      //msxz+=(1.0-mdx)*mdy*soxz[cellY+1][cellX];
      //msxz+=mdx*(1.0-mdy)*soxz[cellY][cellX+1];
      //msxz+=mdx*mdy*soxz[cellY+1][cellX+1];
      //PetscScalar msyz=0.0;
      //msyz+=(1.0-mdx)*(1.0-mdy)*soyz[cellY][cellX];
      //msyz+=(1.0-mdx)*mdy*soyz[cellY+1][cellX];
      //msyz+=mdx*(1.0-mdy)*soyz[cellY][cellX+1];
      //msyz+=mdx*mdy*soyz[cellY+1][cellX+1];
      
      /* mds** */
      PetscScalar mdsxy=0.0;
      mdsxy+=(1.0-mdx)*(1.0-mdy)*dsxy[cellY][cellX];
      mdsxy+=(1.0-mdx)*mdy*dsxy[cellY+1][cellX];
      mdsxy+=mdx*(1.0-mdy)*dsxy[cellY][cellX+1];
      mdsxy+=mdx*mdy*dsxy[cellY+1][cellX+1];
      
      /*       PetscScalar mdsxz=0.0; */
      /*       mdsxz+=(1.0-mdx)*(1.0-mdy)*dsxz[cellY][cellX]; */
      /*       mdsxz+=(1.0-mdx)*mdy*dsxz[cellY+1][cellX]; */
      /*       mdsxz+=mdx*(1.0-mdy)*dsxz[cellY][cellX+1]; */
      /*       mdsxz+=mdx*mdy*dsxz[cellY+1][cellX+1]; */
      /*       PetscScalar mdsyz=0.0; */
      /*       mdsyz+=(1.0-mdx)*(1.0-mdy)*dsyz[cellY][cellX]; */
      /*       mdsyz+=(1.0-mdx)*mdy*dsyz[cellY+1][cellX]; */
      /*       mdsyz+=mdx*(1.0-mdy)*dsyz[cellY][cellX+1]; */
      /*       mdsyz+=mdx*mdy*dsyz[cellY+1][cellX+1]; */
      
      /* use gerya's 'method 1' using nodal stress, nodal stress changes, and constit. rel.*/
#ifndef TEXTURE
      /* if this marker has damage, calculate a new damage rate and density-change rate using last (stresses+dstress) and pressure interpolated from nodes */

      markers[m].e.T11=msxx/2.0/markers[m].eta + mdsxx/2.0/dt/markers[m].mu;
      markers[m].e.T22=msyy/2.0/markers[m].eta + mdsyy/2.0/dt/markers[m].mu;
      //markers[m].e.T33=mszz/2.0/markers[m].eta + mdszz/2.0/dt/markers[m].mu;
      markers[m].e.T12=msxy/2.0/markers[m].eta + mdsxy/2.0/dt/markers[m].mu;
      //markers[m].e.T13=msxz/2.0/markers[m].eta + mdsxz/2.0/dt/markers[m].mu;
      //markers[m].e.T23=msyz/2.0/markers[m].eta + mdsyz/2.0/dt/markers[m].mu;
      markers[m].strainRateResidual = sqrt( pow(markers[m].e.T11 - mexx , 2.0) + pow( markers[m].e.T12 - mexy, 2.0));
      if( isnan(markers[m].strainRateResidual) ){
	printf("marker %d has nan strain rate residual, eta=%e\n",m,markers[m].eta);
	fflush(stdout);  	  
	abort();
      }
#else
      /* use nodal stresses to re-compute viscoplasticity tensor (constitutive eqn) */
      /* NO ELASTICITY!! */

      /* re-calculate M using constitutive equation */
      /* backup marker stress state */
      Tensor33s moldstress = markers[m].s;
      /* assign marker nodal stresses */
      markers[m].s.T11 = msxx;
      markers[m].s.T22 = -msxx;/* note no stress changees - there is no elasticity in the viscous case */
      markers[m].s.T12 = msxy;

      ierr=formViscoplasticMNewtonMarker( &(markers[m]), materials, options);CHKERRQ(ierr);
      PetscScalar (*Ms)[6] = markers[m].texture.M; /* symmetric viscoplasticity tensor */
      /* restore old nodal stresses */
      markers[m].s = moldstress;
      msxz = markers[m].s.T13;
      msyz = markers[m].s.T23;
      //      mszz = markers[m].s.T33;

      markers[m].e.T11 = Ms[1-1][1-1]*msxx - Ms[1-1][2-1]*msxx + Ms[1-1][3-1]*mszz + 2.0*(Ms[1-1][4-1]*msyz + Ms[1-1][5-1]*msxz + Ms[1-1][6-1]*msxy);
      markers[m].e.T22 = Ms[2-1][1-1]*msxx - Ms[2-1][2-1]*msxx + Ms[2-1][3-1]*mszz + 2.0*(Ms[2-1][4-1]*msyz + Ms[2-1][5-1]*msxz + Ms[2-1][6-1]*msxy);
      markers[m].e.T33 = Ms[3-1][1-1]*msxx - Ms[3-1][2-1]*msxx + Ms[3-1][3-1]*mszz + 2.0*(Ms[3-1][4-1]*msyz + Ms[3-1][5-1]*msxz + Ms[3-1][6-1]*msxy);     
      markers[m].e.T13 = Ms[4-1][1-1]*msxx - Ms[4-1][2-1]*msxx + Ms[4-1][3-1]*mszz + 2.0*(Ms[4-1][4-1]*msyz + Ms[4-1][5-1]*msxz + Ms[4-1][6-1]*msxy);
      markers[m].e.T23 = Ms[5-1][1-1]*msxx - Ms[5-1][2-1]*msxx + Ms[5-1][3-1]*mszz + 2.0*(Ms[5-1][4-1]*msyz + Ms[5-1][5-1]*msxz + Ms[5-1][6-1]*msxy);
      markers[m].e.T12 = Ms[6-1][1-1]*msxx - Ms[6-1][2-1]*msxx + Ms[6-1][3-1]*mszz + 2.0*(Ms[6-1][4-1]*msyz + Ms[6-1][5-1]*msxz + Ms[6-1][6-1]*msxy);
      /* calculate residual between macroscopic strain rate and marker strain rate calculated using assumed CR */
      markers[m].strainRateResidual = sqrt( pow(markers[m].e.T11 - mexx , 2.0) + pow( markers[m].e.T12 - mexy, 2.0));
#endif                  
    } else {
      /* marker is out of bounds*/
    }
  }
  /* now restore local work vectors for all of these fields*/
  ierr=DMDAVecRestoreArray( grid->da,dsxxl,&dsxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsyyl,&dsyy);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,dszzl,&dszz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsxyl,&dsxy);CHKERRQ(ierr);
  //  ierr=DMDAVecRestoreArray( grid->da,dsxzl,&dsxz);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,dsyzl,&dsyz);CHKERRQ(ierr);
  /* stress*/
  ierr=DMDAVecRestoreArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,sozzl,&sozz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,soxzl,&soxz);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,soyzl,&soyz);CHKERRQ(ierr);
  /*strain rate*/
  ierr=DMDAVecRestoreArray( grid->da,exxl,&exx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,eyyl,&eyy);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,ezzl,&ezz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,exyl,&exy);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,exzl,&exz);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,eyzl,&eyz);CHKERRQ(ierr);
  /*vorticity*/
  ierr=DMDAVecRestoreArray( grid->da,wxyl,&wxy);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,wxzl,&wxz);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,wyzl,&wyz);CHKERRQ(ierr);
  /* pressure*/
  ierr=DMDAVecRestoreArray( grid->da,pl,&p);CHKERRQ(ierr);

  /* destroy all of the local vectors*/
  /* stress changes*/
  ierr = VecDestroy(&dsxxl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsyyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&dszzl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsxyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&dsxzl);CHKERRQ(ierr);
  //ierr=VecDestroy(&dsyzl);CHKERRQ(ierr);
  /*stress*/
  ierr=VecDestroy(&soxxl);CHKERRQ(ierr);
  ierr=VecDestroy(&soyyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&sozzl);CHKERRQ(ierr);
  ierr=VecDestroy(&soxyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&soxzl);CHKERRQ(ierr);
  //ierr=VecDestroy(&soyzl);CHKERRQ(ierr);
  /*strain rate*/
  ierr=VecDestroy(&exxl);CHKERRQ(ierr);
  ierr=VecDestroy(&eyyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&ezzl);CHKERRQ(ierr);
  ierr=VecDestroy(&exyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&exzl);CHKERRQ(ierr);
  //ierr=VecDestroy(&eyzl);CHKERRQ(ierr);
  /*vorticity*/
  ierr=VecDestroy(&wxyl);CHKERRQ(ierr);
  //ierr=VecDestroy(&wxzl);CHKERRQ(ierr);
  //ierr=VecDestroy(&wyzl);CHKERRQ(ierr);
  /*pressure*/
  ierr=VecDestroy(&pl);CHKERRQ(ierr);
  
  PetscLogStagePop();
  PetscFunctionReturn(ierr);
}
