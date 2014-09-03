#include "fdcode.h"
#include "updateMarkerStress.h"

/* subroutine to update the marker strain rate and pressure fields*/
PetscErrorCode updateMarkerStress( GridData *grid, NodalFields *nodalFields,MarkerSet *markerset, Materials *materials){
  PetscInt m;
/*   PetscScalar LX = grid->LX; */
/*   PetscScalar LY=grid->LY; */
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  PetscErrorCode ierr;

  Marker *markers = markerset->markers;
  //  PetscScalar *p = nodalFields->p;
/*   PetscScalar *dsxx = nodalFields->dsxx; */
/*   PetscScalar *dsyy = nodalFields->dsyy; */
/*   PetscScalar *dszz = nodalFields->dszz; */
/*   PetscScalar *dsxy = nodalFields->dsxy; */
/*   PetscScalar *dsxz = nodalFields->dsxz; */
/*   PetscScalar *dsyz = nodalFields->dsyz; */
  Vec dsxxl, dsyyl, dsxyl;
  PetscScalar **dsxx, **dsyy, **dsxy;

  /* dsxx */
  ierr=DMCreateLocalVector(grid->da,&dsxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&dsyyl);CHKERRQ(ierr);
  //  ierr=VecDuplicate(dsxxl,&dszzl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxl,&dsxyl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&dsxzl);CHKERRQ(ierr);
  //ierr=VecDuplicate(dsxxl,&dsyzl);CHKERRQ(ierr);
  /* populate the local arrays with global values*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsxx,INSERT_VALUES,dsxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsxx,INSERT_VALUES,dsxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsyy,INSERT_VALUES,dsyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsyy,INSERT_VALUES,dsyyl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dszz,INSERT_VALUES,dszzl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dszz,INSERT_VALUES,dszzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsxy,INSERT_VALUES,dsxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsxy,INSERT_VALUES,dsxyl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsxz,INSERT_VALUES,dsxzl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsxz,INSERT_VALUES,dsxzl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalBegin(grid->da,nodalFields->dsyz,INSERT_VALUES,dsyzl);CHKERRQ(ierr);
  //ierr=DMGlobalToLocalEnd(grid->da,nodalFields->dsyz,INSERT_VALUES,dsyzl);CHKERRQ(ierr);
  /* get arrays */
  ierr=DMDAVecGetArray( grid->da,dsxxl,&dsxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsyyl,&dsyy);CHKERRQ(ierr);
  //ierr=DMDAVecGetArray( grid->da,dszzl,&dszz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsxyl,&dsxy);CHKERRQ(ierr);
  //ierr=DMDAVecGetArray( grid->da,dsxzl,&dsxz);CHKERRQ(ierr);
  //ierr=DMDAVecGetArray( grid->da,dsyzl,&dsyz);CHKERRQ(ierr);
  // PetscScalar *soxx = nodalFields->soxx;
  //PetscScalar *soxy = nodalFields->soxy;
  //PetscScalar *exx = nodalFields->edotxx;
  //PetscScalar *exy = nodalFields->edotxy;
  //PetscScalar *w = nodalFields->w;
  /*   PetscScalar *markerY = markers->Y; */
  /*   PetscScalar *markerX = markers->X; */
  
  PetscFunctionBegin;
  for(m=0;m<markerset->nMark;m++){
    /* is marker inside domain?*/
    if( markers[m].cellX != -1){
      /* marker is in bounds. */
      PetscInt cellX = markers[m].cellX;
      PetscInt cellY = markers[m].cellY;
      PetscInt cellXMax = grid->xperiodic ? NX : NX-2;
      /*       cellI = floor((markerY[m]-dy/2.0)/dy)+1; */
      if(cellX < cellXMax && markers[m].X >= grid->xc[cellX+1]) cellX++;
      if(cellY < NY-1 && markers[m].Y >= grid->yc[cellY+1]) cellY++;
      
      /*       cellJ = floor((markerX[m]-dx/2.0)/dx)+1; */
      /* first do fields interpolated on cell centers*/

      /* make sure that marker is assigned to a sensible cell. cellI,J=0 -> ghost*/
      if( cellY < 1){ cellY=1;} else if(cellY>NY-2){ cellY=NY-2;}/* do not let a marker be assigned to one of the ghost cells */
      cellXMax = grid->xperiodic ? NX-1 : NX-2;
      if( cellX < 1){ cellX=1;} else if(cellX>cellXMax){ cellX=cellXMax;}
      /*pressure*/
      PetscScalar mdx = (markers[m].X - grid->xc[cellX])/(grid->xc[cellX+1]-grid->xc[cellX]);
      PetscScalar mdy = (markers[m].Y - grid->yc[cellY])/(grid->yc[cellY+1]-grid->yc[cellY]);

      //      PetscInt idxc = cellI+cellJ*NY;

      PetscScalar mdsxx=0.0;
      mdsxx+=(1.0-mdx)*(1.0-mdy)*dsxx[cellY][cellX];
      mdsxx+=(1.0-mdx)*mdy*      dsxx[cellY+1][cellX];
      mdsxx+=mdx*(1.0-mdy)*      dsxx[cellY][cellX+1];
      mdsxx+=mdx*mdy*            dsxx[cellY+1][cellX+1];
      PetscScalar mdsyy=0.0;
      mdsyy+=(1.0-mdx)*(1.0-mdy)*dsyy[cellY][cellX];
      mdsyy+=(1.0-mdx)*mdy*dsyy[cellY+1][cellX];
      mdsyy+=mdx*(1.0-mdy)*dsyy[cellY][cellX+1];
      mdsyy+=mdx*mdy*dsyy[cellY+1][cellX+1];
#ifndef TEXTURE
      //PetscScalar mdszz=0.0;
      //mdszz+=(1.0-mdx)*(1.0-mdy)*dszz[cellY][cellX];
      //mdszz+=(1.0-mdx)*mdy*dszz[cellY+1][cellX];
      //mdszz+=mdx*(1.0-mdy)*dszz[cellY][cellX+1];
      //mdszz+=mdx*mdy*dszz[cellY+1][cellX+1];
#endif
      /* apply new stress*/
      markers[m].s.T11 += mdsxx;
      markers[m].s.T22 += mdsyy;
#ifndef TEXTURE
      //markers[m].s.T33 += mdszz;
#endif
      /* interpolate quantities from basic nodes*/
      cellX=markers[m].cellX;
      cellY=markers[m].cellY;
      cellXMax = grid->xperiodic ? NX-1 : NX-2;
      if( cellY < 0){ cellY=0;} else if(cellY>NY-2){ cellY=NY-2;}/* do not let a marker be assigned to one of the boundary cells*/
      if( cellX < 0){ cellX=0;} else if(cellX>cellXMax){ cellX=cellXMax;}
      mdx= (markers[m].X-grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
      mdy =(markers[m].Y-grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
      
      PetscScalar mdsxy=0.0;
      mdsxy+=(1.0-mdx)*(1.0-mdy)*dsxy[cellY][cellX];
      mdsxy+=(1.0-mdx)*mdy*dsxy[cellY+1][cellX];
      mdsxy+=mdx*(1.0-mdy)*dsxy[cellY][cellX+1];
      mdsxy+=mdx*mdy*dsxy[cellY+1][cellX+1];
      markers[m].s.T12 += mdsxy;
#ifndef TEXTURE
      //PetscScalar mdsxz=0.0;
      //mdsxz+=(1.0-mdx)*(1.0-mdy)*dsxz[cellY][cellX];
      //mdsxz+=(1.0-mdx)*mdy*dsxz[cellY+1][cellX];
      //mdsxz+=mdx*(1.0-mdy)*dsxz[cellY][cellX+1];
      //mdsxz+=mdx*mdy*dsxz[cellY+1][cellX+1];
      //markers[m].s.T13 += mdsxz;
      //PetscScalar mdsyz=0.0;
      //mdsyz+=(1.0-mdx)*(1.0-mdy)*dsyz[cellY][cellX];
      //mdsyz+=(1.0-mdx)*mdy*dsyz[cellY+1][cellX];
      //mdsyz+=mdx*(1.0-mdy)*dsyz[cellY][cellX+1];
      //mdsyz+=mdx*mdy*dsyz[cellY+1][cellX+1];
      //markers[m].s.T23 += mdsyz;
#endif
      if(isnan(mdsxy)){
	PetscInt x,y,m1,n;
	ierr=DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m1,&n,PETSC_NULL);CHKERRQ(ierr);
	printf("marker %d has nan mdsxy, (cellX,cellY)=(%d,%d),(x,y,m,n)=(%d,%d,%d,%d)\n",m,cellX,cellY,x,y,m1,n);
	fflush(stdout);
	/* 	abort(); */
      }
    }/* end if marker in bounds*/
  }/* end loop over markers*/
  
  ierr=DMDAVecRestoreArray( grid->da,dsxxl,&dsxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsyyl,&dsyy);CHKERRQ(ierr);
  //  ierr=DMDAVecRestoreArray( grid->da,dszzl,&dszz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsxyl,&dsxy);CHKERRQ(ierr);
  //  ierr=DMDAVecRestoreArray( grid->da,dsxzl,&dsxz);CHKERRQ(ierr);
  //ierr=DMDAVecRestoreArray( grid->da,dsyzl,&dsyz);CHKERRQ(ierr);
  
  ierr=VecDestroy(& dsxxl );CHKERRQ(ierr);
  ierr=VecDestroy(& dsyyl );CHKERRQ(ierr);
  //ierr=VecDestroy(& dszzl );CHKERRQ(ierr);
  ierr=VecDestroy(& dsxyl );CHKERRQ(ierr);
  //  ierr=VecDestroy(& dsxzl );CHKERRQ(ierr);
  //  ierr=VecDestroy(& dsyzl );CHKERRQ(ierr);


  PetscFunctionReturn(0);

}

