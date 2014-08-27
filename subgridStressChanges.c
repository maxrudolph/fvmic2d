#include "fdcode.h"
#include "subgridStressChanges.h"

PetscErrorCode subgridStressChanges(GridData *grid, NodalFields *nodalFields, MarkerSet *markerset, Materials *materials, PetscScalar dt, PetscScalar dsubgrid){
  PetscInt NX = grid->NX;
  PetscInt NY = grid->NY;
  Marker *markers = markerset->markers;

  PetscErrorCode ierr;
  /*   ierr = PetscMalloc */
  PetscFunctionBegin;
  
  /* allocate ds??n. these should exist both as global and local vecs with the DA structure. */
  Vec dsxxng, dsyyng, dszzng, dsxyng, dsxzng, dsyzng, wtetang, wtetasg;
  Vec soxxl, soyyl, sozzl, soxyl, soxzl, soyzl;
  Vec dsxxnl, dsyynl, dszznl, dsxynl, dsxznl, dsyznl, wtetanl, wtetasl;
  PetscScalar **soxx, **soyy, **sozz, **soxy, **soxz, **soyz, **wtetas, **wtetan;
  PetscScalar **dsxxn, **dsyyn, **dszzn, **dsxyn, **dsxzn, **dsyzn;
  ierr=DMCreateGlobalVector(grid->da,&dsxxng); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxng,&dsyyng); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxng,&dszzng); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxng,&dsxyng); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxng,&dsxzng); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxng,&dsyzng); CHKERRQ(ierr);
  /* create local vectors*/
  ierr=DMCreateLocalVector(grid->da,&dsxxnl); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&dsyynl); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&dszznl); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&dsxynl); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&dsxznl); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&dsyznl); CHKERRQ(ierr);
  /* nodal weights*/
  ierr=VecDuplicate(dsxxng, &wtetang); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl, &wtetanl); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxng, &wtetasg); CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl, &wtetasl); CHKERRQ(ierr);

  /*get stress - this code copied and pasted from updateMarkerStrainPressure*/
  /*stress*/
  ierr=VecDuplicate(dsxxnl,&soxxl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&soyyl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&sozzl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&soxyl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&soxzl);CHKERRQ(ierr);
  ierr=VecDuplicate(dsxxnl,&soyzl);CHKERRQ(ierr);
  /*stress*/
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxx,INSERT_VALUES,soxxl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyy,INSERT_VALUES,soyyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->sozz,INSERT_VALUES,sozzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->sozz,INSERT_VALUES,sozzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxy,INSERT_VALUES,soxyl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soxz,INSERT_VALUES,soxzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soxz,INSERT_VALUES,soxzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalBegin(grid->da,nodalFields->soyz,INSERT_VALUES,soyzl);CHKERRQ(ierr);
  ierr=DMGlobalToLocalEnd(grid->da,nodalFields->soyz,INSERT_VALUES,soyzl);CHKERRQ(ierr);
  /* stress*/
  ierr=DMDAVecGetArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,sozzl,&sozz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soxzl,&soxz);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,soyzl,&soyz);CHKERRQ(ierr);
  /* nodal stress changes */
  ierr=DMDAVecGetArray( grid->da,dsxxnl,&dsxxn);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsyynl,&dsyyn);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dszznl,&dszzn);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsxynl,&dsxyn);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsxznl,&dsxzn);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,dsyznl,&dsyzn);CHKERRQ(ierr);
  /* weights */
  ierr=DMDAVecGetArray( grid->da,wtetanl,&wtetan);CHKERRQ(ierr);
  ierr=DMDAVecGetArray( grid->da,wtetasl,&wtetas);CHKERRQ(ierr);

  PetscInt m;
  for(m=0;m<markerset->nMark;m++){
    if( markers[m].cellX != -1){
      PetscScalar sdm=markers[m].eta / markers[m].mu;
      /*materials->materialMu[(PetscInt) markers[m].Mat];*/
      PetscScalar sdif=-dsubgrid*dt/sdm;
      if(sdif < -30.0) sdif = -30.0;
      sdif=1.0-exp(sdif);

      /* calculate internal basic cell that I belong to*/
      PetscInt cellY = markers[m].cellY;
      PetscInt cellX = markers[m].cellX;
 
      if(cellX<0) cellX=0;
      if(cellX>NX-2) cellX=NX-2;
      if(cellY<0) cellY=0;
      if(cellY>NY-2) cellY=NY-2;
      PetscScalar mdx=(markers[m].X - grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
      PetscScalar mdy=(markers[m].Y - grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
      PetscScalar sxym=0.0;
      /*       PetscInt idxc = cellI + cellJ*NY; */
      sxym+=(1.0-mdx)*(1.0-mdy)*soxy[cellY][cellX];
      sxym+=(1.0-mdx)*mdy*soxy[cellY+1][cellX];
      sxym+=mdx*(1.0-mdy)*soxy[cellY][cellX+1];
      sxym+=mdx*mdy*soxy[cellY+1][cellX+1];
      PetscScalar dsxym=sxym - markers[m].s.T12;//xy;
      /* sxz */
      PetscScalar sxzm=0.0;
      sxzm+=(1.0-mdx)*(1.0-mdy)*soxz[cellY][cellX];
      sxzm+=(1.0-mdx)*mdy*soxz[cellY+1][cellX];
      sxzm+=mdx*(1.0-mdy)*soxz[cellY][cellX+1];
      sxzm+=mdx*mdy*soxz[cellY+1][cellX+1];

      PetscScalar dsxzm=sxzm - markers[m].s.T13;//xz;
      /* syz */
      PetscScalar syzm=0.0;
      syzm+=(1.0-mdx)*(1.0-mdy)*soyz[cellY][cellX];
      syzm+=(1.0-mdx)*mdy*soyz[cellY+1][cellX];
      syzm+=mdx*(1.0-mdy)*soyz[cellY][cellX+1];
      syzm+=mdx*mdy*soyz[cellY+1][cellX+1];
      PetscScalar dsyzm=syzm - markers[m].s.T23;//yz;
      
      dsxym*=sdif;
      dsxzm*=sdif;
      dsyzm*=sdif;
      markers[m].s.T12 += dsxym;
      markers[m].s.T13 += dsxzm;
      markers[m].s.T23 += dsyzm;
      
      
      /* now interpolate change in stress to surrounding nodes*/
      cellY = markers[m].cellY;
      cellX = markers[m].cellX;
      mdx=(markers[m].X - grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
      mdy=(markers[m].Y - grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
      //      idxc = cellI + cellJ*NY;
      if(mdx<=0.5 && mdy<=0.5){
/* 	dsxyn[idxc]=dsxyn[idxc]+(1.0-mdx)*(1.0-mdy)*dsxym; */
	dsxyn[cellY][cellX] = dsxyn[cellY][cellX] + (1.0-mdx)*(1.0-mdy)*dsxym;
/* 	dsxzn[idxc]=dsxzn[idxc]+(1.0-mdx)*(1.0-mdy)*dsxzm; */
	dsxzn[cellY][cellX] = dsxzn[cellY][cellX] + (1.0-mdx)*(1.0-mdy)*dsxzm;
/* 	dsyzn[idxc]=dsyzn[idxc]+(1.0-mdx)*(1.0-mdy)*dsyzm; */
	dsyzn[cellY][cellX] = dsyzn[cellY][cellX] + (1.0-mdx)*(1.0-mdy)*dsyzm;
/* 	wtetas[idxc]=wtetas[idxc]+(1.0-mdx)*(1.0-mdy); */
	wtetas[cellY][cellX] = wtetas[cellY][cellX]+(1.0-mdx)*(1.0-mdy);
      }
      if(mdx<=0.5 && mdy>=0.5){
	/* 	dsxyn[idxc+1]=dsxyn[idxc+1]+(1.0-mdx)*mdy*dsxym; */
	dsxyn[cellY+1][cellX] = dsxyn[cellY+1][cellX] + (1.0-mdx)*(mdy)*dsxym;
	dsxzn[cellY+1][cellX] = dsxzn[cellY+1][cellX] + (1.0-mdx)*(mdy)*dsxzm;
	dsyzn[cellY+1][cellX] = dsyzn[cellY+1][cellX] + (1.0-mdx)*(mdy)*dsyzm;
	wtetas[cellY+1][cellX] = wtetas[cellY+1][cellX]+(1.0-mdx)*(mdy);

/* 	dsxzn[idxc+1]=dsxzn[idxc+1]+(1.0-mdx)*mdy*dsxzm; */
/* 	dsyzn[idxc+1]=dsyzn[idxc+1]+(1.0-mdx)*mdy*dsyzm; */
/* 	wtetas[idxc+1]=wtetas[idxc+1]+(1.0-mdx)*mdy; */
      }
      if(mdx>=0.5 && mdy<=0.5){
	dsxyn[cellY][cellX+1] = dsxyn[cellY][cellX+1] + (mdx)*(1.0-mdy)*dsxym;
	dsxzn[cellY][cellX+1] = dsxzn[cellY][cellX+1] + (mdx)*(1.0-mdy)*dsxzm;
	dsyzn[cellY][cellX+1] = dsyzn[cellY][cellX+1] + (mdx)*(1.0-mdy)*dsyzm;
	wtetas[cellY][cellX+1] = wtetas[cellY][cellX+1]+(mdx)*(1.0-mdy);

/* 	dsxyn[idxc+NY]=dsxyn[idxc+NY]+mdx*(1.0-mdy)*dsxym; */
/* 	dsxzn[idxc+NY]=dsxzn[idxc+NY]+mdx*(1.0-mdy)*dsxzm; */
/* 	dsyzn[idxc+NY]=dsyzn[idxc+NY]+mdx*(1.0-mdy)*dsyzm; */
/* 	wtetas[idxc+NY]=wtetas[idxc+NY]+mdx*(1.0-mdy); */
      }
      if(mdx>=0.5 && mdy>=0.5){
	dsxyn[cellY+1][cellX+1] = dsxyn[cellY+1][cellX+1] + (mdx)*(mdy)*dsxym;
	dsxzn[cellY+1][cellX+1] = dsxzn[cellY+1][cellX+1] + (mdx)*(mdy)*dsxzm;
	dsyzn[cellY+1][cellX+1] = dsyzn[cellY+1][cellX+1] + (mdx)*(mdy)*dsyzm;
	wtetas[cellY+1][cellX+1] = wtetas[cellY+1][cellX+1]+(mdx)*(mdy);
/* 	dsxyn[idxc+NY+1]=dsxyn[idxc+NY+1]+mdx*mdy*dsxym; */
/* 	dsxzn[idxc+NY+1]=dsxzn[idxc+NY+1]+mdx*mdy*dsxzm; */
/* 	dsyzn[idxc+NY+1]=dsyzn[idxc+NY+1]+mdx*mdy*dsyzm; */
/* 	wtetas[idxc+NY+1]=wtetas[idxc+NY+1]+mdx*mdy; */
      }
     

      /* calculate pressure cell that I belong to*/
      cellX = markers[m].cellX;
      cellY = markers[m].cellY;
      /*       cellI = floor((markerY[m]-dy/2.0)/dy)+1; */
      if(cellX < NX-1 && markers[m].X >= grid->xc[cellX+1]) cellX++;
      if(cellY < NY-1 && markers[m].Y >= grid->yc[cellY+1]) cellY++;
      
      /* make sure that marker is assigned to a sensible cell. cellI,J=0 -> ghost*/
      if( cellY < 1){ cellY=1;} else if(cellY>NY-2){ cellY=NY-2;}/* do not let a marker be assigned to one of the ghost cells */
      if( cellX < 1){ cellX=1;} else if(cellX>NX-2){ cellX=NX-2;}
      /*pressure*/
      mdx = (markers[m].X - grid->xc[cellX])/(grid->xc[cellX+1]-grid->xc[cellX]);
      mdy = (markers[m].Y - grid->yc[cellY])/(grid->yc[cellY+1]-grid->yc[cellY]);
    
      //      idxc = cellI + cellJ*NY;
      PetscScalar sxxm=0.0;
      sxxm+=(1.0-mdx)*(1.0-mdy)*soxx[cellY][cellX];
      sxxm+=(1.0-mdx)*mdy*soxx[cellY+1][cellX];
      sxxm+=mdx*(1.0-mdy)*soxx[cellY][cellX+1];
      sxxm+=mdx*mdy*soxx[cellY+1][cellX+1];
      PetscScalar dsxxm = sxxm-markers[m].s.T11;//xx;
      dsxxm*=sdif;
      markers[m].s.T11+=dsxxm;
      /*syy*/
      PetscScalar syym=0.0;
      syym+=(1.0-mdx)*(1.0-mdy)*soyy[cellY][cellX];
      syym+=(1.0-mdx)*mdy*soyy[cellY+1][cellX];
      syym+=mdx*(1.0-mdy)*soyy[cellY][cellX+1];
      syym+=mdx*mdy*soyy[cellY+1][cellX+1];

      PetscScalar dsyym = syym-markers[m].s.T22;//yy;
      dsyym*=sdif;
      markers[m].s.T22+=dsyym;
      /*szz*/
      PetscScalar szzm=0.0;
      szzm+=(1.0-mdx)*(1.0-mdy)*sozz[cellY][cellX];
      szzm+=(1.0-mdx)*mdy*sozz[cellY+1][cellX];
      szzm+=mdx*(1.0-mdy)*sozz[cellY][cellX+1];
      szzm+=mdx*mdy*sozz[cellY+1][cellX+1];

      PetscScalar dszzm = szzm-markers[m].s.T33;//zz;
      dszzm*=sdif;
      markers[m].s.T33+=dszzm;

      /* now calculate contribution of each marker's subgrid change to the pressure node*/
      cellX = markers[m].cellX;
      cellY = markers[m].cellY;

      /*       idxc=cellI + cellJ*NY; */
      mdy = (markers[m].Y-grid->y[cellY])/(grid->y[cellY+1]-grid->y[cellY]);
      mdx = (markers[m].X-grid->x[cellX])/(grid->x[cellX+1]-grid->x[cellX]);
           
      dsxxn[cellY+1][cellX+1]+=(1.0-fabs(0.5-mdx))*(1.0-fabs(0.5-mdy))*dsxxm;
      dsyyn[cellY+1][cellX+1]+=(1.0-fabs(0.5-mdx))*(1.0-fabs(0.5-mdy))*dsyym;
      dszzn[cellY+1][cellX+1]+=(1.0-fabs(0.5-mdx))*(1.0-fabs(0.5-mdy))*dszzm;
      wtetan[cellY+1][cellX+1]+=(1.0-fabs(0.5-mdx))*(1.0-fabs(0.5-mdy));

    }/* end if marker in domain */
  }/*end loop over markers*/

  /* normalize the stress changes at nodes*/
  /* restore work vectors*/
  ierr=DMDAVecRestoreArray( grid->da,soxxl,&soxx);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,soyyl,&soyy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,sozzl,&sozz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,soxyl,&soxy);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,soxzl,&soxz);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,soyzl,&soyz);CHKERRQ(ierr);
  /* nodal stress changes */
  ierr=DMDAVecRestoreArray( grid->da,dsxxnl,&dsxxn);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsyynl,&dsyyn);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dszznl,&dszzn);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsxynl,&dsxyn);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsxznl,&dsxzn);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,dsyznl,&dsyzn);CHKERRQ(ierr);
  /* weights */
  ierr=DMDAVecRestoreArray( grid->da,wtetanl,&wtetan);CHKERRQ(ierr);
  ierr=DMDAVecRestoreArray( grid->da,wtetasl,&wtetas);CHKERRQ(ierr);

  /* scatter dsxyn and weights to global vectors */
  ierr = DMLocalToGlobalBegin(grid->da,wtetanl,ADD_VALUES,wtetang); CHKERRQ(ierr); /* these always assume ADD_VALUES*/
  ierr = DMLocalToGlobalEnd(grid->da,wtetanl,ADD_VALUES,wtetang); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,wtetasl,ADD_VALUES,wtetasg); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,wtetasl,ADD_VALUES,wtetasg); CHKERRQ(ierr); 
  /* ds??n */
  ierr = DMLocalToGlobalBegin(grid->da,dsxxnl,ADD_VALUES,dsxxng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,dsxxnl,ADD_VALUES,dsxxng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,dsyynl,ADD_VALUES,dsyyng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,dsyynl,ADD_VALUES,dsyyng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,dszznl,ADD_VALUES,dszzng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,dszznl,ADD_VALUES,dszzng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,dsxynl,ADD_VALUES,dsxyng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,dsxynl,ADD_VALUES,dsxyng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,dsxznl,ADD_VALUES,dsxzng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,dsxznl,ADD_VALUES,dsxzng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalBegin(grid->da,dsyznl,ADD_VALUES,dsyzng); CHKERRQ(ierr); 
  ierr = DMLocalToGlobalEnd(grid->da,dsyznl,ADD_VALUES,dsyzng); CHKERRQ(ierr); 

  ierr=VecPointwiseDivide( dsxyng,dsxyng,wtetasg);CHKERRQ(ierr);
  ierr=VecPointwiseDivide( dsxzng,dsxzng,wtetasg);CHKERRQ(ierr);
  ierr=VecPointwiseDivide( dsyzng,dsyzng,wtetasg);CHKERRQ(ierr);
  ierr=VecPointwiseDivide( dsxxng,dsxxng,wtetang);CHKERRQ(ierr);
  ierr=VecPointwiseDivide( dsyyng,dsyyng,wtetang);CHKERRQ(ierr);
  ierr=VecPointwiseDivide( dszzng,dszzng,wtetang);CHKERRQ(ierr);

  ierr=VecAXPY( nodalFields->dsxx, -1.0, dsxxng); CHKERRQ(ierr);
  ierr=VecAXPY( nodalFields->dsyy, -1.0, dsyyng); CHKERRQ(ierr);
  ierr=VecAXPY( nodalFields->dszz, -1.0, dszzng); CHKERRQ(ierr);
  ierr=VecAXPY( nodalFields->dsxy, -1.0, dsxyng); CHKERRQ(ierr);
  ierr=VecAXPY( nodalFields->dsxz, -1.0, dsxzng); CHKERRQ(ierr);
  ierr=VecAXPY( nodalFields->dsyz, -1.0, dsyzng); CHKERRQ(ierr);

  /* global stress changes*/
  ierr=VecDestroy(&dsxxng);CHKERRQ(ierr);
  ierr=VecDestroy(&dsyyng);CHKERRQ(ierr);
  ierr=VecDestroy(&dszzng);CHKERRQ(ierr);
  ierr=VecDestroy(&dsxyng);CHKERRQ(ierr);
  ierr=VecDestroy(&dsxzng);CHKERRQ(ierr);
  ierr=VecDestroy(&dsyzng);CHKERRQ(ierr);
  /* local stress changes*/
  ierr=VecDestroy(&dsxxnl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsyynl);CHKERRQ(ierr);
  ierr=VecDestroy(&dszznl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsxynl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsxznl);CHKERRQ(ierr);
  ierr=VecDestroy(&dsyznl);CHKERRQ(ierr);
  /* stress tensor */
  ierr=VecDestroy(&soxxl);CHKERRQ(ierr);
  ierr=VecDestroy(&soyyl);CHKERRQ(ierr);
  ierr=VecDestroy(&sozzl);CHKERRQ(ierr);
  ierr=VecDestroy(&soxyl);CHKERRQ(ierr);
  ierr=VecDestroy(&soxzl);CHKERRQ(ierr);
  ierr=VecDestroy(&soyzl);CHKERRQ(ierr);
  /* weights*/
  ierr=VecDestroy(&wtetanl);CHKERRQ(ierr);
  ierr=VecDestroy(&wtetang);CHKERRQ(ierr);
  ierr=VecDestroy(&wtetasl);CHKERRQ(ierr);
  ierr=VecDestroy(&wtetasg);CHKERRQ(ierr);

  PetscFunctionReturn(ierr);
}
