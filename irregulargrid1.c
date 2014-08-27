#include "fdcode.h"
#include "grid.h"
#include "math.h"

/* set up a regular grid.*/

PetscErrorCode initializeIrregularGridConstantInnerOuter(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY, Options *options){

  PetscScalar hcontrast = options->gridRatio;
  PetscErrorCode ierr=0;
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc((NX+4)*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc((NY+4)*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);
  grid->x += 2; /* increment pointer by 2 elements */
  grid->xc += 2;
  grid->y += 2;
  grid->yc += 2;


  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;
  grid->xperiodic = 0;
  grid->yperiodic = 0;

  PetscScalar f=1.2407234;
  PetscInt ix,jy;
  
  /* x gridlines */
  /* check to see if DA is periodic */
  if( options->mechBCLeft.type[0] == 2){
    if( NX%2){
      printf("Need even NX for periodic X boundary conditions\n");
      abort();
    }
    grid->xperiodic = 1;
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX ++;
  }else{
    if (!NX%2 || !NY%2 ){
      printf("Need odd number of gridlines for the mesh you've selected\n");
      abort();
    }
  }

  PetscScalar L = LX/2.0;
  PetscInt NE = (NX-1)/2;
  PetscInt Nr=NE/hcontrast;/* number of elements in refined region */
  PetscInt Nc=NE-Nr-5;
  PetscScalar a = L/(10.0+Nr+Nc*hcontrast);

  grid->x[0] = 0.0;
  for(ix=1;ix<=Nr;ix++){
    grid->x[ix] = grid->x[ix-1] + a;
  }
  for(ix=Nr+1;ix<=Nr+5;ix++){
    grid->x[ix] = grid->x[ix-1] + a*pow(f,ix-Nr);
  }
  for(ix=Nr+6;ix<=Nc+Nr+6;ix++){
    grid->x[ix] = grid->x[ix-1] + hcontrast*a;
  }
  /* mirror other half of grid */
  grid->x[NE] = LX/2;
  for(ix=(NX+1)/2;ix<NX-1;ix++){
    grid->x[ix] = LX - grid->x[NX-1 - ix];
  }
  grid->x[NX-1] = LX;
  
  /* y-gridlines */
  L = LY/2.0;
  NE = (NY-1)/2;
  Nr=NE/hcontrast;/* number of elements in refined region */
  Nc=NE-Nr-5;
  a = L/(10.0+Nr+Nc*hcontrast);

  grid->y[0] = 0.0;
  for(jy=1;jy<=Nr;jy++){
    grid->y[jy] = grid->y[jy-1] + a;
  }
  for(jy=Nr+1;jy<=Nr+5;jy++){
    grid->y[jy] = grid->y[jy-1] + a*pow(f,jy-Nr);
  }
  for(jy=Nr+6;jy<=Nc+Nr+6;jy++){
    grid->y[jy] = grid->y[jy-1] + hcontrast*a;
  }
  /* mirror other half of grid */
  grid->y[NE] = LY/2;
  for(jy=(NY+1)/2;jy<NY-1;jy++){
    grid->y[jy] = LY - grid->y[NY-1 - jy];
  }
  grid->y[NY-1] = LY;
  
  /* cell center locations*/
  for(jy=1;jy<NY;jy++){
    grid->yc[jy] = (grid->y[jy-1]+grid->y[jy])/2.0;
  }
  grid->yc[0] = - grid->yc[1];
  grid->yc[NY] = 2.0*grid->LY - grid->yc[NY-1];

  for(ix=1;ix<NX;ix++){
    grid->xc[ix] = (grid->x[ix-1]+grid->x[ix])/2.0;
  }
  grid->xc[0] = -grid->xc[1];
  grid->xc[NX] = 2.0*grid->LX - grid->xc[NX-1];


  /* check to see if DA is periodic */
  if( grid->xperiodic){
    /* periodic in x-direction */
    /* trick mesher in to thinking that there is actually one more node in x-dimension than actually present in the DA */
    NX --;
    grid->x[NX] = LX;
    grid->x[NX+1] = LX+grid->x[1];
    grid->x[-1] = grid->x[NX-1]-LX;
    grid->x[-2] = grid->x[NX-2]-LX;
    grid->xc[NX] = (LX + grid->x[NX-1])/2.0;
    grid->xc[NX+1] = grid->xc[1]+LX;
    grid->xc[0] = 0.0-(LX-grid->xc[NX]);
    grid->xc[-1] = 0.0-(LX-grid->xc[NX-1]);
  }else{
  }  


  PetscFunctionReturn(ierr);
}

PetscErrorCode initializeIrregularGridConstantIncreaseFixedRatio(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY, PetscScalar hcontrast){
  PetscErrorCode ierr=0;
  /* hcontrast is ratio of interior cell dimension to boundary cell dimension */
  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);

  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;


  if(!(NX % 2) || !(NY%2) ){/* N is even. this grid will not work */
    printf("number of gridlines must be odd\n"); fflush(stdout);
    abort();
  }
  
  PetscScalar L=LX;
  PetscScalar NN = (PetscScalar) NX;
  PetscScalar h0 = -(pow(hcontrast,1/(-3 + NN))*(-1 + pow(hcontrast,2/(-3 + NN)))*L)/(2.*(pow(hcontrast,1/(-3 + NN)) - pow(hcontrast,NN/(-3 + NN))));
  PetscInt NE = (NN-1)/2;
  PetscScalar a = pow(hcontrast,1.0/((NN-1)/2-1));
  
  PetscInt ix,jy;
  /* x gridlines */
  grid->x[0] = 0.0;
  for(ix=1;ix<(NX-1)/2;ix++){
    grid->x[ix] = grid->x[ix-1] + pow(a,(PetscScalar)(ix-1))*h0;
  }
  grid->x[(NX-1)/2] = LX/2;
  for(ix=(NX+1)/2;ix<NX-1;ix++){
    grid->x[ix] = LX - grid->x[NX-1 - ix];
  }
  grid->x[NX-1] = LX;
  
  /* y gridlines*/
  L = LY;
  NN = (PetscScalar) NY;
  h0 = -(pow(hcontrast,1/(-3 + NN))*(-1 + pow(hcontrast,2/(-3 + NN)))*L)/(2.*(pow(hcontrast,1/(-3 + NN)) - pow(hcontrast,NN/(-3 + NN))));
  NE = (NN-1)/2;
  a = pow(hcontrast,1.0/((NN-1)/2-1));

  grid->y[0] = 0.0;
  for(jy=1;jy<(NY-1)/2;jy++){
    grid->y[jy] = grid->y[jy-1] + pow(a,(PetscScalar) (jy-1))*h0;
  }
  grid->y[(NY-1)/2] = LY/2;
  for(jy=(NY+1)/2;jy<NY-1;jy++){
    grid->y[jy] = LY - grid->y[NY-1 - jy];
  }
  grid->y[NY-1] = LY;

  /* cell center locations*/
  for(jy=1;jy<NY;jy++){
    grid->yc[jy] = (grid->y[jy-1]+grid->y[jy])/2.0;
  }
  for(ix=1;ix<NX;ix++){
    grid->xc[ix] = (grid->x[ix-1]+grid->x[ix])/2.0;
  }
  PetscFunctionReturn(ierr);

}



PetscErrorCode initializeIrregularGridConstantIncrease(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY,PetscScalar xfactor, PetscScalar yfactor){
  PetscErrorCode ierr=0;

  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);

  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;

  PetscInt ix,jy;
  /* basic node gridlines*/

  if(!(NX % 2) || !(NY%2)){/* N is even. this grid will not work */
    printf("N cannot be even for this mesh!"); fflush(stdout); 
    abort();
  }
  
  grid->x[0] = 0.0;
  PetscScalar h0=LX/2*(xfactor-1)/(-1+pow(xfactor,1+(NX-1)/2));
  for(ix=1;ix<((NX-1)/2);ix++){
    grid->x[ix]=grid->x[ix-1]+pow(xfactor,ix)*h0;
  }
  grid->x[(NX-1)/2] = LX/2.0;
  grid->x[NX-1] = LX;
  for(ix=NX-2;ix>(NX-1)/2;ix--){
    grid->x[ix]=grid->x[ix+1]-pow(xfactor,(NX-1-ix))*h0;
  }

  /* y direction */
  grid->y[0] = 0.0;
  h0=LY/2*(yfactor-1)/(-1+pow(yfactor,1+(NY-1)/2));
  for(jy=1;jy<((NY-1)/2);jy++){
    grid->y[jy]=grid->y[jy-1]+pow(yfactor,jy)*h0;
  }
  grid->y[(NY-1)/2] = LY/2.0;
  grid->y[NY-1] = LY;
  for(jy=NY-2;jy>(NY-1)/2;jy--){
    grid->y[jy]=grid->y[jy+1]-pow(yfactor,(NY-1-jy))*h0;
  }
   
  /* cell center locations*/
  for(jy=1;jy<NY;jy++){
    grid->yc[jy] = (grid->y[jy-1]+grid->y[jy])/2.0;
  }
  for(ix=1;ix<NX;ix++){
    grid->xc[ix] = (grid->x[ix-1]+grid->x[ix])/2.0;
  }
  PetscFunctionReturn(ierr);
}




PetscErrorCode initializeIrregularGrid2Region(GridData *grid, PetscScalar LX, PetscScalar LY, PetscInt NX, PetscInt NY){
  /* This is a grid with constant spacing in a boundary layer region and constant spacing in the interior region */
  PetscErrorCode ierr=0;
  //  PetscScalar dx,dy,dx2,dy2,dxdy;

  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->y);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->x);CHKERRQ(ierr);
  ierr = PetscMalloc(NX*sizeof(PetscScalar),&grid->xc);CHKERRQ(ierr);
  ierr = PetscMalloc(NY*sizeof(PetscScalar),&grid->yc);CHKERRQ(ierr);

  /* put half of resolution in box boundaries */
  PetscScalar fhalf = 1.0/5.0; /* fraction of domain into which half of resolution will go */

  grid->NX = NX;
  grid->NY = NY;
  grid->LX=LX;
  grid->LY=LY;

  PetscInt ix,jy;
  /* basic node gridlines*/
  PetscInt nbx = (PetscInt) (NX/4.0);/* number of nodes in refined boundary region*/
  /*   PetscInt nintx = NX-2*nbx; */

  PetscScalar dx = fhalf/2*LX/(nbx-1);
  grid->x[0] = 0.0;
  for(ix=1;ix<nbx;ix++){
    grid->x[ix] = grid->x[ix-1]+dx;
  }
  dx = (1.0-fhalf)*LX/(NX-2*nbx+1);
  for(ix=nbx;ix<=NX-nbx;ix++){
    grid->x[ix] = grid->x[ix-1] + dx;
  }
  dx = fhalf/2*LX/(nbx-1);
  for(ix=NX-nbx+1;ix<NX;ix++){
    grid->x[ix] = grid->x[ix-1] + dx;
  }
  grid->x[NX-1] = LX;
  /* y-direction */
  PetscInt nby = (PetscInt) (NY/4.0);/* number of nodes in refined boundary region*/
  PetscScalar dy = fhalf/2*LY/(nby-1);
  grid->y[0] = 0.0;
  for(jy=1;jy<nby;jy++){
    grid->y[jy] = grid->y[jy-1]+dy;
  }
  dy = (1.0-fhalf)*LY/(NY-2*nby+1);
  for(jy=nby;jy<=NY-nby;jy++){
    grid->y[jy] = grid->y[jy-1] + dy;
  }
  dy = fhalf/2*LY/(nby-1);
  for(jy=NY-nby+1;jy<NY;jy++){
    grid->y[jy] = grid->y[jy-1] + dy;
  }
  grid->y[NY-1] = LY;

    /* cell center locations*/
  for(jy=1;jy<NY;jy++){
    grid->yc[jy] = (grid->y[jy-1]+grid->y[jy])/2.0;
  }
  for(ix=1;ix<NX;ix++){
    grid->xc[ix] = (grid->x[ix-1]+grid->x[ix])/2.0;
  }
  PetscFunctionReturn(ierr);
}

PetscErrorCode destroyGrid(GridData *grid){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  grid->x -= 2; grid->y -=2; grid->xc -= 2; grid->yc -= 2;
  ierr = PetscFree(grid->x);CHKERRQ(ierr);
  ierr = PetscFree(grid->y);CHKERRQ(ierr);
  ierr = PetscFree(grid->xc);CHKERRQ(ierr);
  ierr = PetscFree(grid->yc);CHKERRQ(ierr);
  PetscFunctionReturn(ierr);
}

/* routine to set locations of basic nodes in DA for spatial fields*/

