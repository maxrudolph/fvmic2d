#include "fdcode.h"
#include "benchmarkInitialConditions.h"
#include "viscosity.h"
#include "math.h"

PetscErrorCode initialConditionsConvectionBarr( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* This is the initial condition from Amy Barr's PhD thesis, equation 2.9 */
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar D = grid->LY;
  const PetscScalar lam = 2.0*grid->LX; /* perturbation wavelength */
  printf("Perturbation wavelength %e\n",lam);
  printf("Tperturb %e\n",options->Tperturb);
  for(m=0;m<markerset->nMark;m++){
    PetscScalar TTop = options->thermalBCTop.value[0];
    PetscScalar TBtm = options->thermalBCBottom.value[0];
    /* note that D essentially falls out of cos term. This is because Barr fixes xmax = lam/(2*D) */	   
    if( grid->xperiodic){/* for periodic grid, enforce smoothness across lateral boundary */
      markers[m].T = TTop + (TBtm-TTop)*markers[m].Y/grid->LY + options->Tperturb * sin(2.0*M_PI*(markers[m].X/grid->LX))*sin( markers[m].Y*M_PI/D);
    }else{
      markers[m].T = TTop + (TBtm-TTop)*markers[m].Y/grid->LY + options->Tperturb * cos(2.0*M_PI*D/lam* (markers[m].X/grid->LX)*lam/2.0/D)*sin( markers[m].Y*M_PI/D);
    }
    markers[m].Mat = 0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}

PetscErrorCode initialConditionsVanKeken( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* This is the initial condition from Amy Barr's PhD thesis, equation 2.9 */
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar D = grid->LY;
  const PetscScalar lam = 2.0*grid->LX; /* perturbation wavelength */
  printf("Perturbation wavelength %e\n",lam);
  printf("Tperturb %e\n",options->Tperturb);
  for(m=0;m<markerset->nMark;m++){
    PetscScalar TTop = options->thermalBCTop.value[0];
    PetscScalar TBtm = options->thermalBCBottom.value[0];
    /* note that D essentially falls out of cos term. This is because Barr fixes xmax = lam/(2*D) */	   
    if( grid->xperiodic){/* for periodic grid, enforce smoothness across lateral boundary */
      markers[m].T = TTop + (TBtm-TTop)*markers[m].Y/grid->LY + options->Tperturb * sin(2.0*M_PI*(markers[m].X/grid->LX))*sin( markers[m].Y*M_PI/D);
    }else{
      markers[m].T = TTop + (TBtm-TTop)*markers[m].Y/grid->LY + options->Tperturb * cos(2.0*M_PI*D/lam* (markers[m].X/grid->LX)*lam/2.0/D)*sin( markers[m].Y*M_PI/D);
    }
    markers[m].Mat = 0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}


PetscErrorCode initialConditionsConvectionTuran( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* This is an initial condition for Turan et al. 2010 J. Non-New. Fl. Mech. */
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar D = grid->LY;
  const PetscScalar lam = 2.0*grid->LX; /* perturbation wavelength */
  printf("Perturbation wavelength %e\n",lam);
  printf("Tperturb %e\n",options->Tperturb);
  for(m=0;m<markerset->nMark;m++){
    PetscScalar Tc = options->thermalBCRight.value[0];
    PetscScalar Th = options->thermalBCLeft.value[0];
    /* not that D essentially falls out of cos term. This is because Barr fixes xmax = lam/(2*D) */	   
    markers[m].T = Th + (Tc-Th)*markers[m].X/grid->LX + options->Tperturb * cos(2.0*M_PI*D/lam* (markers[m].X/grid->LX)*lam/2.0/D)*sin( markers[m].Y*M_PI/D);    
    markers[m].Mat = 0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}

PetscErrorCode initialConditionsMagmaMixing( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* This is the initial condition from Amy Barr's PhD thesis, equation 2.9 */
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar D = grid->LY;
  const PetscScalar lam = 2.0*grid->LX; /* perturbation wavelength */
  const PetscScalar dfrac = 4.0/5.0;    /* depth to interface */
  const PetscScalar dperturb = 0.01;    /* initial y-fractional perturbation amplitude */

  printf("Perturbation wavelength %e\n",lam);
  printf("Tperturb %e\n",options->Tperturb);
  for(m=0;m<markerset->nMark;m++){
    PetscScalar TTop = options->thermalBCTop.value[0];
    PetscScalar TBtm = options->thermalBCBottom.value[0];
    markers[m].T = TTop + (TBtm-TTop)*markers[m].Y/grid->LY;

    /* make bottom composition 0 and top composition 1 */
    PetscScalar ybond = grid->LY * (dfrac + dperturb*cos(2.0*M_PI/lam*markers[m].X/grid->LX));
    if( markers[m].Y > ybond ){
      markers[m].Mat = 1;
    }else{
      markers[m].Mat = 0;
    }
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}


PetscErrorCode initialConditionsMagmaMixingBA( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* Initialize domain with 'dike' in center of fixed width */
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar dw = 5.0; /* dike width */
  const PetscScalar T0 = 1125;/* temperature of material 1 */
  const PetscScalar T1 = 950;/* temperature of material 2 */

  const PetscScalar D = grid->LY;

  for(m=0;m<markerset->nMark;m++){

    /* 0 = INTRUSION, 1 = HOST */
    if( markers[m].X > (grid->LX/2.0 - dw/2.0) && markers[m].X < (grid->LX/2.0 + dw/2.0) ){
      markers[m].Mat = 0;
      markers[m].T = T0;
    }else{
      markers[m].Mat = 1;
      markers[m].T = T1;
    }
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}


PetscErrorCode initialConditionsConvectionGerya( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    PetscScalar TTop = options->thermalBCTop.value[0];
    PetscScalar TBtm = options->thermalBCBottom.value[0];
    
    PetscScalar xwt = markers[m].X/grid->LX;
    PetscScalar ywt = markers[m].Y/grid->LY;
    PetscScalar ymid = 0.5;
    if(ywt>ymid){/* this temperature distribution is from gerya's chapter 16 benchmark for convection */
      ywt=(1-ymid)+ymid*pow(fabs(ywt-ymid)/(1-ymid),10);
    }else{
      ywt=(1-ymid)-(1-ymid)*pow(fabs(ywt-ymid)/(ymid),10);
    }
    markers[m].T = (TTop+(TBtm-TTop)*pow(ywt,2+xwt) + TBtm+(TTop-TBtm)*pow(1-ywt,3-xwt))/2;
    
    markers[m].Mat = 0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}

PetscErrorCode initialConditionsOOPCouette( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    PetscScalar xfrac = markers[m].X/grid->LX;
    markers[m].Mat = 0;
    markers[m].T = options->thermalBCLeft.value[0] + (options->thermalBCRight.value[0]-options->thermalBCLeft.value[0])*xfrac;
    markers[m].rho = markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat];
    markers[m].rhodot = 0.0;
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
  }
  PetscFunctionReturn(ierr);
}

PetscErrorCode initialConditionsBending( MarkerSet *markerset , Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    /* slab bending initial condition */
    if( markers[m].X < grid->LX*4.0/5.0 && markers[m].Y > grid->LX/5.0 && markers[m].Y < grid->LX*4.0/5.0 ){
      /* slab */
      markers[m].Mat = 0;
      markers[m].T = 1.0;
    }else{
      /* not slab */
      markers[m].Mat = 1;
      markers[m].T = 1.0;
    }
    markers[m].rhodot = 0.0;
    //markers[m].rho = materials.materialRho[(PetscInt) markers[m].Mat];
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    markers[m].eta = materials->materialEta[(PetscInt) markers[m].Mat];
    markers[m].mu = materials->materialMu[(PetscInt) markers[m].Mat];
  }
  
  PetscFunctionReturn(ierr);
}

PetscErrorCode initialConditionsDilationTest( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  PetscScalar TTop = options->thermalBCTop.value[0];
  PetscScalar TBtm = options->thermalBCBottom.value[0];

  const PetscScalar w1 = 1000.0; /* width of dilating material */

  for(m=0;m<markerset->nMark;m++){
    if( markers[m].Y <= grid->LY/5.0 ){
      markers[m].Mat = 2;
    }else if( fabs( markers[m].X - grid->LX/2.0 ) < ( w1  - w1*( markers[m].Y - grid->LY/5.0 )/(4.0/5.0*grid->LY)) ){
      markers[m].Mat = 1;
    }else{
      markers[m].Mat = 0;
    }
    markers[m].T = 100.0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );    
  }

  PetscFunctionReturn(ierr);
}

PetscErrorCode initialConditionsDilationTest2( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* This adds a lower layer */
  /* sticky air (material 2) 1/5 of domain */
  /* material 0+1 3/5 of domain */
  /* material 3   1/5 of domain */
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  PetscScalar TTop = options->thermalBCTop.value[0];
  PetscScalar TBtm = options->thermalBCBottom.value[0];

  const PetscScalar w1 = 250.0; /* width of dilating material */

  for(m=0;m<markerset->nMark;m++){
    if( markers[m].Y <= grid->LY/5.0 ){ /* sticky air */
      markers[m].Mat = 2;
      markers[m].T = TTop;
    } else if(markers[m].Y > 4.0*grid->LY/5.0){ /* basement */
      markers[m].Mat = 3;
      /*}else if( fabs( markers[m].X - grid->LX/2.0 ) < ( w1  - w1*( markers[m].Y - grid->LY/5.0 )/(3.0/5.0*grid->LY)) ){*/
    }else if( fabs(markers[m].X - grid->LX/2.0) < w1 && markers[m].Y > (grid->LY/5.0+0.0) && markers[m].Y < 4.0*grid->LY/5.0-0.0){
      /* dilating region */     
      markers[m].Mat = 1;     
    }else{ 
      markers[m].Mat = 0;
    }
    if(markers[m].Mat != 2){  
      markers[m].T = TBtm + (grid->LY - markers[m].Y)/(4.0/5.0*grid->LY)*(TTop-TBtm);
    }
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );    
  }

  PetscFunctionReturn(ierr);
}
PetscErrorCode initialConditionsInclusion( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  /* This is a test from Deubelbeiss and Kaus 2008 */
  /* Circular inclusion, radius 0.1 embedded in a domain with width and height 2 */
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  PetscScalar TTop = options->thermalBCTop.value[0];
  PetscScalar TBtm = options->thermalBCBottom.value[0];
  const PetscScalar r = 0.1;
  const PetscScalar xc = grid->LX/2.0;
  const PetscScalar yc = grid->LY/2.0;

  for(m=0;m<markerset->nMark;m++){
    if( pow(markers[m].X - xc,2.0)+pow(markers[m].Y - yc,2.0) < r*r ){
      markers[m].Mat = 1;
    }else{
      /* material 0 */
      markers[m].Mat = 0;
    }      
    markers[m].T = 100.0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );    
  }
  PetscFunctionReturn(ierr);
}


PetscErrorCode initialConditionsSandbox( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr =0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  PetscScalar TTop = options->thermalBCTop.value[0];
  PetscScalar TBtm = options->thermalBCBottom.value[0];
  
  for(m=0;m<markerset->nMark;m++){
    /* These initial conditions are taken directly from Gerya's Sandbox_shortening_ratio.m for consistency */
    /* 0 = air, 1 = sand1, 2=sand2, 3=microbeads, 4=wall */
    
    //    % Layered structure in the sand
    //    m2=double(int16(MY(ym,xm)/(ysize/24)-0.5))+1;
    markers[m].Mat = (((PetscInt) (markers[m].Y/(grid->LY/24))) % 2) +1; /* should return either 0 or 1, in alternating horizontal bands */
        
    //% Microbeads
    /*     if (markers[m].Y<-0.005 && markers[m].Y>ysize-0.010) */
    if( markers[m].Y < grid->LY-0.005 && markers[m].Y > grid->LY-0.010){
      markers[m].Mat = 3;
    }
    //% Mobile Wall
    if (markers[m].X>grid->LX - 0.01 && markers[m].X<grid->LX-0.005 && markers[m].Y>0.005 && markers[m].Y<grid->LY-0.002){
      markers[m].Mat = 4;
    }
    //% Air
    if (markers[m].X<=grid->LX-0.11 && markers[m].Y<=grid->LY-0.035){
      markers[m].Mat = 0;
    }
     if (markers[m].X>=grid->LX-0.005 && markers[m].Y<=grid->LY-0.002){     
       markers[m].Mat = 0; 
     } 
    if (markers[m].Y<=0.005){
      markers[m].Mat = 0;
    }
    if(markers[m].X>grid->LX-0.11 && markers[m].X<=grid->LX-0.01 && markers[m].Y<=grid->LY-0.035-(markers[m].X-(grid->LX-0.11))*tan(10.0/180.0*M_PI)){
      markers[m].Mat = 0;
    }
    //% Initial temperature profile
    PetscScalar ywt=markers[m].Y/grid->LY;
    markers[m].T=TTop+(TBtm-TTop)*ywt;
    
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    markers[m].eta = materials->materialEta[(PetscInt) markers[m].Mat];
    markers[m].mu = materials->materialMu[(PetscInt) markers[m].Mat];
  }
  PetscFunctionReturn(ierr);
}

PetscErrorCode sandboxMobileWallAndTimesteps( MarkerSet *markerset, GridData *grid, Materials *materials, Options *options, PetscInt iTime, PetscScalar elapsedTime){
  PetscErrorCode ierr = 0;
  PetscFunctionBegin;
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    if (markers[m].X>grid->LX - 0.01 && markers[m].X<grid->LX-0.005 && markers[m].Y>0.005 && markers[m].Y<grid->LY-0.002){
      markers[m].Mat = 4;
      markers[m].Eii = 0.0;
      markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
      markers[m].eta = materials->materialEta[(PetscInt) markers[m].Mat];
      markers[m].mu = materials->materialMu[(PetscInt) markers[m].Mat];
    }
    
    if((markers[m].Mat ==0 || markers[m].Mat==4) && markers[m].Y >= grid->LY-0.002){/* maintain weak (sand) layer */
      markers[m].Mat=2;
    }
    
  }
  if( iTime == 10 || iTime == 20 /*|| iTime == 30 || iTime == 40*/){
    options->dtMax *= sqrt(10.0);
    options->displacementStepLimit *= sqrt(10.0);
    if( options->displacementStepLimit > 0.5 ) options->displacementStepLimit = 0.5;
  }
  if( elapsedTime*2.5/100.0/3600.0 > 0.02 ){
    options->nTime = iTime+1;
    options->saveInterval = 1;
  }

  PetscFunctionReturn(ierr);
}
