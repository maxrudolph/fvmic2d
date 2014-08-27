#include "fdcode.h"
#include "benchmarkInitialConditions.h"
#include "viscosity.h"
#include "math.h"


PetscErrorCode initialConditionsVanKeken( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar slab_angle = options->slabAngle;
  const PetscScalar D = grid->LY;
  const PetscScalar Tlith = 1573;// temperature at base of overriding plates

  for(m=0;m<markerset->nMark;m++){    
    /* note that D essentially falls out of cos term. This is because Barr fixes xmax = lam/(2*D) */	   
    if( grid->xperiodic){/* for periodic grid, enforce smoothness across lateral boundary */
      fprintf(stderr,"Do not use a periodic boundary condition for the subduction benchmark\n");
      abort();
    }else{
      if( in_slab(markers[m].X, markers[m].Y, slab_angle) ){
	//slab IC	
	markers[m].T = slab_inflow_temperature( markers[m].X, markers[m].Y, slab_angle);
      }else if(markers[m].Y < 50000.0){
	//in overriding plate
	markers[m].T = markers[m].Y/50000.0 * (Tlith-273.0) + 273.0;
      }else{
	//wedge temperature
	markers[m].T = Tlith;
      }     
    }
    markers[m].Mat = 0;
    markers[m].rhodot = 0.0;
    markers[m].rho = materials->materialRho[(PetscInt) markers[m].Mat]*(1.0-materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);
    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );
        
  } 
  PetscFunctionReturn(ierr);
}

PetscInt in_plate(PetscScalar x, PetscScalar y, PetscScalar angle){
  const PetscScalar lith_depth = 50000;
  if( !in_slab(x,y,angle) && y < lith_depth ){
    return 1;
  }else{
    return 0;
  }
}

PetscScalar slab_depth(PetscScalar x, PetscScalar angle){
  // calculates depth of slab interface given x and angle
  return x*tan(angle);
}

PetscInt in_slab(PetscScalar x, PetscScalar y, PetscScalar angle ){
  if( y > slab_depth(x,angle) ){
    return 1;
  }else{
    return 0;
  }
}

PetscScalar slab_inflow_temperature(PetscScalar x, PetscScalar y, PetscScalar angle){
  // calculates temperature in the slab based on half-space cooling model
  const PetscScalar Ts = 273;
  const PetscScalar T0 = 1573;
  const PetscScalar kappa = 0.7272e-6;//Van Keken et al. table 1
  const PetscScalar t50 = 1.5778463e15;//50 Myr in seconds
  return Ts + (T0-Ts) * erf( (y-slab_depth(x,angle))/(2.0*sqrt(kappa*t50)) );

}
