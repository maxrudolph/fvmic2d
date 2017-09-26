#include "fdcode.h"
#include "benchmarkInitialConditions.h"
#include "viscosity.h"
#include "math.h"
#include "markers.h"


PetscErrorCode initialConditionsVanKeken( MarkerSet *markerset, Options *options, Materials *materials, GridData *grid){
  PetscErrorCode ierr=0;
  PetscInt m;
  Marker *markers = markerset->markers;
  const PetscScalar slab_angle = options->slabAngle;

  for(m=0;m<markerset->nMark;m++){    
    if( grid->xperiodic){/* for periodic grid, enforce smoothness across lateral boundary */
      fprintf(stderr,"Do not use a periodic boundary condition for the subduction benchmark\n");
      abort();
    }else{
      if( in_slab(markers[m].X, markers[m].Y, options) ){
	//slab IC	
	markers[m].T = slab_inflow_temperature( markers[m].X, markers[m].Y, slab_angle, options->slabAge);
      }else if(in_plate( markers[m].X, markers[m].Y, options )){
	//in overriding plate
	//	markers[m].T = markers[m].Y/50000.0 * (Tlith-273.0) + 273.0;
	markers[m].T = plate_geotherm(markers[m].X, markers[m].Y , options);
      }else{
	//wedge temperature
	markers[m].T = mantle_temperature();
      }     
    }
    markers[m].Mat = 0;
    markers[m].rhodot = 0.0;
    markers[m].rho = -materials->materialRho[(PetscInt) markers[m].Mat]*(materials->materialAlpha[(PetscInt) markers[m].Mat]*markers[m].T);    
    updateMarkerViscosity( &markers[m], options, materials, 0.0 );    
  } 
  /* Set up boundary conditions to match prescribed slab velocity in advection routines */
  options->mechBCLeft.type[0] = 0;
  options->mechBCLeft.type[1] = 0;
  PetscScalar svx = options->slabVelocity * cos( slab_angle );
  PetscScalar svy = options->slabVelocity * sin( slab_angle );
  options->mechBCLeft.value[0] = svx;
  options->mechBCLeft.value[1] = svy;
  options->mechBCTop.type[0] = 0;
  options->mechBCTop.value[0] = 0.0;
  options->mechBCTop.type[1] = 0;
  options->mechBCTop.value[1] = 0.0;


  PetscFunctionReturn(ierr);
}

Marker new_slab_marker(PetscScalar x, PetscScalar y, Options *options, Materials *materials){
  Marker m;
  resetMarker( &m );
  m.X = x;
  m.Y = y;
  m.Mat = 0;
  m.T = slab_inflow_temperature( x, y, options->slabAngle ,options->slabAge);
  updateMarkerViscosity( &m, options, materials, 0.0);
  return m;
}

Marker new_inflow_marker(PetscScalar x, PetscScalar y, Options *options, Materials *materials){
  Marker m;
  resetMarker( &m );
  m.X = x;
  m.Y = y;
  m.Mat = 0;
  m.T = mantle_temperature();
  updateMarkerViscosity( &m, options, materials, 0.0);
  return m;
}


PetscScalar mantle_temperature(){
  const PetscScalar t = 1573.0;
  return t;
}

PetscScalar plate_depth(PetscScalar x, Options *options){
  const double plate_thickness = options->overridingPlateThickness;

  const double root_thickness         = options->rootThickness;// 10 km
  const double root_center_location   = options->rootCenter;// 120 km
  const double root_width             = options->rootWidth;//40 km
  /* introduce root with functional form sin^2( (x-center)/width) ) */
  if( fabs( x-root_center_location ) < root_width ){   
    return plate_thickness + root_thickness * pow( cos( M_PI*(x-root_center_location)/root_width/2.0 ) , 2.0 );
  }else{
    return plate_thickness;
  }
}

PetscScalar plate_geotherm(PetscScalar x, PetscScalar y, Options *options){
  PetscScalar T = y/plate_depth(x,options) * (mantle_temperature() - 273.0) + 273.0;
  return T;
}


PetscInt in_plate(PetscScalar x, PetscScalar y, Options *options){
  const PetscScalar lith_depth = plate_depth(x,options);
  if( !in_slab(x,y,options) && y <= lith_depth ){
    return 1;
  }else{
    return 0;
  }
}

PetscScalar slab_depth(PetscScalar x, PetscScalar angle){
  // calculates depth of slab interface given x and angle
  return x*tan(angle);
}

PetscScalar slab_x(PetscScalar y, PetscScalar angle){
  // calculates horizontal position of slab interface given depth (y)
  return y/tan(angle);
}

PetscInt in_slab(PetscScalar x, PetscScalar y, Options *options ){
  const PetscScalar angle = options ->slabAngle;
  if( y >= slab_depth(x,angle) ){
    return 1;
  }else{
    return 0;
  }
}

PetscScalar slab_inflow_temperature(PetscScalar x, PetscScalar y, PetscScalar angle, PetscScalar age){
  // calculates temperature in the slab based on half-space cooling model
  const PetscScalar Ts = 273;
  const PetscScalar T0 = 1623;
  const PetscScalar kappa = 0.7272e-6;//Van Keken et al. table 1
  const PetscScalar s_in_year = 3.15576e7;
  const PetscScalar slab_age_s = age * 1.0e6 * s_in_year;
  return Ts + (T0-Ts) * erf( (y-slab_depth(x,angle))/(2.0*sqrt(kappa*slab_age_s)) );

}

PetscInt in_either( PetscScalar x, PetscScalar y, Options *options){
  if( in_slab( x, y, options) || in_plate(x,y,options) ){
    return 1;
  }else{
    return 0;
  }
}
