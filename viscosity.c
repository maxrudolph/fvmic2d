#include "fdcode.h"
#include "viscosity.h"
#include <math.h>
#ifdef COMPBA
#include "phaseBA.h"
#endif
/* This file contains subroutines that implement various rheologies on the markers */

PetscScalar goldsbyKohlstedt(PetscScalar , PetscScalar , PetscScalar , PetscScalar );
PetscScalar calculateCompositionViscosityBA( Materials *, Marker *);

PetscErrorCode computeRheology( Marker *marker, Options *options, Materials *materials, PetscScalar s, PetscScalar *eta, PetscScalar *mu){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  PetscScalar thiseta, thismu;
  const PetscInt mymat = (PetscInt) marker -> Mat;
  thismu = materials -> materialMu[ mymat ];

  if( materials->hasEtaT[ mymat] ==  1){
    /* temperature dependent */
    /* this is a typical diffusion creep rheology. Can be made equivalent, for instance, to Han and Showman 2010 */
    thiseta = materials->materialEta[ mymat ]*exp( materials->QonR[mymat] * (1.0/marker->T - 1.0/materials->Tref[mymat]) );	
  }else if(materials->hasEtaT[mymat]==2){
    /* For T+S shear heating couette flow benchmark */
    PetscScalar T0 = materials->Tref[mymat];
    thiseta = materials->materialEta[mymat]*exp(materials->QonR[mymat]/T0*(1.0-(marker->T-T0)/T0));
  }else if(materials->hasEtaT[mymat] ==3 ){
    /* viscosity linear in temperature - useful for some analytic solution benchmarks */
    thiseta = materials->materialEta[mymat]*marker->T;

  }else if(materials->hasEtaT[mymat] == 5){
    /* Form for blankenback benchmarks */
    thiseta = materials->materialEta[mymat]*exp(-materials->QonR[mymat]*marker->T);

  }else if(materials->hasEtaT[mymat] == 0){
    /* no temperature dependence */
    thiseta = materials->materialEta[ mymat ];

  }else if(materials->hasEtaT[mymat] == 6){
    /* viscous layer */ 
    const PetscScalar depth = marker->Y;
    const PetscScalar distance = marker->X;
    PetscScalar slabAngle = options->slabAngle;
    PetscScalar decoupleThickness = options->decoupleThickness;
    PetscScalar minDecoupleDepth = options->minDecoupleDepth;
    PetscScalar maxDecoupleDepth = options->maxDecoupleDepth;
    PetscScalar decoupleEtaScaler = options->decoupleEtaScaler;
    thiseta = materials->materialEta[ mymat ];
    if( distance > depth/tan(slabAngle) ){
      if( distance < (depth/tan(slabAngle) + decoupleThickness) ){
	if( depth > minDecoupleDepth && depth < maxDecoupleDepth ){
	  thiseta *= decoupleEtaScaler;
	}
      }
    }
    /* }else if(materials->hasEtaT[mymat] == 7){
       /* temperature dependent */
    /* thiseta = materials->materialEta[mymat]*exp(materials->QonR[mymat] * (materials->Tref[mymat] - marker->T)); */

  }else{
    printf("Error! Viscosity law undefined\n");
    abort();
  }  
  eta[0] = thiseta;
  mu[0] = thismu;
  PetscFunctionReturn(ierr);
}


void updateMarkerViscosity( Marker *marker, Options *options, Materials *materials, PetscScalar s ){
  
  computeRheology( marker, options, materials, s, &(marker->eta), &(marker->mu) );
   

}

void updateViscosity( MarkerSet *markerset, Options *options, Materials *materials){
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    if( markerset->markers[m].cellX != - 1){
      /* use stress second invariant from marker */
      PetscScalar sii = sqrt(0.5*(markers[m].s.T11*markers[m].s.T11 + markers[m].s.T22*markers[m].s.T22 ) + markers[m].s.T12*markers[m].s.T12 );

      updateMarkerViscosity( &markerset->markers[m], options, materials, sii);
    }
  }
}

void limitViscosity( MarkerSet *markerset, Options *options, Materials *materials){
  PetscInt m;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    limitMarkerViscosity( &markers[m], options, materials);
  }
}

void limitMarkerViscosity( Marker *marker, Options *options, Materials *materials){
  if( marker->eta > options->etamax ){
    marker->eta = options->etamax;
  }else if( marker->eta < options->etamin ){
    marker->eta = options->etamin;
  }
}

PetscScalar goldsbyKohlstedt(PetscScalar s, PetscScalar T, PetscScalar d, PetscScalar P){
  /* returns strain rate predicted by goldsby and kohlstedt 2001 in units of s^-1*/
  PetscScalar etot=0.0; /* total strain-rate predicted from constitutive equation */
  const PetscScalar R = 8.3144621; /* ideal gas constant, J/mol/K */
  PetscScalar etaeff = 0.0;/* effective visocity */
  /* diffusion creep */
  PetscScalar ediff;  
  if( P < 1.0 ) P = 1.0;/* prevent spurious pressures from affecting viscosity */
  {
    /* quantities from table 6: Diffusion Creep Parameters */
    const PetscScalar b = 4.52e-10; /* burgers vector */
    const PetscScalar Vm = 1.97e-5; /* m^3 */
    const PetscScalar Qv = 59.4*1000.0;/* Activation energy for volume diffusion, J/mol */
    const PetscScalar delta = 2.0*b; /* grain boundary width, m */
    const PetscScalar Qb = 49*1000.0;
    const PetscScalar D0v = 9.10e-4; /* m^2/s */
    const PetscScalar D0b = 5.8e-4; /* m^2/s - this is lower of two upper bounde presented in text */
    /* diffusivities */
    const PetscScalar Dv = D0v*exp(-Qv/(R*T));
    const PetscScalar Db = D0b*exp(-Qb/(R*T));
    /* note that all units above are SI, so strain rate will be in 1/s for stress in Pa */
    ediff = 42.0*Vm/(R*T*d*d)*(Dv + M_PI*delta/d*Db );
  }
  PetscScalar egbs;
  if( s < 1e0 ){/* if stress is really small (i.e. for first timestep), we cannot compute dislocation creep or gbs strain rates */
    etaeff = 0.5/ediff;
    return(etaeff);
  }
  
  ediff *= s;
  {
    PetscScalar A,n,Q;
    if( T < 255 ){
      A=3.9e-3; /* MPa^-1.8 m^1.4 s^-1 */
      n=1.8;
      Q=49*1000.0; /* J/mol */
    }else{
      A = 3.0e26; /* MPa^-1.8 m^1.4 S^-1 */
      n = 1.8;
      Q = 192*1000.0;
    }
    const PetscScalar V = -13e-6; /* m^3/mol, activation volume from text */
    const PetscScalar p = 1.4;/* grain size dependence */
    egbs = A*pow(s/1e6,n)/pow(d,p)*exp( (-Q + P*V)/(R*T) );
  }
  PetscScalar edisl;
  {
    PetscScalar A,n,Q;
    if( T < 258 ){
      A=4.0e5; /* MPa^-4.0 s^-1 */
      n=4.0;
      Q=60*1000.0; /* J/mol */
    }else{
      A = 6.0e28; /* MPa^-4.0 S^-1 */
      n = 4.0;
      Q = 180*1000.0; /* Note that there is a mistake in G&K table 5. Should be 180 kJ/mol per text */
    }
    const PetscScalar V = -13e-6; /* m^3/mol, activation volume from text */
    edisl = A*pow(s/1e6,n)*exp( (-Q + P*V)/(R*T) );
  } 
  PetscScalar ebs;
  {
    const PetscScalar A = 5.5e7; /* MPa^-2.4 s^-1 */
    const PetscScalar n = 2.4;
    const PetscScalar Q = 60.0*1000.0; /* J/mol */
    const PetscScalar V = -13e-6; /* m^3/mol, activation volume from text */
    ebs = A*pow(s/1e6,n)*exp( (-Q + P*V)/(R*T) );
  }
  etot = ediff + 1.0/(1.0/ebs + 1.0/egbs) + edisl;
  
  etaeff = s/2.0/etot;

  return(etaeff);
}
