
#include "fdcode.h"
#include "updateDamageViscosity.h"
#include "viscosity.h"

PetscScalar greatestPrincipalStress( PetscScalar, PetscScalar , PetscScalar, PetscScalar , PetscScalar , PetscScalar );

PetscScalar getDamageDensity( Marker *marker, Materials *materials, Options *options ){
  /* calculate density for a marker based on damage */
  const PetscInt mymat = marker->Mat;
  PetscScalar rho0 =  materials->materialRho[mymat]*(1-materials->materialAlpha[mymat]*marker->T);
  PetscScalar rho=rho0;
  if( materials->hasDamage[mymat] && materials->hasDilation[mymat] ){
    /* calculate damage component */
    if( options->dilationModel == 1){
      rho = pow(1-marker->D,materials->damagem[mymat]) * rho0;
    }else if( options->dilationModel == 2){
      PetscScalar D = marker->D;
      PetscScalar Dmax = (1.0-options->Dcrit)*(1.0-pow(1.0-options->maxPorosity,1.0/materials->damagem[mymat])) + options->Dcrit;
      if( D < options->Dcrit ){
	rho = rho0;
      }else if(D >= Dmax){
	rho = (1.0-options->maxPorosity)*rho0;
      }else{
	rho = pow( 1.0-(D-options->Dcrit)/(1.0-options->Dcrit), materials->damagem[mymat])*rho0;
      }      
    }else{
      printf(" ERROR: Damage function not implemented \n");
    }
  }
  return rho;
}

PetscErrorCode updateMarkerDamageRate( Marker *marker, Tensor33s *s, Materials *materials, Options *options){
  /* This routine updates the marker damage rate and viscosity, without changing the marker damage. Timestep is used because we do calculate new damage based on old damage and damage rate when computing viscosity */
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  
  if(materials->hasDamage[marker->Mat] && marker->cellX != -1){    
    /* take stresses from markers and use to evolve damage variable*/
    PetscScalar sxx = s->T11; //xx;
    PetscScalar syy = s->T22; //yy;
    PetscScalar szz = s->T33; //zz;
    PetscScalar sxy = s->T12; //xy;
    PetscScalar sxz = s->T13; //xz;
    PetscScalar syz = s->T23; //yz;
    /* calculate greatest principal stress*/
    PetscScalar gps = greatestPrincipalStress(sxx,syy,szz,sxy,sxz,syz);
    PetscInt mymat = (PetscInt) marker->Mat;
    /* compute damage*/
    /*     PetscScalar Is = sxx+syy+szz; */
    PetscScalar IIs = sxx*syy + syy*szz + szz*sxx - sxy*sxy - syz*syz - sxz*sxz;
    PetscScalar alpha =  materials->hayhurstAlpha[(PetscInt)marker->Mat];
    PetscScalar beta =  materials->hayhurstBeta[(PetscInt)marker->Mat];
    /* g.p.s. is negative */
    if(gps > 0.0){ gps = 0.0;} else {gps = -gps;} /* use negative of gps - this is extensional=positive now*/
    if( IIs < 0.0) IIs = 0.0;
    PetscScalar stressMeasure = 1.0/(1.0-marker->D)*(alpha*gps + beta*sqrt(3.0*IIs)+(1.0-alpha-beta)*(3.0*marker->p));
    if( stressMeasure < 0.0) stressMeasure = 0.0;
        
    if(marker->D < 0.0) marker->D = 0.0;
    
    PetscScalar ddot = materials->damageB[mymat] * pow(stressMeasure, materials->damager[mymat]) * pow((1.0-marker->D),-materials->damagek[mymat]);

    if(marker->p > 0.0){/* dilation with neg pressure is already in hayhurst criterion*/
      ddot -= materials->damageAlpha3[mymat]*(1.0/(1.0-marker->D))*marker->p/marker->eta;
    }
    
    marker->Ddot = ddot;
    /* calculate rate of density change */
    if( materials->hasDamage[mymat] && materials->hasDilation[mymat] ){
      if( options->dilationModel == 1 ){
	/* this expression is for d\rho/dt = m*(1.0-D)^(m-1)*-(dD/dt)*rho_0, i.e. rho= rho_0 * (1.0-D)^m */
	PetscScalar rho0 = materials->materialRho[mymat]*(1-materials->materialAlpha[mymat]*marker->T);/* portion of density due to temperature*/
	PetscScalar rhodot;
	PetscScalar Dmax = 1.0-pow(1.0-options->maxPorosity,1.0/materials->damagem[mymat]);
	
	if( marker->D > Dmax ){
	  rhodot = 0.0;
	}else{
	  rhodot =  materials->damagem[mymat]*pow(1.0-marker->D,materials->damagem[mymat]-1.0)*-marker->Ddot * rho0;
	}
	marker->rhodot = rhodot;
	
      } else if(options->dilationModel == 2){
	/* second formulation */	      
	PetscScalar Dmax = (1.0-options->Dcrit)*(1.0-pow(1.0-options->maxPorosity,1.0/materials->damagem[mymat])) + options->Dcrit;
	
	if( marker->D < options->Dcrit || marker->D >= Dmax ){
	  /* no density changes until D gets very close to 1 */
	  //marker->rhodot = -materials->materialRho[mymat]*materials->materialAlpha[mymat]*marker->Tdot;
	  marker->rhodot = 0.0;
	} else {
	  /* this dilation model has rho = rho0 for D < options->Dcrit */
	  /* rho = (1-D')^m *rho0 where D'=(D-Dcrit)/(1-Dcrit) */
	  /* rhodot = m*(1-D')^(m-1)*ddot/(1-Dcrit)*rho0 */
	  PetscScalar rho0 = materials->materialRho[mymat]*(1-materials->materialAlpha[mymat]*marker->T);
	  marker -> rhodot = -materials->damagem[mymat] * pow( 1.0-(marker->D-options->Dcrit)/(1.0-options->Dcrit), materials->damagem[mymat]-1.0 )*ddot/(1.0-options->Dcrit)*rho0;	 
	}
      }else{/* end if then on dilation Model */
	marker->rhodot = 0.0;
      }
    }
    
  } else {/* marker is out of bounds or does not have damage turned on */    
    marker->Ddot = 0.0;
    marker->rhodot = 0.0;
  }      
  PetscFunctionReturn(ierr);
}


PetscErrorCode updateDamageRate(  GridData *grid,MarkerSet *markerset, Materials *materials, Options *options){
  PetscFunctionBegin;    
  PetscInt m;
  Marker *markers = markerset->markers;  
  for( m=0;m<markerset->nMark;m++){
    updateMarkerDamageRate( &(markers[m]), &(markers[m].s), materials, options);
  }
  PetscFunctionReturn(0); 
}

void updateMarkerDamage( Marker *marker, Materials *materials, Options *options, PetscScalar dt){
  const PetscScalar maxdamage=1.0-1e-6;
  PetscInt mymat = marker->Mat;
 
  if( materials->hasDamage[ mymat ] ){
    /* calculate new damage using damage rate */
    marker->D += marker->Ddot*dt;
    if( marker->D > maxdamage ) marker->D = maxdamage;
  }
}

/* returns the most extensional (negative) p.s. */
PetscScalar greatestPrincipalStress( PetscScalar sxx, PetscScalar syy, PetscScalar szz, PetscScalar sxy, PetscScalar sxz, PetscScalar syz){  
  /* calculate the principal stresses and return the greatest of the three */
  const PetscScalar pi = 3.141592653589793;
  /* evaluate Is, IIs, IIIs (invariants) */
  PetscScalar Is = sxx+syy+szz;
  PetscScalar IIs = sxx*syy + syy*szz + szz*sxx - sxy*sxy - syz*syz - sxz*sxz;
  PetscScalar IIIs = sxx*(syy*szz-syz*syz) - sxy*(sxy*szz - sxz*syz) + sxz*(sxy*syz - sxz*syy);
  
  /* idea: solve for roots using S^3 - Is S^2 + IIs S - IIIs = 0 , where S is a principal stress*/
  PetscScalar a=1.0;
  PetscScalar b=-Is;
  PetscScalar c= IIs;
  PetscScalar d=-IIIs;
  
  PetscScalar tol=1e-15;

  PetscScalar p;
  PetscScalar q;
  if( fabs(b) < tol){
    p = c;
    q = d;
  }else{
    p  = (3.0*a*c-b*b)/(3.0*a*a);
    q=  (2.0*b*b*b - 9.0*a*b*c + 27.0*a*a*d)/(27.0*a*a*a);
  }
  /* now solve monic trinomial t^3 + p*t + q ==0 */
  PetscScalar t1,t2,t3, nt;

  if(fabs(p) < tol){
    t1 = pow(  -q, 1.0/3.0);
    nt = 1;
  } else {
    /* obtain solutions using trig formulae */
    PetscScalar acosterm = acos( 3.0*q/2.0/p*sqrt(-3.0/p) )/3.0;
    t1 = 2.0*sqrt(-p/3.0) * cos(acosterm);
    t2 = 2.0*sqrt(-p/3.0) * cos(acosterm-2*pi/3);
    t3 = 2.0*sqrt(-p/3.0) * cos(acosterm-4*pi/3);
    nt = 3;
  }
  /* convert from t back to s: s=t-b/(3a) */
  PetscScalar s1,s2,s3;
  s1 = t1-b/3.0/a;
  if(nt>1){
    s2 = t2-b/3.0/a;
    s3 = t3-b/3.0/a;
  } else {
    s2 = 0.0;
    s3 = 0.0;
  }
/*   printf("principal stresses: %e %e %e\n",s1,s2,s3); */
  
  PetscScalar val=s1;
  if(nt>1){
  if(s2 < s1) val = s2;
  if(s3 < val) val = s3;
  }
  return val;

}
