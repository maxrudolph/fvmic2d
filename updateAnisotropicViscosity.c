
#include "fdcode.h"
#include "updateDamageViscosity.h"
PetscScalar greatestPrincipalStress( PetscScalar, PetscScalar , PetscScalar, PetscScalar , PetscScalar , PetscScalar );


PetscErrorCode updateAnisotropicViscosity(GridData *grid,Markers *markers, Materials *materials, Options *options, PetscScalar dt){
  PetscScalar maxdamage;
  /* because damage-density relationship is rho/rho_0 = (1-D)^damagem, we can compute maxdamage on the fly */
  /* maxdamage = 1-options->maxPorosity^(1/m) */
  
  PetscScalar LX = grid->LX;
  PetscScalar LY = grid->LY;
  
  PetscInt m;
  
  PetscFunctionBegin;
  for(m=0;m<markers->nMark;m++){
    if(materials->hasDamage[(PetscInt)markers->Mat[m]] &&  markers->X[m] >0 && markers->Y[m] > 0 && markers->X[m] < LX && markers->Y[m] <LY){
      /* take stresses from markers and use to evolve damage variable*/
      PetscScalar sxx = markers->sxx[m];
      PetscScalar syy = markers->syy[m];
      PetscScalar szz = markers->szz[m];
      PetscScalar sxy = markers->sxy[m];
      PetscScalar sxz = markers->sxz[m];
      PetscScalar syz = markers->syz[m];
      /* calculate greatest principal stress*/
      PetscScalar gps = greatestPrincipalStress(sxx,syy,szz,sxy,sxz,syz);
      PetscInt mymat = (PetscInt) markers->Mat[m];
      /* compute damage*/
      PetscScalar Is = sxx+syy+szz;
      PetscScalar IIs = sxx*syy + syy*szz + szz*sxx - sxy*sxy - syz*syz - sxz*sxz;
      PetscScalar alpha =  materials->hayhurstAlpha[(PetscInt)markers->Mat[m]];
      PetscScalar beta =  materials->hayhurstBeta[(PetscInt)markers->Mat[m]];
      /* g.p.s. is negative */
      if(gps > 0.0){ gps = 0.0;} else {gps = -gps;} /* use negative of gps - this is extensional=positive now*/
      if( IIs < 0.0) IIs = 0.0;
      PetscScalar stressMeasure = alpha*gps + beta*sqrt(3.0*IIs)+(1.0-alpha-beta)*(3.0*markers->p[m]);
      if( stressMeasure < 0.0) stressMeasure = 0.0;
      
      maxdamage = 1.0 - pow(options->maxPorosity, 1.0/materials->damagem[mymat]);
      
      if(markers->D[m] > maxdamage) markers->D[m] = maxdamage;
      if(markers->D[m] < 0.0) markers->D[m] = 0.0;
      
      PetscScalar ddot = materials->damageB[mymat] * pow(stressMeasure, materials->damager[mymat]) * pow((1.0-markers->D[m]),-materials->damagek[mymat]);
      if(markers->p[m] > 0.0){/* dilation with neg pressure is already in hayhurst criterion*/
/* 	ddot -= materials->damageAlpha3[mymat]*markers->D[m]*markers->p[m]/markers->eta[m]; */
	ddot -= materials->damageAlpha3[mymat]*markers->p[m]/markers->eta[m];
      }
      /* compute new viscosity*/
      if((markers->D[m] >= maxdamage && ddot > 0.0) || (markers->D[m] <= 0.0 && ddot <= 0.0)) ddot=0.0;
      

      markers->Ddot[m] = ddot;
      
    } else {/* marker is out of bounds or does not have damage turned on */
      PetscInt mymat = (PetscInt) markers->Mat[m];
      markers->Ddot[m] = 0.0;

    }
  }  
  PetscFunctionReturn(0); 
}

PetscErrorCode updateDamageViscosity( GridData *grid,Markers *markers, Materials *materials, Options *options, PetscScalar dt){
  PetscScalar maxdamage;
  /* because damage-density relationship is rho/rho_0 = (1-D)^damagem, we can compute maxdamage on the fly */
  /* maxdamage = 1-options->maxPorosity^(1/m) */
  PetscScalar LX = grid->LX;
  PetscScalar LY = grid->LY;
  PetscInt m;
  PetscFunctionBegin;
  for(m=0;m<markers->nMark;m++){
    if(materials->hasDamage[(PetscInt)markers->Mat[m]] &&  markers->X[m] >0 && markers->Y[m] > 0 && markers->X[m] < LX && markers->Y[m] <LY){
      /* take stresses from markers and use to evolve damage variable*/
      PetscInt mymat = (PetscInt) markers->Mat[m];
      maxdamage = 1.0 - pow(options->maxPorosity, 1.0/materials->damagem[mymat]);
      /* compute new viscosity*/
      PetscScalar ddot = markers->Ddot[m];
      markers->D[m] += ddot*dt;      
      /* temperature dependent viscosity*/
      if( materials->hasEtaT[ mymat] ){
	markers -> eta[m] = materials->materialEta[ mymat ]*exp( materials->QonR[mymat] * (1.0/markers->T[m] - 1.0/materials->Tref[mymat]) );/* Nimmo et al. Icarus 2003, Appendix*/	
      }else{
	markers -> eta[m] = materials->materialEta[ mymat ];
      }      
      if(markers->D[m] > maxdamage) markers->D[m] = maxdamage;
      if(markers->D[m] < 0.0) markers->D[m] = 0.0;
      markers->eta[m] *= (1.0-markers->D[m]); /* this viscosity should be temperature-corrected already and should not account for damage*/      
    } else {/* marker is out of bounds or does not have damage turned on */
      PetscInt mymat = (PetscInt) markers->Mat[m];
      markers->Ddot[m] = 0.0;
      //      markers->eta[m] = materials->materialEta[(PetscInt) markers->Mat[m]];
      /* update temperature-dependent viscosity*/
      if( materials->hasEtaT[ mymat] ){
	markers -> eta[m] = materials->materialEta[ mymat ]*exp( materials->QonR[mymat] * (1.0/markers->T[m] - 1.0/materials->Tref[mymat]) );	
      }else{
	markers -> eta[m] = materials->materialEta[ mymat ];
      }            
      /* no need to update damage or viscosity*/
    }
    /* enforce viscosity limits*/
    if(markers->eta[m] > options->etamax) markers->eta[m] = options->etamax;
    if(markers->eta[m] < options->etamin) markers->eta[m] = options->etamin;
  }  
  PetscFunctionReturn(0);
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
