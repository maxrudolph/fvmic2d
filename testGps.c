



#include "fdcode.h"
//#include <stdlib>
PetscScalar greatestPrincipalStress( PetscScalar, PetscScalar , PetscScalar, PetscScalar , PetscScalar , PetscScalar );
int main(){

  PetscScalar sxx = 2.0;
  PetscScalar sxy = 3.5;
  PetscScalar sxz = 0.1;
  PetscScalar syz = 0.1;
  PetscScalar syy = -9.0;
  PetscScalar szz = 7.0;


  
  greatestPrincipalStress( sxx,syy,szz,sxy,sxz,syz);


  return(0);
}

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
  PetscInt ncp=0;/* number of critical points (where derivative of cubic is zero*/
  
  PetscScalar tol=1e-15;
  //  PetscScalar delta = 18.0*a*b*c*d - 4.0*b*b*b*d + b*b*c*c - 4.0*a*c*c*c - 27.0*a*a*d*d;
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
  printf("principal stresses: %e %e %e\n",s1,s2,s3);
  
  PetscScalar val=s1;
  if(nt>1){
  if(s2 > s1) val = s2;
  if(s3 > val) val = s3;
  }
  return val;

}
