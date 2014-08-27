#include "fdcode.h"
#include <math.h>
#include "updateBCs.h"

/* this is a subroutine to allow for time-varying boundary conditions*/

const PetscScalar pi=3.14159265358979323846;

PetscErrorCode updateBoundaryConditions( Options *options, BoundaryValues *bv , PetscScalar currentTime){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  /* x component */

  /* first just assign options->mechBC* to bv->mechBC* */
  bv->mechBCLeft   = options->mechBCLeft;
  bv->mechBCRight  = options->mechBCRight;
  bv->mechBCTop    = options->mechBCTop;
  bv->mechBCBottom = options->mechBCBottom;

  if( options->oscillatoryX){
    bv->mechBCLeft.value[0]  *= sin( 2*pi* currentTime/options->oscillationPeriod);
    bv->mechBCRight.value[0] *= sin( 2*pi* currentTime/options->oscillationPeriod);
  }
  
  /* z component */
  if( options->oscillatoryZ){
    bv->mechBCLeft.value[2]  *= sin( 2*pi* currentTime/options->oscillationPeriod);
    bv->mechBCRight.value[2] *= sin( 2*pi* currentTime/options->oscillationPeriod);
    
  }

  //printf("vbz = %e, %e\n",bv->mechBCLeft.value[2], bv->mechBCRight.value[2]);




  PetscFunctionReturn(ierr);
}
