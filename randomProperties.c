#include "fdcode.h"
#include "randomProperties.h"

/* This function assigns random values to several parameters for material 0*/
/* the random number generator needs to be created externally so that it can be re-used in successive runs*/

PetscErrorCode getRandomProperties( Materials *materials, Options *options, PetscRandom *r){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  /* which parameters will be varied?*/
  /* vbz */
  /* damageAlpha3 - annealing pre-factor */
  /* damagem - exponent in damage-porosity relationship*/
  /* damageB - accumulation pre-factor */
  
  /* initialize parameters with random values */
  /* vbz*/
  PetscScalar vbz_min = -10;//1e-11;
  PetscScalar vbz_max = -6;//1e-8;
  ierr=PetscRandomSetInterval(*r, vbz_min, vbz_max);CHKERRQ(ierr);
  PetscRandomGetValue(*r, &options->mechBCLeft.value[2]);
  options->mechBCRight.value[2] = pow(10,options->mechBCLeft.value[2]);
  options->mechBCLeft.value[2] = -options->mechBCRight.value[2];
  
  /* vbx */
  PetscScalar vbx_min = -11; 
  PetscScalar vbx_max = -6;
  ierr=PetscRandomSetInterval(*r, vbx_min, vbx_max);CHKERRQ(ierr); 
  PetscRandomGetValue(*r, &options->mechBCLeft.value[0]); 
  options->mechBCRight.value[0] = pow(10,options->mechBCLeft.value[0]); 
  options->mechBCLeft.value[0] = -options->mechBCRight.value[0];
  
  
  /* alpha3 */
  PetscScalar alpha3_min = -5;
  PetscScalar alpha3_max = 3;
  ierr=PetscRandomSetInterval(*r, alpha3_min, alpha3_max);CHKERRQ(ierr);
  PetscRandomGetValue(*r, &materials->damageAlpha3[0]);
  materials->damageAlpha3[0] = pow(10,materials->damageAlpha3[0]);
  //  materials->damageAlpha3[0] = 0.0; 

  /* m */
  PetscScalar m_min = 0.0;
  PetscScalar m_max = 3.0;
  ierr=PetscRandomSetInterval(*r, m_min, m_max);CHKERRQ(ierr);
  PetscRandomGetValue(*r, &materials->damagem[0] );
  materials->damagem[0] = pow(10,materials->damagem[0]);

  /* B */
  PetscScalar B_min = -30;//1e-30;
  PetscScalar B_max = -8;//1e-20;
  ierr=PetscRandomSetInterval(*r, B_min, B_max); CHKERRQ(ierr);
  PetscRandomGetValue(*r, &materials->damageB[0]);
  materials->damageB[0] = pow(10,materials->damageB[0]);

  /* distribute options to all nodes */
  ierr= MPI_Bcast( options, sizeof(Options), MPI_BYTE, 0,PETSC_COMM_WORLD);

  PetscFunctionReturn(ierr);
}
