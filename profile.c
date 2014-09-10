#include "fdcode.h"
#include "profile.h"

static PetscLogStage *stages;

void initializeLogging( ){
  /* log stages */

  PetscMalloc( LOG_N_STAGES*sizeof( PetscLogStage), &stages );

  PetscLogStageRegister("Main",&stages[0]);
  PetscLogStageRegister("Mech Assem",&stages[1]);
  PetscLogStageRegister("Check Yield",&stages[2]);
  PetscLogStageRegister("Linear Solve",&stages[3]);
  PetscLogStageRegister("Exch Mark",&stages[4]);
  PetscLogStageRegister("Proj M to N",&stages[5]);
  PetscLogStageRegister("Node StrStra",&stages[6]);
  PetscLogStageRegister("Subgrid Stres",&stages[7]);
  PetscLogStageRegister("Therm Assem",&stages[8]);
  PetscLogStageRegister("Therm Solve",&stages[9]);
  PetscLogStageRegister("Sub. Temp.",&stages[10]); 
  PetscLogStageRegister("I/O",&stages[11]); 
  PetscLogStageRegister("Mark Str Press",&stages[12]);
  PetscLogStageRegister("Advection",&stages[13]);  
  PetscLogStageRegister("Add Markers",&stages[14]);
  PetscLogStageRegister("Find Cells",&stages[15]);
  PetscLogStagePush(stages[0]);

  //  return( stages );


}

PetscErrorCode setLogStage( LogStage ls ){
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = PetscLogStagePush( stages[ls] );
  PetscFunctionReturn(ierr);
}

PetscErrorCode finalizeLogging( ){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr=PetscFree( stages );
  PetscFunctionReturn(ierr);
}
