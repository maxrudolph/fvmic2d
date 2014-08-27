/* Profiling Information */

#define LOG_N_STAGES 15

typedef enum {
  LOG_MAIN = 0,
  LOG_MECH_ASM = 1,
  LOG_YIELD = 2,
  LOG_SOLVE = 3,
  LOG_MARK_EXCH = 4,
  LOG_PROJECT_M_TO_N = 5,
  LOG_NODAL_STR = 6,
  LOG_SUBGRID_STRESS = 7,
  LOG_THERM_ASM = 8,
  LOG_THERM_SOLVE = 9,
  LOG_SUBGRID_TEMP = 10,
  LOG_IO = 11,
  LOG_MARK_STR_PRESS = 12,
  LOG_MARK_MOVE = 13
} LogStage;


PetscLogStage * initializeLogging();
PetscErrorCode setLogStage( LogStage );
PetscErrorCode finalizeLogging( );
