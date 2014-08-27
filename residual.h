
PetscScalar getGlobalStrainRateResidual( PetscInt );
void setGlobalStrainRateResidual( PetscInt , PetscScalar );
PetscErrorCode updateGlobalStrainRateResidual( GridData *, MarkerSet *, PetscInt  );
PetscInt checkConvergence( PetscScalar , PetscScalar , PetscInt );
