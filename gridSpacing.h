void gridSpacingUniform( PetscScalar *, PetscScalar , PetscScalar, PetscInt );
void makePeriodic( PetscScalar *, PetscScalar *, PetscScalar , PetscInt );
void getCellCenters( PetscScalar *, PetscScalar *, PetscScalar , PetscInt );
void gridSpacingConstantInnerOuter( PetscScalar *, PetscScalar , PetscInt , PetscScalar);
void smoothGrid( PetscScalar *, PetscScalar , PetscInt );
void gridSpacingRefinedCenter( PetscScalar *, PetscScalar, PetscInt, PetscScalar );
void circShift( PetscScalar *, PetscScalar, PetscInt , PetscInt );
void gridSpacingRamp( PetscScalar *x, PetscScalar xmin, PetscScalar xmax, PetscInt NX, PetscScalar hcontrast);
