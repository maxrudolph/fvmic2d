PetscErrorCode initializeRegularGrid(GridData *, PetscScalar , PetscScalar , PetscInt , PetscInt,Options *);

PetscErrorCode initializeIrregularGrid2Region(GridData *, PetscScalar , PetscScalar , PetscInt , PetscInt );
PetscErrorCode initializeIrregularGridConstantIncrease(GridData *, PetscScalar, PetscScalar, PetscInt , PetscInt ,PetscScalar , PetscScalar );
PetscErrorCode initializeIrregularGridConstantIncreaseFixedRatio(GridData *, PetscScalar, PetscScalar, PetscInt , PetscInt ,PetscScalar );

PetscErrorCode destroyGrid(GridData *);
PetscErrorCode initializeIrregularGridConstantInnerOuter(GridData *, PetscScalar , PetscScalar , PetscInt , PetscInt  , Options *);
PetscErrorCode initializeIrregularGridFault(GridData *, PetscScalar , PetscScalar , PetscInt , PetscInt  , Options *);
PetscErrorCode initializeIrregularGridFaultFS(GridData *, PetscScalar, PetscScalar , PetscInt, PetscInt, Options *, PetscScalar);
