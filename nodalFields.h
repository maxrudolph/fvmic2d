PetscErrorCode initializeNodalFields( NodalFields *, GridData *, Options *);
PetscErrorCode finalizeNodalFields( NodalFields *);
PetscErrorCode saveNodalFieldsASCIIMatlab( NodalFields *, GridData *, PetscInt,PetscInt,PetscScalar,PetscScalar);
PetscErrorCode resetNodalFields( NodalFields *, GridData *,Options *);
PetscErrorCode destroyNodalFields( NodalFields *, GridData *);
