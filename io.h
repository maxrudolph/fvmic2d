PetscErrorCode saveRunInfo( Options *, Materials *, PetscInt );
PetscErrorCode saveMarkersBinary( MarkerSet *, PetscInt,PetscInt,PetscScalar);
PetscErrorCode saveGriddedMarkersBinary( MarkerSet *, GridData *, PetscInt, PetscInt, PetscInt,PetscInt,PetscScalar);
PetscErrorCode saveMarkerFieldS(MarkerSet *, PetscScalar *, PetscInt *, FILE *);
PetscErrorCode saveMarkerFieldI(MarkerSet *, PetscInt *, PetscInt *, FILE *);

PetscErrorCode saveNodalFields( NodalFields *, GridData *,PetscInt,PetscInt, PetscScalar , PetscScalar, PetscInt);
PetscErrorCode saveTextureBinary(MarkerSet *, Materials *, PetscInt, PetscInt);
PetscErrorCode saveGrid( GridData *);
