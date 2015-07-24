PetscErrorCode restartFromMarkers( Problem * , PetscInt , PetscInt , PetscScalar *);
PetscErrorCode restartFromMarkersNoTexture( MarkerSet *,GridData *,Materials *, Options *, PetscInt , PetscInt , PetscScalar *);
PetscErrorCode  loadTextureBinary(MarkerSet *, Materials *,PetscInt *, PetscInt, PetscInt);

PetscErrorCode loadMarkerFieldI( PetscInt *, PetscInt *, FILE *);
PetscErrorCode loadMarkerFieldS( PetscScalar *, PetscInt *, FILE *);
