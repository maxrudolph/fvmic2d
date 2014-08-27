void updateMarkerViscosity( Marker *, Options *, Materials * ,PetscScalar);
void updateViscosity( MarkerSet *, Options *, Materials *);
void limitViscosity( MarkerSet *, Options *, Materials *);
void limitMarkerViscosity( Marker *, Options *, Materials *);
PetscErrorCode computeRheology( Marker *, Options *, Materials *, PetscScalar , PetscScalar *, PetscScalar *);
