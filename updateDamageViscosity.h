//PetscErrorCode updateDamageViscosity( GridData *,MarkerSet *, Materials *, Options *, PetscScalar );
PetscErrorCode updateDamageRate( GridData *,MarkerSet *, Materials *, Options *);
PetscErrorCode updateMarkerDamageRate( Marker *, Tensor33s *, Materials *, Options *);
void updateMarkerDamage( Marker *, Materials *, Options *, PetscScalar);
//void markerDamageViscosityNoUpdate( Marker *, Materials *, Options *, PetscScalar);
PetscScalar getDamageDensity( Marker *, Materials *, Options * );
