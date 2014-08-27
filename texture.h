/* this is a data structure for storing information about texture*/
#ifndef NT
#define NT 8
#endif

PetscErrorCode formViscoplasticMMarker(Marker *,  Options *,Tensor33s *, PetscInt );
PetscErrorCode formViscoplasticM(MarkerSet *, Materials *, Options *, PetscInt);
PetscErrorCode invertMtoN(TextureInfo *);
PetscErrorCode formViscoplasticMNewton(MarkerSet *, Materials *,Options *);

PetscErrorCode formViscoplasticMNewtonMarker(Marker *, Materials *, Options *);

PetscErrorCode updateFabric(MarkerSet *, Materials *, PetscScalar );
PetscScalar getFabricDT(MarkerSet *);
PetscErrorCode initializeIsotropicFabric(MarkerSet *);
PetscErrorCode limitAnisotropicViscosity( MarkerSet *, Materials *, Options *);

