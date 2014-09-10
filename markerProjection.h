

typedef enum {
  ARITHMETIC,
  ARITHMETIC_LOCAL,
  HARMONIC,
  GEOMETRIC
} PROJECTION_METHOD;

PetscErrorCode projectMarkerFieldToNodes(GridData *, MarkerSet *, PetscScalar *, PetscScalar *, Vec, const PetscInt, const PetscInt, const PROJECTION_METHOD);
PetscErrorCode projectMarkersNodesAll2( MarkerSet *, GridData *, NodalFields *, Materials *, Options *);
