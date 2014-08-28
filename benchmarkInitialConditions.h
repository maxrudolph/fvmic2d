PetscInt in_slab( PetscScalar, PetscScalar, PetscScalar);
PetscInt in_plate( PetscScalar, PetscScalar, PetscScalar);
PetscScalar slab_inflow_temperature(PetscScalar, PetscScalar, PetscScalar);
PetscErrorCode initialConditionsVanKeken( MarkerSet *, Options *, Materials *, GridData *);
PetscScalar slab_depth(PetscScalar, PetscScalar);
PetscScalar slab_x(PetscScalar, PetscScalar);


