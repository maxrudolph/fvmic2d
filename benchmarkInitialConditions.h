PetscScalar plate_depth(PetscScalar);
PetscInt in_slab( PetscScalar, PetscScalar, PetscScalar);
PetscInt in_plate( PetscScalar, PetscScalar, PetscScalar);
PetscInt in_either( PetscScalar, PetscScalar, PetscScalar);
PetscScalar slab_inflow_temperature(PetscScalar, PetscScalar, PetscScalar);
PetscErrorCode initialConditionsVanKeken( MarkerSet *, Options *, Materials *, GridData *);
PetscScalar slab_depth(PetscScalar, PetscScalar);
PetscScalar slab_x(PetscScalar, PetscScalar);
PetscScalar plate_geotherm( PetscScalar);


PetscScalar mantle_temperature();
Marker new_slab_marker( PetscScalar, PetscScalar, Options *, Materials * );
Marker new_inflow_marker( PetscScalar, PetscScalar, Options *, Materials * );
