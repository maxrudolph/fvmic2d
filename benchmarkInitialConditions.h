#ifndef _BENCHMARK_IC_H
#define _BENCHMARK_IC_H
PetscScalar plate_depth(PetscScalar, Options *);
PetscInt in_slab( PetscScalar, PetscScalar, Options *);
PetscInt in_plate( PetscScalar, PetscScalar, Options *);
PetscInt in_either( PetscScalar, PetscScalar, Options *);
PetscScalar slab_inflow_temperature(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
PetscErrorCode initialConditionsVanKeken( MarkerSet *, Options *, Materials *, GridData *);
PetscScalar slab_depth(PetscScalar, PetscScalar);
PetscScalar slab_x(PetscScalar, PetscScalar);
PetscScalar plate_geotherm(PetscScalar, PetscScalar, Options *);


PetscScalar mantle_temperature();
Marker new_slab_marker( PetscScalar, PetscScalar, Options *, Materials * );
Marker new_inflow_marker( PetscScalar, PetscScalar, Options *, Materials * );

#endif
