#ifndef _grid_generator_h
#define _grid_generator_h

PetscErrorCode initializeRegularGrid(GridData *, Options *);
PetscErrorCode initializeSubductionGrid(GridData *, Options *);
PetscErrorCode initializeIrregularGridConstantInnerOuter(GridData *, Options *);
PetscErrorCode destroyGrid(GridData *);
PetscErrorCode initializeGrid(GridData *, Options *);

#endif
