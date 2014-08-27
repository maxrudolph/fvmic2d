PetscErrorCode initialConditionsBending( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsSandbox( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode sandboxMobileWallAndTimesteps( MarkerSet *, GridData *, Materials *, Options *, PetscInt , PetscScalar);
PetscErrorCode initialConditionsOOPCouette( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsConvectionGerya( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsConvectionBarr( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsMagmaMixing( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsMagmaMixingBA( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsConvectionTuran( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsDilationTest( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsDilationTest2( MarkerSet *, Options *, Materials *, GridData *);
PetscErrorCode initialConditionsInclusion(  MarkerSet *, Options *, Materials *, GridData *);
