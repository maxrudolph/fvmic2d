

#ifdef _PHASE_BA_H

#else

typedef struct{
  /* oxides, H = H20, C = CO2 */
  PetscScalar Si, Ti, Al, Fe, Mn, P, Mg, Ca, Na, K, H, F;
} Norm;

typedef struct{
  Norm melt;
  PetscScalar phi; /* solid fraction */
  PetscScalar phiplag; /* solid fraction plagioclase */
  PetscScalar rho;
} Composition;

#define NPMAX 1
#define NTMAX 20

typedef struct{
  PetscInt nT;/* number of temperatures */
  PetscInt nP;/* number of pressures */
  PetscScalar Tc[NTMAX]; /* Temperatures */
  PetscScalar PPa[NPMAX]; /* Pressure */
  //solids SiO2 TiO2 Al2O3 FeO MnO MgO CaO Na2O K2O H2O Plag Dens(noH2O) Dens Ves
  PetscScalar cdata[NTMAX][15];
} CompositionLookupTableBA;

void printComposition( Composition * );
void lookupCompositionBA( CompositionLookupTableBA *, PetscScalar , Composition *);
CompositionLookupTableBA loadCompositionLookupTableBAFromFile( char * );

PetscScalar giordanoViscosity( Composition *, PetscScalar );
PetscScalar einsteinRoscoe( PetscScalar , PetscScalar );

#define _PHASE_BA_H

#endif



