#include<stdio.h>
#include<stdlib.h>
#include "petscksp.h"
#include "petscdmda.h"
#ifdef COMPBA
#include "phaseBA.h"
#endif

//#define DMCreateMatrix DMGetMatrix /* uncomment for backwards compatibility for petsc-3.2 */

#define MAXMAT 8 /* maximum number of materials allowed*/
/* all structure definitions go here*/
typedef struct {
  DM da, vda; /*da is for temperature and z-velocity, vda is for vx,vy,p, tda is for texture information (3x3 symmetric matrices)*/
#ifdef TEXTURE
  DM tda;
#endif
  PetscScalar LX, LY; /* x,y dimensions*/
  PetscInt NX, NY; /* number of nodes*/
  PetscInt mp,np; /* number of processors in y,x directions*/
  PetscScalar *x, *y, *xc, *yc;

  PetscInt xperiodic, yperiodic; /* flags that can be used globally */
} GridData;


typedef struct{
  PetscScalar T11, T22, T33, T23, T13, T12;
} Tensor33s;

typedef struct{
  PetscScalar T11, T12, T21, T22;
} Tensor22;

typedef struct{/* 6x6 matrix entries for storing full viscoplasticity tensor - only needed in texture subroutine*/
  PetscScalar T[6][6];
} Tensor66s;


typedef struct {
  /* New Organization: Fields defined on xy nodes */
  
  /* fields defined on xz nodes */

  /* fields defined on yz nodes */

  /* fields defined at cell centers */

  /* fields defined at gridline intersections */

  Vec rho, rhodot, etaS, etaN; /* mechanical properties */
  Vec etavx, etavy; /* mechanical properties stored at mid-side nodes for problems with out-of-plane motion */
  Vec thisT, lastT, kThermal, Cp; /*thermal*/
  PetscScalar boundaryT; /* boundary condition value*/

  Vec edotxx, edotxy;          /* nodal strain rate*/
  Vec edotyy, edotzz, edotxz, edotyz;

  Vec soxx, soxy;/* nodal stresses*/
  Vec soyy, sozz, soxz, soyz; /* additional stress components - needed for 2.5d simulations*/
  Vec dsxx, dsyy, dszz, dsxy, dsxz, dsyz;/* stress changes*/
  Vec muS, muN;  /*elastic mu*/
  Vec muvx, muvy; /* mid-side properties */
  Vec vx,vy,vz,p;
  Vec eii,sii,wxy,wxz,wyz;/* stress and strain-rate 2nd invariants, vorticity*/
  Vec ha;         /* adiabatic heating*/
  Vec nodalHeating; /* shear + adiabatic heating */
  Vec VPTensorC; /*viscoplasticity tensor, cell-centered*/
  Vec VPTensorB; /*viscoplasticity tensor, at basic nodes*/
  Vec strainRateResidual;
} NodalFields;


typedef struct {
/* material properties. these are used in conjunction with material id numbers on the markers*/
  PetscInt nMaterials;
  PetscScalar materialEta[MAXMAT];// = {1.0e20, 1.0e22};
  PetscScalar materialRho[MAXMAT];// = {3200.0, 3300.0};
  PetscScalar materialkThermal[MAXMAT];// = {3.0,10.0};
  PetscScalar materialCp[MAXMAT];// = {1000.0,1100.0};
  PetscScalar materialAlpha[MAXMAT];/* thermal expansivity K^-1*/
  PetscScalar materialMu[MAXMAT];// elastic mu;
  PetscScalar materialFriction[MAXMAT];
  PetscScalar materialCohesion[MAXMAT];
  PetscInt hasPlasticity[MAXMAT];
  PetscScalar gamma[MAXMAT][2];/* initial strain */
  PetscScalar C[MAXMAT][2];/* initial, final cohesion */
  PetscScalar F[MAXMAT][2];/* initial, final friction */
  PetscScalar binghamYieldStress[MAXMAT]; /* yield stress for solid to fluid transition in bingham fluid */
  
  PetscInt hasDamage[MAXMAT];
  PetscInt hasDilation[MAXMAT];

  PetscScalar hayhurstAlpha[MAXMAT];
  PetscScalar hayhurstBeta[MAXMAT];
  PetscScalar damagem[MAXMAT];
  PetscScalar damagek[MAXMAT];
  PetscScalar damager[MAXMAT];
  PetscScalar damageB[MAXMAT];
  PetscScalar damageAlpha3[MAXMAT];

  PetscInt hasEtaT[MAXMAT];
  PetscScalar QonR[MAXMAT];/* activation energy / R */
  PetscScalar Tref[MAXMAT];/* temperature to which flow law is referenced*/
#ifdef COMPBA
  CompositionLookupTableBA compositionLookupTable[MAXMAT];
#endif
} Materials;

typedef struct {/* an individual marker */
  PetscScalar X;
  PetscScalar Y;
  PetscScalar Z;
  PetscScalar VX;
  PetscScalar VY;
  PetscScalar VZ;
  PetscInt Mat;
  PetscScalar T, Tdot;
  /*PetscScalar *sxx, *syy, *szz, *sxy, *sxz, *syz;*/
  Tensor33s s;
  PetscScalar eta; /* effective viscosity */
  PetscScalar mu;  /* elastic modulus */
  PetscScalar rho; /* density */
  PetscScalar rhodot; /* drho/dt*/
  PetscScalar p;/* pressure */
  /*PetscScalar *exx, *eyy, *ezz, *exy, *exz, *eyz;*/ /* strain rate*/
  Tensor33s e;
  /*   PetscScalar *Exx, *Exy, *Eyy, *Exz, *Eyz; */ /* total strain*/
  Tensor33s E;
  PetscScalar Eii; /* total strain second invariant */
  PetscScalar wxy, wxz, wyz; /* vorticity vector */
  PetscScalar D, Ddot;/* scalar damage, damage rate*/
  PetscInt cellX, cellY, cpu;/* I,J cell (using upper right basic node) of each marker, cpu that I belong to*/
  PetscScalar strainRateResidual;
  PetscInt isYielding;
  /* texture stuff:*/
#ifdef TEXTURE
  Tensor33s emm;/* strain rate predicted using micro-macro model */
  TextureInfo texture;
  PetscScalar residual;
#endif
} Marker;

typedef struct{ /* data structure containing information about all markers on a CPU */
  PetscInt nMark; /* global number of markers*/
  PetscInt maxMark; /* maximum number of local markers allowed*/
  /*  PetscInt nMarkLocal; */ /*not used yet */
  PetscInt nOut; /* number of out-of-domain markers*/
  Marker *markers;

} MarkerSet;


typedef struct {
  PetscInt    type[3];   /* x,y,z BC type */
  /* 0 = prescribed velocity, 1 = free slip, 2 = periodic */
  PetscScalar value[3];  /* x,y,z BC Value */
} BCInfo;

typedef struct{
  BCInfo mechBCLeft;
  BCInfo mechBCRight;
  BCInfo mechBCTop;
  BCInfo mechBCBottom;
} BoundaryValues;


typedef struct {
  PetscInt NX,NY;/* number of nodes in x,y,directions*/
  PetscScalar LX,LY;/* model dimensions (x,y)*/
  PetscInt NMX;/* number of markers per cell, x*/
  PetscInt NMY;/* markers per cell y*/
  PetscScalar gridRatio; /* contrast between largest and smallest cell size */
  PetscInt nTime;
  PetscInt restartStep;
  PetscScalar totalTime;
  PetscInt maxNumPlasticity;
  PetscInt plasticDelay;
  PetscScalar maxPorosity;
  PetscInt dilationModel;
  PetscScalar Dcrit;

  PetscScalar Tperturb; /* magnitude of initial temperature perturbation */
  PetscInt shearHeating, adiabaticHeating;/* toggle these quantities on/off */

  PetscScalar gy;
  PetscScalar gx;
  PetscScalar gz;
  /* PetscScalar vbx,vby,vbz; */
  BCInfo mechBCLeft, mechBCRight, mechBCTop, mechBCBottom;
  BCInfo thermalBCLeft, thermalBCRight, thermalBCTop, thermalBCBottom;

  /* Variables related to oscillatory boundary conditions*/
  PetscInt oscillatoryX, oscillatoryY, oscillatoryZ;/* flags for oscillating boundary conditions. not all are implemented */
  PetscScalar oscillationPeriod;/* period of oscillation when using oscillating boundary conditions */
  /* options related to solution and time step limits */
  PetscScalar fractionalEtamin; /* fractional minimum viscosity */
  PetscScalar etamin;/* minimum allowed viscosity*/
  PetscScalar etamax;/* maximum allowed viscosity */
  PetscScalar dtMax;/*maximum timestep in seconds*/
  PetscScalar displacementStepLimit; /* maximum fractional cell size that a marker may move in one timestep */
  PetscScalar maxTempChange;
  PetscInt saveInterval;
  PetscScalar subgridStressDiffusivity;
  PetscScalar subgridTemperatureDiffusivity;
  PetscInt nMonte;
  PetscScalar maxMarkFactor;
  PetscInt minMarkers, maxMarkers;
  PetscInt doMonte;/* run monte-carlo simulations */
  PetscInt textureDevelopment; /* allow texture to change */
  PetscScalar grainSize;
  PetscScalar slabAngle;
  PetscScalar slabVelocity;
} Options;

typedef struct{/* a structure holding everything related to solving a linear system */
  Mat lhs;
  Mat pc;
  Vec rhs;
  Vec solution;
  MatNullSpace ns;
  char label[80];
} LinearSystem;

typedef struct{/* a structure that holds all of the global variables */
  Options     options;
  Materials   materials;
  GridData    grid;
  NodalFields nodal_fields;
  MarkerSet   markerset;
  /* mechanical system */
  LinearSystem mech_system;

  /* out-of-plane system */
  LinearSystem z_system;
  /* thermal system */
  LinearSystem thermal_system;

  /* thermo-mechanical coupling */
  Vec          nodal_heating;
} Problem;


/* define ordering of fields in each block of global matrix */
#define DOF_U 0
#define DOF_V 1
#define DOF_W 0
#define DOF_P 2

/* define output file formats */
#define OUTPUT_NODAL_FIELDS_ASCII_MATLAB 1
#define OUTPUT_NODAL_FIELDS_BINARY 0
#define OUTPUT_NODAL_FIELDS_BINARY_MATLAB 2
#define OUTPUT_NODAL_FIELDS_ASCII_VTK 3
#define OUTPUT_NODAL_FIELDS_BINARY_VTK 4

/* define multiple level of debug statements */
/* level 4 should be things you never want to see, level 1 things that might be very useful */

#ifdef debug
#define verbose
#endif

#ifdef verbose
#define debug4
#endif

#ifdef debug4
#define verbose
#define debug3
#define debug2
#define debug1
#endif

#ifdef debug3
#define debug2
#define debug1
#endif

#ifdef debug2
#define debug1
#endif

/* Error codes for Petsc error handling */
#define ERR_MICROMACRO_NAN PETSC_ERR_MAX_VALUE+1
#define FIX_PRESSURE PETSC_TRUE

/* end define debug statements */

