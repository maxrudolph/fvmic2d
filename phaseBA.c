/* This file contains an implementation of phase assemblages on markers */
#include "petscksp.h"
#include "fdcode.h"
#include "phaseBA.h"


/* Table 1, Giordano et al EPSL 2008 */
const PetscScalar b[10] = { 159.6, -173.3, 72.1, 75.7, -39.0, -84.1, 141.5, -2.43, -0.91, 17.6 };
const PetscScalar c[7] = {2.75, 15.7, 8.3, 10.2, -12.3, -99.5, 0.30};
const PetscScalar A = -4.55;

PetscInt bisect(PetscScalar *, PetscScalar, PetscInt);
PetscScalar calculateCompositionViscosityBA( Materials *, Marker *);

PetscScalar giordanoViscosity( Composition *comp, PetscScalar T ){/* melt viscosity */
  /* T is temperature in Kelvin, comp is composition structure */
  /* log(eta) = A + B/(T(K)-C)*/
  
  Norm *melt = &comp->melt; 

  PetscScalar V = melt->H + melt->F;
  PetscScalar TA = melt->Ti + melt->Al;
  PetscScalar FM = melt->Fe + melt->Mn + melt->Mg;
  PetscScalar NK = melt->Na + melt->K;

  /* calculate B */
  PetscScalar B = 0.0;
  B += b[0]*( melt->Si + melt->Ti);
  B += b[1]*melt->Al;
  B += b[2]*(melt->Fe + melt->Mn + melt-> P);
  B += b[3]*melt->Mg;
  B += b[4]*melt->Ca;
  B += b[5]*(melt->Na + melt->H + melt->F);
  B += b[6]*(melt->F + melt->H + log(1.0 + melt->H));
  /* cross terms */
  B += b[7] * (melt->Si+melt->Ti) * (FM );    /*b11 */
  B += b[8] * (melt->Si+TA+melt->P) * (NK+melt->H); /*b12 */
  B += b[9] * (melt->Al)*(NK);          /*b13 */
  
  /* calculate C */
  PetscScalar C = 0.0;
  C += c[0]*(melt->Si);
  C += c[1]*TA;
  C += c[2]*FM;
  C += c[3]*melt->Ca;
  C += c[4]*NK;
  C += c[5]*(log(1+V));
  C += c[6]*(melt->Al+FM+melt->Ca-melt->P)*(NK+V);/* actually c11 in Giordano et al. */

  PetscScalar logeta = A + B/(T-C);
#ifdef DEBUG
  printf("A\t\tB\t\tC\n%e\t%e\t%e\nLog eta=%e\n",A,B,C,logeta);
#endif

  return exp(logeta);
}

/* Einstein-Roscoe viscosity to account for presence of solid particles */
PetscScalar einsteinRoscoe( PetscScalar eta0, PetscScalar phisolid){
  /* returns effective viscosity given fluid viscosity and solid volume fraction */
  PetscScalar eta = eta0*pow(1.0-phisolid/0.6,-2.5);
  return eta;
}

PetscScalar calculateCompositionViscosityBA( Materials *materials, Marker *marker ){
  PetscScalar Tk = marker->T;
  Composition c;
  /* get composition using lookup table for this material */
  lookupCompositionBA( &(materials->compositionLookupTable[ marker->Mat ]), Tk-273.15, &c);
  PetscScalar eta = giordanoViscosity( &c, Tk );/* returns viscosity in Pa-s */
  marker->rho = c.rho;
  /* modify viscosity to account for solid fraction */
  eta = einsteinRoscoe( eta, c.phi );
  /* ADD YIELDING BEHAVIOR HERE */

  /* return viscosity */
  return eta;
}


/* interpolate composition from lookup table */
void lookupCompositionBA( CompositionLookupTableBA *clt, PetscScalar T, Composition *c){
  /* interpolate temperatures in lookup table to get composition c at temperature T */
  /* step 1: bisect lookup table temperatures to find temperatures bounding current temperature */
  PetscInt Tmi, Tpi;
  if( T > clt->Tc[ clt->nT-1]){
    Tmi = clt->nT-2;
    Tpi = clt->nT-1;
    T = clt->Tc[ clt->nT-1 ] ;
  }else if(T < clt->Tc[0]){
    Tmi = 0;
    Tpi = 1;
    T = clt->Tc[ 0 ];
  }else{
    Tmi = bisect( clt->Tc, T, clt->nT );
    Tpi = Tmi+1;
  }
  /*   printf("chosen T=%e, T below = %e, T above = %e\n",T,clt->Tc[Tmi],clt->Tc[Tpi]); */
  Norm *melt = &(c->melt);
  PetscScalar x1 = 1.0-(T-clt->Tc[ Tmi])/(clt->Tc[Tpi]-clt->Tc[Tmi]);
  PetscScalar x2 = 1.0-x1;
  
  PetscScalar *c1 = &(clt->cdata[Tmi][0]);
  PetscScalar *c2 = &(clt->cdata[Tpi][0]);
  c->phi = x1*c1[0] + x2*c2[0];
  //solids 
  //SiO2 
  melt->Si = x1*c1[1] + x2*c2[1]; 
  //TiO2 
  melt->Ti = x1*c1[2] + x2*c2[2];
  //Al2O3 
  melt->Al = x1*c1[3] + x2*c2[3];
  //FeO 
  melt->Fe = x1*c1[4] + x2*c2[4];
  //MnO 
  melt->Mn = x1*c1[5] + x2*c2[5];
  //MgO 
  melt->Mg = x1*c1[6] + x2*c2[6];
  //CaO 
  melt->Ca = x1*c1[7] + x2*c2[7];
  //Na2O 
  melt->Na = x1*c1[8] + x2*c2[8];
  //K2O 
  melt->K = x1*c1[9] + x2*c2[9];
  //H2O 
  melt->H = x1*c1[10] + x2*c2[10];
  //Plag 11
  c->phiplag = x1*c1[11] + x2*c2[11];
  //Dens(noH2O) 12 
  //Dens 13
  c->rho = x1*c1[13] + x2*c2[13];
  //Ves 14
  melt->P = 0.0;
  melt->F = 0.0;
}

/* load a composition lookup table from disk */
CompositionLookupTableBA loadCompositionLookupTableBAFromFile( char * filename ){
  CompositionLookupTableBA lt;
  FILE *fh;
  fh = fopen(filename,"r");
  if( fh == NULL ){
    printf("invalid composition file\n");
    abort();
  }
  /* get number of temperatures */
  fscanf(fh,"%d\n",&(lt.nT));
/*   printf("reading %d temperatures\n",lt.nT); */
  /* read temperatures */
  PetscInt iT;
  for(iT=0;iT < lt.nT-1; iT++){
    fscanf(fh,"%lf,",&lt.Tc[iT]);
  }
  fscanf(fh,"%lf\n",&lt.Tc[lt.nT-1]);
/*   for(iT=0;iT<lt.nT;iT++){ */
/*     printf("%lf, ",lt.Tc[iT]); */
/*   } */
/*   printf("\n"); */

  /* read melt composition and other information */
  const PetscInt nc = 15;
  for(iT=0;iT < lt.nT; iT++){
    PetscInt ic;
    for(ic=0;ic<nc;ic++){
      if(ic<nc-1){
	fscanf(fh,"%lf\t",&(lt.cdata[iT][ic]));
      }else{
	fscanf(fh,"%lf\n",&(lt.cdata[iT][ic]));
      }
      /* convert oxides into mol. % instead of mol. fraction */
      if( ic>0 && ic<=10){
	lt.cdata[iT][ic] *= 100.0;
      }
    }


    /* convert crystallinity into volume fraction rather than volume percent */
    lt.cdata[iT][0] /= 100.0;
  }
/*   for(iT=0;iT< lt.nT; iT++){ */
/*     PetscInt ic; */
/*     for(ic=0;ic<nc;ic++){ */
/*       printf("%lf\t",lt.cdata[iT][ic]); */
/*     } */
/*     printf("\n"); */
/*   } */
/*   fclose(fh); */

  return lt;
}


PetscInt bisect(PetscScalar *X, PetscScalar x, PetscInt N){
  PetscInt idxmin=0;
  PetscInt idxmax= N - 1;
  while( (idxmax - idxmin) >1){
    PetscInt idxc = (idxmax+idxmin)/2;
    if(idxc == idxmin){idxc++;}else if(idxc == idxmax){idxc--;}
    if( x >= X[idxc]){ 
      idxmin=idxc;
    }else{
      idxmax=idxc;
    }
  }
  return idxmin;
}

void printComposition( Composition *c ){
  Norm *m=&(c->melt);

  printf("Si = %lf\n",m->Si);
  printf("Ti = %lf\n",m->Ti);
  printf("Al = %lf\n",m->Al);
  printf("Fe = %lf\n",m->Fe);
  printf("Mn = %lf\n",m->Mn);
  printf("P  = %lf\n",m->P);
  printf("Mg = %lf\n",m->Mg);
  printf("Ca = %lf\n",m->Ca);
  printf("Na = %lf\n",m->Na);
  printf("K  = %lf\n",m->K);
  printf("H  = %lf\n",m->H);
  printf("F  = %lf\n",m->F);
}
