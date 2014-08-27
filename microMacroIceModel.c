#include "fdcode.h"
#include "texture.h"
#include "math.h"

#ifndef NT
#define NT 8 /* */
#endif

#define ANG_MAX M_PI/10.0
//#ifdef verbose
//#undef verbose
//#endif
//#define verbose

int isbad(PetscScalar a){
  if( isnan(a) || isinf(a)){
    return 1;
  } else{
    return 0;
  }
}

/* SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
int dgesv_(int *, int *, double *, int *, int *, double *, int *, int *);
/*      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
     $                   WORK, IWORK, INFO )*/
//int dgesvx_( char *, char *, int *, int *, double *, int *, int *, double *, 
void solve33(PetscScalar *,PetscScalar *);
void print33s(Tensor33s *);
PetscScalar max(PetscScalar, PetscScalar);

PetscInt idxmod(PetscInt);

const PetscScalar zeta = 6.0;/* closure coefficients for Thorsteinsson 2002 model*/
const PetscScalar xi = 1.0;

/* Three types of NNI are used, (z, x) = (1,0), (6,1) and (1,1), which will be called the no- NNI, mild-NNI, and full-NNI cases, respectively. */

const PetscInt nt=NT;
const PetscInt NS=3;/* number of slip systems - note that NS=3 activates only basal slip*/
/* crystallographic tensors */
/* assume that c-axis is collinear with z-axis */
/*const PetscScalar bs[3][3] = {{1.0,0.0,0.0},{-0.5,0.8660254037844386,0.0},{-0.5,-0.8660254037844386,0.0}};*/ /* burgers vectors*/

const PetscScalar bs[12][3] = {{1.,0.,0},{-0.5,0.8660254037844386,0.},{-0.5,-0.8660254037844386,0.},{1.,0.,0}, \
			       {-0.5000000000000001,0.8660254037844388,0.},{-0.5000000000000001,-0.8660254037844388,0.},\
			       {0.261317712810192,0.45261555550494503,-0.8525562807736508},\
			       {-0.261317712810192,0.45261555550494503,-0.8525562807736508},\
			       {-0.522635425620384,-2.7755575615628914e-17,-0.8525562807736508},\
			       {-0.261317712810192,-0.45261555550494503,-0.8525562807736508},\
			       {0.261317712810192,-0.45261555550494503,-0.8525562807736508},\
			       {0.522635425620384,2.7755575615628914e-17,-0.8525562807736508}};


/* Schmid tensors */
/*const PetscScalar rs[3][3][3] = {{{0., 0., 0.5}, {0., 0., 0.}, {0.5, 0., 0.}}, {{0., 0., -0.25}, {0., \
  0., 0.433013}, {-0.25, 0.433013, 0.}}, {{0., 0., -0.25}, {0., 0., -0.433013}, {-0.25, -0.433013, 0.}}};*/

const PetscScalar rs[12][3][3] = {{{0,0,1.},{0,0,0.},{0,0,0}},{{0,0,-0.5},{0,0,0.8660254037844386},{0,0,0.}},
   {{0,0,-0.5},{0,0,-0.8660254037844386},{0,0,0.}},{{0.,1.,0},{0.,0.,0},{0,0,0}},
   {{0.4330127018922194,0.25000000000000006,0.},{-0.7500000000000001,-0.4330127018922194,0.},{0.,0.,0.}},
   {{-0.4330127018922194,0.25000000000000006,0.},{-0.7500000000000001,0.4330127018922194,0.},{0.,0.,0.}},
   {{1.3531351129322665e-7,2.3436987651041358e-7,0.0002254998146740371},
    {2.3436987651041366e-7,4.059405338796799e-7,0.00039057713611279814},
    {-0.0006000599170431337,-0.001039334263904273,-0.9999991781553751}},
   {{1.353135112932266e-7,-2.3436987651041366e-7,-0.00022549981467403712},
    {-2.3436987651041353e-7,4.0594053387968e-7,0.00039057713611279814},
    {0.0006000599170431334,-0.0010393342639042733,-0.9999991781553751}},
   {{5.412540451729065e-7,-4.889751300619499e-23,-0.0004509996293480742},
    {3.2529316592229105e-23,-2.938735877055719e-39,-2.7105034792589583e-20},
    {0.0012001198340862672,-1.0842020622216183e-19,-0.9999991781553751}},
   {{1.3531351129322665e-7,2.3436987651041358e-7,-0.0002254998146740371},
    {2.3436987651041366e-7,4.059405338796799e-7,-0.00039057713611279814},
    {0.0006000599170431337,0.001039334263904273,-0.9999991781553751}},
   {{1.353135112932266e-7,-2.3436987651041366e-7,0.00022549981467403712},
    {-2.3436987651041353e-7,4.0594053387968e-7,-0.00039057713611279814},
    {-0.0006000599170431334,0.0010393342639042733,-0.9999991781553751}},
   {{5.412540451729065e-7,-4.889751300619499e-23,0.0004509996293480742},
    {3.2529316592229105e-23,-2.938735877055719e-39,2.7105034792589583e-20},
	{-0.0012001198340862672,1.0842020622216183e-19,-0.9999991781553751}}};

/* CBS is the coefficient for basal slip*/
#define CBS 33.5e6
const PetscScalar gd0s[12] = {1,1,1,1,1,1,1,1,1,1,1,1};/* reference strain rate*/
PetscScalar tau0s[12] = {CBS,CBS,CBS,CBS,CBS,CBS,CBS,CBS,CBS,CBS,CBS,CBS};/* reference shear stress*/
#define A01  7e-6
const PetscScalar A0[12] = {7.35e-4,7.35e-4,7.35e-4,A01,A01,A01,A01,A01,A01,A01,A01,A01 };/* rheological constant A */
const PetscScalar n1s[12] = {2,2,2,3,3,3,3,3,3,3,3,3};/*power law exponent*/
const PetscScalar QonRs[12] = {7216,7216,7216,7000,7000,7000,7000,7000,7000,7000,7000,7000};

/* settingss related to newton iteration */
const PetscScalar newtonAbsTol = 1e-40; /* determined by looking at small residuals */

/* These subroutines evaluate the micro-macro response of an ice polycrystal according to Thorsteinsson 2002 */ 

PetscScalar sT[NT][NT][NT];
PetscScalar softnessE[NT][NT][NT];

PetscScalar max(PetscScalar X,PetscScalar Y){
  if(X>Y) return X;
  else return Y;
}

PetscErrorCode formViscoplasticMNewton(MarkerSet *markerset, Materials *materials, Options *options){
  PetscFunctionBegin;
  PetscErrorCode ierr=0;
  PetscInt m;
  for(m=0;m<markerset->nMark;m++){
#ifdef debug1
    printf("[%d] Now on marker %d/%d\n",rank,m,markerset->nMark);
#endif
    ierr=formViscoplasticMNewtonMarker( &(markerset->markers[m]), materials, options);CHKERRQ(ierr);    
  }/* end loop over markers */
  PetscFunctionReturn(ierr);
}

PetscErrorCode formViscoplasticMNewtonMarker(Marker *marker, Materials *materials, Options *options){
  PetscFunctionBegin;
  PetscErrorCode ierr=0;
  PetscInt nNewton = 100;

  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  /*   for(m=0;m<markerset->nMark;m++){ */
  /*     marker = &(markerset->markers[m]); */
  /* check marker in bounds? */
  if( marker->cellX != -1 ){
    /* check marker has texture */
    if( materials->hasTexture[(PetscInt) marker->Mat ]==1 ){
      Tensor33s D0;
      PetscInt iNewton;

      for(iNewton=0;iNewton<nNewton;iNewton++){
	PetscScalar lastres;
	const PetscScalar newtonRelTol = 1e-6;
	const PetscScalar Kres=1.0;
	static PetscScalar res;

	marker->s.T33=0.0;
	marker->s.T23=0.0;
	marker->s.T13=0.0;

	formViscoplasticMMarker( marker,options, &D0, 0);
	/* compute residual initially - stop if newton iteration not required */
	res=sqrt(D0.T33*D0.T33*Kres*Kres + D0.T13*D0.T13*Kres*Kres +D0.T23*D0.T23*Kres*Kres);
	if( res < newtonAbsTol ) break;
#ifdef verbose
	printf("res0 = %e\n",res);
#endif

	/* these second invariants are left over from an old version of this routine in which the residual strain rates were normalized by the second derivative. This turned out to be unnecessary. */
	PetscScalar D0ii = 1.0;//sqrt( D0.T11*D0.T11 + D0.T22*D0.T22 + D0.T33*D0.T33 + 2*D0.T12*D0.T12 + 2*D0.T13*D0.T13 + 2*D0.T23*D0.T23);
#ifdef verbose
	printf("Dzz, Dxz, Dyz = %e, %e, %e\n",D0.T33,D0.T13,D0.T23);
#endif
	/* compute derivatives numerically */
	
	Tensor33s D1,D2;/* strain rate tensors for computation of derivatives*/
	PetscScalar dsig;
	/* see Numerical Recipes chapter 5.7 */
	const PetscScalar dsig1 = 1e8;/* value to guess for dsig when sig~=0 */
	/* dsig1 = 1e-4,1 causes NaNs - experiment with this until NaNs go away */
	const PetscScalar sf = 1e-6;/* scale factor equal to assumed cube root of relative error in stain-rate evaluation function */
	Marker sp = *marker;
	Marker sm = *marker;
	/* choose correct derivative step size */
	dsig = max( sf*fabs(sp.s.T33) , sf*dsig1 );/* guess at correct step size to take*/
	PetscScalar tmp; /* declared volatile to prevent compiler optimizations that result in a reduction of accuracy*/
	PetscScalar dsigp,dsigm;
      
	/* **** 33 component **** */

	tmp = sp.s.T33 + dsig; dsigp = tmp - sp.s.T33; sp.s.T33 = tmp;
	tmp = sm.s.T33 - dsig; dsigm = sm.s.T33 - tmp; sm.s.T33 = tmp;
	dsig = dsigp + dsigm;
#ifdef verbose
	printf("dsigp = %e, dsigm = %e, dsig = %e\n",dsigp,dsigm,dsig);
#endif
	formViscoplasticMMarker(&sp,options,&D1,0);
	formViscoplasticMMarker(&sm,options,&D2,0);
#ifdef verbose
	printf("Dplus =\n");
	print33s( &D1 );
#endif
	PetscScalar D1ii = 1.0;//sqrt( D1.T11*D1.T11 + D1.T22*D1.T22 + D1.T33*D1.T33 + 2*D1.T12*D1.T12 + 2*D1.T13*D1.T13 + 2*D1.T23*D1.T23);
	PetscScalar D2ii = 1.0;//sqrt( D2.T11*D2.T11 + D2.T22*D2.T22 + D2.T33*D2.T33 + 2*D2.T12*D2.T12 + 2*D2.T13*D2.T13 + 2*D2.T23*D2.T23);
	/*       D1ii = 1.0; D2ii = 1.0; */
	PetscScalar J11,J12,J13,J21,J22,J23,J31,J32,J33;
	J11 = (D1.T33/D1ii-D2.T33/D2ii)/(dsig)*Kres;
	J12 = (D1.T23/D1ii-D2.T23/D2ii)/(dsig)*Kres;
	J13 = (D1.T13/D1ii-D2.T13/D2ii)/(dsig)*Kres;

	/* **** 23 component **** */
	sp = *marker; sm = *marker;
	/* choose correct derivative step size */
	dsig = max( sf*fabs(sp.s.T23) , sf*dsig1 );/* guess at correct step size to take*/
	tmp = sp.s.T23 + dsig; dsigp = tmp - sp.s.T23; sp.s.T23 = tmp;
	tmp = sm.s.T23 - dsig; dsigm = sm.s.T23 - tmp; sm.s.T23 = tmp;
	dsig = dsigp + dsigm;
	formViscoplasticMMarker(&sp,options,&D1,0);
	formViscoplasticMMarker(&sm,options,&D2,0);

	//D1ii = sqrt( D1.T11*D1.T11 + D1.T22*D1.T22 + D1.T33*D1.T33 + 2*D1.T12*D1.T12 + 2*D1.T13*D1.T13 + 2*D1.T23*D1.T23);
	//D2ii = sqrt( D2.T11*D2.T11 + D2.T22*D2.T22 + D2.T33*D2.T33 + 2*D2.T12*D2.T12 + 2*D2.T13*D2.T13 + 2*D2.T23*D2.T23);
	//D1ii = 1.0; D2ii = 1.0;
	J21 = (D1.T33/D1ii-D2.T33/D2ii)/(dsig)*Kres;
	J22 = (D1.T23/D1ii-D2.T23/D2ii)/(dsig)*Kres;
	J23 = (D1.T13/D1ii-D2.T13/D2ii)/(dsig)*Kres;

	/* **** 13 component **** */
	sp = *marker; sm = *marker;
	dsig = max( sf*fabs(sp.s.T13) , sf*dsig1 );/* guess at correct step size to take*/
	tmp = sp.s.T13 + dsig; dsigp = tmp - sp.s.T13; sp.s.T13 = tmp;
	tmp = sm.s.T13 - dsig; dsigm = sm.s.T13 - tmp; sm.s.T13 = tmp;
	dsig = dsigp + dsigm;
	formViscoplasticMMarker(&sp,options,&D1,0);
	formViscoplasticMMarker(&sm,options,&D2,0);

	//D1ii = sqrt( D1.T11*D1.T11 + D1.T22*D1.T22 + D1.T33*D1.T33 + 2*D1.T12*D1.T12 + 2*D1.T13*D1.T13 + 2*D1.T23*D1.T23);
	//D2ii = sqrt( D2.T11*D2.T11 + D2.T22*D2.T22 + D2.T33*D2.T33 + 2*D2.T12*D2.T12 + 2*D2.T13*D2.T13 + 2*D2.T23*D2.T23);
	//D1ii = 1.0; D2ii = 1.0;
	J31 = (D1.T33/D1ii-D2.T33/D2ii)/(dsig)*Kres;
	J32 = (D1.T23/D1ii-D2.T23/D2ii)/(dsig)*Kres;
	J33 = (D1.T13/D1ii-D2.T13/D2ii)/(dsig)*Kres;

		
	PetscScalar s1=J11, s2=J22, s3=J33;
	if( fabs(s1) < 1.0 ) s1 = 1.0;
	if( fabs(s2) < 1.0 ) s2 = 1.0;
	if( fabs(s3) < 1.0 ) s3 = 1.0;

#ifdef verbose
	printf("Scale factors s1,s2,s3 = %e,%e,%e\n",s1,s2,s3);
#endif
	/* 	PetscScalar J[9] = {J11/s1,J21/s2,J31/s3,J12/s1,J22/s2,J32/s3,J13/s1,J23/s2,J33/s3}; */ /* matrix - must be in column-major order for LAPACK*/
	PetscScalar J[9] = {J11/s1,J12/s1,J13/s1,J21/s2,J22/s2,J23/s2,J31/s3,J32/s3,J33/s3};  /* matrix - must be in column-major order*/
#ifdef verbose
	printf("Jnum=[\n");
	printf("%e\t%e\t%e\n",J[0],J[1],J[2]);
	printf("%e\t%e\t%e\n",J[3],J[4],J[5]);
	printf("%e\t%e\t%e]\n",J[6],J[7],J[8]);
#endif
	/* 	PetscInt piv[3]; */
	PetscScalar B[3] = {-D0.T33/D0ii/s1*Kres,-D0.T23/D0ii/s2*Kres,-D0.T13/D0ii/s3*Kres};/* negative of residual */
#ifdef verbose
	printf("RHS norm = %e\n",sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));
#endif
	/* scale right hand side*/
#ifdef verbose
	printf("RHS =[ %e %e %e ]\n",B[0],B[1],B[2]);
#endif
	/* 	PetscInt info; */
	/* 	PetscInt status; */
	/* 	PetscInt three = 3; */
	/* 	PetscInt one = 1; */
	/* code to use LAPACK to solve linear system: */
	/* SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO ) */
	//	status = dgesv_(&three,&one,&J[0],&three,&piv[0],&B[0],&three,&info);
	solve33( J,B); /* use my analytic expression "solver" */
	/* scale solution */
#ifdef verbose
	printf("unscaled B = (%e, %e, %e)\n",B[0],B[1],B[2]);
#endif
	//B[0] *= s1;
	//B[1] *= s2;
	//B[2] *= s3;
	/* 	PetscInt i; */
	/* 	for(i=1;i<=3;i++){ */
	/* 	  if(piv[i] != i){ */
	/* 	    swap rows in B */
	/* 	    PetscScalar val2 = B[piv[i]-1]; */
	/* 	    B[piv[i]-1] = B[i-1]; */
	/* 	    B[i-1] = val2; */
	/* 	  } */
	
	/* 	} */
	/* 	printf("status = %d, info = %d\n",status,info); *//* dump status information about LAPACK solve*/
	/* 	printf("pivots: %d %d %d\n",piv[0],piv[1],piv[2]); */

	const PetscInt nline = 1;
	PetscInt iline;
	PetscScalar amax=1.0;
	PetscScalar amin=0.0;
	PetscScalar resmin;
	PetscScalar mina;
	/*       PetscScalar res; */
	PetscScalar alpha;

	for(iline=0;iline<nline;iline++){
	  alpha = -(amax-amin)/nline*iline + amax;
	  Marker m1= *marker;
	  m1.s.T33 += B[0]*alpha;
	  m1.s.T23 += B[1]*alpha;
	  m1.s.T13 += B[2]*alpha;
	  formViscoplasticMMarker(&m1,options,&D1,0);/* calculate strain rate */
	  /*calculate residual*/
	  PetscScalar D1ii = sqrt( D1.T11*D1.T11 + D1.T22*D1.T22 + D1.T33*D1.T33 + 2*D1.T12*D1.T12 + 2*D1.T13*D1.T13 + 2*D1.T23*D1.T23);
	  D1ii = 1.0;
	  D1.T33 /= D1ii;
	  D1.T13 /= D1ii;
	  D1.T23 /= D1ii;
	  res=sqrt(D1.T33*D1.T33*Kres*Kres + D1.T13*D1.T13*Kres*Kres +D1.T23*D1.T23*Kres*Kres);
#ifdef verbose
	  printf("res%d = %e\n",iline,res);
#endif
	  if(iline==0){
	    resmin=res;
	    mina=alpha;
	  }
	  if(res<resmin){ resmin=res; mina=alpha;} 
	}/* end line search loop*/
	if( mina == 0.0 ){
#ifdef verbose
	  printf("Couldn't improve residual. iNewton = %d\n",iNewton);
#endif
	  break;/* stop newton iteration - we cannot improve*/	
	} else if( iNewton > 1 && (fabs(res-lastres)/lastres < newtonRelTol ) ){
#ifdef verbose
	  printf("relative tolerance met\n");
#endif
	  break;

	} else if( res < newtonAbsTol){
#ifdef verbose
	  printf("absolute tolerance met\n");
#endif
	  break;
	}

	marker->residual = res;

	marker->s.T33 += B[0]*mina;/* note that here we update the stress in the actual marker, not a copy*/
	marker->s.T23 += B[1]*mina;
	marker->s.T13 += B[2]*mina;		
#ifdef verbose
	printf("adjusting stresses: dzsz, dSyz,dSxz = (%e,%e,%e)\n",B[0]*mina,B[1]*mina,B[2]*mina);
#endif

#ifdef verbose      
	printf("D=\n");print33s(&D0);
#endif      
	lastres = res;
      } /* end newton raphson loop*/
	/* dump stress state and strain rate */
      formViscoplasticMMarker(marker,options,&D0,1);
#ifdef verbose    
      printf("final stress state:\n");
      print33s(&(marker->s));
      printf("final strain rate:\n");
      print33s( &D0 );
#endif
      /* invert marker viscoplasticity tensor*/
      invertMtoN( &(marker->texture) );
    } else {/* end if material hasTexture */
      /* stuff to do if this marker has no texture */
      
      PetscInt i,j;
      for(i=0;i<6;i++){
	for(j=0;j<6;j++){
	  marker->texture.M[i][j] = 0.0;
	}
      }
      PetscScalar thiseta;
      PetscInt thismat = (PetscInt) marker->Mat;
      if( materials->hasEtaT[thismat]){
	thiseta = materials->materialEta[thismat]*exp(-materials->QonR[thismat]*marker->T);
      }else{
	thiseta = materials->materialEta[thismat];
      }
	
      marker->texture.M[0][0] = 1.0/(2*thiseta);
      marker->texture.M[1][1] = 1.0/(2*thiseta);
      marker->texture.M[2][2] = 1.0/(2*thiseta);
      marker->texture.M[3][3] = 1.0/(4*thiseta);
      marker->texture.M[4][4] = 1.0/(4*thiseta);
      marker->texture.M[5][5] = 1.0/(4*thiseta);
      /* special case to test anisotropic viscosity */
      if( materials->hasTexture[ marker->Mat ] ==2 ){
	/* 00 relates dxx = m00 * sxx */
	marker->texture.M[0][0] *= 0.5;
      }
      invertMtoN( &(marker->texture) );
      
    }
      
  }/* end if in bounds*/
    

  
  PetscFunctionReturn(ierr);
}/* end function */


PetscErrorCode formViscoplasticM(MarkerSet *markerset, Materials *materials, Options *options, PetscInt updateSpin){
  /* loop over markers in markerset */
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  PetscInt m;
  Tensor33s D;

  for(m=0;m<markerset->nMark;m++){
    if( markerset->markers[m].cellX != -1 ){/* check to see whether marker is in bounds*/
      if( materials->hasTexture[(PetscInt) markerset->markers[m].Mat ] == 1){
	formViscoplasticMMarker( &(markerset->markers[m]), options, &D, updateSpin);
	invertMtoN( &(markerset->markers[m].texture) );
#ifdef verbose
	printf("final stress state:\n");
	print33s(&(markerset->markers[m].s));
	printf("N tensor (2d):\n");
	PetscInt i,j;
	for(i=0;i<2;i++){
	  for(j=0;j<2;j++){
	    printf("%e\t",markerset->markers[m].texture.N[i][j]);
	  }
	  printf("\n");
	}
#endif
      } else {/* end if hasTexture */
	PetscInt i,j;
	for(i=0;i<6;i++){
	  for(j=0;j<6;j++){
	    markerset->markers[m].texture.M[i][j] = 0.0;
	  }
	}
	PetscScalar thiseta;
	PetscInt thismat = (PetscInt) markerset->markers[m].Mat;
	if( materials->hasEtaT[thismat]){
	  thiseta = materials->materialEta[thismat]*exp(-materials->QonR[thismat]*markerset->markers[m].T);
	}else{
	  thiseta = materials->materialEta[thismat];
	}

	markerset->markers[m].texture.M[0][0] = 1/(2*thiseta);
	markerset->markers[m].texture.M[1][1] = 1/(2*thiseta);
	markerset->markers[m].texture.M[2][2] = 1/(2*thiseta);
	markerset->markers[m].texture.M[3][3] = 1/(4*thiseta);
	markerset->markers[m].texture.M[4][4] = 1/(4*thiseta);
	markerset->markers[m].texture.M[5][5] = 1/(4*thiseta);

	if( materials->hasTexture[ markerset->markers[m].Mat ]==2 ){
	  /* 00 relates dxx = m00 * sxx */
	  markerset->markers[m].texture.M[0][0] *= 0.5;
	}	
	invertMtoN( &(markerset->markers[m].texture) );
	/* material does not have texture - just assume that the material is isoviscous */
	//	markerset->markers[m].texture.N[0][0] = 2.0*materials->materialEta[(PetscInt) markerset->markers[m].Mat ];
	//	markerset->markers[m].texture.N[0][1] = 0.0;
	//	markerset->markers[m].texture.N[1][0] = 0.0;
	//	markerset->markers[m].texture.N[1][1] = 2.0*materials->materialEta[(PetscInt) markerset->markers[m].Mat ];
      }
    }/* end if in bounds */
  }/* end loop over markers*/

  PetscFunctionReturn(ierr);
}
				 



PetscScalar getFabricDT(MarkerSet *markerset){
  const PetscScalar pi = M_PI;
  const PetscScalar ang_max = ANG_MAX;
    /* compute maximum timestep to allow rotation through ang_max */
  PetscScalar dtmax = 1e99;
  PetscInt i,j,k,m;
  PetscScalar thetadotmax=0.0; 
  PetscScalar phidotmax=0.0;
  for(m=0;m<markerset->nMark;m++){
    Marker *M = &(markerset->markers[m]);
    if(M->cellX != -1){
      for(i=0;i<NT;i++){
	for(j=0;j<NT;j++){
	  for(k=0;k<NT;k++){
	    if( fabs((PetscScalar) M->texture.cthetadot[i][j][k]) > thetadotmax) thetadotmax =(PetscScalar) fabs((PetscScalar)M->texture.cthetadot[i][j][k]);
	    if( fabs((PetscScalar) M->texture.cphidot[i][j][k]) > phidotmax) phidotmax = (PetscScalar)fabs((PetscScalar)M->texture.cphidot[i][j][k]);
	  }
	}
      }
    }
  }
  {
    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    printf("[%d] max thetadot, phidot = %e,%e\n",rank,thetadotmax,phidotmax);
  }
  /* do an MPI_Allreduce on thetadotmax */
  {
    PetscScalar tdm = thetadotmax;
    PetscScalar pdm = phidotmax;
    PetscErrorCode ierr;
    ierr=MPI_Allreduce(&tdm,&thetadotmax,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr=MPI_Allreduce(&pdm,&phidotmax,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); CHKERRQ(ierr);
  }


  if(thetadotmax > phidotmax){
    return ang_max/thetadotmax;
  }  else{
    return ang_max/phidotmax;
  }
}

PetscErrorCode initializeIsotropicFabric(MarkerSet *markerset){
  PetscFunctionBegin;
  const PetscScalar pi = 3.1415926535;  
  PetscInt m,i,j,k;
  Marker *markers = markerset->markers;
  for(m=0;m<markerset->nMark;m++){
    for(i=0;i<NT;i++){
      for(j=0;j<NT;j++){
	for(k=0;k<NT;k++){
	  markers[m].texture.ctheta[i][j][k] = 2.0*pi*drand48();
	  PetscScalar z = -1.0 + 2.0*drand48();
	  /* cos phi = z/1.0 */
	  markers[m].texture.cphi[i][j][k] = acos(z);
	}
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode updateFabric(MarkerSet *markerset, Materials *materials, PetscScalar dt){
  PetscFunctionBegin;
  const PetscScalar pi = M_PI;  
  const PetscScalar ang_max = ANG_MAX; /* maximum angle through which c-axis may rotate */
  PetscInt m,i,j,k;
  for(m=0;m<markerset->nMark;m++){
    Marker *M = &(markerset->markers[m]);
    if( materials->hasTexture[(PetscInt) markerset->markers[m].Mat ] ){
      for(i=0;i<NT;i++){
	for(j=0;j<NT;j++){
	  for(k=0;k<NT;k++){
	    if( isbad(M->texture.ctheta[i][j][k]) || isbad(M->texture.cphi[i][j][k]) || isbad(M->texture.cthetadot[i][j][k]) || isbad(M->texture.cphidot[i][j][k] ) ){
	      printf("NAN found before updating texture %e, %e, rates = %e,%e \n",M->texture.ctheta[i][j][k],M->texture.cphi[i][j][k],M->texture.cthetadot[i][j][k],M->texture.cphidot[i][j][k]); fflush(stdout);
	    }
	    PetscScalar dtheta =  M->texture.cthetadot[i][j][k]*dt;
	    PetscScalar dphi =  M->texture.cphidot[i][j][k]*dt;
	    if( (fabs(dtheta) > ang_max) || (fabs(dphi) > ang_max)){
	      if( fabs(dtheta) > fabs(dphi) ){
		dphi *= fabs(ang_max/dtheta);
		dtheta *= fabs(ang_max/dtheta);
	      }else{/* dphi > dtheta */
		dtheta *= fabs(ang_max/dphi);
		dphi *= fabs(ang_max/dphi);
	      }
	    }
	    M->texture.ctheta[i][j][k] += dtheta;
	    M->texture.cphi[i][j][k] += dphi;
	  }
	}
      }
      /* check to make sure that theta and phi stay in their appropriate ranges */
      for(i=0;i<NT;i++){
	for(j=0;j<NT;j++){
	  for(k=0;k<NT;k++){
	    while( M->texture.ctheta[i][j][k] < 0.0 ){
	      M->texture.ctheta[i][j][k] += 2.0*pi;
	    }
	    
	    if(M->texture.ctheta[i][j][k] > 2*pi ) M->texture.ctheta[i][j][k] = fmod(M->texture.ctheta[i][j][k],2*pi);
	  
	    if( M->texture.cphi[i][j][k] > pi){ 
	      M->texture.cphi[i][j][k] = 2*pi - M->texture.cphi[i][j][k];
	      M->texture.ctheta[i][j][k] += pi;
	      M->texture.ctheta[i][j][k] = fmod(M->texture.ctheta[i][j][k],2*pi);
	    }
	    if(  M->texture.cphi[i][j][k] < 0.0){
	      M->texture.cphi[i][j][k] = -M->texture.cphi[i][j][k];
	      M->texture.ctheta[i][j][k] += pi;
	      M->texture.ctheta[i][j][k] = fmod(M->texture.ctheta[i][j][k],2*pi);	
	    }
	    if( M->texture.cphi[i][j][k] < 1e-6) M->texture.cphi[i][j][k] = 1e-6;/* keep this a small positive number to prevent divides by zero */
	    if( M_PI-M->texture.cphi[i][j][k] < 1e-6) M->texture.cphi[i][j][k] = M_PI-1e-6;

	    if( isbad(M->texture.ctheta[i][j][k]) || isbad(M->texture.cphi[i][j][k])){
	      printf("NAN found while updating texture %e, %e, rates = %e,%e \n",M->texture.ctheta[i][j][k],M->texture.cphi[i][j][k],M->texture.cthetadot[i][j][k],M->texture.cphidot[i][j][k]); fflush(stdout);
	    }
	  }
	}
      }
    } /* end if marker has texture */


  }
  PetscFunctionReturn(0);
}

PetscErrorCode formViscoplasticMMarker(Marker *marker, Options *options, Tensor33s *D, PetscInt updateSpin){
  /* last argument is whether or not to update the spin for each crystal - this should only be done once the solution has converged */
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  
  const PetscScalar z6xi = 1/(zeta+6*xi);
  const PetscScalar NT3 = (PetscScalar) nt*nt*nt;
  
  TextureInfo *tex = &(marker->texture);
  
  /* Define macroscopic stress tensor components*/
  PetscScalar S11 = marker->s.T11;
  PetscScalar S22 = marker->s.T22;
  PetscScalar S33 = marker->s.T33;
  PetscScalar S23 = marker->s.T23;
  PetscScalar S13 = marker->s.T13;
  PetscScalar S12 = marker->s.T12;
  
  PetscScalar trS = (S11+S22+S33)/3;
  S11 -= trS;
  S22 -= trS;
  S33 -= trS;
  
  /* BEGIN NEWTON-RAPHSON ITERATION OUTERMOST LOOP */
  PetscInt nNewton = 1;
  PetscInt iNewton;
  for(iNewton=0;iNewton<nNewton;iNewton++){
    /*       printf("iNewton = %d\n",iNewton); */
    /* compute effective stress*/
    PetscScalar sii = sqrt( S11*S11 +S22*S22 + S33*S33 + 2*(S12*S12+S23*S23+S13*S13) );    
    /* Define Temperature and Rheological Parameters*/
    PetscScalar T = marker->T;
    PetscScalar d = options->grainSize;     /* define grain size */
    /* loop over slip systems and compute temperature-corrected tau0 */
    {
      PetscInt s;
      for(s=0;s<NS;s++){
	tau0s[s] = pow( gd0s[s]/A0[s]/exp(-QonRs[s]/T) ,1.0/n1s[s]);
      }
    }
      
    /* form script-tau for each crystal in aggregate */
    PetscInt i,j,k;
    for(i=0;i<nt;i++){
      for(j=0;j<nt;j++){
	for(k=0;k<nt;k++){
	    
	  /* sum over slip systems*/
	  PetscScalar CT = (PetscScalar) tex->ctheta[i][j][k];
	  PetscScalar CP = (PetscScalar) tex->cphi[i][j][k];
	  /* temporary variables created by mathematica */
	  PetscScalar v00 = sin(CT);
	  PetscScalar v01 = cos(CT);
	  PetscScalar v02 = cos(CP);
	  PetscScalar v06 = pow(v00,2);
	  PetscScalar v07 = pow(v01,2);
	  PetscScalar v08 = -1 + v02;
	  PetscScalar v09 = pow(sin(CP/2.),2);
	  PetscScalar v10 = sqrt(pow(sin(CP),2));
	    
	  PetscInt s;
	  PetscScalar sT1 = 0.0;
	  PetscScalar sT2 = 0.0;
	  PetscScalar sT3 = 0.0;
	    
	  //	  sT[i][j][k] = 0.0;
	  for(s=0;s<NS;s++){
	    PetscScalar b1 = bs[s][0];
	    PetscScalar b2 = bs[s][1];
	    PetscScalar b3 = bs[s][2];
	    PetscScalar r11 = rs[s][0][0];
	    PetscScalar r12 = rs[s][0][1];
	    PetscScalar r13 = rs[s][0][2];
	    PetscScalar r22 = rs[s][1][1];
	    PetscScalar r23 = rs[s][1][2];
	    /* PetscScalar r33 = rs[s][2][2]; */
	      
	    /* mathematica temporary variables*/
	    PetscScalar v03 = r22*S11;
	    PetscScalar v04 = r12*S12;
	    PetscScalar v05 = r11*S11;
	    PetscScalar v11 = -v03 + 2*v04 + v05 + 2*(v03 - v04)*v06*v09 - 2*(v04 \
									      + v05)*v07*v09 - (r23*S11 - r13*S12)*v00*v10 + v01*(-2*(r11 + \
																      r22)*S12*v00*v09 + (r13*S11 + r23*S12)*v10);
	    sT1 += (b1 + b1*v07*v08 + v01*(b2*v00*v08 + b3*v10))*v11;
	    sT2 += (b2 + b2*v06*v08 - 2*b1*v00*v01*v09 +	\
		    b3*v00*v10)*v11;
	    sT3 += (b3*v02 - (b2*v00 + b1*v01)*v10)*v11;

	    if( isbad(sT1) ||  isbad(sT2) ||  isbad(sT3)){
	      printf("NAN found computing sT. %d,%d,%d, CT=%e,CP=%e, sT123=%e,%e,%e\n",i,j,k,CT,CP,sT1,sT2,sT3); fflush(stdout);
	    }

	      
	  }/* end summation over slip systems*/
	  /* take 2-norm of ST1,ST2,ST3, assign to sT*/ 
	  sT[i][j][k] = sqrt(sT1*sT1+sT2*sT2+sT3*sT3);
	  /* 	  printf("sT[%d][%d][%d][%d] = %e\n",i,j,k,s,sT[i][j][k]);  */
	    
	}
      }
    }
    for(i=0;i<nt;i++){
      for(j=0;j<nt;j++){
	for(k=0;k<nt;k++){
	  PetscScalar iterm = sT[idxmod(i+1)][j][k] + sT[idxmod(i-1)][j][k] + sT[i][idxmod(j+1)][k] + sT[i][idxmod(j-1)][k] + sT[i][j][idxmod(k+1)] + sT[i][j][idxmod(k-1)];
	  /* force softnessE = 0 if sT is really small so that softnessE does not blow up due to divide by zero or small number*/
	    
	  softnessE[i][j][k] = sT[i][j][k] > 1e-8 ? z6xi*(zeta+xi*fabs(iterm/sT[i][j][k])) : 0.0;
	  /*  	  printf("softnessE[%d,%d,%d],%e\t iterm=%e, st[i,j,k]=%e\n",i,j,k,softnessE[i][j][k],iterm,sT[i][j][k]);  */
	  if(isbad(iterm) || isbad(softnessE[i][j][k]) || isbad(sT[i][j][k])) {
	    //printf("NAN DEBUG:%d,%d,%d %e %e %e %e %e %e\n",i,j,k,sT[idxmod(i+1)][j][k] , sT[idxmod(i-1)][j][k] , sT[i][idxmod(j+1)][k] ,sT[i][idxmod(j-1)][k],sT[i][j][idxmod(k+1)], sT[i][j][idxmod(k-1)]);
	    printf("NAN DEBUG:%d,%d,%d,iterm=%e sT=%e, softnessE=%e\n",i,j,k,iterm,softnessE[i][j][k],sT[i][j][k]);
	    
	    printf("Stress tensor:\n");
	    printf("%e %e %e %e %e %e\n",S11,S22,S33,S13,S23,S12); fflush(stdout);

	    SETERRQ(PETSC_COMM_SELF,ERR_MICROMACRO_NAN,"NAN found in microMacro step\n"); fflush(stdout);
	    
	  }
	}
      }
    }
      
    PetscScalar (*Ms)[6] = marker->texture.M;
    PetscScalar Mwc[6][6], Mwtot[6][6];
    for(i=0;i<6;i++){
      for(j=0;j<6;j++){
	Ms[i][j] = 0.0;
	Mwc[i][j] = 0.0;
	Mwtot[i][j] = 0.0;
      }
    }
    /* calculate M tensor by summation */        
    /* loop over i,j,k */
    for(i=0;i<NT;i++){
      for(j=0;j<NT;j++){
	for(k=0;k<NT;k++){
	  PetscScalar CT = tex->ctheta[i][j][k];
	  PetscScalar CP = tex->cphi[i][j][k];
	  PetscInt s;
	  /* loop over slip systems */
	  for(s=0;s<NS;s++){
	    /* BEGIN NEW CODE*/
	    PetscScalar v000 = sin(CT);
	    PetscScalar v001 = cos(CT);
	    PetscScalar v002 = cos(CP);
	    PetscScalar v003 = 2*CT;
	    PetscScalar v004 = -rs[s][-1 + 2][-1 + 2];
	    PetscScalar v005 = -S33;
	    PetscScalar v006 = rs[s][-1 + 3][-1 + 3]*S33;
	    PetscScalar v007 = rs[s][-1 + 2][-1 + 2]*S22;
	    PetscScalar v008 = rs[s][-1 + 1][-1 + 3]*S13;
	    PetscScalar v009 = rs[s][-1 + 1][-1 + 2]*S12;
	    PetscScalar v010 = rs[s][-1 + 2][-1 + 3]*S23;
	    PetscScalar v011 = rs[s][-1 + 1][-1 + 1]*S11;
	    PetscScalar v012 = 1/tau0s[s];
	    PetscScalar v013 = -rs[s][-1 + 2][-1 + 3];
	    PetscScalar v014 = sin(v003);
	    PetscScalar v015 = rs[s][-1 + 2][-1 + 2]*v000;
	    PetscScalar v016 = rs[s][-1 + 1][-1 + 2]*v001;
	    PetscScalar v017 = cos(v003);
	    PetscScalar v018 = -1 + v002;
	    PetscScalar v019 = pow(v001,2);
	    PetscScalar v020 = pow(v000,2);
	    PetscScalar v021 = 2*v010;
	    PetscScalar v022 = 2*v008;
	    PetscScalar v023 = rs[s][-1 + 2][-1 + 3]*v000;
	    PetscScalar v024 = rs[s][-1 + 1][-1 + 3]*v001;
	    PetscScalar v025 = rs[s][-1 + 1][-1 + 2]*v000;
	    PetscScalar v026 = rs[s][-1 + 1][-1 + 1]*v001;
	    PetscScalar v027 = rs[s][-1 + 2][-1 + 3]*v001;
	    PetscScalar v028 = rs[s][-1 + 2][-1 + 3]*v017;
	    PetscScalar v029 = rs[s][-1 + 2][-1 + 2]*v017;
	    PetscScalar v030 = (rs[s][-1 + 1][-1 + 1] + rs[s][-1 + 2][-1 + \
								      2])*v000;
	    PetscScalar v031 = rs[s][-1 + 2][-1 + 3]*v014;
	    PetscScalar v032 = rs[s][-1 + 1][-1 + 3]*v017;
	    PetscScalar v033 = rs[s][-1 + 1][-1 + 3]*v014;
	    PetscScalar v034 = v013*v014;
	    PetscScalar v035 = -v032;
	    PetscScalar v036 = -2*v000*v024;
	    PetscScalar v037 = -(rs[s][-1 + 1][-1 + 2]*v014);
	    PetscScalar v038 = -v033;
	    PetscScalar v039 = -2*v000*v016;
	    PetscScalar v040 = v013*v017;
	    PetscScalar v041 = sqrt(pow(sin(CP),2));
	    PetscScalar v042 = pow(sin(CP/2.),2);
	    PetscScalar v043 = v015 + v016;
	    PetscScalar v044 = v023 + v024;
	    PetscScalar v045 = rs[s][-1 + 2][-1 + 3]*v041;
	    PetscScalar v046 = rs[s][-1 + 1][-1 + 3]*v041;
	    PetscScalar v047 = rs[s][-1 + 1][-1 + 2]*v018*v019;
	    PetscScalar v048 = rs[s][-1 + 3][-1 + 3]*v041;
	    PetscScalar v049 = v000*v048;
	    PetscScalar v050 = v023*v041;
	    PetscScalar v051 = 2*v001*v048;
	    PetscScalar v052 = 2*v049;
	    PetscScalar v053 = 2*v050;
	    PetscScalar v054 = 2*v015*v041;
	    PetscScalar v055 = 2*v016*v041;
	    PetscScalar v056 = v018*v025 + v046;
	    PetscScalar v057 = rs[s][-1 + 3][-1 + 3]*v002 - v041*v044;
	    PetscScalar v058 = rs[s][-1 + 1][-1 + 2] + v001*(v015*v018 + v045) + \
	      v047;
	    PetscScalar v059 = rs[s][-1 + 1][-1 + 1] + rs[s][-1 + 1][-1 + \
								     1]*v018*v019 + v001*v056;
	    PetscScalar v060 = rs[s][-1 + 1][-1 + 3] + rs[s][-1 + 1][-1 + \
								     3]*v018*v019 + v001*(v018*v023 + v048);
	    PetscScalar v061 = rs[s][-1 + 2][-1 + 3] + rs[s][-1 + 2][-1 + \
								     3]*v018*v020 + v036*v042 + v049;
	    PetscScalar v062 = rs[s][-1 + 2][-1 + 2] + rs[s][-1 + 2][-1 + \
								     2]*v018*v020 + v039*v042 + v050;
	    PetscScalar v063 = rs[s][-1 + 2][-1 + 3] + v028 + v036 + \
	      2*v000*v002*v044 + v052;
	    PetscScalar v064 = rs[s][-1 + 2][-1 + 2] + v029 + v039 + \
	      2*v000*v002*v043 + v053;
	    PetscScalar v065 = rs[s][-1 + 1][-1 + 2] + v004*v014 - rs[s][-1 + \
									 1][-1 + 2]*v017 + 2*v027*v041 + 2*v001*v002*v043;
	    PetscScalar v066 = rs[s][-1 + 1][-1 + 3] + v034 + v035 + \
	      2*v001*v002*v044 + v051;
	    PetscScalar v067 = rs[s][-1 + 1][-1 + 1] - rs[s][-1 + 1][-1 + 1]*v017 \
	      + 2*v001*v002*(v025 + v026) + v037 + 2*v024*v041;
	    PetscScalar v068 = rs[s][-1 + 2][-1 + 2] + rs[s][-1 + 1][-1 + \
								     2]*v002*v014 + 2*rs[s][-1 + 2][-1 + 2]*v002*v020 + v029 + v037 + \
	      v053;
	    PetscScalar v069 = rs[s][-1 + 2][-1 + 3] + 2*rs[s][-1 + 2][-1 + \
								       3]*v002*v020 + v028 + v002*v033 + v038 + v052;
	    PetscScalar v070 = rs[s][-1 + 1][-1 + 2] + v002*(rs[s][-1 + 1][-1 + \
									   2] + v001*v030) - v001*(v030 + v013*v041) + v000*v046;
	    PetscScalar v071 = rs[s][-1 + 1][-1 + 3] - 2*rs[s][-1 + 1][-1 + \
								       3]*v002*v020 + v002*v031 + v034 + v035 + 2*(rs[s][-1 + 1][-1 + 1] + \
														   rs[s][-1 + 3][-1 + 3])*v001*v041 + 2*v025*v041;
	    PetscScalar v072 = v001*(2*v000*(rs[s][-1 + 1][-1 + 1] + v004)*v042 + \
				     v045) + v047 - v000*v056;
	    PetscScalar v073 = rs[s][-1 + 2][-1 + 3] - 2*v001*v002*(-(rs[s][-1 + \
									    1][-1 + 3]*v000) + v027) + v028 + v038 + v052 + v054 + v055;
	    PetscScalar v074 = rs[s][-1 + 1][-1 + 3] + v002*(3*rs[s][-1 + 1][-1 + \
									     3] + v031 + v032) + v034 + v035 - 2*v025*v041 - 2*v026*v041 + v051;
	    PetscScalar v075 = rs[s][-1 + 2][-1 + 3] + v028 + v038 + \
	      v002*(3*rs[s][-1 + 2][-1 + 3] + v033 + v040) - 2*v015*v041 - \
	      2*v016*v041 + v052;
	    PetscScalar v076 = v013 + v033 + v002*(-3*rs[s][-1 + 2][-1 + 3] + \
						   v028 + v038) + v040 - 2*v049 + v054 + v055;
	    PetscScalar v077 = pow(fabs(v012*(v006 + v007 + 2*v009 + v011 + v021 \
					      + v022 + v000*(rs[s][-1 + 1][-1 + 3]*S12 - rs[s][-1 + 1][-1 + 2]*S13 \
							     + rs[s][-1 + 2][-1 + 3]*S22 + rs[s][-1 + 3][-1 + 3]*S23 + S23*v004 + \
							     rs[s][-1 + 2][-1 + 3]*v005)*v041 - 2*v020*(v006 + v007 + v008 + v009 \
													+ v021)*v042 - 2*v019*(v006 + v009 + v010 + v011 + v022)*v042 + \
					      v001*((rs[s][-1 + 2][-1 + 3]*S12 - rs[s][-1 + 1][-1 + 1]*S13 + \
						     rs[s][-1 + 3][-1 + 3]*S13 - rs[s][-1 + 1][-1 + 2]*S23 + rs[s][-1 + \
														   1][-1 + 3]*(S11 + v005))*v041 - 2*(rs[s][-1 + 1][-1 + 1]*S12 + \
																		      rs[s][-1 + 2][-1 + 2]*S12 + rs[s][-1 + 2][-1 + 3]*S13 + rs[s][-1 + \
																										    1][-1 + 2]*(S11 + S22) + rs[s][-1 + 1][-1 + 3]*S23)*v000*v042))),-1 + \
				   n1s[s]);
	    PetscScalar v078 = gd0s[s]*softnessE[i][j][k]*v012*v057*v059*v077;
	    PetscScalar v079 = gd0s[s]*softnessE[i][j][k]*v012*v057*v062*v077;
	    PetscScalar v080 = gd0s[s]*softnessE[i][j][k]*v012*v059*v062*v077;
	    Ms[0][0] += gd0s[s]*softnessE[i][j][k]*v012*pow(v059,2)*v077;
	    Ms[0][1] += v080;
	    Ms[0][2] += v078;
	    Ms[0][3] += gd0s[s]*softnessE[i][j][k]*v012*v059*v061*v077;
	    Ms[0][4] += gd0s[s]*softnessE[i][j][k]*v012*v059*v060*v077;
	    Ms[0][5] += gd0s[s]*softnessE[i][j][k]*v012*v058*v059*v077;
	    Ms[1][0] += v080;
	    Ms[1][1] += gd0s[s]*softnessE[i][j][k]*v012*pow(v062,2)*v077;
	    Ms[1][2] += v079;
	    Ms[1][3] += gd0s[s]*softnessE[i][j][k]*v012*v061*v062*v077;
	    Ms[1][4] += gd0s[s]*softnessE[i][j][k]*v012*v060*v062*v077;
	    Ms[1][5] += gd0s[s]*softnessE[i][j][k]*v012*v058*v062*v077;
	    Ms[2][0] += v078;
	    Ms[2][1] += v079;
	    Ms[2][2] += gd0s[s]*softnessE[i][j][k]*v012*pow(v057,2)*v077;
	    Ms[2][3] += gd0s[s]*softnessE[i][j][k]*v012*v057*v061*v077;
	    Ms[2][4] += (gd0s[s]*softnessE[i][j][k]*v012*v057*v066*v077)/2.;
	    Ms[2][5] += gd0s[s]*softnessE[i][j][k]*v012*v057*v058*v077;
	    Ms[3][0] += (gd0s[s]*softnessE[i][j][k]*v012*v067*v075*v077)/8.;
	    Ms[3][1] += -(gd0s[s]*softnessE[i][j][k]*v012*v068*v076*v077)/8.;
	    Ms[3][2] += (gd0s[s]*softnessE[i][j][k]*v012*v057*v075*v077)/4.;
	    Ms[3][3] += -(gd0s[s]*softnessE[i][j][k]*v012*v069*v076*v077)/8.;
	    Ms[3][4] += (gd0s[s]*softnessE[i][j][k]*v012*v066*v075*v077)/8.;
	    Ms[3][5] += -(gd0s[s]*softnessE[i][j][k]*v012*v065*v076*v077)/8.;
	    Ms[4][0] += (gd0s[s]*softnessE[i][j][k]*v012*v067*v074*v077)/8.;
	    Ms[4][1] += (gd0s[s]*softnessE[i][j][k]*v012*v064*v074*v077)/8.;
	    Ms[4][2] += (gd0s[s]*softnessE[i][j][k]*v012*v057*v074*v077)/4.;
	    Ms[4][3] += (gd0s[s]*softnessE[i][j][k]*v012*v069*v074*v077)/8.;
	    Ms[4][4] += (gd0s[s]*softnessE[i][j][k]*v012*v066*v074*v077)/8.;
	    Ms[4][5] += (gd0s[s]*softnessE[i][j][k]*v012*v065*v074*v077)/8.;
	    Ms[5][0] += (gd0s[s]*softnessE[i][j][k]*v012*v067*v070*v077)/4.;
	    Ms[5][1] += (gd0s[s]*softnessE[i][j][k]*v012*v068*v070*v077)/4.;
	    Ms[5][2] += (gd0s[s]*softnessE[i][j][k]*v012*v057*v070*v077)/2.;
	    Ms[5][3] += (gd0s[s]*softnessE[i][j][k]*v012*v063*v070*v077)/4.;
	    Ms[5][4] += (gd0s[s]*softnessE[i][j][k]*v012*v066*v070*v077)/4.;
	    Ms[5][5] += (gd0s[s]*softnessE[i][j][k]*v012*v065*v070*v077)/4.;
	    if(updateSpin){
	      Mwc[0][0] = 0.0;
	      Mwc[0][1] = 0.0;
	      Mwc[0][2] = 0.0;
	      Mwc[0][3] = 0.0;
	      Mwc[0][4] = 0.0;
	      Mwc[0][5] = 0.0;
	      Mwc[1][0] = 0.0;
	      Mwc[1][1] = 0.0;
	      Mwc[1][2] = 0.0;
	      Mwc[1][3] = 0.0;
	      Mwc[1][4] = 0.0;
	      Mwc[1][5] = 0.0;
	      Mwc[2][0] = 0.0;
	      Mwc[2][1] = 0.0;
	      Mwc[2][2] = 0.0;
	      Mwc[2][3] = 0.0;
	      Mwc[2][4] = 0.0;
	      Mwc[2][5] = 0.0;
	      Mwc[3][0] = (gd0s[s]*softnessE[i][j][k]*v012*v067*v073*v077)/8.;
	      Mwc[3][1] = (gd0s[s]*softnessE[i][j][k]*v012*v064*v073*v077)/8.;
	      Mwc[3][2] = (gd0s[s]*softnessE[i][j][k]*v012*v057*v073*v077)/4.;
	      Mwc[3][3] = (gd0s[s]*softnessE[i][j][k]*v012*v069*v073*v077)/8.;
	      Mwc[3][4] = (gd0s[s]*softnessE[i][j][k]*v012*v066*v073*v077)/8.;
	      Mwc[3][5] = (gd0s[s]*softnessE[i][j][k]*v012*v065*v073*v077)/8.;
	      Mwc[4][0] = (gd0s[s]*softnessE[i][j][k]*v012*v067*v071*v077)/8.;
	      Mwc[4][1] = (gd0s[s]*softnessE[i][j][k]*v012*v068*v071*v077)/8.;
	      Mwc[4][2] = (gd0s[s]*softnessE[i][j][k]*v012*v057*v071*v077)/4.;
	      Mwc[4][3] = (gd0s[s]*softnessE[i][j][k]*v012*v069*v071*v077)/8.;
	      Mwc[4][4] = (gd0s[s]*softnessE[i][j][k]*v012*v066*v071*v077)/8.;
	      Mwc[4][5] = (gd0s[s]*softnessE[i][j][k]*v012*v065*v071*v077)/8.;
	      Mwc[5][0] = (gd0s[s]*softnessE[i][j][k]*v012*v067*v072*v077)/4.;
	      Mwc[5][1] = (gd0s[s]*softnessE[i][j][k]*v012*v068*v072*v077)/4.;
	      Mwc[5][2] = (gd0s[s]*softnessE[i][j][k]*v012*v057*v072*v077)/2.;
	      Mwc[5][3] = (gd0s[s]*softnessE[i][j][k]*v012*v063*v072*v077)/4.;
	      Mwc[5][4] = (gd0s[s]*softnessE[i][j][k]*v012*v066*v072*v077)/4.;
	      Mwc[5][5] = (gd0s[s]*softnessE[i][j][k]*v012*v065*v072*v077)/4.;
	      PetscInt l,m;
	      for(l=0;l<6;l++){
		for(m=0;m<6;m++){ /* update sum of rotation tensor */
		  Mwtot[l][m] += Mwc[l][m];
		}
	      }
	    }
	    /*END MACHINE GENERATED CODE*/
	    
	    if(updateSpin){
	      /* compute nx, ny, nz */
	      PetscScalar nx = cos(CT) * sin(CP);
	      PetscScalar ny = sin(CT) * sin(CP);
	      PetscScalar nz = cos(CP);
	      /* compute ndot = Mw * n */
	      /* 	      PetscScalar Wxx = Mwc[1-1][1-1]*S11 + Mwc[1-1][2-1]*S22 + Mwc[1-1][3-1]*S33 + 2.0*(Mwc[1-1][4-1]*S23 + Mwc[1-1][5-1]*S13 + Mwc[1-1][6-1]*S12); */
	      /* 	      PetscScalar Wyy = Mwc[2-1][1-1]*S11 + Mwc[2-1][2-1]*S22 + Mwc[2-1][3-1]*S33 + 2.0*(Mwc[2-1][4-1]*S23 + Mwc[2-1][5-1]*S13 + Mwc[2-1][6-1]*S12); */
	      /* 	      PetscScalar Wzz = Mwc[3-1][1-1]*S11 + Mwc[3-1][2-1]*S22 + Mwc[3-1][3-1]*S33 + 2.0*(Mwc[3-1][4-1]*S23 + Mwc[3-1][5-1]*S13 + Mwc[3-1][6-1]*S12); */
	      const PetscScalar Wxx = 0.0;
	      const PetscScalar Wyy = 0.0;
	      const PetscScalar Wzz = 0.0;
	      PetscScalar Wxz = Mwc[4-1][1-1]*S11 + Mwc[4-1][2-1]*S22 + Mwc[4-1][3-1]*S33 + 2.0*(Mwc[4-1][4-1]*S23 + Mwc[4-1][5-1]*S13 + Mwc[4-1][6-1]*S12);
	      PetscScalar Wyz = Mwc[5-1][1-1]*S11 + Mwc[5-1][2-1]*S22 + Mwc[5-1][3-1]*S33 + 2.0*(Mwc[5-1][4-1]*S23 + Mwc[5-1][5-1]*S13 + Mwc[5-1][6-1]*S12);
	      PetscScalar Wxy = Mwc[6-1][1-1]*S11 + Mwc[6-1][2-1]*S22 + Mwc[6-1][3-1]*S33 + 2.0*(Mwc[6-1][4-1]*S23 + Mwc[6-1][5-1]*S13 + Mwc[6-1][6-1]*S12);
#ifdef verbose
	      printf("Wxx Wyy Wxx Wxz Wyz Wxz\n%e %e %e %e %e %e\n",Wxx,Wyy,Wzz,Wxz,Wyz,Wxy);
#endif
	      /* Wxy is skew symmetric - Wij = 0.5*(dv_i/dx_j - dv_j/dx_i) */
	      PetscScalar nxdot =  Wxx*nx + Wxy*ny + Wxz*nz;
	      PetscScalar nydot = -Wxy*nx + Wyy*ny + Wyz*nz;
	      PetscScalar nzdot = -Wxz*nx - Wyz*ny + Wzz*nz;

	      /* compute thetadot, phidot given ndot */
	      PetscScalar thetadot = -(ny/(nx*nx+ny*ny))*nxdot + nx/(nx*nx+ny*ny)*nydot;
	      if( isbad(thetadot) ){
		printf("thetadot is infinite, nx = %e, ny=%e\n",nx,ny);
	      }
	      PetscScalar phidot = -1/sqrt(1.0 - nz*nz)*nzdot;
	      tex->cthetadot[i][j][k] = thetadot;
	      tex->cphidot[i][j][k] = phidot;

	      if(isbad(tex->cphidot[i][j][k])){
		printf("NAN found line 964 phidot=%e, nz = %e, nzdot=%e, CT=%e, CP=%e\n",phidot,nz,nzdot,CT,CP);
		tex->cphidot[i][j][k] = 0.0;
	      }
	    }/* end update spin */

	  }/* end loop over s*/
	}/* end loop over k*/
      }/* end loop over j*/
    }/* end loop over i*/

    for(i=0;i<6;i++){
      for(j=0;j<6;j++){
	Ms[i][j] /= NT3;
      }
    }
    if(updateSpin){
      for(i=0;i<6;i++){
	for(j=0;j<6;j++){
	  Mwtot[i][j] /= NT3;
	}
      }
    }

    /* Compute grain-boundary sliding strain-rate */
    
    {
      const PetscScalar Agbs = 6.2e-14; /*m^1.4 Pa^-1.8 s^-1 */
      const PetscScalar pgbs = 1.4; /* grain size exponent */
      const PetscScalar ngbs = 1.8; /* creep power law exponent */
      const PetscScalar QonRgbs = -5893.32; /* -49 kJ/mol / R */
      PetscScalar Mgbs = Agbs/pow(d,pgbs)*pow(sii,ngbs-1)*exp(QonRgbs/T);
#ifdef verbose
      printf("Mgbs=%e\n",Mgbs);
#endif
      /* compute inhibition factor alpha */
      /* micro-macro strain-rate - ASSUMES PLANE STRAIN CONDITION MET */
      //	PetscScalar d11mm = M[1][1] * S11 - M[1][2]*S11 + M[1][3]*S12;
      PetscScalar d11mm = Ms[1-1][1-1]*S11 + Ms[1-1][2-1]*S22 + Ms[1-1][3-1]*S33 + 2.0*(Ms[1-1][4-1]*S23 + Ms[1-1][5-1]*S13 + Ms[1-1][6-1]*S12);
      PetscScalar d12mm = Ms[6-1][1-1]*S11 + Ms[6-1][2-1]*S22 + Ms[6-1][3-1]*S33 + 2.0*(Ms[6-1][4-1]*S23 + Ms[6-1][5-1]*S13 + Ms[6-1][6-1]*S12);
      //	PetscScalar d12mm = M[1][3] * S11 - M[2][3]*S11 + M[3][3]*S12;
      PetscScalar dmm = sqrt( d11mm*d11mm + d12mm*d12mm );
      /* gbs strain rate */
      PetscScalar d11gbs = Mgbs*S11;
      PetscScalar d12gbs = Mgbs*S12;
      PetscScalar dgbs = sqrt( d11gbs*d11gbs + d12gbs*d12gbs );
      /* prevent expression from blowing up when dmm ~= 0 (low stress) */
      PetscScalar alpha = dgbs/dmm;
      if( dmm < 1e-99 ) alpha = 0.0;
#ifdef verbose
      printf("gbs alpha = %e\n",alpha);
#endif
      /* scale M tensor */
      alpha /= (1+alpha);
      for(i=0;i<6;i++){
	for(j=0;j<6;j++){
	  Ms[i][j] *= alpha;
	}
      }
      if(updateSpin){
	/* loop over crystals */
	for(i=0;i<NT;i++){
	  for(j=0;j<NT;j++){
	    for(k=0;k<NT;k++){
	      tex->cthetadot[i][j][k] *= alpha;
	      if( isbad(tex->cthetadot[i][j][k]) ){
		printf("thetadot is infinite after mult by alpha, alpha = %e\n",alpha);
	      }
	      tex->cphidot[i][j][k] *= alpha;
	      if( isbad(tex->cphidot[i][j][k]) ){
		printf("NAN DEBUG: phidot is nan after alpha step: alpha=%e\n",alpha);fflush(stdout);
	      }	      
	    }
	  }
	}
      }
    }/* end gbs term */


    { /* compute crystal rotation rates, accounting also for macroscopic spin, which should be precomputed and stored in marker[m].w?? */
      if(updateSpin){
	for(i=0;i<NT;i++){
	  for(j=0;j<NT;j++){
	    for(k=0;k<NT;k++){	  
	      PetscScalar CT = tex->ctheta[i][j][k];
	      PetscScalar CP = tex->cphi[i][j][k];    
	      /* compute nx, ny, nz */
	      PetscScalar nx = cos(CT) * sin(CP);
	      PetscScalar ny = sin(CT) * sin(CP);
	      PetscScalar nz = cos(CP);
	      
	      /* W_total_c = W_bulk - W_micromacro_total + W_micromacro_c */
	      const PetscScalar Wxx = 0.0;
	      const PetscScalar Wyy = 0.0;
	      const PetscScalar Wzz = 0.0;
	      PetscScalar Wxz = 0.0; //marker->wxz;
	      PetscScalar Wyz = 0.0; //marker->wyz;
	      PetscScalar Wxy = 0.0; //marker->wxy;     
	      /* this is Wmicromacro total */
	      /* 	      Wxx -= Mwtot[1-1][1-1]*S11 + Mwtot[1-1][2-1]*S22 + Mwtot[1-1][3-1]*S33 + 2.0*(Mwtot[1-1][4-1]*S23 + Mwtot[1-1][5-1]*S13 + Mwtot[1-1][6-1]*S12); */
	      /* 	      Wyy -= Mwtot[2-1][1-1]*S11 + Mwtot[2-1][2-1]*S22 + Mwtot[2-1][3-1]*S33 + 2.0*(Mwtot[2-1][4-1]*S23 + Mwtot[2-1][5-1]*S13 + Mwtot[2-1][6-1]*S12); */
	      /* 	      Wzz -= Mwtot[3-1][1-1]*S11 + Mwtot[3-1][2-1]*S22 + Mwtot[3-1][3-1]*S33 + 2.0*(Mwtot[3-1][4-1]*S23 + Mwtot[3-1][5-1]*S13 + Mwtot[3-1][6-1]*S12); */
	      Wxz -= Mwtot[4-1][1-1]*S11 + Mwtot[4-1][2-1]*S22 + Mwtot[4-1][3-1]*S33 + 2.0*(Mwtot[4-1][4-1]*S23 + Mwtot[4-1][5-1]*S13 + Mwtot[4-1][6-1]*S12);
	      Wyz -= Mwtot[5-1][1-1]*S11 + Mwtot[5-1][2-1]*S22 + Mwtot[5-1][3-1]*S33 + 2.0*(Mwtot[5-1][4-1]*S23 + Mwtot[5-1][5-1]*S13 + Mwtot[5-1][6-1]*S12);
	      Wxy -= Mwtot[6-1][1-1]*S11 + Mwtot[6-1][2-1]*S22 + Mwtot[6-1][3-1]*S33 + 2.0*(Mwtot[6-1][4-1]*S23 + Mwtot[6-1][5-1]*S13 + Mwtot[6-1][6-1]*S12);      
	      
	      /* calculate macroscopic-micromacro_total contribution to ndot */
	      PetscScalar nxdot =  Wxx*nx + Wxy*ny + Wxz*nz;
	      PetscScalar nydot = -Wxy*nx + Wyy*ny + Wyz*nz;
	      PetscScalar nzdot = -Wxz*nx - Wyz*ny + Wzz*nz;
	      /* thetadot and phidot are linear in nxdot, so we can just add the macroscopic-micromacro_total term to the single crystal term */
	      PetscScalar thetadot = -(ny/(nx*nx+ny*ny))*nxdot + nx/(nx*nx+ny*ny)*nydot;
	      if( isbad(thetadot) ){
		printf("thetadot is infinite, nx = %e, ny=%e\n",nx,ny);
	      }
	      PetscScalar phidot = -1/sqrt(1.0 - nz*nz)*nzdot;

	      tex->cthetadot[i][j][k] += thetadot;
	      tex->cphidot[i][j][k] += phidot;
	      if( isbad(tex->cphidot[i][j][k]) ){
		printf("NAN DEBUG: phidot is nan: nz=%e,nzdot=%e my contribution to phidot = %e\n",nz,nzdot,phidot);fflush(stdout);
	      }
	    }
	  }
	}
      }
    }


    /* Add diffusion creep rate */
    {
#ifdef DC_BARRSTILLMAN
      const PetscScalar QonRdiff = -7144.15;/* Q/R, Kelvin*/
      const PetscScalar Adiff = 3.019e-8; /*42/3*Vm*D_{0,v}/R*/
      PetscScalar Mdiff = Adiff/(d*d)/options->thermalBCBottom.value[0] * exp(QonRdiff/T);
      //	printf("Mdiff=%e\n",Mdiff);

#else
      /* diffusion creep from Goldsby and Kohlstedt 2001 */
      /* quantities from table 6: Diffusion Creep Parameters */
      const PetscScalar R = 8.3144621; /* ideal gas constant, J/mol/K */
      const PetscScalar b = 4.52e-10; /* burgers vector */
      const PetscScalar Vm = 1.97e-5; /* m^3 */
      const PetscScalar Qv = 59.4*1000.0;/* Activation energy for volume diffusion, J/mol */
      const PetscScalar delta = 2.0*b; /* grain boundary width, m */
      const PetscScalar Qb = 49*1000.0;
      const PetscScalar D0v = 9.10e-4; /* m^2/s */
      const PetscScalar D0b = 5.8e-4; /* m^2/s - this is lower of two upper bounde presented in text */
      /* diffusivities */
      const PetscScalar Dv = D0v*exp(-Qv/(R*T));
      const PetscScalar Db = D0b*exp(-Qb/(R*T));
      /* note that all units above are SI, so strain rate will be in 1/s for stress in Pa */
      PetscScalar Mdiff = 42.0*Vm/(R*T*d*d)*(Dv + M_PI*delta/d*Db );/* this is strain rate / stress */      
#endif
      Ms[0][0] += Mdiff;
      Ms[1][1] += Mdiff;
      Ms[2][2] += Mdiff;
      Ms[3][3] += Mdiff;
      Ms[4][4] += Mdiff;
      Ms[5][5] += Mdiff;
    }

    /* compute in- and out- of plane components of strain-rate tensor */
    PetscScalar Dxx = Ms[1-1][1-1]*S11 + Ms[1-1][2-1]*S22 + Ms[1-1][3-1]*S33 + 2.0*(Ms[1-1][4-1]*S23 + Ms[1-1][5-1]*S13 + Ms[1-1][6-1]*S12);
    PetscScalar Dyy = Ms[2-1][1-1]*S11 + Ms[2-1][2-1]*S22 + Ms[2-1][3-1]*S33 + 2.0*(Ms[2-1][4-1]*S23 + Ms[2-1][5-1]*S13 + Ms[2-1][6-1]*S12);
    PetscScalar Dzz = Ms[3-1][1-1]*S11 + Ms[3-1][2-1]*S22 + Ms[3-1][3-1]*S33 + 2.0*(Ms[3-1][4-1]*S23 + Ms[3-1][5-1]*S13 + Ms[3-1][6-1]*S12);
      
    PetscScalar Dxz = Ms[4-1][1-1]*S11 + Ms[4-1][2-1]*S22 + Ms[4-1][3-1]*S33 + 2.0*(Ms[4-1][4-1]*S23 + Ms[4-1][5-1]*S13 + Ms[4-1][6-1]*S12);
    PetscScalar Dyz = Ms[5-1][1-1]*S11 + Ms[5-1][2-1]*S22 + Ms[5-1][3-1]*S33 + 2.0*(Ms[5-1][4-1]*S23 + Ms[5-1][5-1]*S13 + Ms[5-1][6-1]*S12);
    PetscScalar Dxy = Ms[6-1][1-1]*S11 + Ms[6-1][2-1]*S22 + Ms[6-1][3-1]*S33 + 2.0*(Ms[6-1][4-1]*S23 + Ms[6-1][5-1]*S13 + Ms[6-1][6-1]*S12);
    /* PetscScalar Dii = sqrt( 0.5*(Dxx*Dxx+Dyy*Dyy)+Dxy*Dxy);  */
      
    D->T11 = Dxx;
    D->T22 = Dyy;
    D->T33 = Dzz;
    D->T13 = Dxz;
    D->T23 = Dyz;
    D->T12 = Dxy;

#ifdef verbose
    /* print out M */
    printf("M = \n");
    for(i=0;i<6;i++){
      for(j=0;j<6;j++){
	printf("%e\t",Ms[i][j]);
      }
      printf("\n");
    }
#endif

  }/* END OUTER NEWTON LOOP*/

    


  /*}*//* end loop over markers*/
  PetscFunctionReturn(ierr);
}

PetscInt idxmod(PetscInt i){
  switch (i)
    {
    case -1:
      return NT-1;
    case NT:
      return 0;
    default:
      return i;
    }
}

void solve33(PetscScalar *a,PetscScalar *b){
  /* solve 3x3 linear system using analytic expressions */
  PetscScalar v000 = b[2];
  PetscScalar v001 = a[4];
  PetscScalar v002 = a[0];
  PetscScalar v003 = a[3];
  PetscScalar v004 = a[1];
  PetscScalar v005 = b[1];
  PetscScalar v006 = a[7];
  PetscScalar v007 = a[6];
  PetscScalar v008 = b[0];
  PetscScalar v009 = a[8];
  PetscScalar v010 = a[5];
  PetscScalar v011 = a[2];
  PetscScalar v012 = 1/(-(v001*v002*v009) + v003*v004*v009 +		\
			v002*v006*v010 - v004*v007*v010 - v003*v006*v011 + v001*v007*v011);

  b[0] = (v004*v005*v009 - v001*v008*v009 - v000*v004*v010 +	\
	  v006*v008*v010 + v000*v001*v011 - v005*v006*v011)*v012;
  b[1] = ((v000*v002 - v007*v008)*v010 + v003*(v008*v009 - v000*v011) + \
	  v005*(-(v002*v009) + v007*v011))/(v002*(-(v001*v009) + v006*v010) + \
					    v004*(v003*v009 - v007*v010) + (-(v003*v006) + v001*v007)*v011);
  b[2] = (-(v000*v001*v002) + v000*v003*v004 + v002*v005*v006 - \
	  v004*v005*v007 - v003*v006*v008 + v001*v007*v008)*v012;

}

PetscErrorCode invertMtoN(TextureInfo *tex){
  /* this subroutine analytically computes a portion of the matrix inverse of the viscoplasticity tensor M */

  PetscErrorCode ierr=0;
  
  PetscScalar (*M)[6] = tex->M;
  PetscScalar (*N)[2] = tex->N;

  
  /* Begin machine-generated code */

PetscScalar v000000 = M[3-1][2-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000001 = M[3-1][1-1]*M[4-1][5-1]*M[5-1][2-1];
PetscScalar v000002 = M[3-1][5-1]*M[4-1][2-1]*M[5-1][1-1];
PetscScalar v000003 = M[3-1][1-1]*M[4-1][2-1]*M[5-1][3-1];
PetscScalar v000004 = M[3-1][3-1]*M[4-1][1-1]*M[5-1][2-1];
PetscScalar v000005 = M[3-1][2-1]*M[4-1][3-1]*M[5-1][1-1];
PetscScalar v000006 = M[3-1][1-1]*M[4-1][2-1]*M[5-1][4-1];
PetscScalar v000007 = M[3-1][4-1]*M[4-1][1-1]*M[5-1][2-1];
PetscScalar v000008 = M[3-1][2-1]*M[4-1][4-1]*M[5-1][1-1];
PetscScalar v000009 = M[3-1][2-1]*M[4-1][1-1]*M[5-1][3-1];
PetscScalar v000010 = M[3-1][1-1]*M[4-1][3-1]*M[5-1][2-1];
PetscScalar v000011 = M[3-1][3-1]*M[4-1][2-1]*M[5-1][1-1];
PetscScalar v000012 = M[3-1][2-1]*M[4-1][1-1]*M[5-1][6-1];
PetscScalar v000013 = M[3-1][1-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000014 = M[3-1][6-1]*M[4-1][2-1]*M[5-1][1-1];
PetscScalar v000015 = M[3-1][1-1]*M[4-1][2-1]*M[5-1][5-1];
PetscScalar v000016 = M[3-1][5-1]*M[4-1][1-1]*M[5-1][2-1];
PetscScalar v000017 = M[3-1][2-1]*M[4-1][5-1]*M[5-1][1-1];
PetscScalar v000018 = M[3-1][2-1]*M[4-1][1-1]*M[5-1][4-1];
PetscScalar v000019 = M[3-1][1-1]*M[4-1][4-1]*M[5-1][2-1];
PetscScalar v000020 = M[3-1][4-1]*M[4-1][2-1]*M[5-1][1-1];
PetscScalar v000021 = M[3-1][3-1]*M[4-1][1-1]*M[5-1][6-1];
PetscScalar v000022 = M[3-1][1-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000023 = M[3-1][6-1]*M[4-1][3-1]*M[5-1][1-1];
PetscScalar v000024 = M[3-1][1-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000025 = M[3-1][4-1]*M[4-1][1-1]*M[5-1][3-1];
PetscScalar v000026 = M[3-1][3-1]*M[4-1][4-1]*M[5-1][1-1];
PetscScalar v000027 = M[3-1][1-1]*M[4-1][3-1]*M[5-1][5-1];
PetscScalar v000028 = M[3-1][5-1]*M[4-1][1-1]*M[5-1][3-1];
PetscScalar v000029 = M[3-1][3-1]*M[4-1][5-1]*M[5-1][1-1];
PetscScalar v000030 = M[3-1][3-1]*M[4-1][1-1]*M[5-1][4-1];
PetscScalar v000031 = M[3-1][1-1]*M[4-1][4-1]*M[5-1][3-1];
PetscScalar v000032 = M[3-1][4-1]*M[4-1][3-1]*M[5-1][1-1];
PetscScalar v000033 = M[3-1][3-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000034 = M[3-1][2-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000035 = M[3-1][6-1]*M[4-1][3-1]*M[5-1][2-1];
PetscScalar v000036 = M[3-1][2-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000037 = M[3-1][4-1]*M[4-1][2-1]*M[5-1][3-1];
PetscScalar v000038 = M[3-1][3-1]*M[4-1][4-1]*M[5-1][2-1];
PetscScalar v000039 = M[3-1][2-1]*M[4-1][3-1]*M[5-1][5-1];
PetscScalar v000040 = M[3-1][5-1]*M[4-1][2-1]*M[5-1][3-1];
PetscScalar v000041 = M[3-1][3-1]*M[4-1][5-1]*M[5-1][2-1];
PetscScalar v000042 = M[3-1][3-1]*M[4-1][2-1]*M[5-1][4-1];
PetscScalar v000043 = M[3-1][2-1]*M[4-1][4-1]*M[5-1][3-1];
PetscScalar v000044 = M[3-1][4-1]*M[4-1][3-1]*M[5-1][2-1];
PetscScalar v000045 = M[4-1][1-1]*M[5-1][6-1]*M[6-1][2-1];
PetscScalar v000046 = M[4-1][6-1]*M[5-1][2-1]*M[6-1][1-1];
PetscScalar v000047 = M[4-1][1-1]*M[5-1][2-1]*M[6-1][3-1];
PetscScalar v000048 = M[4-1][3-1]*M[5-1][1-1]*M[6-1][2-1];
PetscScalar v000049 = M[4-1][2-1]*M[5-1][3-1]*M[6-1][1-1];
PetscScalar v000050 = M[4-1][1-1]*M[5-1][2-1]*M[6-1][4-1];
PetscScalar v000051 = M[4-1][4-1]*M[5-1][1-1]*M[6-1][2-1];
PetscScalar v000052 = M[4-1][2-1]*M[5-1][4-1]*M[6-1][1-1];
PetscScalar v000053 = M[4-1][2-1]*M[5-1][1-1]*M[6-1][3-1];
PetscScalar v000054 = M[4-1][1-1]*M[5-1][3-1]*M[6-1][2-1];
PetscScalar v000055 = M[4-1][3-1]*M[5-1][2-1]*M[6-1][1-1];
PetscScalar v000056 = M[4-1][1-1]*M[5-1][2-1]*M[6-1][5-1];
PetscScalar v000057 = M[4-1][5-1]*M[5-1][1-1]*M[6-1][2-1];
PetscScalar v000058 = M[4-1][2-1]*M[5-1][5-1]*M[6-1][1-1];
PetscScalar v000059 = M[4-1][2-1]*M[5-1][1-1]*M[6-1][4-1];
PetscScalar v000060 = M[4-1][1-1]*M[5-1][4-1]*M[6-1][2-1];
PetscScalar v000061 = M[4-1][4-1]*M[5-1][2-1]*M[6-1][1-1];
PetscScalar v000062 = M[4-1][2-1]*M[5-1][3-1]*M[6-1][5-1];
PetscScalar v000063 = M[4-1][5-1]*M[5-1][2-1]*M[6-1][3-1];
PetscScalar v000064 = M[4-1][3-1]*M[5-1][5-1]*M[6-1][2-1];
PetscScalar v000065 = M[4-1][3-1]*M[5-1][2-1]*M[6-1][4-1];
PetscScalar v000066 = M[4-1][2-1]*M[5-1][4-1]*M[6-1][3-1];
PetscScalar v000067 = M[4-1][4-1]*M[5-1][3-1]*M[6-1][2-1];
PetscScalar v000068 = -v000015;
PetscScalar v000069 = -v000016;
PetscScalar v000070 = -v000017;
PetscScalar v000071 = M[2-1][2-1]*v000027;
PetscScalar v000072 = M[2-1][1-1]*M[3-1][3-1]*M[4-1][2-1]*M[5-1][5-1];
PetscScalar v000073 = -v000009;
PetscScalar v000074 = -v000010;
PetscScalar v000075 = -v000011;
PetscScalar v000076 = M[2-1][1-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000077 = M[2-1][2-1]*v000028;
PetscScalar v000078 = M[2-1][1-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][2-1];
PetscScalar v000079 = M[2-1][2-1]*v000029;
PetscScalar v000080 = -v000018;
PetscScalar v000081 = -v000019;
PetscScalar v000082 = -v000020;
PetscScalar v000083 = M[2-1][1-1]*v000036;
PetscScalar v000084 = M[2-1][2-1]*v000030;
PetscScalar v000085 = -v000003;
PetscScalar v000086 = -v000004;
PetscScalar v000087 = -v000005;
PetscScalar v000088 = M[2-1][2-1]*v000031;
PetscScalar v000089 = M[2-1][1-1]*v000037;
PetscScalar v000090 = M[2-1][1-1]*v000038;
PetscScalar v000091 = M[2-1][2-1]*v000032;
PetscScalar v000092 = -(M[3-1][1-1]*M[4-1][2-1]*M[5-1][6-1]);
PetscScalar v000093 = -(M[3-1][2-1]*M[4-1][6-1]*M[5-1][1-1]);
PetscScalar v000094 = M[2-1][2-1]*M[3-1][1-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000095 = M[2-1][1-1]*v000033;
PetscScalar v000096 = M[2-1][1-1]*v000034;
PetscScalar v000097 = M[2-1][2-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][3-1];
PetscScalar v000098 = M[2-1][1-1]*v000035;
PetscScalar v000099 = M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000100 = -v000000;
PetscScalar v000101 = -v000001;
PetscScalar v000102 = -v000002;
PetscScalar v000103 = M[2-1][1-1]*v000039;
PetscScalar v000104 = M[2-1][2-1]*M[3-1][3-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000105 = M[2-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000106 = M[2-1][1-1]*v000040;
PetscScalar v000107 = M[2-1][1-1]*v000041;
PetscScalar v000108 = M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][1-1];
PetscScalar v000109 = M[2-1][2-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000110 = M[2-1][1-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000111 = M[2-1][1-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000112 = M[2-1][2-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][4-1];
PetscScalar v000113 = M[2-1][1-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1];
PetscScalar v000114 = M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000115 = M[2-1][1-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000116 = M[2-1][2-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000117 = -v000006;
PetscScalar v000118 = -v000007;
PetscScalar v000119 = -v000008;
PetscScalar v000120 = M[2-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000121 = M[2-1][1-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][4-1];
PetscScalar v000122 = M[2-1][1-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1];
PetscScalar v000123 = M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][1-1];
PetscScalar v000124 = -(M[3-1][1-1]*M[4-1][3-1]*M[5-1][6-1]);
PetscScalar v000125 = -(M[3-1][6-1]*M[4-1][1-1]*M[5-1][3-1]);
PetscScalar v000126 = -(M[3-1][3-1]*M[4-1][6-1]*M[5-1][1-1]);
PetscScalar v000127 = M[2-1][3-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000128 = M[2-1][1-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000129 = -v000030;
PetscScalar v000130 = -v000031;
PetscScalar v000131 = -v000032;
PetscScalar v000132 = M[2-1][1-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000133 = M[2-1][3-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][4-1];
PetscScalar v000134 = M[2-1][1-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1];
PetscScalar v000135 = M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000136 = -(M[3-1][3-1]*M[4-1][1-1]*M[5-1][5-1]);
PetscScalar v000137 = -(M[3-1][1-1]*M[4-1][5-1]*M[5-1][3-1]);
PetscScalar v000138 = -(M[3-1][5-1]*M[4-1][3-1]*M[5-1][1-1]);
PetscScalar v000139 = M[2-1][1-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000140 = M[2-1][3-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000141 = -v000024;
PetscScalar v000142 = -v000025;
PetscScalar v000143 = -v000026;
PetscScalar v000144 = M[2-1][3-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000145 = M[2-1][1-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000146 = M[2-1][1-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000147 = M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][1-1];
PetscScalar v000148 = -(M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1]);
PetscScalar v000149 = -(M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1]);
PetscScalar v000150 = -(M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1]);
PetscScalar v000151 = M[2-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000152 = M[2-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000153 = -v000042;
PetscScalar v000154 = -v000043;
PetscScalar v000155 = -v000044;
PetscScalar v000156 = M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000157 = M[2-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1];
PetscScalar v000158 = M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1];
PetscScalar v000159 = M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000160 = -(M[3-1][3-1]*M[4-1][2-1]*M[5-1][5-1]);
PetscScalar v000161 = -(M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1]);
PetscScalar v000162 = -(M[3-1][5-1]*M[4-1][3-1]*M[5-1][2-1]);
PetscScalar v000163 = M[2-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000164 = M[2-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1];
PetscScalar v000165 = -v000036;
PetscScalar v000166 = -v000037;
PetscScalar v000167 = -v000038;
PetscScalar v000168 = M[2-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000169 = M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000170 = M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000171 = M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1];
PetscScalar v000172 = -(M[4-1][6-1]*M[5-1][1-1]*M[6-1][2-1]);
PetscScalar v000173 = -(M[4-1][2-1]*M[5-1][6-1]*M[6-1][1-1]);
PetscScalar v000174 = -v000055;
PetscScalar v000175 = M[6-1][2-1]*v000022;
PetscScalar v000176 = -v000060;
PetscScalar v000177 = -v000061;
PetscScalar v000178 = -v000047;
PetscScalar v000179 = -v000048;
PetscScalar v000180 = -v000049;
PetscScalar v000181 = -(M[4-1][1-1]*M[5-1][5-1]*M[6-1][2-1]);
PetscScalar v000182 = -(M[4-1][5-1]*M[5-1][2-1]*M[6-1][1-1]);
PetscScalar v000183 = M[6-1][5-1]*v000003;
PetscScalar v000184 = M[6-1][5-1]*v000005;
PetscScalar v000185 = M[6-1][3-1]*v000000;
PetscScalar v000186 = M[6-1][3-1]*v000001;
PetscScalar v000187 = M[6-1][2-1]*v000027;
PetscScalar v000188 = M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][1-1];
PetscScalar v000189 = -v000050;
PetscScalar v000190 = -v000051;
PetscScalar v000191 = -v000052;
PetscScalar v000192 = -(M[4-1][2-1]*M[5-1][5-1]*M[6-1][3-1]);
PetscScalar v000193 = -(M[4-1][5-1]*M[5-1][3-1]*M[6-1][2-1]);
PetscScalar v000194 = M[6-1][5-1]*v000036;
PetscScalar v000195 = M[6-1][5-1]*v000038;
PetscScalar v000196 = -(M[4-1][4-1]*M[5-1][2-1]*M[6-1][3-1]);
PetscScalar v000197 = -(M[4-1][3-1]*M[5-1][4-1]*M[6-1][2-1]);
PetscScalar v000198 = M[3-1][3-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][4-1];
PetscScalar v000199 = M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][4-1];
PetscScalar v000200 = M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][3-1];
PetscScalar v000201 = M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][2-1];
PetscScalar v000202 = -v000103;
PetscScalar v000203 = -v000104;
PetscScalar v000204 = -v000105;
PetscScalar v000205 = -v000106;
PetscScalar v000206 = -v000107;
PetscScalar v000207 = -v000108;
PetscScalar v000208 = M[1-1][1-1]*v000163;
PetscScalar v000209 = M[1-1][3-1]*v000115;
PetscScalar v000210 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000211 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1];
PetscScalar v000212 = M[1-1][1-1]*M[2-1][4-1]*v000039;
PetscScalar v000213 = M[1-1][1-1]*v000164;
PetscScalar v000214 = M[1-1][3-1]*M[2-1][4-1]*v000015;
PetscScalar v000215 = M[1-1][3-1]*v000116;
PetscScalar v000216 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000217 = M[2-1][2-1]*v000141;
PetscScalar v000218 = M[2-1][1-1]*v000153;
PetscScalar v000219 = M[2-1][1-1]*v000154;
PetscScalar v000220 = M[2-1][2-1]*v000142;
PetscScalar v000221 = M[2-1][1-1]*v000155;
PetscScalar v000222 = M[2-1][2-1]*v000143;
PetscScalar v000223 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000224 = M[1-1][1-1]*v000168;
PetscScalar v000225 = M[1-1][3-1]*v000120;
PetscScalar v000226 = M[1-1][1-1]*v000169;
PetscScalar v000227 = M[1-1][2-1]*M[2-1][5-1]*v000024;
PetscScalar v000228 = M[1-1][3-1]*v000121;
PetscScalar v000229 = M[1-1][1-1]*M[2-1][5-1]*v000042;
PetscScalar v000230 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1];
PetscScalar v000231 = M[1-1][3-1]*M[2-1][5-1]*v000018;
PetscScalar v000232 = M[1-1][1-1]*v000170;
PetscScalar v000233 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000234 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1];
PetscScalar v000235 = M[1-1][1-1]*M[2-1][5-1]*v000043;
PetscScalar v000236 = M[1-1][1-1]*M[2-1][4-1]*v000040;
PetscScalar v000237 = M[1-1][2-1]*M[2-1][5-1]*v000025;
PetscScalar v000238 = M[1-1][3-1]*v000122;
PetscScalar v000239 = M[1-1][1-1]*M[2-1][4-1]*v000041;
PetscScalar v000240 = M[1-1][1-1]*v000171;
PetscScalar v000241 = M[1-1][3-1]*M[2-1][5-1]*v000019;
PetscScalar v000242 = M[1-1][1-1]*M[2-1][5-1]*v000044;
PetscScalar v000243 = M[1-1][3-1]*M[2-1][4-1]*v000016;
PetscScalar v000244 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1];
PetscScalar v000245 = M[1-1][3-1]*M[2-1][4-1]*v000017;
PetscScalar v000246 = M[1-1][3-1]*v000123;
PetscScalar v000247 = M[1-1][2-1]*M[2-1][5-1]*v000026;
PetscScalar v000248 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][1-1];
PetscScalar v000249 = M[1-1][3-1]*M[2-1][5-1]*v000020;
PetscScalar v000250 = M[2-1][1-1]*v000148;
PetscScalar v000251 = -(M[2-1][2-1]*v000021);
PetscScalar v000252 = -(M[2-1][2-1]*v000022);
PetscScalar v000253 = M[2-1][1-1]*v000149;
PetscScalar v000254 = M[2-1][1-1]*v000150;
PetscScalar v000255 = -(M[2-1][2-1]*v000023);
PetscScalar v000256 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000257 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000258 = M[1-1][2-1]*v000127;
PetscScalar v000259 = M[1-1][2-1]*v000128;
PetscScalar v000260 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000261 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1];
PetscScalar v000262 = M[1-1][2-1]*M[2-1][4-1]*v000021;
PetscScalar v000263 = M[1-1][2-1]*v000132;
PetscScalar v000264 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000265 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000266 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000267 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1];
PetscScalar v000268 = M[1-1][2-1]*v000133;
PetscScalar v000269 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000270 = M[1-1][2-1]*M[2-1][4-1]*v000022;
PetscScalar v000271 = M[1-1][2-1]*v000134;
PetscScalar v000272 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000273 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1];
PetscScalar v000274 = M[1-1][2-1]*v000135;
PetscScalar v000275 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1];
PetscScalar v000276 = M[1-1][2-1]*M[2-1][4-1]*v000023;
PetscScalar v000277 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000278 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000279 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000280 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000281 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000282 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000283 = M[1-1][3-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000284 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1];
PetscScalar v000285 = M[1-1][2-1]*M[2-1][5-1]*v000021;
PetscScalar v000286 = -v000071;
PetscScalar v000287 = -v000072;
PetscScalar v000288 = -v000076;
PetscScalar v000289 = -v000077;
PetscScalar v000290 = -v000079;
PetscScalar v000291 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000292 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000293 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000294 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1];
PetscScalar v000295 = M[1-1][2-1]*M[2-1][6-1]*v000027;
PetscScalar v000296 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1];
PetscScalar v000297 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000298 = M[1-1][3-1]*M[2-1][6-1]*v000000;
PetscScalar v000299 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000300 = M[1-1][2-1]*M[2-1][5-1]*v000022;
PetscScalar v000301 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000302 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1];
PetscScalar v000303 = M[1-1][2-1]*M[2-1][6-1]*v000028;
PetscScalar v000304 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000305 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000306 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1];
PetscScalar v000307 = M[1-1][3-1]*M[2-1][6-1]*v000001;
PetscScalar v000308 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000309 = M[1-1][3-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000310 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1];
PetscScalar v000311 = M[1-1][2-1]*M[2-1][6-1]*v000029;
PetscScalar v000312 = M[1-1][2-1]*M[2-1][5-1]*v000023;
PetscScalar v000313 = M[1-1][3-1]*M[2-1][6-1]*v000002;
PetscScalar v000314 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000315 = M[1-1][4-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000316 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000317 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000318 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000319 = M[1-1][1-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000320 = M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1];
PetscScalar v000321 = M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1];
PetscScalar v000322 = -(M[2-1][2-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1]);
PetscScalar v000323 = -(M[2-1][1-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1]);
PetscScalar v000324 = -(M[2-1][1-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]);
PetscScalar v000325 = -(M[2-1][2-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1]);
PetscScalar v000326 = -(M[2-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1]);
PetscScalar v000327 = -(M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1]);
PetscScalar v000328 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000329 = M[1-1][1-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000330 = M[1-1][4-1]*M[2-1][2-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000331 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000332 = M[1-1][2-1]*M[2-1][6-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000333 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000334 = M[1-1][1-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000335 = M[1-1][2-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000336 = M[1-1][2-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000337 = M[1-1][1-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000338 = M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1];
PetscScalar v000339 = M[1-1][4-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000340 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000341 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000342 = M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1];
PetscScalar v000343 = M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1];
PetscScalar v000344 = M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1];
PetscScalar v000345 = -(M[2-1][1-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]);
PetscScalar v000346 = -(M[2-1][3-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1]);
PetscScalar v000347 = -(M[2-1][3-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1]);
PetscScalar v000348 = -(M[2-1][1-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]);
PetscScalar v000349 = -(M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]);
PetscScalar v000350 = -(M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1]);
PetscScalar v000351 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000352 = M[1-1][4-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000353 = M[1-1][3-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000354 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000355 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000356 = M[1-1][1-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000357 = M[1-1][4-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000358 = -(M[2-1][3-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1]);
PetscScalar v000359 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000360 = M[1-1][1-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000361 = M[1-1][4-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000362 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000363 = M[1-1][3-1]*M[2-1][6-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000364 = M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000365 = M[1-1][4-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][1-1]*M[5-1][5-1];
PetscScalar v000366 = M[1-1][1-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000367 = M[1-1][3-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000368 = M[1-1][3-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000369 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000370 = M[1-1][4-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000371 = M[1-1][1-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000372 = M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000373 = M[1-1][4-1]*M[2-1][6-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000374 = M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000375 = M[1-1][4-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][1-1];
PetscScalar v000376 = -(M[2-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]);
PetscScalar v000377 = -(M[2-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1]);
PetscScalar v000378 = -(M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1]);
PetscScalar v000379 = -(M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]);
PetscScalar v000380 = -(M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]);
PetscScalar v000381 = -(M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1]);
PetscScalar v000382 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000383 = M[1-1][4-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000384 = M[1-1][3-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1];
PetscScalar v000385 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000386 = M[1-1][2-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1];
PetscScalar v000387 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1];
PetscScalar v000388 = M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000389 = M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1];
PetscScalar v000390 = -(M[2-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1]);
PetscScalar v000391 = -(M[2-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1]);
PetscScalar v000392 = -(M[2-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1]);
PetscScalar v000393 = -(M[2-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][4-1]);
PetscScalar v000394 = -(M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1]);
PetscScalar v000395 = -(M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1]);
PetscScalar v000396 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000397 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000398 = M[1-1][4-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1];
PetscScalar v000399 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000400 = M[1-1][3-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1];
PetscScalar v000401 = M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1];
PetscScalar v000402 = M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1];
PetscScalar v000403 = M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1];
PetscScalar v000404 = M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000405 = M[1-1][3-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1];
PetscScalar v000406 = M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000407 = M[1-1][2-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1];
PetscScalar v000408 = M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1];
PetscScalar v000409 = M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][4-1];
PetscScalar v000410 = M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000411 = M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1];
PetscScalar v000412 = M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1];
PetscScalar v000413 = M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1];
PetscScalar v000414 = M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1];
PetscScalar v000415 = M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1];
PetscScalar v000416 = M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1];
PetscScalar v000417 = M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1];
PetscScalar v000418 = M[6-1][5-1]*v000073;
PetscScalar v000419 = M[6-1][5-1]*v000074;
PetscScalar v000420 = M[6-1][3-1]*v000068;
PetscScalar v000421 = M[6-1][3-1]*v000070;
PetscScalar v000422 = M[6-1][2-1]*v000137;
PetscScalar v000423 = -(M[6-1][1-1]*v000039);
PetscScalar v000424 = M[2-1][2-1]*M[6-1][5-1]*v000024;
PetscScalar v000425 = M[2-1][1-1]*M[6-1][5-1]*v000042;
PetscScalar v000426 = M[2-1][3-1]*M[6-1][5-1]*v000018;
PetscScalar v000427 = M[2-1][1-1]*M[6-1][5-1]*v000043;
PetscScalar v000428 = M[2-1][2-1]*M[6-1][5-1]*v000025;
PetscScalar v000429 = M[2-1][3-1]*M[6-1][5-1]*v000019;
PetscScalar v000430 = M[2-1][1-1]*M[6-1][5-1]*v000044;
PetscScalar v000431 = M[2-1][2-1]*M[6-1][5-1]*v000026;
PetscScalar v000432 = M[2-1][3-1]*M[6-1][5-1]*v000020;
PetscScalar v000433 = M[6-1][4-1]*v000103;
PetscScalar v000434 = M[2-1][3-1]*M[6-1][4-1]*v000015;
PetscScalar v000435 = M[6-1][4-1]*v000104;
PetscScalar v000436 = M[6-1][4-1]*v000105;
PetscScalar v000437 = M[6-1][4-1]*v000106;
PetscScalar v000438 = M[6-1][4-1]*v000107;
PetscScalar v000439 = M[2-1][3-1]*M[6-1][4-1]*v000016;
PetscScalar v000440 = M[2-1][3-1]*M[6-1][4-1]*v000017;
PetscScalar v000441 = M[6-1][4-1]*v000108;
PetscScalar v000442 = M[2-1][2-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][3-1];
PetscScalar v000443 = M[2-1][1-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][3-1];
PetscScalar v000444 = M[2-1][1-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1];
PetscScalar v000445 = M[2-1][2-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1]*M[6-1][3-1];
PetscScalar v000446 = M[2-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][3-1];
PetscScalar v000447 = M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1]*M[6-1][3-1];
PetscScalar v000448 = M[6-1][2-1]*v000139;
PetscScalar v000449 = M[6-1][2-1]*v000140;
PetscScalar v000450 = M[6-1][2-1]*v000144;
PetscScalar v000451 = M[6-1][2-1]*v000145;
PetscScalar v000452 = M[6-1][2-1]*v000146;
PetscScalar v000453 = M[2-1][3-1]*M[3-1][5-1]*v000051;
PetscScalar v000454 = M[2-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][1-1];
PetscScalar v000455 = M[2-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1]*M[6-1][1-1];
PetscScalar v000456 = M[2-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][1-1];
PetscScalar v000457 = M[2-1][3-1]*M[3-1][5-1]*v000052;
PetscScalar v000458 = M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1]*M[6-1][1-1];
PetscScalar v000459 = M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1]*M[6-1][1-1];
PetscScalar v000460 = M[6-1][5-1]*v000153;
PetscScalar v000461 = M[6-1][5-1]*v000154;
PetscScalar v000462 = -(M[6-1][4-1]*v000039);
PetscScalar v000463 = -(M[6-1][4-1]*v000041);
PetscScalar v000464 = -(M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][2-1]);
PetscScalar v000465 = -(M[1-1][2-1]*v000139);
PetscScalar v000466 = M[1-1][1-1]*v000390;
PetscScalar v000467 = M[1-1][3-1]*v000322;
PetscScalar v000468 = M[1-1][1-1]*v000391;
PetscScalar v000469 = -(M[1-1][2-1]*M[2-1][4-1]*v000027);
PetscScalar v000470 = M[1-1][3-1]*v000323;
PetscScalar v000471 = M[1-1][1-1]*M[2-1][4-1]*v000160;
PetscScalar v000472 = -(M[1-1][2-1]*v000140);
PetscScalar v000473 = M[1-1][3-1]*M[2-1][4-1]*v000100;
PetscScalar v000474 = M[1-1][1-1]*v000392;
PetscScalar v000475 = M[1-1][3-1]*v000324;
PetscScalar v000476 = -(M[1-1][2-1]*v000144);
PetscScalar v000477 = -(M[1-1][2-1]*v000145);
PetscScalar v000478 = M[1-1][1-1]*M[2-1][5-1]*v000165;
PetscScalar v000479 = M[1-1][1-1]*v000393;
PetscScalar v000480 = M[1-1][3-1]*M[2-1][5-1]*v000117;
PetscScalar v000481 = M[1-1][3-1]*v000325;
PetscScalar v000482 = M[1-1][2-1]*M[2-1][5-1]*v000129;
PetscScalar v000483 = -(M[1-1][2-1]*v000146);
PetscScalar v000484 = M[1-1][1-1]*M[2-1][4-1]*v000161;
PetscScalar v000485 = M[1-1][1-1]*v000394;
PetscScalar v000486 = M[1-1][2-1]*M[2-1][5-1]*v000130;
PetscScalar v000487 = M[1-1][1-1]*M[2-1][5-1]*v000166;
PetscScalar v000488 = -(M[1-1][2-1]*M[2-1][4-1]*v000028);
PetscScalar v000489 = M[1-1][1-1]*v000395;
PetscScalar v000490 = M[1-1][3-1]*M[2-1][4-1]*v000101;
PetscScalar v000491 = M[1-1][3-1]*v000326;
PetscScalar v000492 = M[1-1][1-1]*M[2-1][5-1]*v000167;
PetscScalar v000493 = M[1-1][1-1]*M[2-1][4-1]*v000162;
PetscScalar v000494 = M[1-1][3-1]*M[2-1][5-1]*v000118;
PetscScalar v000495 = M[1-1][3-1]*v000327;
PetscScalar v000496 = -(M[1-1][2-1]*M[2-1][4-1]*v000029);
PetscScalar v000497 = -(M[1-1][2-1]*v000147);
PetscScalar v000498 = M[1-1][3-1]*M[2-1][5-1]*v000119;
PetscScalar v000499 = M[1-1][2-1]*M[2-1][5-1]*v000131;
PetscScalar v000500 = M[1-1][3-1]*M[2-1][4-1]*v000102;
PetscScalar v000501 = -(M[1-1][2-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]);
PetscScalar v000502 = -(M[1-1][1-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]);
PetscScalar v000503 = -(M[1-1][3-1]*M[2-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1]);
PetscScalar v000504 = -(M[1-1][1-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]);
PetscScalar v000505 = M[1-1][2-1]*M[2-1][5-1]*v000124;
PetscScalar v000506 = -(M[1-1][3-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]);
PetscScalar v000507 = -(M[1-1][1-1]*M[2-1][5-1]*v000033);
PetscScalar v000508 = -(M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1]);
PetscScalar v000509 = -(M[1-1][3-1]*M[2-1][5-1]*v000012);
PetscScalar v000510 = -(M[1-1][1-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]);
PetscScalar v000511 = -(M[1-1][3-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]);
PetscScalar v000512 = -(M[1-1][2-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1]);
PetscScalar v000513 = -(M[1-1][2-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1]);
PetscScalar v000514 = -(M[1-1][1-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1]);
PetscScalar v000515 = M[1-1][3-1]*M[2-1][6-1]*v000068;
PetscScalar v000516 = -(M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1]);
PetscScalar v000517 = M[1-1][2-1]*M[2-1][6-1]*v000136;
PetscScalar v000518 = -(M[1-1][2-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]);
PetscScalar v000519 = -(M[1-1][1-1]*M[2-1][5-1]*v000034);
PetscScalar v000520 = -(M[1-1][1-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1]);
PetscScalar v000521 = M[1-1][2-1]*M[2-1][6-1]*v000137;
PetscScalar v000522 = -(M[1-1][1-1]*M[2-1][6-1]*v000040);
PetscScalar v000523 = M[1-1][2-1]*M[2-1][5-1]*v000125;
PetscScalar v000524 = -(M[1-1][1-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1]);
PetscScalar v000525 = -(M[1-1][3-1]*M[2-1][5-1]*v000013);
PetscScalar v000526 = -(M[1-1][3-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1]);
PetscScalar v000527 = -(M[1-1][1-1]*M[2-1][5-1]*v000035);
PetscScalar v000528 = M[1-1][3-1]*M[2-1][6-1]*v000069;
PetscScalar v000529 = -(M[1-1][3-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1]);
PetscScalar v000530 = M[1-1][2-1]*M[2-1][5-1]*v000126;
PetscScalar v000531 = -(M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1]);
PetscScalar v000532 = M[1-1][3-1]*M[2-1][6-1]*v000070;
PetscScalar v000533 = M[1-1][2-1]*M[2-1][6-1]*v000138;
PetscScalar v000534 = -(M[1-1][3-1]*M[2-1][5-1]*v000014);
PetscScalar v000535 = -(M[1-1][3-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]);
PetscScalar v000536 = -(M[1-1][1-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]);
PetscScalar v000537 = -(M[1-1][4-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1]);
PetscScalar v000538 = -(M[1-1][1-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]);
PetscScalar v000539 = -(M[1-1][3-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][6-1]);
PetscScalar v000540 = -(M[1-1][4-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]);
PetscScalar v000541 = -(M[1-1][1-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1]);
PetscScalar v000542 = -(M[1-1][4-1]*M[2-1][5-1]*v000021);
PetscScalar v000543 = -(M[1-1][1-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]);
PetscScalar v000544 = -(M[1-1][4-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]);
PetscScalar v000545 = -(M[1-1][3-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1]);
PetscScalar v000546 = -(M[1-1][3-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]);
PetscScalar v000547 = -(M[1-1][3-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]);
PetscScalar v000548 = -(M[1-1][1-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1]);
PetscScalar v000549 = -(M[1-1][1-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]);
PetscScalar v000550 = -(M[1-1][3-1]*M[2-1][6-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1]);
PetscScalar v000551 = -(M[1-1][1-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]);
PetscScalar v000552 = -(M[1-1][4-1]*M[2-1][6-1]*v000028);
PetscScalar v000553 = -(M[1-1][4-1]*M[2-1][6-1]*v000029);
PetscScalar v000554 = -(M[1-1][4-1]*M[2-1][5-1]*v000023);
PetscScalar v000555 = -(M[6-1][5-1]*v000083);
PetscScalar v000556 = M[2-1][3-1]*M[6-1][5-1]*v000117;
PetscScalar v000557 = -(M[6-1][5-1]*v000084);
PetscScalar v000558 = -(M[6-1][5-1]*v000088);
PetscScalar v000559 = -(M[6-1][5-1]*v000089);
PetscScalar v000560 = -(M[6-1][5-1]*v000090);
PetscScalar v000561 = M[2-1][3-1]*M[6-1][5-1]*v000118;
PetscScalar v000562 = M[2-1][3-1]*M[6-1][5-1]*v000119;
PetscScalar v000563 = -(M[6-1][5-1]*v000091);
PetscScalar v000564 = M[6-1][4-1]*v000286;
PetscScalar v000565 = M[6-1][4-1]*v000287;
PetscScalar v000566 = M[2-1][3-1]*M[6-1][4-1]*v000100;
PetscScalar v000567 = M[6-1][4-1]*v000288;
PetscScalar v000568 = M[6-1][4-1]*v000289;
PetscScalar v000569 = M[2-1][3-1]*M[6-1][4-1]*v000101;
PetscScalar v000570 = -(M[2-1][1-1]*M[3-1][5-1]*v000065);
PetscScalar v000571 = M[6-1][4-1]*v000290;
PetscScalar v000572 = M[2-1][3-1]*M[6-1][4-1]*v000102;
PetscScalar v000573 = -(M[6-1][3-1]*v000115);
PetscScalar v000574 = -(M[6-1][3-1]*v000116);
PetscScalar v000575 = -(M[6-1][3-1]*v000120);
PetscScalar v000576 = -(M[2-1][1-1]*M[3-1][5-1]*v000066);
PetscScalar v000577 = -(M[2-1][1-1]*M[3-1][4-1]*v000063);
PetscScalar v000578 = -(M[6-1][3-1]*v000123);
PetscScalar v000579 = M[6-1][2-1]*v000358;
PetscScalar v000580 = -(M[2-1][1-1]*M[3-1][4-1]*v000064);
PetscScalar v000581 = -(M[2-1][1-1]*v000201);
PetscScalar v000582 = M[2-1][3-1]*M[3-1][5-1]*v000176;
PetscScalar v000583 = -(M[2-1][1-1]*M[3-1][5-1]*v000067);
PetscScalar v000584 = -(M[2-1][3-1]*M[3-1][4-1]*v000057);
PetscScalar v000585 = -(M[6-1][1-1]*v000163);
PetscScalar v000586 = -(M[2-1][3-1]*M[3-1][4-1]*v000058);
PetscScalar v000587 = -(M[6-1][1-1]*v000168);
PetscScalar v000588 = -(M[6-1][1-1]*v000169);
PetscScalar v000589 = -(M[6-1][1-1]*v000170);
PetscScalar v000590 = M[2-1][3-1]*M[3-1][5-1]*v000177;
PetscScalar v000591 = v000003 + v000004 + v000005 + v000073 + v000074 + v000075;
PetscScalar v000592 = v000006 + v000007 + v000008 + v000080 + v000081 + v000082;
PetscScalar v000593 = v000009 + v000010 + v000011 + v000085 + v000086 + v000087;
PetscScalar v000594 = -(M[3-1][6-1]*M[4-1][1-1]*M[5-1][2-1]) + v000012 + v000013 + v000014 + v000092 + v000093;
PetscScalar v000595 = v000015 + v000016 + v000017 + v000100 + v000101 + v000102;
PetscScalar v000596 = v000036 + v000037 + v000038 + v000153 + v000154 + v000155;
PetscScalar v000597 = v000042 + v000043 + v000044 + v000165 + v000166 + v000167;
PetscScalar v000598 = v000047 + v000048 + v000049 - v000053 - v000054 + v000174;
PetscScalar v000599 = v000053 + v000054 + v000055 + v000178 + v000179 + v000180;
PetscScalar v000600 = v000083 + v000084 + v000088 + v000089 + v000090 + v000091 + v000217 + v000218 + v000219 + v000220 + v000221 + v000222 + M[2-1][3-1]*v000592 + M[2-1][4-1]*v000593;
PetscScalar v000601 = v000094 + v000095 + v000096 + v000097 + v000098 + v000099 + v000250 + v000251 + v000252 + v000253 + v000254 + v000255 + M[2-1][6-1]*v000591 + M[2-1][3-1]*v000594;
PetscScalar v000602 = M[6-1][3-1]*v000006 + M[6-1][3-1]*v000008 + M[6-1][4-1]*v000009 + M[6-1][4-1]*v000010 + M[6-1][2-1]*v000031 + M[6-1][1-1]*v000036 + M[6-1][3-1]*v000080 + M[6-1][3-1]*v000081 + M[6-1][4-1]*v000085 + M[6-1][4-1]*v000087 + M[6-1][2-1]*v000141 + M[6-1][1-1]*v000154 + M[3-1][3-1]*(v000059 + v000060 + v000061 + v000189 + v000190 + v000191) + M[3-1][4-1]*v000598;
PetscScalar v000603 = -(M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1]) + M[3-1][4-1]*(-(M[4-1][3-1]*M[5-1][2-1]*M[6-1][5-1]) + v000062 + v000063 + v000064 + v000192 + v000193) + v000194 + v000195 + M[3-1][5-1]*(-(M[4-1][2-1]*M[5-1][3-1]*M[6-1][4-1]) + v000065 + v000066 + v000067 + v000196 + v000197) + v000198 + v000199 + v000200 + v000201 + v000460 + v000461 + v000462 + v000463 + v000464;
PetscScalar v000604 = M[2-1][5-1]*v000602;
PetscScalar v000605 = v000208 + v000209 + v000210 + v000211 + v000212 + v000213 + v000214 + v000215 + v000216 + v000223 + v000224 + v000225 + v000226 + v000227 + v000228 + v000229 + v000230 + v000231 + v000232 + v000233 + v000234 + v000235 + v000236 + v000237 + v000238 + v000239 + v000240 + v000241 + v000242 + v000243 + v000244 + v000245 + v000246 + v000247 + v000248 + v000249 + v000465 + v000466 + v000467 + v000468 + v000469 + v000470 + v000471 + v000472 + v000473 + v000474 + v000475 + v000476 + v000477 + v000478 + v000479 + v000480 + v000481 + v000482 + v000483 + v000484 + v000485 + v000486 + v000487 + v000488 + v000489 + v000490 + v000491 + v000492 + v000493 + v000494 + v000495 + v000496 + v000497 + v000498 + v000499 + v000500 + M[1-1][4-1]*(M[2-1][3-1]*(v000000 + v000001 + v000002 + v000068 + v000069 + v000070) + v000071 + v000072 + v000076 + v000077 + v000078 + v000079 + v000202 + v000203 + v000204 + v000205 + v000206 + v000207 + M[2-1][5-1]*v000591) + M[1-1][5-1]*v000600;
PetscScalar v000606 = 1/(M[6-1][2-1]*(M[1-1][4-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][1-1] - M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][1-1] + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1] + M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1] - M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1] - M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][1-1] - M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][3-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1] - M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1] - M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1] - M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1] + M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1] - M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][4-1] - M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1] - M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1] - M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1] + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1] - M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1] - M[1-1][1-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1] + M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1] - M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1] - M[1-1][4-1]*M[2-1][5-1]*v000022 - M[1-1][4-1]*M[2-1][6-1]*v000027 + M[1-1][5-1]*(M[2-1][4-1]*(v000021 + v000022 + v000023 + v000124 + v000125 + v000126) + v000127 + v000128 + M[2-1][6-1]*(v000024 + v000025 + v000026 + v000129 + v000130 + v000131) + v000132 + v000133 + v000134 + v000135 + v000345 + v000346 + v000347 + v000348 + v000349 + v000350) + v000351 + v000352 + v000353 + v000354 + v000355 + v000356 + v000357 + M[1-1][6-1]*(-(M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1]) - M[2-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1] - M[2-1][3-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1] - M[2-1][1-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1] - M[2-1][1-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1] + M[2-1][4-1]*(v000027 + v000028 + v000029 + v000136 + v000137 + v000138) + v000139 + v000140 + M[2-1][5-1]*(v000030 + v000031 + v000032 + v000141 + v000142 + v000143) + v000144 + v000145 + v000146 + v000147 + v000358) + v000359 + v000360 + v000361 + v000362 + v000363 + v000364 + v000365 + v000366 + v000367 + v000368 + v000369 + v000370 + v000371 + v000372 + v000373 + v000374 + v000375 + v000535 + v000536 + v000537 + v000538 + v000539 + v000540 + v000541 + v000542 + v000543 + v000544 + v000545 + v000546 + v000547 + v000548 + v000549 + v000550 + v000551 + v000552 + v000553 + v000554) - M[6-1][3-1]*(-(M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][1-1]) - M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][1-1] - M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][1-1] - M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][2-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1] - M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1] - M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1] - M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1] - M[1-1][1-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1] - M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][4-1] - M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][4-1] + M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1] - M[1-1][2-1]*M[2-1][6-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1] - M[1-1][1-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1] - M[1-1][1-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1] - M[1-1][2-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1] - M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1] - M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1] + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1] - M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1] - M[1-1][1-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1] - M[1-1][2-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1] - M[1-1][2-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1] - M[1-1][4-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1] - M[1-1][1-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1] - M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][2-1]*M[5-1][6-1] - M[1-1][1-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1] - M[1-1][4-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1] - M[1-1][2-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][6-1] - M[1-1][1-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1] - M[1-1][4-1]*M[2-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1] - M[1-1][1-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1] - M[1-1][2-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][6-1]*v000000 + M[1-1][4-1]*M[2-1][6-1]*v000001 + M[1-1][4-1]*M[2-1][6-1]*v000002 - M[1-1][4-1]*M[2-1][5-1]*v000012 - M[1-1][4-1]*M[2-1][5-1]*v000013 - M[1-1][4-1]*M[2-1][5-1]*v000014 + M[1-1][4-1]*M[2-1][6-1]*v000068 + M[1-1][4-1]*M[2-1][6-1]*v000069 + M[1-1][4-1]*M[2-1][6-1]*v000070 + v000314 + v000315 + v000316 + v000317 + v000318 + v000319 + v000320 + v000321 + v000328 + v000329 + v000330 + v000331 + v000332 + v000333 + v000334 + v000335 + v000336 + v000337 + v000338 + v000339 + v000340 + v000341 + v000342 + v000343 + v000344 + M[1-1][5-1]*(-(M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1]) - M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1] - M[2-1][1-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1] - M[2-1][2-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1] - M[2-1][2-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1] - M[2-1][1-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1] + v000109 + v000110 + v000111 + v000112 + v000113 + v000114 + M[2-1][6-1]*v000592 + M[2-1][4-1]*v000594) + M[1-1][6-1]*(v000115 + v000116 + M[2-1][5-1]*(v000018 + v000019 + v000020 + v000117 + v000118 + v000119) + v000120 + v000121 + v000122 + v000123 + v000322 + v000323 + v000324 + v000325 + v000326 + v000327 + M[2-1][4-1]*v000595)) - M[6-1][1-1]*(M[1-1][4-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][2-1] - M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1] - M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1] - M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1] - M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1] - M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1] + M[1-1][4-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1] - M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1] - M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1] - M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1] - M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1] - M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1] - M[1-1][3-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1] - M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1] - M[1-1][2-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1] - M[1-1][3-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1] + M[1-1][4-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][2-1]*M[5-1][5-1] - M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1] - M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1] - M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1] - M[1-1][2-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1] - M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1] - M[1-1][3-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1] - M[1-1][4-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1] - M[1-1][2-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1] - M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1] - M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1] - M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1] - M[1-1][3-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1] - M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1] - M[1-1][4-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1] - M[1-1][2-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1] - M[1-1][3-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1] - M[1-1][4-1]*M[2-1][5-1]*v000033 - M[1-1][4-1]*M[2-1][5-1]*v000034 - M[1-1][4-1]*M[2-1][5-1]*v000035 - M[1-1][4-1]*M[2-1][6-1]*v000039 - M[1-1][4-1]*M[2-1][6-1]*v000040 - M[1-1][4-1]*M[2-1][6-1]*v000041 + v000382 + v000383 + v000384 + v000385 + v000386 + v000387 + v000388 + v000389 + v000396 + v000397 + v000398 + v000399 + v000400 + v000401 + v000402 + v000403 + v000404 + v000405 + v000406 + v000407 + v000408 + v000409 + v000410 + v000411 + v000412 + v000413 + v000414 + v000415 + v000416 + v000417 + M[1-1][5-1]*(M[2-1][4-1]*(v000033 + v000034 + v000035 + v000148 + v000149 + v000150) + v000151 + v000152 + v000156 + v000157 + v000158 + v000159 + v000376 + v000377 + v000378 + v000379 + v000380 + v000381 + M[2-1][6-1]*v000596) + M[1-1][6-1]*(M[2-1][4-1]*(v000039 + v000040 + v000041 + v000160 + v000161 + v000162) + v000163 + v000164 + v000168 + v000169 + v000170 + v000171 + v000390 + v000391 + v000392 + v000393 + v000394 + v000395 + M[2-1][5-1]*v000597)) - M[6-1][5-1]*(M[1-1][3-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][1-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][2-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][2-1]*M[5-1][6-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1] - M[1-1][3-1]*M[2-1][4-1]*v000012 - M[1-1][3-1]*M[2-1][4-1]*v000013 - M[1-1][3-1]*M[2-1][4-1]*v000014 + M[1-1][3-1]*M[2-1][6-1]*v000018 + M[1-1][3-1]*M[2-1][6-1]*v000019 + M[1-1][3-1]*M[2-1][6-1]*v000020 + M[1-1][2-1]*M[2-1][6-1]*v000024 + M[1-1][2-1]*M[2-1][6-1]*v000025 + M[1-1][2-1]*M[2-1][6-1]*v000026 - M[1-1][1-1]*M[2-1][4-1]*v000033 - M[1-1][1-1]*M[2-1][4-1]*v000034 - M[1-1][1-1]*M[2-1][4-1]*v000035 + M[1-1][1-1]*M[2-1][6-1]*v000042 + M[1-1][1-1]*M[2-1][6-1]*v000043 + M[1-1][1-1]*M[2-1][6-1]*v000044 - M[1-1][3-1]*v000109 - M[1-1][3-1]*v000110 - M[1-1][3-1]*v000111 - M[1-1][3-1]*v000112 - M[1-1][3-1]*v000113 - M[1-1][3-1]*v000114 + M[1-1][3-1]*M[2-1][6-1]*v000117 + M[1-1][3-1]*M[2-1][6-1]*v000118 + M[1-1][3-1]*M[2-1][6-1]*v000119 + M[1-1][2-1]*M[2-1][4-1]*v000124 + M[1-1][2-1]*M[2-1][4-1]*v000125 + M[1-1][2-1]*M[2-1][4-1]*v000126 + M[1-1][2-1]*M[2-1][6-1]*v000129 + M[1-1][2-1]*M[2-1][6-1]*v000130 + M[1-1][2-1]*M[2-1][6-1]*v000131 - M[1-1][1-1]*v000151 - M[1-1][1-1]*v000152 - M[1-1][1-1]*v000156 - M[1-1][1-1]*v000157 - M[1-1][1-1]*v000158 - M[1-1][1-1]*v000159 + M[1-1][1-1]*M[2-1][6-1]*v000165 + M[1-1][1-1]*M[2-1][6-1]*v000166 + M[1-1][1-1]*M[2-1][6-1]*v000167 + v000256 + v000257 + v000258 + v000259 + v000260 + v000261 + v000262 + v000263 + v000264 + v000265 + v000266 + v000267 + v000268 + v000269 + v000270 + v000271 + v000272 + v000273 + v000274 + v000275 + v000276 + M[1-1][2-1]*v000345 + M[1-1][2-1]*v000346 + M[1-1][2-1]*v000347 + M[1-1][2-1]*v000348 + M[1-1][2-1]*v000349 + M[1-1][2-1]*v000350 + M[1-1][6-1]*v000600 + M[1-1][4-1]*v000601) + M[6-1][4-1]*(M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][2-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][2-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1] + M[1-1][1-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][2-1]*M[5-1][5-1] - M[1-1][1-1]*M[2-1][6-1]*v000039 - M[1-1][1-1]*M[2-1][6-1]*v000041 + v000277 + v000278 + v000279 + v000280 + v000281 + v000282 + v000283 + v000284 + v000285 + v000291 + v000292 + v000293 + v000294 + v000295 + v000296 + v000297 + v000298 + v000299 + v000300 + v000301 + v000302 + v000303 + v000304 + v000305 + v000306 + v000307 + v000308 + v000309 + v000310 + v000311 + v000312 + v000313 + v000501 + v000502 + v000503 + v000504 + v000505 + v000506 + v000507 + v000508 + v000509 + v000510 + v000511 + v000512 + v000513 + v000514 + v000515 + v000516 + v000517 + v000518 + v000519 + v000520 + v000521 + v000522 + v000523 + v000524 + v000525 + v000526 + v000527 + v000528 + v000529 + v000530 + v000531 + v000532 + v000533 + v000534 + M[1-1][6-1]*(-v000078 + v000103 + v000104 + v000105 + v000106 + v000107 + v000108 + v000286 + v000287 + v000288 + v000289 + v000290 + M[2-1][5-1]*v000593 + M[2-1][3-1]*v000595) + M[1-1][5-1]*v000601) + M[6-1][6-1]*v000605);
PetscScalar v000607 = 1/(M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][1-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][3-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][1-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][1-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][1-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][1-1]*M[6-1][3-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][1-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][4-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][4-1]*M[2-1][1-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][1-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1]*M[6-1][5-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][1-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][1-1]*M[6-1][5-1] + M[1-1][4-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][1-1]*M[5-1][3-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][1-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][1-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][1-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][4-1]*M[2-1][1-1]*M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][2-1]*M[2-1][1-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][4-1]*M[2-1][3-1]*M[6-1][6-1]*v000000 + M[1-1][4-1]*M[2-1][3-1]*M[6-1][6-1]*v000001 + M[1-1][4-1]*M[2-1][3-1]*M[6-1][6-1]*v000002 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][6-1]*v000003 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][6-1]*v000004 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][6-1]*v000005 + M[1-1][3-1]*M[2-1][6-1]*M[6-1][5-1]*v000006 + M[1-1][3-1]*M[2-1][6-1]*M[6-1][5-1]*v000007 + M[1-1][3-1]*M[2-1][6-1]*M[6-1][5-1]*v000008 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][5-1]*v000009 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][5-1]*v000010 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][5-1]*v000011 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][3-1]*v000012 - M[1-1][4-1]*M[2-1][3-1]*M[6-1][5-1]*v000012 + M[1-1][3-1]*M[2-1][4-1]*M[6-1][5-1]*v000012 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][3-1]*v000013 - M[1-1][4-1]*M[2-1][3-1]*M[6-1][5-1]*v000013 + M[1-1][3-1]*M[2-1][4-1]*M[6-1][5-1]*v000013 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][3-1]*v000014 - M[1-1][4-1]*M[2-1][3-1]*M[6-1][5-1]*v000014 + M[1-1][3-1]*M[2-1][4-1]*M[6-1][5-1]*v000014 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][3-1]*v000015 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][3-1]*v000016 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][3-1]*v000017 + M[1-1][4-1]*M[2-1][2-1]*M[6-1][5-1]*v000021 + M[1-1][4-1]*M[2-1][2-1]*M[6-1][5-1]*v000022 + M[1-1][4-1]*M[2-1][2-1]*M[6-1][5-1]*v000023 + M[1-1][2-1]*M[2-1][6-1]*M[6-1][5-1]*v000030 + M[1-1][2-1]*M[2-1][6-1]*M[6-1][5-1]*v000031 + M[1-1][2-1]*M[2-1][6-1]*M[6-1][5-1]*v000032 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][1-1]*v000033 + M[1-1][1-1]*M[2-1][4-1]*M[6-1][5-1]*v000033 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][1-1]*v000034 + M[1-1][1-1]*M[2-1][4-1]*M[6-1][5-1]*v000034 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][1-1]*v000035 + M[1-1][1-1]*M[2-1][4-1]*M[6-1][5-1]*v000035 + M[1-1][1-1]*M[2-1][6-1]*M[6-1][5-1]*v000037 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][1-1]*v000039 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][1-1]*v000040 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][1-1]*v000041 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*v000045 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*v000045 - M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*v000045 - M[1-1][4-1]*M[2-1][5-1]*M[3-1][3-1]*v000046 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*v000046 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*v000046 - M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*v000046 + M[1-1][4-1]*M[2-1][6-1]*M[3-1][5-1]*v000048 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*v000050 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*v000051 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*v000052 + M[1-1][4-1]*M[2-1][5-1]*M[3-1][6-1]*v000054 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*v000056 - M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*v000056 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*v000057 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*v000057 - M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*v000057 - M[1-1][4-1]*M[2-1][6-1]*M[3-1][3-1]*v000058 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*v000058 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*v000058 - M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*v000058 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*v000060 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*v000061 + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*v000062 - M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*v000062 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*v000063 + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*v000063 - M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*v000063 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*v000064 + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*v000064 - M[1-1][1-1]*M[2-1][4-1]*M[3-1][6-1]*v000064 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*v000065 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*v000066 - M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*v000066 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*v000067 - M[1-1][1-1]*M[2-1][5-1]*M[3-1][6-1]*v000067 + M[1-1][4-1]*M[2-1][3-1]*M[6-1][6-1]*v000068 + M[1-1][4-1]*M[2-1][3-1]*M[6-1][6-1]*v000069 + M[1-1][4-1]*M[2-1][3-1]*M[6-1][6-1]*v000070 + M[1-1][4-1]*M[6-1][6-1]*v000071 + M[1-1][4-1]*M[6-1][6-1]*v000072 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][6-1]*v000073 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][6-1]*v000074 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][6-1]*v000075 + M[1-1][4-1]*M[6-1][6-1]*v000076 + M[1-1][4-1]*M[6-1][6-1]*v000077 + M[1-1][4-1]*M[6-1][6-1]*v000078 + M[1-1][4-1]*M[6-1][6-1]*v000079 + M[1-1][3-1]*M[2-1][6-1]*M[6-1][5-1]*v000080 + M[1-1][3-1]*M[2-1][6-1]*M[6-1][5-1]*v000081 + M[1-1][3-1]*M[2-1][6-1]*M[6-1][5-1]*v000082 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][5-1]*v000085 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][5-1]*v000086 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][5-1]*v000087 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][3-1]*v000092 + M[1-1][3-1]*M[2-1][4-1]*M[6-1][5-1]*v000092 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][3-1]*v000093 + M[1-1][3-1]*M[2-1][4-1]*M[6-1][5-1]*v000093 - M[1-1][4-1]*M[6-1][5-1]*v000094 - M[1-1][4-1]*M[6-1][5-1]*v000095 - M[1-1][4-1]*M[6-1][5-1]*v000096 - M[1-1][4-1]*M[6-1][5-1]*v000097 - M[1-1][4-1]*M[6-1][5-1]*v000098 - M[1-1][4-1]*M[6-1][5-1]*v000099 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][3-1]*v000100 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][3-1]*v000101 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][3-1]*v000102 + M[1-1][3-1]*M[6-1][5-1]*v000109 + M[1-1][3-1]*M[6-1][5-1]*v000110 + M[1-1][3-1]*M[6-1][5-1]*v000111 + M[1-1][3-1]*M[6-1][5-1]*v000112 + M[1-1][3-1]*M[6-1][5-1]*v000113 + M[1-1][3-1]*M[6-1][5-1]*v000114 + M[1-1][2-1]*M[2-1][6-1]*M[6-1][5-1]*v000141 + M[1-1][2-1]*M[2-1][6-1]*M[6-1][5-1]*v000142 + M[1-1][2-1]*M[2-1][6-1]*M[6-1][5-1]*v000143 + M[1-1][4-1]*M[2-1][5-1]*M[6-1][1-1]*v000148 + M[1-1][1-1]*M[2-1][4-1]*M[6-1][5-1]*v000148 + M[1-1][1-1]*M[2-1][4-1]*M[6-1][5-1]*v000150 + M[1-1][1-1]*M[6-1][5-1]*v000151 + M[1-1][1-1]*M[6-1][5-1]*v000152 + M[1-1][1-1]*M[2-1][6-1]*M[6-1][5-1]*v000155 + M[1-1][1-1]*M[6-1][5-1]*v000156 + M[1-1][1-1]*M[6-1][5-1]*v000157 + M[1-1][1-1]*M[6-1][5-1]*v000158 + M[1-1][1-1]*M[6-1][5-1]*v000159 + M[1-1][4-1]*M[2-1][6-1]*M[6-1][1-1]*v000161 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*v000172 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*v000172 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*v000173 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*v000173 + M[1-1][4-1]*M[2-1][6-1]*M[3-1][5-1]*v000174 - M[1-1][4-1]*M[2-1][5-1]*v000175 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*v000176 + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*v000177 + M[1-1][4-1]*M[2-1][5-1]*M[3-1][6-1]*v000178 + M[1-1][4-1]*M[2-1][5-1]*M[3-1][6-1]*v000180 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*v000181 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*v000181 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*v000182 + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*v000182 - M[1-1][4-1]*M[2-1][6-1]*v000187 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*v000190 + M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*v000191 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*v000192 + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*v000192 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][4-1]*v000193 + M[1-1][4-1]*M[2-1][1-1]*M[3-1][6-1]*v000193 + M[1-1][1-1]*M[2-1][6-1]*v000194 + M[1-1][1-1]*M[2-1][6-1]*v000195 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*v000196 + M[1-1][1-1]*M[2-1][6-1]*M[3-1][5-1]*v000197 + M[1-1][1-1]*M[2-1][6-1]*v000198 + M[1-1][1-1]*M[2-1][6-1]*v000199 + M[1-1][1-1]*M[2-1][6-1]*v000200 + M[1-1][1-1]*M[2-1][6-1]*v000201 + M[1-1][4-1]*M[6-1][6-1]*v000202 + M[1-1][4-1]*M[6-1][6-1]*v000203 + M[1-1][4-1]*M[6-1][6-1]*v000204 + M[1-1][4-1]*M[6-1][6-1]*v000205 + M[1-1][4-1]*M[6-1][6-1]*v000206 + M[1-1][4-1]*M[6-1][6-1]*v000207 + M[6-1][6-1]*v000208 + M[6-1][6-1]*v000209 + M[6-1][6-1]*v000210 + M[6-1][6-1]*v000211 + M[6-1][6-1]*v000212 + M[6-1][6-1]*v000213 + M[6-1][6-1]*v000214 + M[6-1][6-1]*v000215 + M[6-1][6-1]*v000216 + M[6-1][6-1]*v000223 + M[6-1][6-1]*v000224 + M[6-1][6-1]*v000225 + M[6-1][6-1]*v000226 + M[6-1][6-1]*v000227 + M[6-1][6-1]*v000228 + M[6-1][6-1]*v000229 + M[6-1][6-1]*v000230 + M[6-1][6-1]*v000231 + M[6-1][6-1]*v000232 + M[6-1][6-1]*v000233 + M[6-1][6-1]*v000234 + M[6-1][6-1]*v000235 + M[6-1][6-1]*v000236 + M[6-1][6-1]*v000237 + M[6-1][6-1]*v000238 + M[6-1][6-1]*v000239 + M[6-1][6-1]*v000240 + M[6-1][6-1]*v000241 + M[6-1][6-1]*v000242 + M[6-1][6-1]*v000243 + M[6-1][6-1]*v000244 + M[6-1][6-1]*v000245 + M[6-1][6-1]*v000246 + M[6-1][6-1]*v000247 + M[6-1][6-1]*v000248 + M[6-1][6-1]*v000249 - M[6-1][5-1]*v000256 - M[6-1][5-1]*v000257 - M[6-1][5-1]*v000258 - M[6-1][5-1]*v000259 - M[6-1][5-1]*v000260 - M[6-1][5-1]*v000261 - M[6-1][5-1]*v000262 - M[6-1][5-1]*v000263 - M[6-1][5-1]*v000264 - M[6-1][5-1]*v000265 - M[6-1][5-1]*v000266 - M[6-1][5-1]*v000267 - M[6-1][5-1]*v000268 - M[6-1][5-1]*v000269 - M[6-1][5-1]*v000270 - M[6-1][5-1]*v000271 - M[6-1][5-1]*v000272 - M[6-1][5-1]*v000273 - M[6-1][5-1]*v000274 - M[6-1][5-1]*v000275 - M[6-1][5-1]*v000276 + M[6-1][4-1]*v000277 + M[6-1][4-1]*v000278 + M[6-1][4-1]*v000279 + M[6-1][4-1]*v000280 + M[6-1][4-1]*v000281 + M[6-1][4-1]*v000282 + M[6-1][4-1]*v000283 + M[6-1][4-1]*v000284 + M[6-1][4-1]*v000285 + M[6-1][4-1]*v000291 + M[6-1][4-1]*v000292 + M[6-1][4-1]*v000293 + M[6-1][4-1]*v000294 + M[6-1][4-1]*v000295 + M[6-1][4-1]*v000296 + M[6-1][4-1]*v000297 + M[6-1][4-1]*v000298 + M[6-1][4-1]*v000299 + M[6-1][4-1]*v000300 + M[6-1][4-1]*v000301 + M[6-1][4-1]*v000302 + M[6-1][4-1]*v000303 + M[6-1][4-1]*v000304 + M[6-1][4-1]*v000305 + M[6-1][4-1]*v000306 + M[6-1][4-1]*v000307 + M[6-1][4-1]*v000308 + M[6-1][4-1]*v000309 + M[6-1][4-1]*v000310 + M[6-1][4-1]*v000311 + M[6-1][4-1]*v000312 + M[6-1][4-1]*v000313 - M[6-1][3-1]*v000314 - M[6-1][3-1]*v000315 - M[6-1][3-1]*v000316 - M[6-1][3-1]*v000317 - M[6-1][3-1]*v000318 - M[6-1][3-1]*v000319 - M[6-1][3-1]*v000320 - M[6-1][3-1]*v000321 - M[6-1][3-1]*v000328 - M[6-1][3-1]*v000329 - M[6-1][3-1]*v000330 - M[6-1][3-1]*v000331 - M[6-1][3-1]*v000332 - M[6-1][3-1]*v000333 - M[6-1][3-1]*v000334 - M[6-1][3-1]*v000335 - M[6-1][3-1]*v000336 - M[6-1][3-1]*v000337 - M[6-1][3-1]*v000338 - M[6-1][3-1]*v000339 - M[6-1][3-1]*v000340 - M[6-1][3-1]*v000341 - M[6-1][3-1]*v000342 - M[6-1][3-1]*v000343 - M[6-1][3-1]*v000344 + M[6-1][2-1]*v000351 + M[6-1][2-1]*v000352 + M[6-1][2-1]*v000353 + M[6-1][2-1]*v000354 + M[6-1][2-1]*v000355 + M[6-1][2-1]*v000356 + M[6-1][2-1]*v000357 + M[6-1][2-1]*v000359 + M[6-1][2-1]*v000360 + M[6-1][2-1]*v000361 + M[6-1][2-1]*v000362 + M[6-1][2-1]*v000363 + M[6-1][2-1]*v000364 + M[6-1][2-1]*v000365 + M[6-1][2-1]*v000366 + M[6-1][2-1]*v000367 + M[6-1][2-1]*v000368 + M[6-1][2-1]*v000369 + M[6-1][2-1]*v000370 + M[6-1][2-1]*v000371 + M[6-1][2-1]*v000372 + M[6-1][2-1]*v000373 + M[6-1][2-1]*v000374 + M[6-1][2-1]*v000375 - M[6-1][1-1]*v000382 - M[6-1][1-1]*v000383 - M[6-1][1-1]*v000384 - M[6-1][1-1]*v000385 - M[6-1][1-1]*v000386 - M[6-1][1-1]*v000387 - M[6-1][1-1]*v000396 - M[6-1][1-1]*v000397 - M[6-1][1-1]*v000398 - M[6-1][1-1]*v000399 - M[6-1][1-1]*v000400 - M[6-1][1-1]*v000401 - M[6-1][1-1]*v000402 - M[6-1][1-1]*v000404 - M[6-1][1-1]*v000405 - M[6-1][1-1]*v000406 - M[6-1][1-1]*v000407 - M[6-1][1-1]*v000408 - M[6-1][1-1]*v000410 - M[6-1][1-1]*v000411 - M[6-1][1-1]*v000412 - M[6-1][1-1]*v000413 + M[1-1][1-1]*M[2-1][6-1]*v000460 + M[1-1][1-1]*M[2-1][6-1]*v000461 + M[1-1][1-1]*M[2-1][6-1]*v000462 + M[1-1][1-1]*M[2-1][6-1]*v000463 + M[1-1][1-1]*M[2-1][6-1]*v000464 + M[6-1][6-1]*v000465 + M[6-1][6-1]*v000466 + M[6-1][6-1]*v000467 + M[6-1][6-1]*v000468 + M[6-1][6-1]*v000469 + M[6-1][6-1]*v000470 + M[6-1][6-1]*v000471 + M[6-1][6-1]*v000472 + M[6-1][6-1]*v000473 + M[6-1][6-1]*v000474 + M[6-1][6-1]*v000475 + M[6-1][6-1]*v000476 + M[6-1][6-1]*v000477 + M[6-1][6-1]*v000478 + M[6-1][6-1]*v000479 + M[6-1][6-1]*v000480 + M[6-1][6-1]*v000481 + M[6-1][6-1]*v000482 + M[6-1][6-1]*v000483 + M[6-1][6-1]*v000484 + M[6-1][6-1]*v000485 + M[6-1][6-1]*v000486 + M[6-1][6-1]*v000487 + M[6-1][6-1]*v000488 + M[6-1][6-1]*v000489 + M[6-1][6-1]*v000490 + M[6-1][6-1]*v000491 + M[6-1][6-1]*v000492 + M[6-1][6-1]*v000493 + M[6-1][6-1]*v000494 + M[6-1][6-1]*v000495 + M[6-1][6-1]*v000496 + M[6-1][6-1]*v000497 + M[6-1][6-1]*v000498 + M[6-1][6-1]*v000499 + M[6-1][6-1]*v000500 + M[6-1][4-1]*v000501 + M[6-1][4-1]*v000502 + M[6-1][4-1]*v000503 + M[6-1][4-1]*v000504 + M[6-1][4-1]*v000505 + M[6-1][4-1]*v000506 + M[6-1][4-1]*v000507 + M[6-1][4-1]*v000508 + M[6-1][4-1]*v000509 + M[6-1][4-1]*v000510 + M[6-1][4-1]*v000511 + M[6-1][4-1]*v000512 + M[6-1][4-1]*v000513 + M[6-1][4-1]*v000514 + M[6-1][4-1]*v000515 + M[6-1][4-1]*v000516 + M[6-1][4-1]*v000517 + M[6-1][4-1]*v000518 + M[6-1][4-1]*v000519 + M[6-1][4-1]*v000520 + M[6-1][4-1]*v000521 + M[6-1][4-1]*v000522 + M[6-1][4-1]*v000523 + M[6-1][4-1]*v000524 + M[6-1][4-1]*v000525 + M[6-1][4-1]*v000526 + M[6-1][4-1]*v000527 + M[6-1][4-1]*v000528 + M[6-1][4-1]*v000529 + M[6-1][4-1]*v000530 + M[6-1][4-1]*v000531 + M[6-1][4-1]*v000532 + M[6-1][4-1]*v000533 + M[6-1][4-1]*v000534 + M[6-1][2-1]*v000535 + M[6-1][2-1]*v000536 + M[6-1][2-1]*v000537 + M[6-1][2-1]*v000538 + M[6-1][2-1]*v000539 + M[6-1][2-1]*v000540 + M[6-1][2-1]*v000541 + M[6-1][2-1]*v000542 + M[6-1][2-1]*v000543 + M[6-1][2-1]*v000544 + M[6-1][2-1]*v000545 + M[6-1][2-1]*v000546 + M[6-1][2-1]*v000547 + M[6-1][2-1]*v000548 + M[6-1][2-1]*v000549 + M[6-1][2-1]*v000550 + M[6-1][2-1]*v000551 + M[6-1][2-1]*v000552 + M[6-1][2-1]*v000553 + M[6-1][2-1]*v000554 + M[1-1][5-1]*(M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][1-1] + M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][1-1] + M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][1-1] + M[2-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][1-1] + M[2-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][1-1] + M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][1-1]*M[6-1][3-1] + M[2-1][1-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][3-1] + M[2-1][2-1]*M[3-1][1-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] + M[2-1][2-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][6-1]*M[6-1][3-1] + M[2-1][1-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] + M[2-1][3-1]*M[6-1][6-1]*v000006 + M[2-1][3-1]*M[6-1][6-1]*v000007 + M[2-1][3-1]*M[6-1][6-1]*v000008 + M[2-1][3-1]*M[6-1][4-1]*v000012 + M[2-1][3-1]*M[6-1][4-1]*v000013 + M[2-1][3-1]*M[6-1][4-1]*v000014 - M[2-1][3-1]*M[3-1][4-1]*v000045 - M[2-1][3-1]*M[3-1][4-1]*v000046 + M[2-1][3-1]*M[3-1][6-1]*v000060 + M[2-1][3-1]*M[3-1][6-1]*v000061 + M[2-1][1-1]*M[3-1][6-1]*v000066 + M[2-1][1-1]*M[3-1][6-1]*v000067 + M[2-1][3-1]*M[6-1][6-1]*v000080 + M[2-1][3-1]*M[6-1][6-1]*v000081 + M[2-1][3-1]*M[6-1][6-1]*v000082 + M[6-1][6-1]*v000083 + M[6-1][6-1]*v000084 + M[6-1][6-1]*v000088 + M[6-1][6-1]*v000089 + M[6-1][6-1]*v000090 + M[6-1][6-1]*v000091 + M[2-1][3-1]*M[6-1][4-1]*v000092 + M[2-1][3-1]*M[6-1][4-1]*v000093 + M[6-1][4-1]*v000094 + M[6-1][4-1]*v000095 + M[6-1][4-1]*v000096 + M[6-1][4-1]*v000097 + M[6-1][4-1]*v000098 + M[6-1][4-1]*v000099 - M[6-1][3-1]*v000109 - M[6-1][3-1]*v000110 - M[6-1][3-1]*v000111 - M[6-1][3-1]*v000112 - M[6-1][3-1]*v000113 - M[6-1][3-1]*v000114 + M[6-1][2-1]*v000127 + M[6-1][2-1]*v000128 + M[6-1][2-1]*v000132 + M[6-1][2-1]*v000135 - M[6-1][1-1]*v000151 - M[6-1][1-1]*v000152 - M[6-1][1-1]*v000156 - M[6-1][1-1]*v000158 + M[2-1][3-1]*M[3-1][6-1]*v000189 + M[2-1][3-1]*M[3-1][6-1]*v000190 + M[2-1][3-1]*M[3-1][6-1]*v000191 + M[2-1][1-1]*M[3-1][6-1]*v000197 + M[6-1][6-1]*v000217 + M[6-1][6-1]*v000218 + M[6-1][6-1]*v000219 + M[6-1][6-1]*v000220 + M[6-1][6-1]*v000221 + M[6-1][6-1]*v000222 + M[6-1][4-1]*v000250 + M[6-1][4-1]*v000251 + M[6-1][4-1]*v000252 + M[6-1][4-1]*v000253 + M[6-1][4-1]*v000254 + M[6-1][4-1]*v000255 + M[6-1][2-1]*v000345 + M[6-1][2-1]*v000347 + M[6-1][2-1]*v000349 + M[2-1][4-1]*(M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][1-1] + M[3-1][2-1]*M[4-1][6-1]*M[5-1][1-1]*M[6-1][3-1] + M[3-1][1-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] + M[6-1][6-1]*v000009 + M[6-1][6-1]*v000010 - M[6-1][3-1]*v000012 - M[6-1][3-1]*v000013 - M[6-1][1-1]*v000034 + M[6-1][6-1]*v000085 + M[6-1][6-1]*v000087 + M[6-1][2-1]*v000124 + M[3-1][3-1]*(M[4-1][2-1]*M[5-1][1-1]*M[6-1][6-1] - M[4-1][1-1]*M[5-1][2-1]*M[6-1][6-1] + v000045 + v000046 + v000172 + v000173) + v000175 + M[3-1][6-1]*v000598) + M[2-1][6-1]*(M[6-1][4-1]*v000003 + M[6-1][4-1]*v000005 + M[6-1][3-1]*v000018 + M[6-1][3-1]*v000019 + M[6-1][2-1]*v000024 + M[6-1][1-1]*v000043 + M[6-1][4-1]*v000073 + M[6-1][4-1]*v000074 + M[6-1][3-1]*v000117 + M[6-1][3-1]*v000119 + M[6-1][2-1]*v000130 + M[6-1][1-1]*v000165 + M[3-1][3-1]*(v000050 + v000051 + v000052 - v000059 + v000176 + v000177) + M[3-1][4-1]*v000599)) + M[1-1][6-1]*(v000424 + v000425 + v000426 + v000427 + v000428 + v000429 + v000430 + v000431 + v000432 + v000433 + v000434 + v000435 + v000436 + v000437 + v000438 + v000439 + v000440 + v000441 + v000442 + v000443 + v000444 + v000445 + v000446 + v000447 + v000448 + v000449 + v000450 + v000451 + v000452 + v000453 + v000454 + v000455 + v000456 + v000457 + v000458 + v000459 + v000555 + v000556 + v000557 + v000558 + v000559 + v000560 + v000561 + v000562 + v000563 + v000564 + v000565 + v000566 + v000567 + v000568 + v000569 + v000570 + v000571 + v000572 + v000573 + v000574 + v000575 + v000576 + v000577 + v000578 + v000579 + v000580 + v000581 + v000582 + v000583 + v000584 + v000585 + v000586 + v000587 + v000588 + v000589 + v000590 + M[2-1][4-1]*(M[3-1][3-1]*(-(M[4-1][2-1]*M[5-1][1-1]*M[6-1][5-1]) + v000056 + v000057 + v000058 + v000181 + v000182) + v000183 + v000184 + v000185 + v000186 + v000187 + v000188 + v000418 + v000419 + v000420 + v000421 + v000422 + v000423 + M[3-1][5-1]*v000599) + v000604));
N[0][0] = (M[1-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][2-1] + M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][2-1] + M[1-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][2-1] + M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][2-1] - M[1-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][2-1] - M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][2-1] + M[1-1][5-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][2-1] + M[2-1][5-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][2-1] - M[1-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][2-1] - M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][2-1] - M[1-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][2-1] - M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][2-1] + M[1-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][2-1] + M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][2-1] + M[1-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][2-1] + M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][2-1] + M[1-1][4-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][2-1] + M[2-1][4-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][2-1] - M[1-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][2-1] - M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][2-1] - M[1-1][5-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][2-1] - M[2-1][5-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][2-1] + M[1-1][4-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][2-1] + M[2-1][4-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][2-1] + M[1-1][5-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][2-1] + M[2-1][5-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][2-1] - M[1-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][2-1] - M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][2-1] - M[1-1][4-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][2-1] - M[2-1][4-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][2-1] + M[1-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][2-1] + M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][2-1] + M[1-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][3-1] + M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][3-1] - M[1-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][3-1] - M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][3-1] + M[1-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][3-1] + M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1] + M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] + M[2-1][5-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] - M[1-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] - M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][4-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][3-1] + M[2-1][4-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][3-1] - M[1-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][3-1] - M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][3-1] - M[1-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] - M[2-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] + M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][5-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] + M[2-1][5-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] - M[1-1][4-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] - M[2-1][4-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][3-1] - M[1-1][5-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] - M[2-1][5-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] + M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][4-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] + M[2-1][4-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] - M[1-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] - M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][3-1] + M[1-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1]*M[6-1][4-1] + M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1]*M[6-1][4-1] + M[1-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][4-1] + M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][4-1] - M[1-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][4-1] - M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][4-1] + M[1-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1]*M[6-1][4-1] + M[2-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1]*M[6-1][4-1] - M[1-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][4-1] - M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][4-1] + M[1-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][4-1] + M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][4-1] - M[1-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][4-1] - M[2-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][4-1] + M[1-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1]*M[6-1][4-1] + M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1]*M[6-1][4-1] + M[1-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][4-1] + M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][4-1] - M[1-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][4-1] - M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1]*M[6-1][4-1] + M[1-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][4-1] + M[2-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][4-1] + M[1-1][5-1]*M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][4-1] + M[2-1][5-1]*M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][4-1] - M[1-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][4-1] - M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][4-1] - M[1-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][4-1] - M[2-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][4-1] + M[1-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][4-1] + M[2-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1]*M[6-1][4-1] - M[1-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][5-1] + M[1-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1]*M[6-1][5-1] + M[1-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1]*M[6-1][5-1] - M[1-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1]*M[6-1][5-1] + M[1-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1]*M[6-1][5-1] - M[1-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][5-1] - M[1-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][5-1] + M[1-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1]*M[6-1][5-1] - M[1-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][5-1] - M[1-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1]*M[6-1][5-1] + M[1-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][6-1] - M[1-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1]*M[6-1][6-1] - M[1-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1]*M[6-1][6-1] + M[1-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][6-1] - M[1-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][4-1]*M[6-1][6-1] + M[1-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][6-1] + M[1-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][6-1] - M[1-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][6-1] + M[1-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][6-1] - M[1-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1]*M[6-1][6-1] - M[1-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][6-1] + M[1-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][6-1] - M[1-1][5-1]*M[6-1][4-1]*v000033 - M[2-1][5-1]*M[6-1][4-1]*v000033 + M[1-1][4-1]*M[6-1][5-1]*v000033 + M[2-1][4-1]*M[6-1][5-1]*v000033 - M[1-1][5-1]*M[6-1][4-1]*v000034 - M[2-1][5-1]*M[6-1][4-1]*v000034 + M[1-1][4-1]*M[6-1][5-1]*v000034 + M[2-1][4-1]*M[6-1][5-1]*v000034 - M[1-1][5-1]*M[6-1][4-1]*v000035 - M[2-1][5-1]*M[6-1][4-1]*v000035 + M[1-1][4-1]*M[6-1][5-1]*v000035 + M[2-1][4-1]*M[6-1][5-1]*v000035 + M[1-1][4-1]*M[6-1][6-1]*v000039 + M[2-1][4-1]*M[6-1][6-1]*v000039 + M[1-1][4-1]*M[6-1][6-1]*v000040 + M[2-1][4-1]*M[6-1][6-1]*v000040 + M[1-1][4-1]*M[6-1][6-1]*v000041 + M[2-1][4-1]*M[6-1][6-1]*v000041 + M[1-1][5-1]*M[6-1][6-1]*v000042 + M[2-1][5-1]*M[6-1][6-1]*v000042 + M[1-1][5-1]*M[6-1][6-1]*v000043 + M[2-1][5-1]*M[6-1][6-1]*v000043 + M[1-1][5-1]*M[6-1][6-1]*v000044 + M[2-1][5-1]*M[6-1][6-1]*v000044 - M[1-1][4-1]*M[3-1][6-1]*v000062 - M[2-1][4-1]*M[3-1][6-1]*v000062 - M[1-1][4-1]*M[3-1][6-1]*v000063 - M[2-1][4-1]*M[3-1][6-1]*v000063 - M[1-1][4-1]*M[3-1][6-1]*v000064 - M[2-1][4-1]*M[3-1][6-1]*v000064 - M[1-1][5-1]*M[3-1][6-1]*v000066 - M[2-1][5-1]*M[3-1][6-1]*v000066 - M[1-1][5-1]*M[3-1][6-1]*v000067 - M[2-1][5-1]*M[3-1][6-1]*v000067 + M[1-1][4-1]*M[6-1][5-1]*v000148 + M[2-1][4-1]*M[6-1][5-1]*v000148 + M[1-1][4-1]*M[6-1][5-1]*v000150 + M[2-1][4-1]*M[6-1][5-1]*v000150 + M[6-1][5-1]*v000151 + M[6-1][5-1]*v000152 + M[6-1][5-1]*v000156 + M[6-1][5-1]*v000157 + M[6-1][5-1]*v000158 + M[6-1][5-1]*v000159 + M[1-1][4-1]*M[6-1][6-1]*v000160 + M[2-1][4-1]*M[6-1][6-1]*v000160 + M[1-1][4-1]*M[6-1][6-1]*v000161 + M[2-1][4-1]*M[6-1][6-1]*v000161 + M[1-1][4-1]*M[6-1][6-1]*v000162 + M[2-1][4-1]*M[6-1][6-1]*v000162 + M[6-1][6-1]*v000163 + M[6-1][6-1]*v000164 + M[1-1][5-1]*M[6-1][6-1]*v000165 + M[2-1][5-1]*M[6-1][6-1]*v000165 + M[1-1][5-1]*M[6-1][6-1]*v000166 + M[2-1][5-1]*M[6-1][6-1]*v000166 + M[1-1][5-1]*M[6-1][6-1]*v000167 + M[2-1][5-1]*M[6-1][6-1]*v000167 + M[6-1][6-1]*v000168 + M[6-1][6-1]*v000169 + M[6-1][6-1]*v000170 + M[6-1][6-1]*v000171 + M[6-1][5-1]*v000376 + M[6-1][5-1]*v000377 + M[6-1][5-1]*v000378 + M[6-1][5-1]*v000379 + M[6-1][5-1]*v000380 + M[6-1][5-1]*v000381 + M[6-1][6-1]*v000390 + M[6-1][6-1]*v000391 + M[6-1][6-1]*v000392 + M[6-1][6-1]*v000393 + M[6-1][6-1]*v000394 + M[6-1][6-1]*v000395 + M[1-1][6-1]*v000603 + M[2-1][6-1]*v000603)*v000607;
N[0][1] = (M[1-1][3-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][2-1] + M[1-1][3-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][2-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][2-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][3-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][3-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][3-1] + M[1-1][3-1]*M[2-1][5-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][4-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1] + M[1-1][3-1]*M[2-1][6-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][5-1]*M[5-1][4-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][4-1] + M[1-1][3-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][6-1]*M[5-1][4-1] + M[1-1][3-1]*M[2-1][6-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][6-1]*M[4-1][2-1]*M[5-1][5-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][5-1] + M[1-1][2-1]*M[2-1][6-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1] + M[1-1][3-1]*M[2-1][2-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][5-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][5-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][3-1]*M[4-1][6-1]*M[5-1][5-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][5-1] + M[1-1][3-1]*M[2-1][4-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][6-1] + M[1-1][2-1]*M[2-1][5-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][6-1] + M[1-1][3-1]*M[2-1][5-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][6-1] + M[1-1][2-1]*M[2-1][3-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][6-1] + M[1-1][2-1]*M[2-1][4-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][6-1] + M[1-1][3-1]*M[2-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][6-1] + M[1-1][4-1]*M[2-1][5-1]*v000033 + M[1-1][4-1]*M[2-1][5-1]*v000034 + M[1-1][4-1]*M[2-1][5-1]*v000035 + M[1-1][4-1]*M[2-1][6-1]*v000039 + M[1-1][4-1]*M[2-1][6-1]*v000040 + M[1-1][4-1]*M[2-1][6-1]*v000041 + M[1-1][4-1]*M[2-1][5-1]*v000148 + M[1-1][4-1]*M[2-1][5-1]*v000149 + M[1-1][4-1]*M[2-1][5-1]*v000150 + M[1-1][4-1]*M[2-1][6-1]*v000160 + M[1-1][4-1]*M[2-1][6-1]*v000161 + M[1-1][4-1]*M[2-1][6-1]*v000162 - v000382 - v000383 - v000384 - v000385 - v000386 - v000387 - v000388 - v000389 - v000396 - v000397 - v000398 - v000399 - v000400 - v000401 - v000402 - v000403 - v000404 - v000405 - v000406 - v000407 - v000408 - v000409 - v000410 - v000411 - v000412 - v000413 - v000414 - v000415 - v000416 - v000417 + M[1-1][6-1]*(M[2-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1] + M[2-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1] + M[2-1][3-1]*M[3-1][5-1]*M[4-1][2-1]*M[5-1][4-1] + M[2-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1] + M[2-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1] + M[2-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1] + M[2-1][4-1]*(M[3-1][5-1]*M[4-1][3-1]*M[5-1][2-1] + M[3-1][2-1]*M[4-1][5-1]*M[5-1][3-1] + M[3-1][3-1]*M[4-1][2-1]*M[5-1][5-1] - v000039 - v000040 - v000041) - v000163 - v000164 - v000168 - v000169 - v000170 - v000171 + M[2-1][5-1]*v000596) + M[1-1][5-1]*(M[2-1][3-1]*M[3-1][6-1]*M[4-1][4-1]*M[5-1][2-1] + M[2-1][2-1]*M[3-1][4-1]*M[4-1][6-1]*M[5-1][3-1] + M[2-1][2-1]*M[3-1][6-1]*M[4-1][3-1]*M[5-1][4-1] + M[2-1][3-1]*M[3-1][2-1]*M[4-1][6-1]*M[5-1][4-1] + M[2-1][3-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][6-1] + M[2-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][6-1] + M[2-1][4-1]*(M[3-1][3-1]*M[4-1][6-1]*M[5-1][2-1] + M[3-1][6-1]*M[4-1][2-1]*M[5-1][3-1] + M[3-1][2-1]*M[4-1][3-1]*M[5-1][6-1] - v000033 - v000034 - v000035) - v000151 - v000152 - v000156 - v000157 - v000158 - v000159 + M[2-1][6-1]*v000597))*v000606;
N[1][0] = ((M[1-1][3-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][2-1]*M[6-1][1-1] + M[1-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][3-1]*M[6-1][1-1] - M[1-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][1-1] - M[1-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][1-1] - M[1-1][3-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][2-1]*M[3-1][3-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][1-1] + M[1-1][2-1]*M[3-1][4-1]*M[4-1][3-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][3-1]*M[3-1][2-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][1-1] - M[1-1][2-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][1-1] + M[1-1][1-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][2-1] + M[1-1][1-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][4-1]*M[6-1][2-1] + M[1-1][3-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][2-1] + M[1-1][3-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1]*M[6-1][2-1] - M[1-1][3-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][2-1] + M[1-1][1-1]*M[3-1][3-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][2-1] - M[1-1][2-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][1-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][4-1]*M[4-1][5-1]*M[5-1][1-1]*M[6-1][3-1] + M[1-1][1-1]*M[3-1][5-1]*M[4-1][4-1]*M[5-1][2-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][5-1]*M[4-1][1-1]*M[5-1][4-1]*M[6-1][3-1] - M[1-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1] + M[1-1][1-1]*M[3-1][2-1]*M[4-1][5-1]*M[5-1][4-1]*M[6-1][3-1] - M[1-1][2-1]*M[3-1][4-1]*M[4-1][1-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][1-1]*M[3-1][4-1]*M[4-1][2-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][1-1]*M[4-1][4-1]*M[5-1][5-1]*M[6-1][3-1] + M[1-1][2-1]*M[3-1][5-1]*M[4-1][3-1]*M[5-1][1-1]*M[6-1][4-1] + M[1-1][2-1]*M[3-1][1-1]*M[4-1][5-1]*M[5-1][3-1]*M[6-1][4-1] + M[1-1][2-1]*M[3-1][3-1]*M[4-1][1-1]*M[5-1][5-1]*M[6-1][4-1] + M[1-1][4-1]*M[6-1][3-1]*v000002 + M[2-1][4-1]*M[6-1][3-1]*v000002 + M[1-1][4-1]*M[6-1][5-1]*v000004 + M[2-1][4-1]*M[6-1][5-1]*v000004 + M[1-1][3-1]*M[6-1][4-1]*v000015 + M[1-1][3-1]*M[6-1][4-1]*v000016 + M[1-1][3-1]*M[6-1][4-1]*v000017 + M[1-1][3-1]*M[6-1][5-1]*v000018 + M[1-1][3-1]*M[6-1][5-1]*v000019 + M[1-1][3-1]*M[6-1][5-1]*v000020 + M[1-1][2-1]*M[6-1][5-1]*v000024 + M[1-1][2-1]*M[6-1][5-1]*v000025 + M[1-1][2-1]*M[6-1][5-1]*v000026 - M[1-1][2-1]*M[6-1][4-1]*v000027 + M[1-1][4-1]*M[6-1][2-1]*v000028 + M[2-1][4-1]*M[6-1][2-1]*v000028 - M[1-1][2-1]*M[6-1][4-1]*v000028 + M[1-1][4-1]*M[6-1][2-1]*v000029 + M[2-1][4-1]*M[6-1][2-1]*v000029 - M[1-1][2-1]*M[6-1][4-1]*v000029 + M[1-1][1-1]*M[6-1][4-1]*v000039 - M[1-1][4-1]*M[6-1][1-1]*v000040 - M[2-1][4-1]*M[6-1][1-1]*v000040 + M[1-1][1-1]*M[6-1][4-1]*v000040 - M[1-1][4-1]*M[6-1][1-1]*v000041 - M[2-1][4-1]*M[6-1][1-1]*v000041 + M[1-1][1-1]*M[6-1][4-1]*v000041 + M[1-1][1-1]*M[6-1][5-1]*v000042 + M[1-1][1-1]*M[6-1][5-1]*v000043 + M[1-1][1-1]*M[6-1][5-1]*v000044 + M[1-1][3-1]*M[3-1][5-1]*v000051 + M[1-1][3-1]*M[3-1][5-1]*v000052 + M[1-1][4-1]*M[3-1][5-1]*v000055 + M[2-1][4-1]*M[3-1][5-1]*v000055 - M[1-1][3-1]*M[3-1][4-1]*v000057 + M[1-1][4-1]*M[3-1][3-1]*v000058 + M[2-1][4-1]*M[3-1][3-1]*v000058 - M[1-1][3-1]*M[3-1][4-1]*v000058 - M[1-1][1-1]*M[3-1][4-1]*v000063 - M[1-1][1-1]*M[3-1][4-1]*v000064 - M[1-1][1-1]*M[3-1][5-1]*v000065 - M[1-1][1-1]*M[3-1][5-1]*v000066 - M[1-1][1-1]*M[3-1][5-1]*v000067 + M[1-1][4-1]*M[6-1][3-1]*v000069 + M[2-1][4-1]*M[6-1][3-1]*v000069 + M[1-1][4-1]*M[6-1][5-1]*v000075 + M[2-1][4-1]*M[6-1][5-1]*v000075 + M[1-1][3-1]*M[6-1][4-1]*v000100 + M[1-1][3-1]*M[6-1][4-1]*v000101 + M[1-1][3-1]*M[6-1][4-1]*v000102 + M[1-1][3-1]*M[6-1][5-1]*v000117 + M[1-1][3-1]*M[6-1][5-1]*v000118 + M[1-1][3-1]*M[6-1][5-1]*v000119 + M[1-1][2-1]*M[6-1][5-1]*v000129 + M[1-1][2-1]*M[6-1][5-1]*v000130 + M[1-1][2-1]*M[6-1][5-1]*v000131 + M[1-1][4-1]*M[6-1][2-1]*v000136 + M[2-1][4-1]*M[6-1][2-1]*v000136 + M[1-1][1-1]*M[6-1][4-1]*v000160 + M[1-1][1-1]*M[6-1][4-1]*v000161 + M[1-1][1-1]*M[6-1][5-1]*v000165 + M[1-1][1-1]*M[6-1][5-1]*v000166 + M[1-1][1-1]*M[6-1][5-1]*v000167 + M[1-1][3-1]*M[3-1][5-1]*v000176 + M[1-1][3-1]*M[3-1][5-1]*v000177 + M[1-1][4-1]*M[3-1][5-1]*v000179 + M[2-1][4-1]*M[3-1][5-1]*v000179 + M[1-1][4-1]*v000183 + M[2-1][4-1]*v000183 + M[1-1][4-1]*v000184 + M[2-1][4-1]*v000184 + M[1-1][4-1]*v000185 + M[2-1][4-1]*v000185 + M[1-1][4-1]*v000186 + M[2-1][4-1]*v000186 + M[1-1][4-1]*v000187 + M[2-1][4-1]*v000187 + M[1-1][4-1]*v000188 + M[2-1][4-1]*v000188 - M[1-1][1-1]*v000200 - M[1-1][1-1]*v000201 + M[1-1][4-1]*v000418 + M[2-1][4-1]*v000418 + M[1-1][4-1]*v000419 + M[2-1][4-1]*v000419 + M[1-1][4-1]*v000420 + M[2-1][4-1]*v000420 + M[1-1][4-1]*v000421 + M[2-1][4-1]*v000421 + M[1-1][4-1]*v000422 + M[2-1][4-1]*v000422 + M[1-1][4-1]*v000423 + M[2-1][4-1]*v000423 + v000424 + v000425 + v000426 + v000427 + v000428 + v000429 + v000430 + v000431 + v000432 + v000433 + v000434 + v000435 + v000436 + v000437 + v000438 + v000439 + v000440 + v000441 + v000442 + v000443 + v000444 + v000445 + v000446 + v000447 + v000448 + v000449 + v000450 + v000451 + v000452 + v000453 + v000454 + v000455 + v000456 + v000457 + v000458 + v000459 + v000555 + v000556 + v000557 + v000558 + v000559 + v000560 + v000561 + v000562 + v000563 + v000564 + v000565 + v000566 + v000567 + v000568 + v000569 + v000570 + v000571 + v000572 + v000573 + v000574 + v000575 + v000576 + v000577 + v000578 + v000579 + v000580 + v000581 + v000582 + v000583 + v000584 + v000585 + v000586 + v000587 + v000588 + v000589 + v000590 + M[1-1][5-1]*v000602 + v000604)*v000607)/2.;
N[1][1] = (v000605*v000606)/2.;


  /* End machine-generated code */


  PetscFunctionReturn(0);
}

void print33s(Tensor33s *T){
  printf("%e\t%e\t%e\nX\t\t%e\t%e\t\nx\t\tX\t\t%e\n",T->T11,T->T12,T->T13,T->T22,T->T23,T->T33);

}

