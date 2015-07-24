#include "fdcode.h"
#include "options.h"
#include <float.h>
#include <string.h>
#include "io.h"

#ifndef LINEWIDTH
#define LINEWIDTH 128
#endif


/* Improved options handler */

/* maximum number of options */
#define MAX_OPTIONS 200
/* maximum length of string describing option */
#define OPTSTR_LEN  30

/* statically-allocated option database */
static char option_name_db[MAX_OPTIONS][OPTSTR_LEN];
static char option_value_db[MAX_OPTIONS][OPTSTR_LEN];
static option_type option_type_db[MAX_OPTIONS];/* type of option */
static void *option_ptr_db[MAX_OPTIONS];/* pointer to option value */
static int  option_default_db[MAX_OPTIONS];
static int nopt=0;

PetscErrorCode declare_option( const char *pattern, option_type opt, void *option_ptr, const char *default_value ){
  /* This function should declare an option */
  PetscErrorCode ierr=0;
  const int pattern_len = strlen(pattern);
  /* first, check list of options to see whether this option has yet been declared */
  int declared=0;
  int i;
  for( i=0;i<nopt;i++){
    if( !strncmp( pattern, option_name_db[i], pattern_len) ){
      declared=1;
    }
  }
  if( declared && opt != OPTION_SPECIAL1 ){
    /* if so, abort */
    SETERRQ(PETSC_COMM_WORLD,101,"Parameters may not be declared twice!");
    CHKERRQ(ierr);
  }else{
    if( nopt >= MAX_OPTIONS ) SETERRQ(PETSC_COMM_SELF,101,"Too many parameters.");
    /* if not, add it to the list of options */
    strncpy( option_name_db[nopt] , pattern       , strlen(pattern) );
    strncpy( option_value_db[nopt], default_value , strlen(default_value) );
    option_type_db[nopt] = opt;
    option_ptr_db[nopt] = option_ptr;
    option_default_db[nopt] = 1;
    printf("declared option `%s` to '%s'\n",option_name_db[nopt],option_value_db[nopt]);
    /* parse the default value */
    nopt++;
    parse_option( pattern, default_value );
    option_default_db[i] = 1;/* indicate that default option was set */
  }
  return ierr;
}

PetscErrorCode parse_option( const char *key, const char *value ){
  /* assume input is in form of key/value pair */
  /* compare key with each entry in option database */
  int i;
  for(i=0;i<nopt;i++){
    if( !strncmp( key, option_name_db[i], strlen( option_name_db[i]) ) ){
      /* We have a match */
      printf("found match for key %s\n",key);
      option_default_db[i] = 0;
      strcpy( option_value_db[i], value );
      if( option_type_db[i] == OPTION_SCALAR ){
	sscanf( value, "%le", (double *) option_ptr_db[i] );
	return 0;
      }else if( option_type_db[i] == OPTION_INTEGER ){
	sscanf( value, "%d", (int *) option_ptr_db[i] );
	return 0;
      }else if( option_type_db[i] == OPTION_III ){
	sscanf( value, "%d,%d,%d", ((int *) option_ptr_db[i] ), ((int *) option_ptr_db[i])+1, ((int *) option_ptr_db[i])+2 );
	return 0;
      }else if( option_type_db[i] == OPTION_SSS ){
	sscanf( value, "%le,%le,%le", ((PetscScalar *) option_ptr_db[i] ), ((PetscScalar *) option_ptr_db[i])+1, ((PetscScalar *) option_ptr_db[i])+2 );  
	return 0;
      }else if( option_type_db[i] == OPTION_II ){
	sscanf( value, "%d,%d", ((int *) option_ptr_db[i] ), ((int *) option_ptr_db[i])+1 );  
	return 0;
      }else if( option_type_db[i] == OPTION_SS ){
	sscanf( value, "%le,%le", ((PetscScalar *) option_ptr_db[i] ), ((PetscScalar *) option_ptr_db[i])+1);  
	return 0;
      }else if( option_type_db[i] == OPTION_ISPAIR ){
	int tmp1;
	double tmp2;
	sscanf( value, "%d,%le",&tmp1,&tmp2 );
	PetscScalar *val = (PetscScalar *) option_ptr_db[i];
	val[tmp1]=tmp2;
	return 0;
      }else if( option_type_db[i] == OPTION_IIPAIR ){
	int tmp1;
	int tmp2;
	sscanf( value, "%d,%d",&tmp1,&tmp2 );
	int *val = (int *) option_ptr_db[i];
	val[tmp1]=tmp2;
	return 0;
      }else if( option_type_db[i] == OPTION_SPECIAL1 ){
	/* for strain-weakening plasticity, put two values into correct location in materials structure */
	int tmp1,tmp2;
	double  tmp3;
	sscanf(value,"%d,%d,%le",&tmp1,&tmp2,&tmp3);
	double *plastarray = (double *) option_ptr_db[i];
	*(plastarray + 2*tmp1 + tmp2)= tmp3;/* option_prt_db[i] is a pointer to a statically-allocated MAXMAT x 2 array */
	return 0;
      }else if( option_type_db[i] == OPTION_PSO ){
	/* read a string and match with typedef'd enum */
	PsoType *pso = (PsoType *) option_ptr_db[i];
	if( !strncmp(value,"none",4) ){
	  *pso = PSO_NONE;
	  return 0;
	}else if( !strncmp(value,"KinematicSubduction",19 ) ){
	  *pso = PSO_KINEMATIC_WEDGE;
	  return 0;
	}else if( !strncmp(value,"Sandbox",7) ){
	  *pso = PSO_SANDBOX;
	  return 0;
	}
      }else if( option_type_db[i] == OPTION_GRIDTYPE ){
	GridType *g = (GridType *) option_ptr_db[i];
	if( !strncmp(value,"regular",7) ){
	  *g = GRID_REGULAR;
	  return 0;
	}else if( !strncmp(value,"Subduction",10) ){
	  *g = GRID_SUBDUCTION;
	  return 0;
	}else if( !strncmp(value,"ConstantInnerOuter",18) ){
	  *g = GRID_CONSTANTINNEROUTER;
	  return 0; 
	}
      }else{
	SETERRQ(PETSC_COMM_SELF,101,"Unknown Option type");
      }
    }
  }
  fprintf(stderr,"No match for key '%s' value '%s'\n",key,value);
  SETERRQ(PETSC_COMM_SELF,101,"No match for specified option in option db");
  abort();
}

/* csvOptions initializes options from a csv (comma-separated values) file */
PetscErrorCode csvOptions(Options *options, Materials *materials){
  PetscErrorCode ierr=0;
  PetscFunctionBegin;
  FILE *csvFile;
  //  PetscErrorCode ierr;
  char line[LINEWIDTH];
  PetscMPIInt rank,size;
  const char default_file_name[12] = "options.csv";
  char *ifn;
  char csv_file_name[80];
  PetscFunctionBegin;

  /* Get input filename using petsc options interface */
  PetscBool set = PETSC_FALSE;
  ierr=PetscOptionsGetString(PETSC_NULL,"-input_file",csv_file_name,sizeof(csv_file_name),&set);CHKERRQ(ierr);

  if( set ){
    ifn = csv_file_name;
  } else {
    ifn = default_file_name;
  }

  PetscFunctionBegin;
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  if(!rank) printf("reading from %s\n",ifn);

  if(!rank){ 
    /* declare each option */
    declare_option("LX",OPTION_SCALAR, &(options->LX), "1.0");
    declare_option("LY",OPTION_SCALAR, &(options->LY), "1.0");
    declare_option("NX",OPTION_INTEGER, &(options->NX), "21");
    declare_option("NY",OPTION_INTEGER, &(options->NY), "20");
    declare_option("gridRatio",OPTION_SCALAR, &(options->gridRatio), "1.0");
    declare_option("gridType",OPTION_GRIDTYPE, &(options->gridType), "regular");
    declare_option("NMX",OPTION_INTEGER, &(options->NMX), "4");
    declare_option("NMY",OPTION_INTEGER, &(options->NMY), "4");
    declare_option("saveInterval",OPTION_INTEGER,&(options->saveInterval),"100");
    declare_option("maxMarkFactor",OPTION_SCALAR,&(options->maxMarkFactor),"2.0");
    declare_option("minMarkers",OPTION_INTEGER, &(options->minMarkers), "1");
    declare_option("maxMarkers",OPTION_INTEGER, &(options->maxMarkers), "100");
    declare_option("dtMax",OPTION_SCALAR,&(options->dtMax),"1.0e13");
    declare_option("displacementStepLimit",OPTION_SCALAR,&(options->displacementStepLimit),"0.1");
    declare_option("maxTempChange",OPTION_SCALAR,&(options->maxTempChange),"100.0");
    declare_option("nTime",OPTION_INTEGER,&(options->nTime),"1");
    declare_option("restartStep",OPTION_INTEGER,&(options->restartStep),"0");
    declare_option("totalTime",OPTION_SCALAR,&(options->totalTime),"3.15e16");
    declare_option("maxNumPlasticity",OPTION_INTEGER,&(options->maxNumPlasticity),"1");
    declare_option("plasticDelay",OPTION_INTEGER,&(options->plasticDelay),"0");
    declare_option("etamin",OPTION_SCALAR,&(options->etamin),"1.0e14");
    declare_option("fractionalEtamin",OPTION_SCALAR,&(options->fractionalEtamin),"1e-4");
    declare_option("etamax",OPTION_SCALAR,&(options->etamax),"1.0e25");
    declare_option("grainSize",OPTION_SCALAR,&(options->grainSize),"NaN");
    declare_option("Tperturb",OPTION_SCALAR,&(options->Tperturb),"0.0");
    declare_option("shearHeating",OPTION_INTEGER,&(options->shearHeating),"0");
    declare_option("adiabaticHeating",OPTION_INTEGER,&(options->adiabaticHeating),"0");
    declare_option("mechBCLeftType",OPTION_III,options->mechBCLeft.type,"0,0,0");
    declare_option("mechBCLeftValue",OPTION_SSS,options->mechBCLeft.value,"0,0,0");
    declare_option("mechBCRightType",OPTION_III,options->mechBCRight.type,"0,0,0");
    declare_option("mechBCRightValue",OPTION_SSS,options->mechBCRight.value,"0,0,0");
    declare_option("mechBCTopType",OPTION_III,options->mechBCTop.type,"0,0,0");
    declare_option("mechBCTopValue",OPTION_SSS,options->mechBCTop.value,"0,0,0");
    declare_option("mechBCBottomType",OPTION_III,options->mechBCBottom.type,"0,0,0");
    declare_option("mechBCBottomValue",OPTION_SSS,options->mechBCBottom.value,"0,0,0");
    /* thermal BCs */
    declare_option("thermalBCBottomType",OPTION_INTEGER,options->thermalBCBottom.type,"0");
    declare_option("thermalBCBottomValue",OPTION_SCALAR,options->thermalBCBottom.value,"0.0");
    declare_option("thermalBCTopType",OPTION_INTEGER,   options->thermalBCTop.type,"0");
    declare_option("thermalBCTopValue",OPTION_SCALAR,   options->thermalBCTop.value,"0.0");
    declare_option("thermalBCLeftType",OPTION_INTEGER,  options->thermalBCLeft.type,"0");
    declare_option("thermalBCLeftValue",OPTION_SCALAR,  options->thermalBCLeft.value,"0.0");
    declare_option("thermalBCRightType",OPTION_INTEGER, options->thermalBCRight.type,"0");
    declare_option("thermalBCRightValue",OPTION_SCALAR, options->thermalBCRight.value,"0.0");
    
    declare_option("gy",OPTION_SCALAR,&(options->gy),"9.8");
    declare_option("gx",OPTION_SCALAR,&(options->gx),"0.0");
    declare_option("subgridStressDiffusivity",OPTION_SCALAR,&(options->subgridStressDiffusivity),"1.0");
    declare_option("subgridTemperatureDiffusivity",OPTION_SCALAR,&(options->subgridTemperatureDiffusivity),"1.0");
    declare_option("nMaterials",OPTION_INTEGER,&(materials->nMaterials),"1");
    /* material properties */
    declare_option("materialEta",OPTION_ISPAIR,materials->materialEta,"0,1.0e21");
    declare_option("materialRho",OPTION_ISPAIR,materials->materialRho,"0,4000.0");
    declare_option("materialkThermal",OPTION_ISPAIR,materials->materialkThermal,"0,4.0");
    declare_option("materialCp",OPTION_ISPAIR,materials->materialCp,"0,1250.0");
    declare_option("materialAlpha",OPTION_ISPAIR,materials->materialAlpha,"0,0.0");
    declare_option("materialMu",OPTION_ISPAIR,materials->materialMu,"0,1.0e99");
    declare_option("materialCohesion",OPTION_ISPAIR,materials->materialCohesion,"0,1.0e99");
    declare_option("materialFriction",OPTION_ISPAIR,materials->materialFriction,"0,1.0");
    declare_option("hasPlasticity",OPTION_ISPAIR,materials->hasPlasticity,"0,0");
    declare_option("hasEtaT",OPTION_IIPAIR,materials->hasEtaT,"0,0");
    declare_option("QonR",OPTION_ISPAIR,materials->QonR,"0,0.0");
    declare_option("Tref",OPTION_ISPAIR,materials->Tref,"0,273.0");
    /* special parameters for strain-weakening plasticity */
    declare_option("materialC",OPTION_SPECIAL1,materials->C,"0,0,1.0e99");
    declare_option("materialC",OPTION_SPECIAL1,materials->C,"0,1,1.0e99");
    declare_option("materialF",OPTION_SPECIAL1,materials->F,"0,0,1.0e99");
    declare_option("materialF",OPTION_SPECIAL1,materials->F,"0,1,1.0e99");
    declare_option("materialGamma",OPTION_SPECIAL1,materials->gamma,"0,0,1e99"); /* Note - this is not working right now because it will require double-declaration of this parameter*/
    declare_option("materialGamma",OPTION_SPECIAL1,materials->gamma,"0,1,1e98"); 
    /* enable problem-specific choices */
    declare_option("problemSpecificOptions",OPTION_PSO,&options->pso,"none");

    /* subduction - specific stuff */
    declare_option("slabAngle",OPTION_SCALAR, &options->slabAngle,"45");
    declare_option("slabVelocity",OPTION_SCALAR, &(options->slabVelocity),"1.0e-10");
    declare_option("rootThickness",OPTION_SCALAR,&options->rootThickness,"0.0");
    declare_option("rootCenter",OPTION_SCALAR,&options->rootCenter,"1.20e5");
    declare_option("rootWidth",OPTION_SCALAR,&options->rootWidth,"4.0e4");
    declare_option("staticVelocity",OPTION_INTEGER,&options->staticVelocity,"0");

    /* Parse the options file line-by-line */
    csvFile = fopen( ifn, "r");
    while( fgets(line, sizeof(line), csvFile) != NULL ){
      /* parse the line */
      PetscInt idxComma=0;
      /* find the location of the comma */
      if( !strncmp( "#",line,1) ){
	
      } else {
	while( line[++idxComma] != ','){
	  /* do nothing - just increment the comma locator */
	}
	
	{
	  char key[OPTSTR_LEN];
	  char val[OPTSTR_LEN];
	  strncpy( key, line, idxComma );
	  //strcpy( value, line+idxComma+1, 
	  key[ idxComma ] = '\0';/* put in termination character */
	  int ii=0;
	  while( line[idxComma+ii+1] != '\0' && line[idxComma+ii+1] != '\n' && line[idxComma+ii+1] != '#' ){
	    val[ii] = line[idxComma+ii+1];
	    ii++;
	  }
	  val[ii] = '\0';
	  parse_option( key, val );
	}	
      }
    }
      
    fclose( csvFile );
    print_options( PETSC_NULL );
  }
  options->slabAngle = options->slabAngle/180.0*M_PI;
  /* send parameters to all nodes */
  MPI_Bcast( options, sizeof( Options ), MPI_BYTE, 0, PETSC_COMM_WORLD);
  MPI_Bcast( materials, sizeof( Materials ), MPI_BYTE, 0, PETSC_COMM_WORLD);
  /* save the run info to file */
  if( !rank )  ierr = saveRunInfo( options, materials, 0);
  PetscFunctionReturn(ierr);
}


void print_options( FILE *fp ){
  /* loop over entries in the option database */
  int i;
  for(i=0;i<nopt;i++){
    printf("%s,%s\n",option_name_db[i],option_value_db[i]);
  }

}
