#include "fdcode.h"
#define FN_MAX_LEN 160
#define NLAT 9
#define NLON 18


PetscScalar as[NLAT-1][NLON-1];

PetscScalar areaquad(PetscScalar, PetscScalar, PetscScalar, PetscScalar);
PetscScalar textureIndex( PetscScalar *, PetscScalar *, PetscInt , PetscScalar , PetscScalar );

int main(int argc, char *argv[]){
  if(argc != 3){
    printf("usage: ./ProgramName input.file output.file\n");
    return(0);
  }
  
  /* set up lat-lon grid */
  PetscScalar lat[NLAT];
  PetscScalar lon[NLON];
  {/* set up lat-lon grid */    
    PetscInt i;
    for(i=0;i<NLON;i++){
      lon[i] = i * 2.0*M_PI/(NLON-1);
      printf("%e ",lon[i]);
    }printf("\n");
    for(i=0;i<NLAT;i++){
      lat[i] = -M_PI/2.0 + i * M_PI/(NLAT-1);
      printf("%e ",lat[i]);
    }printf("\n");   
  }
  PetscScalar dlat = lat[1]-lat[0];
  PetscScalar dlon = lon[1]-lon[0];
  
  {
    PetscInt ilat,ilon;
    PetscScalar asum = 0.0;
    for(ilat=0;ilat<NLAT-1;ilat++){
      for(ilon=0;ilon<NLON-1;ilon++){
	as[ilat][ilon] = areaquad(lat[ilat],lon[ilon],lat[ilat+1],lon[ilon+1]);
	asum += as[ilat][ilon];
      }
    }
    printf("asum = %e\n",asum);
  }


  /* read the texture information */
  PetscInt ncpu;
  PetscInt nmark;
  PetscInt nt;
  FILE *ifile, *ofile;

  char ifn[FN_MAX_LEN];
  
  sprintf(ifn,"%s.%d",argv[1],0);
  ifile = fopen(ifn,"rb");
  printf("Reading from %s\n",ifn);
  fread(&ncpu,1,sizeof(PetscInt),ifile);
  printf("ncpu = %d\n",ncpu);
  fclose(ifile);
  PetscInt icpu;
  PetscInt nmarkg=0;
  for(icpu=0;icpu<ncpu;icpu++){
    sprintf(ifn,"%s.%d",argv[1],icpu);
    ifile = fopen(ifn,"rb");
    //printf("Reading from %s\n",ifn);
    fread(&ncpu,1,sizeof(PetscInt),ifile);
    //    printf("ncpu = %d\n",ncpu);
    fread(&nmark,1,sizeof(PetscInt),ifile);
    printf("nmark = %d\n",nmark);
    nmarkg += nmark;
    fread(&nt,1,sizeof(PetscInt),ifile); printf("NT = %d\n",nt);
    fclose(ifile);
  }
  printf("TOTAL MARKERS %d\n",nmarkg);
  ofile = fopen(argv[2],"wb");
  fwrite(&nmarkg,sizeof(PetscInt),1,ofile);/* global number of markers */
  
  for(icpu=0;icpu<ncpu;icpu++){
    sprintf(ifn,"%s.%d",argv[1],icpu);
    ifile = fopen(ifn,"rb");
    printf("Reading from %s\n",ifn);
    fread(&ncpu,1,sizeof(PetscInt),ifile);
    printf("ncpu = %d\n",ncpu);
    fread(&nmark,1,sizeof(PetscInt),ifile);
    printf("nmark = %d\n",nmark);
    fread(&nt,1,sizeof(PetscInt),ifile); printf("NT = %d\n",nt);
  

    /* allocate storage for textures*/
    PetscInt nt3 = nt*nt*nt;
    PetscScalar **ctheta;
    PetscScalar **cphi;
    PetscInt imark;
    PetscMalloc( nmark*sizeof(PetscScalar *), &ctheta);
    PetscMalloc( nmark*sizeof(PetscScalar *), &cphi);
    for(imark=0;imark<nmark;imark++){
      PetscMalloc( nt3*sizeof(PetscScalar), &ctheta[imark]);
      PetscMalloc( nt3*sizeof(PetscScalar), &cphi[imark]);
    }
    /* read the texture */
    for(imark=0;imark<nmark;imark++){
      fread( ctheta[imark], nt3, sizeof(PetscScalar), ifile);
      fread( cphi[imark], nt3, sizeof(PetscScalar), ifile);     
    }
    /* calclate texture index */
    PetscScalar *ti;
    PetscMalloc( nmark*sizeof(PetscScalar), &ti);
    for(imark=0;imark<nmark;imark++){
      /*      ti[imark] =*/ 
      ti[imark] = textureIndex( ctheta[imark], cphi[imark], nt3, dlat, dlon );
    }
    fwrite(ti,sizeof(PetscScalar),nmark,ofile);
    PetscFree( ti );

    /* free texture storage*/
    for(imark=0;imark<nmark;imark++){
      PetscFree(ctheta[imark]);
      PetscFree(cphi[imark]);
    }
    PetscFree(ctheta);
    PetscFree(cphi);

    /* close input file */
    fclose(ifile);
  }
  fclose(ofile);
  return(0);
}

/* implements Matlab built-in function areaquad.m for sphere of radius 1*/
PetscScalar areaquad(PetscScalar lat1, PetscScalar lon1, PetscScalar lat2, PetscScalar lon2){
  PetscScalar a;
  a = fabs(lon1-lon2) *  fabs( sin(lat1)-sin(lat2) )/(4.0*M_PI);

  return(a);


}

PetscScalar textureIndex( PetscScalar *clon, PetscScalar *clat, PetscInt nt3, PetscScalar dlat, PetscScalar dlon){
  /* clon == ctheta, clat == cphi */
  /* count number of crystals in each lat-lon box */
  PetscInt nc[NLAT-1][NLON-1];
  {PetscInt i,j;
    for(i=0;i<NLAT-1;i++){
      for(j=0;j<NLON-1;j++){
	nc[i][j] = 0;
      }
    }
  }
  PetscInt i;
  for(i=0;i<nt3;i++){
    PetscInt ilat=floor(clat[i]/dlat);  if(ilat>NLAT-2){ printf("error: lat=%e ilat =%d too large, NLAT = %d\n",clat[i],ilat,NLAT);}
    PetscInt ilon=floor(clon[i]/dlon);  if(ilon>NLON-2){ printf("error: ilon too large\n");}
    nc[ilat][ilon]++;
  }
  {
    PetscInt i,j;
    PetscInt ncmax=0;PetscInt ncmin=9999;
    for(i=0;i<NLAT-1;i++){
      for(j=0;j<NLON-1;j++){
	if(ncmax < nc[i][j]) ncmax = nc[i][j] ;
	if(ncmax > nc[i][j]) ncmin = nc[i][j] ;
      }
    }
    printf("NC min max [%d %d]\n",ncmin,ncmax);
  }
  PetscScalar f[NLAT-1][NLON-1];
  PetscScalar integrand=0.0;
  {PetscInt i,j;
    for(i=0;i<NLAT-1;i++){
      for(j=0;j<NLON-1;j++){
	f[i][j] = ((PetscScalar) nc[i][j])/((PetscScalar) nt3)/as[i][j];
      }
    }
    for(i=0;i<NLAT-1;i++){
      for(j=0;j<NLON-1;j++){
	integrand += f[i][j]*f[i][j]*as[i][j];
      }
    }
  }
  printf("TI=%e\n",integrand);
  return(integrand);

}
