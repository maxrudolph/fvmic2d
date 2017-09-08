#include "fdcode.h"
#include "thermalSystem.h"
#include "benchmarkInitialConditions.h"

PetscErrorCode enforceThermalBCs1( GridData *grid, Options *options, NodalFields *nodalFields){
  PetscErrorCode ierr;
  PetscFunctionBegin;
  
  /* enforce thermal boundary conditions on projection of last timestep's solution*/
  PetscInt x,y,m,n;
  PetscScalar **lastT;
  ierr = DMDAGetCorners(grid->da,&x,&y,PETSC_NULL,&m,&n,PETSC_NULL);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(grid->da,nodalFields->lastT,&lastT);CHKERRQ(ierr);  

  PetscInt ix,jy;
  /* First enforce dirichlet BCs */
  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){
      if(jy == 0){/* TOP */
	lastT[jy][ix] = options->thermalBCTop.value[0];
      } else if(jy == grid->NY-1){              /* BOTTOM */
	// do nothing
      }
      
      if(ix ==0){/* LEFT */
	lastT[jy][ix] = slab_inflow_temperature( 0.0, grid->y[jy], options->slabAngle ,options->slabAge);

      }else if(ix == grid->NX-1){/* RIGHT */
	if( grid->y[jy] <= plate_depth(grid->LX,options) ){
	  lastT[jy][ix] = plate_geotherm( grid->LX, grid->y[jy],options );	  
	} 
      }
    }/* end loop over y*/
  }/* end loop over x */
  

  ierr =  DMDAVecRestoreArray(grid->da,nodalFields->lastT,&lastT);CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}
