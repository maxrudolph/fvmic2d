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
      if(jy == 0){
	if(options->thermalBCTop.type[0] == 0){ /* TOP */
	  lastT[jy][ix] = options->thermalBCTop.value[0];
	}
      } else if(jy == grid->NY-1){              /* BOTTOM */
	if(options->thermalBCBottom.type[0] == 0){
	  //lastT[jy][ix] = options->thermalBCBottom.value[0];
	}
      }
      
      if(!grid->xperiodic && ix ==0){/* LEFT */
	if(options->thermalBCLeft.type[0] == 0){
	  // left boundary gets half-space cooling
	  lastT[jy][ix] = slab_inflow_temperature( 0.0, grid->y[jy], options->slabAngle );
	}   
      }else if(!grid->xperiodic && ix == grid->NX-1){/* RIGHT */
	//if(options->thermalBCRight.type[0] == 0){
	//lastT[jy][ix] = options->thermalBCRight.value[0];
	//}
      }
    }/* end loop over y*/
  }/* end loop over x */

  /* now enforce Neumann BCs (may use Dirichlet info in corners! */
  for(ix=x;ix<x+m;ix++){
    for(jy=y;jy<y+n;jy++){
      if(jy == 0){
	if(options->thermalBCTop.type[0] == 1){
	  lastT[jy][ix] = lastT[jy+1][ix] - options->thermalBCTop.value[0]*(grid->y[jy+1]-grid->y[jy]);
	}
      }else if(jy == grid->NY-1){
	if(options->thermalBCBottom.type[0] == 1){
	  lastT[jy][ix] = lastT[jy-1][ix] + options->thermalBCBottom.value[0]*(grid->y[jy]-grid->y[jy-1]);
	}
      }
      if(!grid->xperiodic && ix == 0){
	if(options->thermalBCLeft.type[0] == 1){
	   lastT[jy][ix] = lastT[jy][ix+1] - options->thermalBCLeft.value[0]*(grid->x[ix+1]-grid->x[ix]);
	 }	
      }else if(!grid->xperiodic && ix == grid->NX-1){
	if(options->thermalBCRight.type[0] == 1){
	  lastT[jy][ix] = lastT[jy][ix-1] + options->thermalBCRight.value[0]*(grid->x[ix]-grid->x[ix-1]);
	}
      }
    }/* end loop over y*/
  }/* end loop over x */

  ierr =  DMDAVecRestoreArray(grid->da,nodalFields->lastT,&lastT);CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}
