/* subroutine to form LHS and RHS for stokes flow problem */
#include "fdcode.h"

#define NADD 12 /* number of entries to add to matrix at a time*/

PetscErrorCode formStokesVarViscSystem(NodalFields *nodalFields, GridData *grid, Mat LHS, Vec RHS, PetscScalar Kbond, PetscScalar Kcont, PetscScalar gy){
  PetscInt i,j;
  PetscInt NX = grid->NX;
  PetscInt NY = grid-> NY;
/*   PetscScalar LX = grid-> LX; */
/*   PetscScalar LY = grid->LY; */
  PetscScalar dx = grid->dx;
  PetscScalar dx2 = grid->dx2;
  PetscScalar dy = grid->dy;
  PetscScalar dy2 = grid->dy2;
  PetscScalar dxdy=grid->dxdy;
  PetscScalar *etaN = nodalFields->etaN;
  PetscScalar *etaS = nodalFields->etaS;
  PetscScalar *rho = nodalFields -> rho;

  PetscErrorCode ierr;

  PetscFunctionBegin;

  for (i=1;i<=NY;i++){
    for(j=1;j<=NX;j++){
      PetscInt idxnode = (j-1)*(NY)+(i-1);
      PetscInt pdof=3*idxnode+0;
      PetscInt vxdof=3*idxnode+1;
      PetscInt vydof=3*idxnode+2;
      
      PetscScalar vals[NADD];
      PetscInt rowidx[NADD];
      PetscInt colidx[NADD];
      
      
      /**              %continuity*/		       
      if( (i == 1 || j == 1) ){
	/*       %pressure ghost nodes. Set equal to zero. */
	ierr = MatSetValue(LHS,pdof,pdof,1.0*Kbond,INSERT_VALUES); CHKERRQ(ierr);
	ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES); CHKERRQ(ierr);
      } else if( (i==2 && j == 2) || (i==NY && j==2)){
	//upper and lower left corner - dp/dx = 0;
	//LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	ierr = MatSetValue(LHS,pdof,pdof, 1.0*Kbond,INSERT_VALUES); CHKERRQ(ierr);
	//      LHS(p(i,j),p(i,j+1)) = -1.0*Kbond;
	ierr = MatSetValue(LHS,pdof,pdof + 3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	//      RHS(p(i,j)) = 0.0;
	ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES); CHKERRQ(ierr);
      } else if( (i==2 && j==NX) || (i==NY && j== NX)) {
	//upper and lower right corners - dp/dx = 0;
	//LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	ierr = MatSetValue(LHS,pdof,pdof,1.0*Kbond,INSERT_VALUES); CHKERRQ(ierr);
	//LHS(p(i,j),p(i,j-1)) = -1.0*Kbond;
	ierr = MatSetValue(LHS,pdof,pdof - 3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	//	RHS(p(i,j)) = 0.0;
	ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES); CHKERRQ(ierr);
	
      } else if( i==2 && j==3){
	//	LHS(p(i,j),p(i,j)) = 1.0*Kbond;
	ierr = MatSetValue(LHS,pdof,pdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	//RHS(p(i,j)) = 0.0;
	ierr = VecSetValue(RHS,pdof,0.0,INSERT_VALUES);CHKERRQ(ierr);
      } else {
	//LHS(p(i,j), vx(i-1,j)) =  Kcont/dx;
	ierr = MatSetValue(LHS,pdof,vxdof - 3,Kcont/dx,INSERT_VALUES);CHKERRQ(ierr);
	//LHS(p(i,j), vx(i-1,j-1) ) = -Kcont/dx;
	ierr = MatSetValue(LHS,pdof,vxdof-3-3*NY,-Kcont/dx,INSERT_VALUES);CHKERRQ(ierr);
	//LHS(p(i,j), vy(i,j-1) ) =  Kcont/dy;
	ierr = MatSetValue(LHS,pdof,vydof-3*NY,Kcont/dy,INSERT_VALUES);CHKERRQ(ierr);
	//LHS(p(i,j), vy(i-1,j-1) ) = -Kcont/dy;
	ierr = MatSetValue(LHS,pdof,vydof-3-3*NY,-Kcont/dy,INSERT_VALUES);CHKERRQ(ierr);
	//      end
      }
      
      /*         %x-momentum */
      if(i == NY){
	/*             %ghost vx nodes. set to zero. */
	/*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; RHS(vx(i,j)) = 0.0; */
	ierr =  MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/* RHS is automatic*/
	
      } else if( (j==1 || j == NX) && i < NY){
	/*             %left boundary. vx=0; */
	/*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; RHS(vx(i,j)) = 0.0; */
	ierr =  MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
      }else if( (i ==1) && j>1 && j<NX ){
	/*             %top boundary - free slip, dvx/y = 0 */
	/*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; */
	ierr=MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             LHS(vx(i,j),vx(i+1,j)) = -1.0*Kbond; */
	ierr=MatSetValue(LHS,vxdof,vxdof+3,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             RHS(vx(i,j)) = 0.0; */
      }else if( i==NY-1 && j>1 && j<NX ) {
	/*             LHS(vx(i,j),vx(i,j)) = 1.0*Kbond; */
	ierr=MatSetValue(LHS,vxdof,vxdof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             LHS(vx(i,j),vx(i-1,j)) = -1.0*Kbond; */
	ierr=MatSetValue(LHS,vxdof,vxdof-3,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             RHS(vx(i,j)) = 0.0; */
      }else{
	
	/*   LHS(vx(i,j),vx(i,j+1)) =  2*etaN(i+1,j+1)/dx2; */
	rowidx[0] = vxdof; colidx[0] = vxdof+3*NY; vals[0] = 2*etaN[idxnode+NY+1]/dx2;
	/*             LHS(vx(i,j),vx(i,j)) =   -2*etaN(i+1,j+1)/dx2- 2*etaN(i+1,j)/dx2 -etaS(i+1,j)/dy2 - etaS(i,j)/dy2; */
	rowidx[1] = vxdof; colidx[1] = vxdof; vals[1] = -2.0*etaN[idxnode+NY+1]/dx2 -2.0*etaN[idxnode+1]/dx2 - etaS[idxnode+1]/dy2 - etaS[idxnode]/dy2;
	/*             LHS(vx(i,j),vx(i,j-1)) = 2*etaN(i+1,j)/dx2; */
	rowidx[2] = vxdof;
	colidx[2] = vxdof-3*NY;
	vals[2] = 2*etaN[idxnode+1]/dx2;
	/*             LHS(vx(i,j),vx(i+1,j)) = etaS(i+1,j)/dy2;  */
	rowidx[3] = vxdof;
	colidx[3] = vxdof+3;
	vals[3] = etaS[idxnode+1]/dy2;
	/*             LHS(vx(i,j),vy(i+1,j)) = etaS(i+1,j)/dxdy; */
	rowidx[4] = vxdof;
	colidx[4] = vydof+3;
	vals[4] = etaS[idxnode+1]/dxdy;
	/*             LHS(vx(i,j),vy(i+1,j-1)) = -etaS(i+1,j)/dxdy; */
	rowidx[5] = vxdof;
	colidx[5] = vydof+3-3*NY;
	vals[5] = -1.0*etaS[idxnode+1]/dxdy;
	/*             LHS(vx(i,j),vx(i-1,j)) = +etaS(i,j)/dy2; */
	rowidx[6] = vxdof;
	colidx[6] = vxdof-3;
	vals[6] = 1.0*etaS[idxnode]/dy2;
	/*             LHS(vx(i,j),vy(i,j)) = -etaS(i,j)/dxdy; */
	rowidx[7] = vxdof;
	colidx[7] = vydof;
	vals[7] = -1.0*etaS[idxnode]/dxdy;
	/*             LHS(vx(i,j),vy(i,j-1)) = etaS(i,j)/dxdy; */
	rowidx[8] = vxdof;
	colidx[8] = vydof-3*NY;
	vals[8] = etaS[idxnode]/dxdy;
	/*             LHS(vx(i,j),p(i+1,j+1)) = - Kcont/dx; */
	rowidx[9] = vxdof;
	colidx[9] = pdof+3*NY+3;
	vals[9] = -Kcont/dx;
	/*             LHS(vx(i,j),p(i+1,j)) =  Kcont/dx; */
	rowidx[10] = vxdof;
	colidx[10] = pdof+3;
	vals[10] = Kcont/dx;
	/*             RHS(vx(i,j)) =  0; */
	ierr = MatSetValues( LHS, 1,&rowidx[0],11,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
	
	
      } 
      
      /*         %y-momentum */
      if( j==NX ) {
	/*             %ghost vy nodes. set to zero. */
	/*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; RHS(vy(i,j)) =0.0; */
	ierr=MatSetValue(LHS,vydof,vydof,Kbond,INSERT_VALUES);CHKERRQ(ierr);
      }else if( i==1 || i == NY) {
	/*             %top and bottom boundary. vy=0 */
	/*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; RHS(vy(i,j)) =0.0; */
	ierr=MatSetValue(LHS,vydof,vydof,Kbond,INSERT_VALUES);CHKERRQ(ierr);
      }else if( j==1 && i>1 && i < NY){
	/*             %left wall - free slip - dvy/dx = 0; */
	/*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; */
	ierr=MatSetValue(LHS,vydof,vydof, Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             LHS(vy(i,j),vy(i,j+1)) = -1.0*Kbond; */
	ierr=MatSetValue(LHS,vydof,vydof+3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             RHS(vy(i,j)) = 0.0; */
      }else if( j == NX-1 && i >1 && i < NY) {
	/*             %right wall - free slip - dvy/dx =0; */
	/*             LHS(vy(i,j),vy(i,j)) = 1.0*Kbond; */
	ierr=MatSetValue(LHS,vydof,vydof,1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             LHS(vy(i,j),vy(i,j-1)) = -1.0*Kbond; */
	ierr=MatSetValue(LHS,vydof,vydof-3*NY,-1.0*Kbond,INSERT_VALUES);CHKERRQ(ierr);
	/*             RHS(vy(i,j)) = 0.0; */
      }else {
	
	/*             LHS(vy(i,j),vy(i+1,j)) = 2*etaN(i+1,j+1)/dy2; */
	rowidx[0] = vydof; colidx[0] = vydof+3; vals[0] = 2*etaN[idxnode+NY+1]/dy2;
	
	/* LHS(vy(i,j),vy(i,j)) = -2*etaN(i+1,j+1)/dy2 - 2*etaN(i,j+1)/dy2 -etaS(i,j+1)/dx2 - etaS(i,j)/dx2 ; */
	rowidx[1] = vydof;
	colidx[1] = vydof;
	vals[1] = -2.0*etaN[idxnode+NY+1]/dy2 -2.0*etaN[idxnode+NY]/dy2 - etaS[idxnode+NY]/dx2 - etaS[idxnode]/dx2;
	/* LHS(vy(i,j),vy(i-1,j)) = 2*etaN(i,j+1)/dy2; */
	rowidx[2] = vydof;
	colidx[2] = vydof-3;
	vals[2] = 2*etaN[idxnode+NY]/dy2;
	/* LHS(vy(i,j),vy(i,j+1)) = etaS(i,j+1)/dx2; */
	rowidx[3] = vydof;
	colidx[3] = vydof+3*NY;
	vals[3] = etaS[idxnode+NY]/dx2;
	/*             LHS(vy(i,j),vx(i,j+1)) = etaS(i,j+1)/dxdy; */
	rowidx[4] = vydof;
	colidx[4] = vxdof+3*NY;
	vals[4] = etaS[idxnode+NY]/dxdy;
	/*             LHS(vy(i,j),vx(i-1,j+1)) = -etaS(i,j+1)/dxdy; */
	rowidx[5] = vydof;
	colidx[5] = vxdof-3+3*NY;
	vals[5] = -1.0*etaS[idxnode+NY]/dxdy;
	/*  LHS(vy(i,j),vy(i,j-1)) = etaS(i,j)/dx2; */
	rowidx[6] = vydof;
	colidx[6] = vydof-3*NY;
	vals[6] = etaS[idxnode]/dx2;
	/*  LHS(vy(i,j),vx(i,j))   = - etaS(i,j)/dxdy; */
	rowidx[7] = vydof;
	colidx[7] = vxdof;
	vals[7] = -1.0*etaS[idxnode]/dxdy;
	/* LHS(vy(i,j),vx(i-1,j)) = etaS(i,j)/dxdy; */
	rowidx[8] = vydof;
	colidx[8] = vxdof-3;
	vals[8] = etaS[idxnode]/dxdy;
	/* LHS(vy(i,j),p(i+1,j+1)) = -Kcont/dy; */
	rowidx[9] = vydof;
	colidx[9] = pdof+3*NY+3;
	vals[9] = -Kcont/dy;
	/* LHS(vy(i,j),p(i,j+1)) = Kcont/dy; */
	rowidx[10] = vydof;
	colidx[10] = pdof+3*NY;
	vals[10] = Kcont/dy;
	/*   RHS(vy(i,j)) = -gy*(RHO(i,j)+RHO(i,j+1))/2; */
	ierr = VecSetValue( RHS, vydof, -gy*(rho[idxnode]+rho[idxnode+NY])/2.0,INSERT_VALUES);CHKERRQ(ierr);
	ierr = MatSetValues( LHS, 1,&rowidx[0],11,&colidx[0],&vals[0],INSERT_VALUES);CHKERRQ(ierr);
      }
    }/* end loop over i (Y)*/
  }/* end loop over j (x) */
  
  /* finalize assembly*/
  ierr = VecAssemblyBegin(RHS);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(RHS);CHKERRQ(ierr);
  
  ierr = MatAssemblyBegin( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd( LHS, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  PetscFunctionReturn(ierr);
}
