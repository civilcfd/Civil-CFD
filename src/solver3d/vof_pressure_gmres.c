/* vof_pressure_gmres.c
 *
 * implement pressure iteration
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <petscksp.h>

#include "vtk.h"
#include "vof.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "track.h"

#include "vof_macros.h"


int vof_pressure_gmres(struct solver_data *solver) {
	PetscInt i,j,k,l,m,n,size;
	PetscInt  row_idx[7];
	PetscInt  nidx, nlmn;
	PetscInt	iter;
	PetscScalar rhs, dt_rho, dt_rhodx2, dt_rhody2, dt_rhodz2, pijk, dpijk, mpeta;
	PetscScalar delp;
	PetscScalar row[7];
	PetscErrorCode ierr;
	Vec x, b;
	Mat A;
  KSP            ksp;         /* linear solver context */
  PC             pc;           /* preconditioner context */
  static int initialize = 0;
  
	if(!initialize) {
		PetscInitialize(NULL, NULL, NULL, NULL);
		initialize = 1;
	}
	
	size = IMAX * JMAX * KMAX;
	
	dt_rho = solver->delt / solver->rho;
	dt_rhodx2 = dt_rho / pow(DELX,2);
	dt_rhody2 = dt_rho / pow(DELY,2);
	dt_rhodz2 = dt_rho / pow(DELZ,2);
	
	/*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,size);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	
	  /*
     Create matrix.  When using MatCreate(), the matrix format can
     be specified at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  //ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,size,size,7,NULL,&A);
  //ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRQ(ierr);
  //ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  //ierr = MatSetUp(A);CHKERRQ(ierr);
  
  /* Assemble matrix */
  
#define emf solver->emf

  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {      
      
        
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
               
        if(N_VOF(i,j,k) != 0) {
          l = i;
          m = j;
          n = k;
          switch(N_VOF(i,j,k)) {
          case east:
            l=i+1;
            break;
          case west:
            l=i-1;
            break;
          case north:
            m=j+1;
            break;
          case south:
            m=j-1;
            break;
          case top:
            n=k+1;
            break;
          case bottom:
            n=k-1;
            break;
          case none:
          default: 
            continue;
          }
          
          nidx = mesh_index(solver->mesh,i,j,k);
          nlmn = mesh_index(solver->mesh,l,m,n);
          
          if(N_VOF(l,m,n) != 0) {
          	ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          	ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          	continue;
          }
          
          mpeta = 1 - PETA(i,j,k);
          pijk  = mpeta * P(l,m,n);
          dpijk = pijk - P(i,j,k);
          
          ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          ierr = MatSetValue(A,nidx,nlmn,-1.0 * mpeta,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecSetValue(b,nidx,dpijk,INSERT_VALUES);CHKERRQ(ierr);
          
          
        } else {
        /* interior non-void cell */
        
        	dpijk = dt_rhodx2 * (AE(i,j,k) + AE(i-1,j,k)) + 
        				  dt_rhody2 * (AN(i,j,k) + AN(i,j-1,k)) +
        				  dt_rhodz2 * (AT(i,j,k) + AT(i,j,k-1));

        	
        	row_idx[0] = mesh_index(solver->mesh,i,j,k);
        	row[0] = dpijk * -1.0;
        	
        	row_idx[1] = mesh_index(solver->mesh,i+1,j,k);
        	row_idx[2] = mesh_index(solver->mesh,i-1,j,k);
        	row[1] = dt_rhodx2 * AE(i,j,k);
        	row[2] = dt_rhodx2 * AE(i-1,j,k);
        
        	row_idx[3] = mesh_index(solver->mesh,i,j+1,k);
        	row_idx[4] = mesh_index(solver->mesh,i,j-1,k);
        	row[3] = dt_rhody2 * AN(i,j,k);
        	row[4] = dt_rhody2 * AN(i,j-1,k);
        	
        	row_idx[5] = mesh_index(solver->mesh,i,j,k+1);
        	row_idx[6] = mesh_index(solver->mesh,i,j,k-1);
        	row[5] = dt_rhodz2 * AT(i,j,k);
        	row[6] = dt_rhodz2 * AT(i,j,k-1);
        	
        	nidx = mesh_index(solver->mesh,i,j,k);
    			ierr   = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
    			
    			rhs  = (AE(i,j,k) * U(i,j,k) - AE(i-1,j,k) * U(i-1,j,k)) * RDX;
    			rhs += (AN(i,j,k) * V(i,j,k) - AN(i,j-1,k) * V(i,j-1,k)) * RDY;
    			rhs += (AT(i,j,k) * W(i,j,k) - AT(i,j,k-1) * W(i,j,k-1)) * RDZ;     
    			   	
        	/* de-foaming */
          if(VOF(i,j,k) < (1-emf)) {
             /* uncomment to not de-foam next to boundaries *
            if(!(FV(i+1,j,k) < emf || FV(i-1,j,k) < emf ||
                 FV(i,j+1,k) < emf || FV(i,j-1,k) < emf ||
                 FV(i,j,k+1) < emf || FV(i,j,k-1) < emf)) {  */
             	rhs += min(solver->epsi * solver->rho, 
                                    0.1 * (1.0 - VOF(i,j,k)) / solver->delt) / 10;
            /* } */
          }
          
    			ierr   = VecSetValue(b,nidx,rhs,INSERT_VALUES);CHKERRQ(ierr);
        }
        
        
      }
    }
  }
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  
    /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-5,solver->epsi,PETSC_DEFAULT,solver->niter);CHKERRQ(ierr);
  
  /*
     Solve linear system
  */
	ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
	/* ierr = KSPMonitorSet(ksp, KSPMonitorDefault, NULL, NULL); Uncomment for verbose residuals */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp, &solver->resimax);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp, &iter);CHKERRQ(ierr);
  solver->iter = iter;

  /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  /* ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); */
  
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {      
      
        
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
        
        n = mesh_index(solver->mesh,i,j,k);
        ierr = VecGetValues(x,1,&n,&delp); CHKERRQ(ierr);
        
        P(i,j,k) += delp;

        if(AE(i,j,k) > emf && i < IMAX-1)
          U(i,j,k)=U(i,j,k) + solver->delt* RDX * delp / (solver->rho /* AE(i,j,k) */);
        if(AE(i-1,j,k) > emf && i > 0)
          U(i-1,j,k)=U(i-1,j,k) - solver->delt* RDX * delp / (solver->rho /* AE(i-1,j,k) */);
        if(AN(i,j,k) > emf && j < JMAX-1)
          V(i,j,k)=V(i,j,k) + solver->delt * RDY * delp / (solver->rho /* AN(i,j,k) */);
        if(AN(i,j-1,k) > emf && j > 0)
          V(i,j-1,k)=V(i,j-1,k) - solver->delt * RDY * delp / (solver->rho /* AN(i,j-1,k) */);
        if(AT(i,j,k) > emf && k < KMAX-1)
          W(i,j,k)=W(i,j,k) + solver->delt * RDZ * delp / (solver->rho /* AT(i,j,k) */);
        if(AT(i,j,k-1) > emf && k > 0)
          W(i,j,k-1)=W(i,j,k-1) - solver->delt * RDZ * delp / (solver->rho /* AT(i,j,k-1) */);        
               
        
  		}
  	}
  }
  ierr = VecDestroy(&x);CHKERRQ(ierr); 
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  
  solver->p_flag = 0;
  return 0;

}

