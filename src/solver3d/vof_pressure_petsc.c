/* vof_pressure_gmres.c
 *
 * implement pressure iteration
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>

#include "vtk.h"
#include "solver.h"
#include "solver_mpi.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "track.h"
#include "vof_boundary.h"
#include "vof_mpi.h"

#include "vof_macros.h"

extern struct mesh_data *mesh_n; /* describes mesh at previous timestep for explicit calcs */

int vof_pressure_gmres_assemble(struct solver_data *solver, Mat A, Vec b);
int vof_pressure_gmres_update(struct solver_data *solver, Mat A, Vec b);

/* these boundary functions handle dirchelet boundaries for velocity */
int vof_pressure_gmres_boundary(struct solver_data *solver, Mat A, Vec b);
int vof_pressure_gmres_write(Mat A, double timestep);
int vof_pressure_gmres_boundary_edges(struct solver_data *solver, Mat A, Vec b);

int vof_pressure_gmres_boundary_edges(struct solver_data *solver, Mat A, Vec b) {
  PetscInt i,j,k,nidx;
	PetscErrorCode ierr;
  int offset;
  if(!solver->rank) offset = 0;
  else offset=1;

  for(j=0; j<JMAX-1; j++) {

    if(ISTART == 0) {
      nidx = mesh_index(solver->mesh,0,j,0);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

      nidx = mesh_index(solver->mesh,0,j,KMAX-1);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
    }

    if(ISTART + IRANGE == IMAX) {
      nidx = mesh_index(solver->mesh,IMAX-1,j,0);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

      nidx = mesh_index(solver->mesh,IMAX-1,j,KMAX-1);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  for(i=offset; i<IRANGE; i++) {
    nidx = mesh_index(solver->mesh,i+ISTART,0,0);
    ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

    nidx = mesh_index(solver->mesh,i+ISTART,JMAX-1,0);
    ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

    nidx = mesh_index(solver->mesh,i+ISTART,0,KMAX-1);
    ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

    nidx = mesh_index(solver->mesh,i+ISTART,JMAX-1,KMAX-1);
    ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
  }

  for(k=0; k<KMAX-1; k++) {

    if(ISTART == 0) {
      nidx = mesh_index(solver->mesh,0,0,k);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

      nidx = mesh_index(solver->mesh,0,JMAX-1,k);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
    }

    if(ISTART + IRANGE == IMAX) {
      nidx = mesh_index(solver->mesh,IMAX-1,0,k);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

      nidx = mesh_index(solver->mesh,IMAX-1,JMAX-1,k);
      ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  return 0;
}

int vof_pressure_gmres_boundary(struct solver_data *solver, Mat A, Vec b) {
  PetscInt i,j,k,nidx,nlmn;
	PetscErrorCode ierr;
	enum special_boundaries sb;
	enum wall_boundaries wb;
  int offset;
  double *vec;
	
  if(!solver->rank) offset = 0;
  else offset=1;
  ierr = VecGetArray(b,&vec); CHKERRQ(ierr);

	for(j=1; j<JMAX-1; j++) {
    for(k=1; k<KMAX-1; k++) {

      /* west boundary */
      if(ISTART == 0) {
        wb = solver->mesh->wb[0];
        sb = vof_boundaries_check_inside_sb(solver, j, k, 0);
        if(sb == fixed_velocity || sb == mass_outflow || 
          (sb == wall && (wb == slip || wb == no_slip) ) ) {   
          
          if(FV(1,j,k) < solver->emf || AE(0,j,k) < solver->emf) {
            nidx = mesh_index(solver->mesh,0,j,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
          }
          
          nidx = mesh_index(solver->mesh,0,j,k);
          nlmn = mesh_index(solver->mesh,1,j,k);
          
          if(N_VOF(1,j,k) != 0) {
            /* explicit zero out since this could change */
            ierr   = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr   = MatSetValue(A,nidx,nlmn,0,INSERT_VALUES);CHKERRQ(ierr);
            vec[mesh_index(solver->mesh,0,j,k)] = 0;
            //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          } else {           
            ierr   = MatSetValue(A,nidx,nidx,1/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
            ierr   = MatSetValue(A,nidx,nlmn,-1.0/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
            //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
            vec[mesh_index(solver->mesh,0,j,k)] = 0;
          }
        } else {
            nidx = mesh_index(solver->mesh,0,j,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
        }
      }
    }
  }
	for(j=1; j<JMAX-1; j++) {
    for(k=1; k<KMAX-1; k++) {

      /* east boundary */
      if(ISTART + IRANGE == IMAX) {
        wb = solver->mesh->wb[1];
        sb = vof_boundaries_check_inside_sb(solver, j, k, 1);
        if(sb == fixed_velocity || sb == mass_outflow || 
          (sb == wall && (wb == slip || wb == no_slip) ) ) {     
          
          if(FV(IRANGE-2,j,k) < solver->emf || AE(IRANGE-2,j,k) < solver->emf) {
            nidx = mesh_index(solver->mesh,IMAX-1,j,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
          } 
          
          nidx = mesh_index(solver->mesh,IMAX-1,j,k);
          nlmn = mesh_index(solver->mesh,IMAX-2,j,k);
          
          if(N_VOF(IRANGE-2,j,k) != 0) {
            /* explicit zero out since this could change */
            ierr   = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr   = MatSetValue(A,nidx,nlmn,0,INSERT_VALUES);CHKERRQ(ierr);
            //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
            vec[mesh_index(solver->mesh,IRANGE-1,j,k)] = 0;
          } else {        	
            ierr   = MatSetValue(A,nidx,nidx,1/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
            ierr   = MatSetValue(A,nidx,nlmn,-1.0/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
            //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
            vec[mesh_index(solver->mesh,IRANGE-1,j,k)] = 0;
          }
        } else {
            nidx = mesh_index(solver->mesh,IMAX-1,j,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
        }      
      }
      
    }
  }
  
  /* boundaries on y-axis */
  for(i=1; i<IRANGE-1; i++) {
    for(k=1; k<KMAX-1; k++) {
         
      /* south boundary */
      wb = solver->mesh->wb[2];
      sb = vof_boundaries_check_inside_sb(solver, i, k, 2);
      if(sb == fixed_velocity || sb == mass_outflow || 
      	 (sb == wall && (wb == slip || wb == no_slip) ) ) {    
      	 
      	if(FV(i,1,k) < solver->emf || AN(i,0,k) < solver->emf) {
            nidx = mesh_index(solver->mesh,i+ISTART,0,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
          }
        
      	nidx = mesh_index(solver->mesh,i+ISTART,0,k);
      	nlmn = mesh_index(solver->mesh,i+ISTART,1,k);
        
        if(N_VOF(i,1,k) != 0) {
          /* explicit zero out since this could change */
          ierr   = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,0,INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,0,k)] = 0;
        } else {           
          ierr   = MatSetValue(A,nidx,nidx,1/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,-1.0/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,0,k)] = 0;
        }
      } else {
            nidx = mesh_index(solver->mesh,i+ISTART,0,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
        }
    }
  }
  for(i=1; i<IRANGE-1; i++) {
    for(k=1; k<KMAX-1; k++) {    
 
      /* north boundary */
      wb = solver->mesh->wb[3];
      sb = vof_boundaries_check_inside_sb(solver, i, k, 3);
      if(sb == fixed_velocity || sb == mass_outflow || 
      	 (sb == wall && (wb == slip || wb == no_slip) ) ) {   
      	 
      	if(FV(i,JMAX-2,k) < solver->emf || AN(i,JMAX-2,k) < solver->emf) {
            nidx = mesh_index(solver->mesh,i+ISTART,JMAX-1,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
          }
        
      	nidx = mesh_index(solver->mesh,i+ISTART,JMAX-1,k);
      	nlmn = mesh_index(solver->mesh,i+ISTART,JMAX-2,k);
        
        if(N_VOF(i,JMAX-2,k) != 0) {
          /* explicit zero out since this could change */
          ierr   = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,0,INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,JMAX-1,k)] = 0;
        } else {        	
          ierr   = MatSetValue(A,nidx,nidx,1/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,-1.0/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,JMAX-1,k)] = 0;
        }
   
      } else {
            nidx = mesh_index(solver->mesh,i+ISTART,JMAX-1,k);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
        }        
    }
  }
  
    
  /* boundaries on z-axis */
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
   	
      /* bottom boundary */
      wb = solver->mesh->wb[4];
      sb = vof_boundaries_check_inside_sb(solver, i, j, 4);
      if(sb == fixed_velocity || sb == mass_outflow || 
      	 (sb == wall && (wb == slip || wb == no_slip) ) ) {  
      	 
      	if(FV(i,j,1) < solver->emf || AT(i,j,0) < solver->emf) {
            nidx = mesh_index(solver->mesh,i+ISTART,j,0);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
          }  
        
      	nidx = mesh_index(solver->mesh,i+ISTART,j,0);
      	nlmn = mesh_index(solver->mesh,i+ISTART,j,1);
        
        if(N_VOF(i,j,1) != 0) {
          /* explicit zero out since this could change */
          ierr   = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,0,INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,0)] = 0;
        } else {           
          ierr   = MatSetValue(A,nidx,nidx,1/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,-1.0/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,0)] = 0;
        }
            
      } else {
            nidx = mesh_index(solver->mesh,i+ISTART,j,0);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
        } 
    }
  }
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {

 
      /* top boundary */
      wb = solver->mesh->wb[5];
      sb = vof_boundaries_check_inside_sb(solver, i, j, 5);
      if(sb == fixed_velocity || sb == mass_outflow || 
      	 (sb == wall && (wb == slip || wb == no_slip) ) ) {    
      	 
      	if(FV(i,j,KMAX-2) < solver->emf || AT(i,j,KMAX-2) < solver->emf) {
            nidx = mesh_index(solver->mesh,i+ISTART,j,KMAX-1);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
          }    
        
      	nidx = mesh_index(solver->mesh,i+ISTART,j,KMAX-1);
      	nlmn = mesh_index(solver->mesh,i+ISTART,j,KMAX-2);
        
        if(N_VOF(i,j,KMAX-2) != 0) {
          /* explicit zero out since this could change */
          ierr   = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,0,INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,KMAX-1)] = 0;
        } else {        	
          ierr   = MatSetValue(A,nidx,nidx,1/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          ierr   = MatSetValue(A,nidx,nlmn,-1.0/(solver->rho * solver->delt),INSERT_VALUES);CHKERRQ(ierr);
          //ierr   = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,KMAX-1)] = 0;
        }
      } else {
            nidx = mesh_index(solver->mesh,i+ISTART,j,KMAX-1);
            ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

            continue;
        }        
         	
   	} 
  }   


  VecRestoreArray(b,&vec);
  
  return 0;
}

/* this function assembles the matrix for the first run */
int vof_pressure_gmres_assemble(struct solver_data *solver, Mat A, Vec b) {
	PetscInt i,j,k,l,m,n,a;
	PetscInt  row_idx[7];
	PetscInt  nidx, ridx;
	PetscScalar rhs, pijk, dpijk;
	PetscScalar r_rhodx2, r_rhody2, r_rhodz2;
	PetscScalar row[7];
	PetscErrorCode ierr;
  double mpeta;
	
	r_rhodx2 = 1/solver->rho * 1/pow(DELX,2);
	r_rhody2 = 1/solver->rho * 1/pow(DELY,2);
	r_rhodz2 = 1/solver->rho * 1/pow(DELZ,2);
 /* Assemble matrix */
 #define emf solver->emf

  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {      
      
        
        if(FV(i,j,k)<emf) {
          nidx = mesh_index(solver->mesh,i+ISTART,j,k);
          ierr = MatSetValue(A,nidx,nidx,1,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);

          continue;
        }
                       
        row_idx[0] = mesh_index(solver->mesh,i+ISTART,j,k);   
        row_idx[1] = mesh_index(solver->mesh,i+1+ISTART,j,k);  
        row_idx[2] = mesh_index(solver->mesh,i-1+ISTART,j,k);  
        row_idx[3] = mesh_index(solver->mesh,i+ISTART,j+1,k);
        row_idx[4] = mesh_index(solver->mesh,i+ISTART,j-1,k);
        row_idx[5] = mesh_index(solver->mesh,i+ISTART,j,k+1);
        row_idx[6] = mesh_index(solver->mesh,i+ISTART,j,k-1);
               
        if(N_VOF(i,j,k) != 0) {
        
          l = i;
          m = j;
          n = k;
          for(a=0;a<7;a++) row[a] = 0;
          
          switch(N_VOF(i,j,k)) {
          case east:
            l=i+1;
            ridx = 1;
            break;
          case west:
            l=i-1;
            ridx = 2;
            break;
          case north:
            m=j+1;
            ridx = 3;
            break;
          case south:
            m=j-1;
            ridx = 4;
            break;
          case top:
            n=k+1;
            ridx = 5;
            break;
          case bottom:
            n=k-1;
            ridx = 6;
            break;
          case none:
          default: 
            row[0] = 1;
            nidx = mesh_index(solver->mesh,i+ISTART,j,k);
            ierr = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
            continue;
          }
          
          nidx = mesh_index(solver->mesh,i+ISTART,j,k);
          
          if(N_VOF(l,m,n) != 0) {
            row[0] = 1 / (solver->rho * solver->delt);
            dpijk = 0 - P(i,j,k);
            dpijk /= (solver->rho * solver->delt);
            ierr = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
            ierr = VecSetValue(b,nidx,dpijk,INSERT_VALUES);CHKERRQ(ierr);
          	continue;
          }
          
          mpeta = 1.0 - surface_interpolate(solver,i,j,k);
          pijk  = mpeta * P(l,m,n);
          dpijk = pijk - P(i,j,k);
          dpijk /= (solver->rho * solver->delt);
          
          row[0] = 1 / (solver->rho * solver->delt);
          row[ridx] = -1.0 * mpeta / (solver->rho * solver->delt);
          
    			ierr = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecSetValue(b,nidx,dpijk,INSERT_VALUES);CHKERRQ(ierr);
          
          
        } else {
        /* interior non-void cell */
        
        	dpijk = r_rhodx2 * (AE(i,j,k) + AE(i-1,j,k)) + 
        				  r_rhody2 * (AN(i,j,k) + AN(i,j-1,k)) +
        				  r_rhodz2 * (AT(i,j,k) + AT(i,j,k-1));

        	
        	row[0] = dpijk * -1.0;
        	
        	row[1] = r_rhodx2 * AE(i,j,k);
        	row[2] = r_rhodx2 * AE(i-1,j,k);
        
        	row[3] = r_rhody2 * AN(i,j,k);
        	row[4] = r_rhody2 * AN(i,j-1,k);
        	
        	row[5] = r_rhodz2 * AT(i,j,k);
        	row[6] = r_rhodz2 * AT(i,j,k-1);
        	
        	nidx = mesh_index(solver->mesh,i+ISTART,j,k);
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
          
          rhs /= solver->delt;
          
    			ierr   = VecSetValue(b,nidx,rhs,INSERT_VALUES);CHKERRQ(ierr);
        }
        
        
      }
    }
  }
  
  return 0;
}

/* update the matrix */
int vof_pressure_gmres_update(struct solver_data *solver, Mat A, Vec b) {
	PetscInt i,j,k,l,m,n,a;
	PetscInt  row_idx[7];
	PetscInt  nidx, ridx;
	PetscScalar rhs, pijk, dpijk, mpeta;
	PetscScalar r_rhodx2, r_rhody2, r_rhodz2;
	PetscScalar row[7];
	PetscErrorCode ierr;
  int offset;
  double *vec;
	
  if(!solver->rank) offset = 0;
  else offset=1;
  ierr = VecGetArray(b,&vec); CHKERRQ(ierr);
	
	r_rhodx2 = 1/solver->rho * 1/pow(DELX,2);
	r_rhody2 = 1/solver->rho * 1/pow(DELY,2);
	r_rhodz2 = 1/solver->rho * 1/pow(DELZ,2);
 /* Assemble matrix */
 #define emf solver->emf
  
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {      
      
        
        if(FV(i,j,k)<emf) continue;
                       
        row_idx[0] = mesh_index(solver->mesh,i+ISTART,j,k);   
        row_idx[1] = mesh_index(solver->mesh,i+1+ISTART,j,k);  
        row_idx[2] = mesh_index(solver->mesh,i-1+ISTART,j,k);  
        row_idx[3] = mesh_index(solver->mesh,i+ISTART,j+1,k);
        row_idx[4] = mesh_index(solver->mesh,i+ISTART,j-1,k);
        row_idx[5] = mesh_index(solver->mesh,i+ISTART,j,k+1);
        row_idx[6] = mesh_index(solver->mesh,i+ISTART,j,k-1);
               
        if(N_VOF(i,j,k) != 0) {
        
          l = i;
          m = j;
          n = k;
          for(a=0;a<7;a++) row[a] = 0;
          
          switch(N_VOF(i,j,k)) {
          case east:
            l=i+1;
            ridx = 1;
            break;
          case west:
            l=i-1;
            ridx = 2;
            break;
          case north:
            m=j+1;
            ridx = 3;
            break;
          case south:
            m=j-1;
            ridx = 4;
            break;
          case top:
            n=k+1;
            ridx = 5;
            break;
          case bottom:
            n=k-1;
            ridx = 6;
            break;
          case none:
          default: 
            row[0] = 1;
            if(N_VOF_N(i,j,k) == N_VOF(i,j,k)) continue;
            nidx = mesh_index(solver->mesh,i+ISTART,j,k);
            ierr = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
            //ierr = VecSetValue(b,nidx,0,INSERT_VALUES);CHKERRQ(ierr);
            vec[mesh_index(solver->mesh,i-offset,j,k)] = 0;
            continue;
          }
          
          nidx = mesh_index(solver->mesh,i+ISTART,j,k);
          
          if(N_VOF(l,m,n) != 0) {
            row[0] = 1 / (solver->rho * solver->delt);
            dpijk = 0 - P(i,j,k);
            dpijk /= (solver->rho * solver->delt);
            ierr = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
            //ierr = VecSetValue(b,nidx,dpijk,INSERT_VALUES);CHKERRQ(ierr);
            vec[mesh_index(solver->mesh,i-offset,j,k)] = dpijk;
          	continue;
          }
          
          mpeta = 1 - surface_interpolate(solver,i,j,k);
          pijk  = mpeta * P(l,m,n);
          dpijk = pijk - P(i,j,k);
          dpijk /= (solver->rho * solver->delt);
          
          row[0] = 1 / (solver->rho * solver->delt);
          row[ridx] = -1.0 * mpeta / (solver->rho * solver->delt);
          
    			ierr = MatSetValues(A,1,&nidx,7,row_idx,row,INSERT_VALUES);CHKERRQ(ierr);
          //ierr = VecSetValue(b,nidx,dpijk,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,k)] = dpijk;
          
          
        }
        else if(N_VOF_N(i,j,k) != 0) {
        /* interior non-void cell */
        
        	dpijk = r_rhodx2 * (AE(i,j,k) + AE(i-1,j,k)) + 
        				  r_rhody2 * (AN(i,j,k) + AN(i,j-1,k)) +
        				  r_rhodz2 * (AT(i,j,k) + AT(i,j,k-1));

        	
        	row[0] = dpijk * -1.0;
        	
        	row[1] = r_rhodx2 * AE(i,j,k);
        	row[2] = r_rhodx2 * AE(i-1,j,k);
        
        	row[3] = r_rhody2 * AN(i,j,k);
        	row[4] = r_rhody2 * AN(i,j-1,k);
        	
        	row[5] = r_rhodz2 * AT(i,j,k);
        	row[6] = r_rhodz2 * AT(i,j,k-1);
        	
        	nidx = mesh_index(solver->mesh,i+ISTART,j,k);
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
          
          rhs /= solver->delt;
          
    			//ierr   = VecSetValue(b,nidx,rhs,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,k)] = rhs;
        }
        else {
    			
        	nidx = mesh_index(solver->mesh,i+ISTART,j,k);
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
          
          rhs /= solver->delt;
          
    			//ierr   = VecSetValue(b,nidx,rhs,INSERT_VALUES);CHKERRQ(ierr);
          vec[mesh_index(solver->mesh,i-offset,j,k)] = rhs;
        }
        
        
      }
    }
  }
  VecRestoreArray(b,&vec);
  
  return 0;
}

int vof_pressure_gmres(struct solver_data *solver) {
	PetscInt i,j,k,n,size;
	PetscInt	iter;
	PetscErrorCode ierr;
	PetscScalar delp;
	static Vec x, b;
	static Mat A;
  static KSP  ksp;         /* linear solver context */
  static PC  pc;           /* preconditioner context */
  static int initialize = 0;
  
	size = IMAX * JMAX * KMAX;
	
	if(!initialize) {
		PetscInitialize(NULL, NULL, NULL, NULL);	
    
    ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,size,size,7,NULL,&A); 
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,size);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
    
    vof_pressure_gmres_assemble(solver, A, b);
		initialize = 1;
	}
  else {
    vof_pressure_gmres_update(solver, A, b);  
  }
  
  vof_pressure_gmres_boundary(solver, A, b); 
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* vof_pressure_gmres_write(A,solver->t); */
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  
  
  /*
     Solve linear system
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-3,solver->epsi/(solver->rho * solver->delt),PETSC_DEFAULT,solver->niter);CHKERRQ(ierr);
  //ierr = KSPSetTolerances(ksp,1e-3,PETSC_DEFAULT,PETSC_DEFAULT,solver->niter);CHKERRQ(ierr);
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
  
  /* ierr = VecDestroy(&x);CHKERRQ(ierr); 
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);  */
  
  solver->p_flag = 0;
  return 0;

}

int vof_pressure_gmres_mpi(struct solver_data *solver) {
	PetscInt i,j,k,size;
	PetscInt	iter;
	PetscErrorCode ierr;
	PetscScalar delp;
	static Vec x, b;
	static Mat A;
  static KSP  ksp;         /* linear solver context */
  static PC  pc;           /* preconditioner context */
  KSP *subksp;
  PC subpc;
  static int initialize = 0;
  int offset = 1;
  int range;
  int Istart, Iend;
  int Cstart, Cend;
  static int nlocal, first;
  int *blks;
  double *results;
  IS diag_zeros;
  
  range = IRANGE-1;
  if(solver->rank < solver->size - 1 && solver->rank > 0) {
    range--;
  }

	if(!initialize) {
    //PetscLogBegin();
	  size = IMAX * JMAX * KMAX;
    
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,range * JMAX * KMAX,range * JMAX * KMAX,size,size);CHKERRQ(ierr);
    ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatSetOption(A,MAT_IGNORE_OFF_PROC_ENTRIES,PETSC_TRUE);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A, 7, PETSC_NULL, 7, PETSC_NULL);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,7,NULL); CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);
    ierr = MatGetOwnershipRangeColumn(A,&Cstart,&Cend);
    ierr = VecCreateMPI(PETSC_COMM_WORLD,range * JMAX * KMAX,size,&x);CHKERRQ(ierr);
    ierr = VecSetOption(x, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD,range * JMAX * KMAX,size,&b);CHKERRQ(ierr);
    ierr = VecSetOption(b, VEC_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);CHKERRQ(ierr);
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    
    vof_pressure_gmres_assemble(solver, A, b);
    vof_pressure_gmres_boundary_edges(solver, A, b); 
    vof_pressure_gmres_boundary(solver, A, b); 
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		initialize = 1;
    
#ifdef DEBUG
    MatFindZeroDiagonals(A, &diag_zeros);
    printf("Checking for non-zero diagonals:\n");
    ISView(diag_zeros, PETSC_VIEWER_STDOUT_WORLD);
#endif


    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    ierr = PCSetType(pc,PCBJACOBI);CHKERRQ(ierr);

    //PCASMSetLocalSubdomains(pc, range, NULL, NULL);
    //PCASMSetOverlap(pc, KMAX);
    
    ierr = PetscMalloc1(range,&blks);CHKERRQ(ierr);
    for (i=0; i<range; i++) blks[i] = JMAX * KMAX;
    ierr = PCBJacobiSetLocalBlocks(pc,range,blks);CHKERRQ(ierr);
    ierr = PetscFree(blks);

    ierr = KSPSetUp(ksp); CHKERRQ(ierr);
    //ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp,.001,solver->epsi/(solver->rho * solver->delt),PETSC_DEFAULT,solver->niter);CHKERRQ(ierr);
    //ierr = PCASMGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
    ierr = PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);CHKERRQ(ierr);
    for (i=0; i<nlocal; i++) {
      ierr = KSPGetPC(subksp[i],&subpc);CHKERRQ(ierr);
      ierr = PCSetType(subpc,PCSOR);CHKERRQ(ierr);
      ierr = KSPSetType(subksp[i],KSPPREONLY);CHKERRQ(ierr);
      ierr = KSPSetTolerances(subksp[i],.001,solver->epsi/(solver->rho * solver->delt),PETSC_DEFAULT,solver->niter);CHKERRQ(ierr);
    }

    printf("Built matrix with range: %d to %d on proc %d\n",Istart, Iend, solver->rank);

	}
  else {
    vof_pressure_gmres_update(solver, A, b);  
    vof_pressure_gmres_boundary(solver, A, b); 
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  

  /* vof_pressure_gmres_write(A,solver->t); */
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  
  
  /*
     Solve linear system
  */

	/* ierr = KSPMonitorSet(ksp, KSPMonitorDefault, NULL, NULL); */ 
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    
  ierr = KSPGetResidualNorm(ksp, &solver->resimax);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp, &iter);CHKERRQ(ierr);
  solver->iter = iter;
  
  if(!solver->rank) offset = 0;
  VecGetArray(x,&results);
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {      
      
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
        
        delp = results[mesh_index(solver->mesh,i-offset,j,k)];        
        P(i,j,k) += delp;

        if(AE(i,j,k) > emf && i+ISTART < IMAX-1)
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
  
  /* boundary edges */
  for(j=1; j<JMAX-1; j++) {
    for(k=1; k<KMAX-1; k++) {      
        
      i = IRANGE - 2;
      solver->mesh->delu_downstream[KMAX * j + k] = 0;
      if(FV(i,j,k)>emf && VOF(i,j,k) > emf) {

        /*n = mesh_index(solver->mesh,i+ISTART,j,k);
        ierr = VecGetValues(x,1,&n,&delp); CHKERRQ(ierr);*/
        delp = results[mesh_index(solver->mesh,i-offset,j,k)]; 
        
        if(AE(i,j,k) > emf && i+ISTART < IMAX-2)
          solver->mesh->delu_downstream[KMAX * j + k] = solver->delt* RDX * delp / (solver->rho /* AE(i,j,k) */);
      }

      i = 1;
      solver->mesh->delu_upstream[KMAX * j + k] = 0;
      if(FV(i,j,k)>emf && VOF(i,j,k) > emf) {

        /*n = mesh_index(solver->mesh,i+ISTART,j,k);
        ierr = VecGetValues(x,1,&n,&delp); CHKERRQ(ierr);*/
        delp = results[mesh_index(solver->mesh,i-offset,j,k)]; 
        
        if(AE(i-1,j,k) > emf && i+ISTART > 1)
          solver->mesh->delu_upstream[KMAX * j + k] = -1 * solver->delt* RDX * delp / (solver->rho /* AE(i-1,j,k) */);

      }     
    }
  }
  VecRestoreArray(x,&results);

  solver_sendrecv_delu(solver);

  for(j=1; j<JMAX-1; j++) {
    for(k=1; k<KMAX-1; k++) {      
      i = IRANGE - 2;
      U(i,j,k) += solver->mesh->delu_downstream[KMAX * j + k];

      i = 0;
      U(i,j,k) += solver->mesh->delu_upstream[KMAX * j + k];
    }
  }
  
  solver->p_flag = 0;
  return 0;

}

int vof_pressure_gmres_write(Mat A, double timestep) 
{
  PetscViewer viewer;
  char filename[256];

  sprintf(filename,"%4.3lfMat",timestep);

  PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
  MatView(A, viewer);
  PetscViewerDestroy(&viewer);

  return 0;
}