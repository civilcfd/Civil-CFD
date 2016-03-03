/* vof_baffles.c
 *
 * implement baffles
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "vtk.h"
#include "vof.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "vof_boundary.h"
#include "vof_baffles.h"

#include "vof_macros.h"

int vof_baffles(struct solver_data *solver) {
  int x;
  struct baffle_data *baffle;
#define emf solver->emf  

  for(x=0; x < 3; x++) {
    for(baffle = solver->mesh->baffles[x]; baffle != NULL; baffle = baffle->next) {    
        switch(baffle->type) {
        case flow:
          baffle_flow(solver, x, baffle->extent_a[0], baffle->extent_a[1], 
                                     baffle->extent_b[0], baffle->extent_b[1], 
                                     &(baffle->value), baffle->pos);
          break;
        case slip:
          baffle_slip(solver, x, baffle->extent_a[0], baffle->extent_a[1], 
                                     baffle->extent_b[0], baffle->extent_b[1], 
                                     baffle->value, baffle->pos);
          break;     
        case k:
          baffle_k(solver, x, baffle->extent_a[0], baffle->extent_a[1], 
                                     baffle->extent_b[0], baffle->extent_b[1], 
                                     baffle->value, baffle->pos);
          break;         
        }
    }
  }

  return 0;
#undef emf
}

int baffle_k(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, long int pos) {

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  double delp;
  double g;
  double sgn = 1.0;
  
  baffle_setup(solver, x, pos, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {   
        switch(x) {
        case 0:
          delp = P(i+1,j,k) - P(i,j,k);
          sgn = delp * sgn / fabs(delp);
          U(i,j,k) = -1.0 * sgn * sqrt((delp * 2) / (solver->rho * value));
          break;
        case 1:
          delp = P(i,j,k) - P(i,j+1,k);
          sgn = delp * sgn / fabs(delp);
          V(i,j,k) = -1.0 * sgn * sqrt((delp * 2) / (solver->rho * value));
          break;
        case 2:
          delp = P(i,j,k+1) - P(i,j,k);
          sgn = delp * sgn / fabs(delp);
          W(i,j,k) = -1.0 * sgn * sqrt((delp * 2) / (solver->rho * value));
          break;
        }          
      }
    }
  }

  return 0;
}

int baffle_slip(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, long int pos) {

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  
  baffle_setup(solver, x, pos, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {   
        switch(x) {
        case 0:
          U(i,j,k) = 0;
          break;
        case 1:
          V(i,j,k) = 0;
          break;
        case 2:
          W(i,j,k) = 0;
          break;
        }        
      }
    }
  }

  return 0;
}

int baffle_flow(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double *value, long int pos) {

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  
  baffle_setup(solver, x, pos, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {   
        switch(x) {
        case 0:
          *value = U(i,j,k) * (VOF(i,j,k) + VOF(i+1,j,k)) / 2 * AE(i,j,k) + *value;
          break;
        case 1:
          *value = W(i,j,k) * (VOF(i,j,k) + VOF(i,j+1,k)) / 2 * AN(i,j,k) + *value;
          break;
        case 2:
          *value = V(i,j,k) * (VOF(i,j,k) + VOF(i,j,k+1)) / 2 * AT(i,j,k) + *value;
          break;
        }        
      }
    }
  }

  return 0;
}

int baffle_setup(struct solver_data *solver, int x, long int pos, 
                           long int *imin, long int *jmin, long int *kmin, 
                           long int *imax, long int *jmax, long int *kmax,
                           double min_1, double min_2, double max_1, double max_2) {
  switch(x) {
  case 0:
    *imin = pos;
    *imax = pos;
    break;
  case 1:
    *jmin = pos;
    *jmax = pos;
    break;
  case 2:
    *kmin = pos;
    *kmax = pos;
    break;
  }

  switch(x) {
  case 0:
    *jmin = min_1;
    *jmax = min(JMAX-2,max_1);
    *kmin = min_2;
    *kmax = min(KMAX-1,max_2);    
    break;
  case 1
    *imin = min_1;
    *imax = min(IMAX-2,max_1);
    *kmin = min_2;
    *kmax = min(KMAX-2,max_2);   
    break;    
  case 2
    *imin = min_1;
    *imax = min(IMAX-2,max_1);
    *jmin = min_2;
    *jmax = min(JMAX-2,max_2);   
    break;   
  }
  return 0;
}