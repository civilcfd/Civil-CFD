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

int vof_baffles_output(struct solver_data *solver) {
  struct baffle_data *baffle;
  int x;
  int count = 0;
  
  for(x=0; x < 3; x++) {
    for(baffle = solver->mesh->baffles[x]; baffle != NULL; baffle = baffle->next) {    
        switch(baffle->type) {
        case flow:
          printf("Flow baffle %d, axis %d position %ld, extent %ld %ld to %ld %ld: %.2lf L/s\n", \
                 count, x, baffle->pos, baffle->extent_a[0], baffle->extent_a[1], \
                 baffle->extent_b[0], baffle->extent_b[1], baffle->value * 1000);
          count++;
        break;
      }
    }
  }

  return 0;
}

int vof_baffles_write(struct solver_data *solver) {
  FILE *fp;
  struct baffle_data *baffle;
  int x, count=0;
  
  if(solver->t < 0.000001) {
    fp = fopen("baffles.csv", "w"); 
    if(fp == NULL) {
      printf("error: cannot open tracking file\n");
      return(1);
    }    
    
    for(x=0; x < 3; x++) {
      for(baffle = solver->mesh->baffles[x]; baffle != NULL; baffle = baffle->next) {    
          switch(baffle->type) {
          case flow:
            if(count > 0) {
              fprintf(fp, ", ");
            }
            fprintf(fp, "%d", count);
            count++;
          break;
        }
      }
    }
    
    if(count>0) fprintf(fp, "\n");
    
  } else {
    fp = fopen("baffles.csv", "a"); 
  }

  if(fp == NULL) {
    printf("error: cannot open tracking file\n");
    return(1);
  }

  count = 0;
  for(x=0; x < 3; x++) {
    for(baffle = solver->mesh->baffles[x]; baffle != NULL; baffle = baffle->next) {    
        switch(baffle->type) {
        case flow:
          if(count > 0) {
            fprintf(fp, ", ");
          }
          fprintf(fp, "%.2lf", \
                 baffle->value * 1000);
          count++;
        break;
      }
    }
  }

  if(count>0) fprintf(fp,"\n");
  
  fclose(fp);

  return 0;
}

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
        case barrier:
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
  double delp, v_prime;
  double sgn = 1.0;
  
  baffle_setup(solver, x, pos, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  if(value < 0.01) {
    printf("Headloss baffle with k-factor < 0.01.  Ignoring.\n");
    return 0;
  }
  
  /* if(solver->p_flag != 0) return 0; */
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {   
        sgn = 1.0;
        switch(x) {
        case 0:
          if(AE(i,j,k) < solver->emf) continue;
          delp = P(i+1,j,k) - P(i,j,k);
          sgn = delp * sgn / fabs(delp);
          if(isnan(sgn)) {
            U(i,j,k) = 0;
            continue;
          }
          U(i,j,k) = sgn * sqrt((delp * 2) / (solver->rho * value));
          break;
        case 1:
          /* TODO: generalize */
          if(AN(i,j,k) < solver->emf) continue;
          if(fabs(V(i,j,k)) < solver->emf) continue;
          
          delp = P(i,j,k) - P(i,j+1,k);
          sgn = delp * sgn / fabs(delp);
          
          if(isnan(sgn)) {
            v_prime = 0;
          } else {
            v_prime = sgn * sqrt((fabs(delp) * 2) / (solver->rho * value));
          }
          
          if(fabs(v_prime) > fabs(V(i,j,k))) continue; /*
          if(fabs(V(i,j,k) - v_prime)/fabs(V(i,j,k)) > 0.1) V(i,j,k) = 0.9 * V(i,j,k); */
          else V(i,j,k) = v_prime;
          
          if(isnan(V(i,j,k)) || isnan(v_prime)) {
            printf("break\n");
          }
          
          break;
        case 2:
          if(AT(i,j,k) < solver->emf) continue;
          delp = P(i,j,k+1) - P(i,j,k);
          sgn = delp * sgn / fabs(delp);
          if(isnan(sgn)) {
            W(i,j,k) = 0;
            continue;
          }
          W(i,j,k) = sgn * sqrt((delp * 2) / (solver->rho * value));
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
  
  *value = 0;
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {   
        switch(x) {
        case 0:
          *value = DELY * DELZ * (U(i,j,k) * ((VOF(i,j,k) + VOF(i+1,j,k)) / 2) * AE(i,j,k)) + *value;
          break;
        case 1:
          *value = DELX * DELZ * (V(i,j,k) * ((VOF(i,j,k) + VOF(i,j+1,k)) / 2) * AN(i,j,k)) + *value;
          break;
        case 2:
          *value = DELX * DELY * (W(i,j,k) * ((VOF(i,j,k) + VOF(i,j,k+1)) / 2) * AT(i,j,k)) + *value;
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
    *kmax = min(KMAX-2,max_2);    
    break;
  case 1:
    *imin = min_1;
    *imax = min(IMAX-2,max_1);
    *kmin = min_2;
    *kmax = min(KMAX-2,max_2);   
    break;    
  case 2:
    *imin = min_1;
    *imax = min(IMAX-2,max_1);
    *jmin = min_2;
    *jmax = min(JMAX-2,max_2);   
    break;   
  }
  return 0;
}