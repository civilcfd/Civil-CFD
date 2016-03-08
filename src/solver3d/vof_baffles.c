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

#ifdef _WIN32
#define M_PI 3.14
#endif

int vof_baffles_output(struct solver_data *solver) {
  struct baffle_data *baffle;
  int x;
  int flow_count = 0;
  int swirl_count = 0;
  
  for(x=0; x < 3; x++) {
    for(baffle = solver->mesh->baffles[x]; baffle != NULL; baffle = baffle->next) {    
        switch(baffle->type) {
        case flow:
          printf("Flow baffle %d : %.2lf L/s | axis %d position %ld  extent %ld %ld to %ld %ld\n", \
                 flow_count, baffle->value * 1000, x, baffle->pos, baffle->extent_a[0], baffle->extent_a[1], \
                 baffle->extent_b[0], baffle->extent_b[1]);
          flow_count++;
        break;
        case swirl_angle:
          printf("Swirl angle baffle %d : %.2lf degrees | axis %d position %ld  extent %ld %ld to %ld %ld\n", \
                 swirl_count, baffle->value, x, baffle->pos, baffle->extent_a[0], baffle->extent_a[1], \
                 baffle->extent_b[0], baffle->extent_b[1]);
          swirl_count++;
        break;
        case v_deviation:
          printf("Velocity deviation baffle %d : %.2lf %% | axis %d position %ld  extent %ld %ld to %ld %ld\n", \
                 swirl_count, baffle->value * 100, x, baffle->pos, baffle->extent_a[0], baffle->extent_a[1], \
                 baffle->extent_b[0], baffle->extent_b[1]);
          swirl_count++;
        break;
      }
    }
  }

  return 0;
}

int vof_baffles_write(struct solver_data *solver) {
  FILE *fp;
  struct baffle_data *baffle;
  int x, count=0, flow_count=0, swirl_count=0, v_count=0;
  
  if(solver->t < 0.000001) {
    fp = fopen("baffles.csv", "w"); 
    if(fp == NULL) {
      printf("error: cannot open baffles.csv for writing\n");
      return(1);
    }    
    
    for(x=0; x < 3; x++) {
      for(baffle = solver->mesh->baffles[x]; baffle != NULL; baffle = baffle->next) {    
          switch(baffle->type) {
          case flow:
            if(count > 0) {
              fprintf(fp, ", ");
            }
            fprintf(fp, "flow %d (L/s)", flow_count);
            flow_count++;
            count++;
            break;
          case swirl_angle:
            if(count > 0) {
              fprintf(fp, ", ");
            }
            fprintf(fp, "swirl angle %d (degrees)", swirl_count);
            swirl_count++;
            count++;
            break;    
          case v_deviation:
            if(count > 0) {
              fprintf(fp, ", ");
            }
            fprintf(fp, "velocity deviation %d (%%)", v_count);
            v_count++;
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
        case swirl_angle:
          if(count > 0) {
            fprintf(fp, ", ");
          }
          fprintf(fp, "%.2lf", \
                 baffle->value);
          count++;
          break;
        case v_deviation:
          if(count > 0) {
            fprintf(fp, ", ");
          }
          fprintf(fp, "%.2lf", \
                 baffle->value * 100);
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
        case swirl_angle:
          baffle_swirl(solver, x, baffle->extent_a[0], baffle->extent_a[1], 
                                     baffle->extent_b[0], baffle->extent_b[1], 
                                     &(baffle->value), baffle->pos);
          break;     
        case v_deviation:
          baffle_velocity_dev(solver, x, baffle->extent_a[0], baffle->extent_a[1], 
                                     baffle->extent_b[0], baffle->extent_b[1], 
                                     &(baffle->value), baffle->pos);
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
          if(fabs(U(i,j,k)) < solver->emf) continue;
          
          delp = P(i,j,k) - P(i+1,j,k);
          
          if(delp < 0.0) sgn = -1.0;
          v_prime = sgn * sqrt((fabs(delp) * 2) / (solver->rho * value));
                    
          if(v_prime * U(i,j,k) < solver->emf * -1.0) v_prime = 0; /* this shouldn't cause a direction change */
          if(fabs(v_prime) > fabs(U(i,j,k)) && v_prime * U(i,j,k) > solver->emf) continue; /* and never speed up */
          else U(i,j,k) = v_prime;

          break;
        case 1:
          if(AN(i,j,k) < solver->emf) continue;
          if(fabs(V(i,j,k)) < solver->emf) continue;
          
          delp = P(i,j,k) - P(i,j+1,k);
          
          if(delp < 0.0) sgn = -1.0;
          v_prime = sgn * sqrt((fabs(delp) * 2) / (solver->rho * value));
                    
          if(v_prime * V(i,j,k) < solver->emf * -1.0) v_prime = 0; /* this shouldn't cause a direction change */
          if(fabs(v_prime) > fabs(V(i,j,k)) && v_prime * V(i,j,k) > solver->emf) continue; /* and never speed up */
          else V(i,j,k) = v_prime;

          break;
        case 2:
          if(AT(i,j,k) < solver->emf) continue;
          if(fabs(W(i,j,k)) < solver->emf) continue;
          
          delp = P(i,j,k) - P(i,j,k+1);
          
          if(delp < 0.0) sgn = -1.0;
          v_prime = sgn * sqrt((fabs(delp) * 2) / (solver->rho * value));
                    
          if(v_prime * W(i,j,k) < solver->emf * -1.0) v_prime = 0; /* this shouldn't cause a direction change */
          if(fabs(v_prime) > fabs(W(i,j,k)) && v_prime * W(i,j,k) > solver->emf) continue; /* and never speed up */
          else W(i,j,k) = v_prime;

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

int baffle_swirl(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double *value, long int pos) {

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  double count;
  double swirl, u_ave, v_ave, w_ave;
  
  baffle_setup(solver, x, pos, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  *value = 0;
  
  imin = max(imin-1,0);
  jmin = max(jmin-1,0);
  kmin = max(kmin-1,0);
  switch(x) {
  case 0:
    imax = max(imax-1,0);
    break;
  case 1:
    jmax = max(jmax-1,0);
    break;
  case 2:
    kmax = max(kmax-1,0);
    break;
  }
  
  swirl = 0;
  count = 0;
  
  for(i=imin; i <= imax-1; i++) {
    for(j=jmin; j <= jmax-1; j++) {
      for(k=kmin; k <= kmax-1; k++) {   
        /* find cell-centered velocity */
        u_ave = (U(i,j,k) + U(i+1,j,k)) / 2;
        v_ave = (V(i,j,k) + V(i,j+1,k)) / 2;
        w_ave = (W(i,j,k) + W(i,j,k+1)) / 2;        
      
        switch(x) {
        case 0:
          if((AE(i,j,k) + AE(i+1,j,k)) < (1-solver->emf)) continue;
          swirl = swirl + atan(sqrt(pow(v_ave,2) + pow(w_ave,2)) / fabs(u_ave)) * 180 / M_PI;
          count = count + 1;
          break;
        case 1:
          if((AN(i,j,k) + AN(i,j+1,k)) < (1-solver->emf)) continue;
          swirl = swirl + atan(sqrt(pow(u_ave,2) + pow(w_ave,2)) / fabs(v_ave)) * 180 / M_PI;
          count = count + 1;
          break;
        case 2:
          if((AT(i,j,k) + AT(i+1,j,k)) < (1-solver->emf)) continue;
          swirl = swirl + atan(sqrt(pow(u_ave,2) + pow(v_ave,2)) / fabs(w_ave)) * 180 / M_PI;
          count = count + 1;
          break;
        }        
      }
    }
  }
  
  *value = swirl / count;
  
  return 0;
}

int baffle_velocity_dev(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double *value, long int pos) {

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  double count;
  double v_ave, v_dev;
  
  baffle_setup(solver, x, pos, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  *value = 0;
    
  v_ave = 0;
  v_dev = 0;
  count = 0;
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {         
      
        switch(x) {
        case 0:
          if(AE(i,j,k) > solver->emf) {
            v_ave += U(i,j,k);
            count = count+1;
          }
          break;
        case 1:
          if(AN(i,j,k) > solver->emf) {
            v_ave += V(i,j,k);
            count = count+1;
          }
          break;
        case 2:
          if(AT(i,j,k) > solver->emf) {
            v_ave += W(i,j,k);
            count = count+1;
          }
          break;
        }        
      }
    }
  }
  
  v_ave = fabs(v_ave / count);
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {         
      
        switch(x) {
        case 0:
          v_dev = max(v_dev, fabs(U(i,j,k) / v_ave));
          break;
        case 1:
          v_dev = max(v_dev, fabs(V(i,j,k) / v_ave));
          break;
        case 2:
          v_dev = max(v_dev, fabs(W(i,j,k) / v_ave));
          break;
        }        
      }
    }
  }
  
  *value = v_dev;
  
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