/*
 * vof_nvof.c: Calculate NVOF and interpolation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <omp.h>
#include <mpi.h>
#include <petscsys.h>
 
#include "vtk.h"
#include "vof_mpi.h"
#include "vof_boundary.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "track.h"
#include "vof_baffles.h"
#include "mesh_mpi.h"
#include "solver_mpi.h"

#include "vof_macros.h"

#define emf solver->emf
#define emf_c solver->emf_c

double surface_interpolate(struct solver_data *solver, long int i, long int j, long int k) {
  double interpolate=0;
  long int l,m,n;

  l = i; m = j; n = k;

  switch (N_VOF(i,j,k))
  {
  case west:
    l=i-1;
    break;

  case east:
    l=i+1;
    break;

  case south:
    m=j-1;
    break;

  case north:
    m=j+1;  
    break;

  case bottom:
    n=k-1;
    break;

  case top:
    n=k+1;
    break;

  case none:
  case empty:
    return 1;

  }

  if(FV(l,m,n) > emf) {
    interpolate = 1 / ((VOF(i,j,k) + VOF(l,m,n) * 0.5));    
  }
  else {
    interpolate = 1 / (0.5 + VOF(i,j,k));
  }

  if ((FV(l,m,n) >= emf) && (N_VOF(l,m,n)!=0.0)) interpolate=1.0;

  interpolate = min(interpolate , 2);
  interpolate = max(interpolate , 0);
            
  return interpolate;
}

int vof_mpi_nvof(struct solver_data *solver) {
  long int i,j,k;
  int norm[6][3] = { {  1, 0, 0 },
                     { -1, 0, 0 },
                     {  0, 1, 0 },
                     {  0,-1, 0 },
                     {  0, 0, 1 },
                     {  0, 0,-1 }};
  int g[6];
  int l,m,n,x,obs;
  double score[6] = { 0,0,0,0,0,0 };
  double mult[3][3] = { { 0.7, 0.9, 0.7 }, 
                        { 0.9, 1.0, 0.9 },
                        { 0.7, 0.9, 0.7 } };
  double top_score = 0;
  double lvof[6];

#define emf solver->emf
#define emf_c solver->emf_c

  g[0] = solver->gx;
  g[1] = solver->gx * -1;
  g[2] = solver->gy;
  g[3] = solver->gy * -1;
  g[4] = solver->gz;
  g[5] = solver->gz * -1;

  for(i=1; i<IRANGE-1; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
        N_VOF(i,j,k) = none;
        if(j==0 || k==0 || j==JMAX-1 || k==KMAX-1 || FV(i,j,k) == 0)
          N_VOF(i,j,k) = 0;
        else if(VOF(i+1,j,k) >= emf && VOF(i,j+1,k) >= emf 
          && VOF(i-1,j,k) >= emf && VOF(i,j-1,k) >= emf
          && VOF(i,j,k+1) >= emf && VOF(i,j,k-1) >= emf) {
          N_VOF(i,j,k) = 0;          
        }
      }
    }
  }

  if(ISTART ==0) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        N_VOF(0,j,k) = 0;
      }
    }
  }
  if(IRANGE+ISTART == IMAX) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        N_VOF(IRANGE-1,j,k) = 0;
      }
    }
  }
 
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        if(FV(i,j,k) == 0) {
          N_VOF(i,j,k) = 0;
          continue;
        }

        if(VOF(i,j,k) < emf) {
          N_VOF(i,j,k) = 8; 
          continue;
        }
  
        if(VOF(i,j,k) > emf_c) {
          N_VOF(i,j,k) = 0;
          continue;
        }


        /* iterate on all sides and check if we are bounded by either an obstacle or fluid for all cells, if so NVOF=0 */
        obs = 0;
        for(n=0; n<6; n++) {
          lvof[n] = VOF(i,j,k);
          score[n] = 0;

          switch(n) {
          case 0:
            if (AE(i,j,k) >emf) lvof[n] = VOF(i+norm[n][0], j+norm[n][1], k+norm[n][2]);
            break; 
          case 1:
            if (AE(i-1,j,k) >emf) lvof[n] = VOF(i+norm[n][0], j+norm[n][1], k+norm[n][2]);
            break; 
          case 2:
            if (AN(i,j,k) >emf) lvof[n] = VOF(i+norm[n][0], j+norm[n][1], k+norm[n][2]);
            break; 
          case 3:
            if (AN(i,j-1,k) >emf) lvof[n] = VOF(i+norm[n][0], j+norm[n][1], k+norm[n][2]);
            break; 
          case 4:
            if (AT(i,j,k) >emf) lvof[n] = VOF(i+norm[n][0], j+norm[n][1], k+norm[n][2]);
            break; 
          case 5:
            if (AT(i,j,k-1) >emf) lvof[n] = VOF(i+norm[n][0], j+norm[n][1], k+norm[n][2]);
            break; 
          default:
            continue;
          }
          
          if(FV(i+norm[n][0], j+norm[n][1], k+norm[n][2]) < emf && g[n] < emf) {
            obs++;
          }
          else if(lvof[n] > emf) {
            obs++;
          }
        }
        if(obs == n) N_VOF(i,j,k) = 0;

        if(N_VOF(i,j,k) == 0) continue;

        top_score = emf;
        for(x=0; x<6; x++) {
          l = norm[x][0]; m = norm[x][1]; n = norm[x][2];

          if(x < 2) {
            for(m=-1; m<=1; m++) {
              for(n=-1; n<=1; n++) {
                score[x] += VOF(l+i,m+j,n+k) * mult[m+1][n+1];
              }
            }

            if(score[x] > top_score && VOF(i+norm[x][0], j+norm[x][1], k+norm[x][2]) > emf) {
              if(FV(i+norm[x][0], j+norm[x][1], k+norm[x][2]) > emf || g[x] > emf) {
                N_VOF(i,j,k) = x+1;
                top_score = score[x];
              }
            }
          }
          else if(x < 4) {
            for(l=-1; l<=1; l++) {
              for(n=-1; n<=1; n++) {
                score[x] += VOF(l+i,m+j,n+k)  * mult[l+1][n+1];
              }
            }

            if(score[x] > top_score && VOF(i+norm[x][0], j+norm[x][1], k+norm[x][2]) > emf) {
              if(FV(i+norm[x][0], j+norm[x][1], k+norm[x][2]) > emf || g[x] > emf) {
                N_VOF(i,j,k) = x+1;
                top_score = score[x];
              }
            }
          }
          else {
            for(l=-1; l<=1; l++) {
              for(m=-1; m<=1; m++) {
                score[x] += VOF(l+i,m+j,n+k)  * mult[l+1][m+1];
              }
            }
            if(score[x] > top_score && VOF(i+norm[x][0], j+norm[x][1], k+norm[x][2]) > emf) {
              if(FV(i+norm[x][0], j+norm[x][1], k+norm[x][2]) > emf || g[x] > emf) {
                N_VOF(i,j,k) = x+1;
                top_score = score[x];
              }
            }
          }

        }
          

               
      }
    }
  }


  return 0;

#undef emf
#undef emf_c
}
