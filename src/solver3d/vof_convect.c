
/* vof_convect.c
 *
 * fluid convection function
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

extern struct mesh_data *mesh_n; /* describes mesh at previous timestep for explicit calcs */

double calc_dVOF(struct solver_data *solver, long int i, long int j, long int k, int x);

double calc_dVOF(struct solver_data *solver, long int i, long int j, long int k, int x) {
  double v, del, A, rb, ra, rd, VOFdm, CF, dVOF;
  long int range, acceptor[3], acceptor_d[3], donor[3], donor_m[3];
  int surf_flag = 0;

  acceptor[0] = i; acceptor[1] = j; acceptor[2] = k;
  donor[0] = i; donor[1] = j; donor[2] = k;
  donor_m[0] = i; donor_m[1] = j; donor_m[2] = k;

  if(FV(i,j,k) < solver->emf) return 0;

  switch(x) {
    case 0:
      v = U(i,j,k);
      del = DELX;
      range = IRANGE;
      rb = AE(i,j,k);
      break;
    case 1:
      v = V(i,j,k);
      del = DELY;
      range = JMAX;
      rb = AN(i,j,k);
      break;
    case 2:
      v = W(i,j,k);
      del = DELZ;
      range = KMAX;
      rb = AT(i,j,k);
      break;
  }
  v *= solver->delt;

  if(fabs(v) > 0.5 * DELX) solver->vof_flag = 1;

  if(v > 0) {
    acceptor[x] += 1;
    donor_m[x] -= 1;
    if(donor_m[x] < 0) donor_m[x] = 0;
    
    switch(x) {
      case 0:
        A = AE(donor_m[0], j, k);
        break;
      case 1:
        A = AN(i, donor_m[1], k);
        break;
      case 2:
        A = AT(i, j, donor_m[2]);
        break;
    }
  }
  else {
    donor[x] += 1;
    donor_m[x] += 2;
    if(donor_m[x] > range - 1) donor_m[x] = range - 1;

    switch(x) {
      case 0:
        A = AE(donor_m[0]-1, j, k);
        break;
      case 1:
        A = AN(i, donor_m[1]-1, k);
        break;
      case 2:
        A = AT(i, j, donor_m[2]-1);
        break;    
    }
  }

  ra = FV(acceptor[0], acceptor[1], acceptor[2]);
  rd = FV(donor[0], donor[1], donor[2]);

  /* check if there is a surface parallel to U, if so set surf_flag */
  switch(N_VOF(donor[0], donor[1], donor[2])) {
    case east:
    case west:
      if(x != 0) surf_flag = 1;
      break;
    case north:
    case south:
      if(x != 1) surf_flag = 1;
      break;
    case top:
    case bottom:
      if(x != 2) surf_flag = 1;
      break;
    case empty:
    case none:
      break;
  }

  acceptor_d[0] = acceptor[0];
  acceptor_d[1] = acceptor[1];
  acceptor_d[2] = acceptor[2];
  
  if(surf_flag) {
    if(VOF_N(acceptor[0], acceptor[1], acceptor[2]) > solver->emf &&
       VOF_N(donor_m[0], donor_m[1], donor_m[2]) > solver->emf) {
      acceptor_d[0] = donor[0];
      acceptor_d[1] = donor[1];
      acceptor_d[2] = donor[2];
    }
  }
  
  VOFdm = max(VOF_N(donor_m[0], donor_m[1], donor_m[2]),
              VOF_N(donor[0], donor[1], donor[2]));
  if(A < solver->emf) VOFdm = 1.0;

  if(rb > solver->emf && ra > solver->emf && rd > solver->emf) {
    CF = max((VOFdm-VOF_N(acceptor_d[0], acceptor_d[1], acceptor_d[2]))*fabs(v) - 
             (VOFdm-VOF_N(donor[0],donor[1], donor[2]))*del, 0.0);
    dVOF  = VOF_N(acceptor_d[0], acceptor_d[1], acceptor_d[2])*fabs(v) + CF;                  
    
  /*  # check - you can't give more than you got
    # also this would mean the timestep is too long */
    dVOF = min(dVOF,VOF_N(donor[0],donor[1], donor[2])*del*rd/rb);

    if(v < 0) dVOF *= -1.0; /* adjust sign */
    
    /* test stability */
    if(fabs(v) / del * (rb/rd) * VOF_N(donor[0],donor[1], donor[2]) > 0.8 * ra) 
      solver->vof_flag = 1;
    
  }
  else dVOF = 0;

  return dVOF; 

}

int vof_mpi_convect(struct solver_data *solver) {
  long int i,j,k;
  double vchg = 0.0;
  double dVOF;

#define emf solver->emf

  solver->vof_flag = 0;
  
  if(solver->t > 0) {
  /* this code only executes after the first timestep */
    for(i=0; i<IRANGE-1; i++) {
      for(j=0; j<JMAX-1; j++) {
        for(k=0; k<KMAX-1; k++) { 

          if(FV(i,j,k) < emf) continue;

          if(FV(i+1,j,k) > emf) {
            dVOF = calc_dVOF(solver, i, j, k, 0);
            VOF(i,j,k) -= dVOF * RDX * AE(i,j,k) / FV(i,j,k);
            VOF(i+1,j,k) += dVOF * RDX * AE(i,j,k) / FV(i+1,j,k);
          }

          if(FV(i,j+1,k) > emf) {
            dVOF = calc_dVOF(solver, i, j, k, 1);
            VOF(i,j,k) -= dVOF * RDY * AN(i,j,k) / FV(i,j,k);
            VOF(i,j+1,k) += dVOF * RDY * AN(i,j,k) / FV(i,j+1,k);
          }

          if(FV(i,j,k+1) > emf) {
            dVOF = calc_dVOF(solver, i, j, k, 2);
            VOF(i,j,k) -= dVOF * RDZ * AT(i,j,k) / FV(i,j,k);
            VOF(i,j,k+1) += dVOF * RDZ * AT(i,j,k) / FV(i,j,k+1);
          }
        }
      }
    }
  } 

#define min_vof solver->min_vof
#define max_vof solver->max_vof

  /* # this code executes on any timestep
  # it calculates how much VOF is being lost or gained in the solution */
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        if(FV(i,j,k) == 0)
          continue;

        vchg = 0;

        if(VOF(i,j,k) <= min_vof || VOF(i,j,k) >= max_vof) {
          if(VOF(i,j,k) <= min_vof) {
            /* # in this case the cell is emtpy */
            vchg=VOF(i,j,k);
            VOF(i,j,k)=0.0;
            P(i,j,k)=0.0;
            
            U(i,j,k) = min(U(i,j,k), 0);
            U(i-1,j,k) = max(U(i-1,j,k), 0);
            V(i,j,k) = min(V(i,j,k), 0);
            V(i,j-1,k) = max(V(i,j-1,k), 0);
            W(i,j,k) = min(W(i,j,k), 0);
            W(i,j,k-1) = max(W(i,j,k-1), 0);
            
            /* this code is intended to stabilize velocities around empty cells */
            if(VOF(i+1,j,k) <= min_vof) U(i,j,k) = 0.0;
            if(VOF(i,j+1,k) <= min_vof) V(i,j,k) = 0.0;
            if(VOF(i,j,k+1) <= min_vof) W(i,j,k) = 0.0;
          }
          else if(VOF(i,j,k) >= max_vof) {
            /*# in this case the cell is full */
            vchg = -(1.0-VOF(i,j,k));
            VOF(i,j,k)=1.0;
          }
        }

        
        solver->vchgt = solver->vchgt + vchg*DELX*DELY*DELZ*FV(i,j,k);
        
        /*# if there is a full cell with an empty cell adjacent
        # that full cell just loses 1.1 * emf of fluid*/
        if(VOF(i,j,k) >= max_vof) {
          if(VOF(i+1,j,k) < min_vof || VOF(i-1,j,k) < min_vof ||
             VOF(i,j+1,k) < min_vof || VOF(i,j-1,k) < min_vof ||
             VOF(i,j,k+1) < min_vof || VOF(i,j,k-1) < min_vof) {

             VOF(i,j,k) = VOF(i,j,k) - 1.1*min_vof;
             vchg=1.1*min_vof;

             solver->vchgt = solver->vchgt +vchg*DELX*DELY*DELZ*FV(i,j,k);
          }
        }
      }
    }
  }

  return 0;
#undef emf
#undef min_vof
#undef max_vof
}
