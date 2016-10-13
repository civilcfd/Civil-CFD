/* vof_velocity.c
 *
 * velocity predictor routines
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <omp.h>
#include <mpi.h>
 
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


int vof_mpi_velocity_upwind(struct solver_data *solver) {
  double vel[3][27];
  double af[3][27];
  double vis[3];
  double Flux, Viscocity, Q_C, Q_W, H_vel, upwind, sum_fv, delp, delv, nu;

  long int i,j,k;
  int n,m,o;

#define dim(i,j,k) i+3*(j+k*3)
  const int pdim[3][3] = { {  2,1,1 }, { 1,2,1 }, { 1,1,2 } };
  const int ndim[3][3] = { {  0,1,1 }, { 1,0,1 }, { 1,1,0 } };
  const int odim[3][3] = { {  1,0,0 }, { 0,1,0 }, { 0,0,1 } };
  int ro_p1, ro_m1, ro_mp1, ro_mm1, ro_nmm1; 
  
  const int ro=dim(1,1,1);

  const double del[3] = { DELX, DELY, DELZ };

  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        U(i,j,k) = 0;
        V(i,j,k) = 0;
        W(i,j,k) = 0;
          
        if (FV(i,j,k) == 0.0) continue;

        for(m=0; m<3; m++) { /* fixed 6/16 from n,m,o */
          for(n=0; n<3; n++) {
            for(o=0; o<3; o++) {
              /* vel[n][i*j*k] and af["]["] define a matrix
               * where n is { U, V, W }
               * and the second dimension represents the values of the scalar in the
               * current cell (i,j,k) and the 8 surrounding cells
               * the current cell is given the location 1,1,1 
               * this caches the data and allows the velocity predictor
               * calcs to be generalized */

              vel[0][dim(m,n,o)] = UN(i-1+m,j-1+n,k-1+o);
              vel[1][dim(m,n,o)] = VN(i-1+m,j-1+n,k-1+o);
              vel[2][dim(m,n,o)] = WN(i-1+m,j-1+n,k-1+o);

              af[0][dim(m,n,o)] = AE(i-1+m,j-1+n,k-1+o);
              af[1][dim(m,n,o)] = AN(i-1+m,j-1+n,k-1+o);
              af[2][dim(m,n,o)] = AT(i-1+m,j-1+n,k-1+o);
            }
          }
        }

        for(n=0; n<3; n++) {
        
          /* ADDED 9/12 to eliminate pointless calcs that mess things up */
          if(af[n][ro] < solver->emf) continue;

          if(VOF(i,j,k) + VOF(i+odim[n][0],j+odim[n][1],k+odim[n][2]) < solver->emf /*
             || (N_VOF(i,j,k) >  7 && N_VOF(i+odim[n][0],j+odim[n][1],k+odim[n][2]) > 0)
             || (N_VOF(i,j,k) >  0 && N_VOF(i+odim[n][0],j+odim[n][1],k+odim[n][2]) > 7) */) { /* added 09/13 */
            switch(n) {
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
            continue; 
          }

          /* ro: represents the position 1,1,1 and is the Relative Origin
           * ro_p1: represents the position of the origin plus 1 in the n dimension
           * ro_m1: represents the position of the origin minus 1 in the n dimension
           */
          ro_p1 = dim(pdim[n][0], pdim[n][1], pdim[n][2]);
          ro_m1 = dim(ndim[n][0], ndim[n][1], ndim[n][2]);

          /* cell centered fluxes 
           * this is the first term of the NS equation
           * for example, if n = 0, then this is: u * du/dx */
          Q_C = 0;

          /* Flux from cell centered to the east/north/top */
          H_vel = (vel[n][ro] * af[n][ro] + 
                   vel[n][ro_p1] * af[n][ro_p1]) / 2;
          
          if (vel[n][ro] >= 0) 
            upwind = vel[n][ro]; 
          else
            upwind = vel[n][ro_p1];

          if(af[n][ro_p1] > 0.01)
            Q_C += 2/del[n] * H_vel * ((upwind - vel[n][ro]));
          
          /* Flux from cell centered to the west/south/bottom */
          H_vel = (vel[n][ro] * af[n][ro] + 
                   vel[n][ro_m1] * af[n][ro_m1]) / 2;
          
          if (vel[n][ro] >= 0) 
            upwind = vel[n][ro_m1]; 
          else
            upwind = vel[n][ro];

          if(af[n][ro_m1] > 0.01)
            Q_C += 2/del[n] * H_vel * ((vel[n][ro] - upwind));
 
          

          /* Viscocity Calculation */
          vis[n] = 0;
          if(af[n][ro_p1] > 0.01) 
            /* the first term is the average of the area fractions
             * the second term is du/dx if n=0 */
            vis[n] += ((af[n][ro]+af[n][ro_p1]) / 2) * (vel[n][ro_p1]-vel[n][ro]);
          if(af[n][ro_m1] > 0.01) 
            vis[n] -= ((af[n][ro]+af[n][ro_m1]) / 2) * (vel[n][ro]-vel[n][ro_m1]);

          /* Wall centered fluxes
           * this is the next two terms in the NS equation
           * for example, if n=0, this is: v du/dy + w du/dz
           */

          Q_W = 0; /* flux from walls */

          for(m=0; m<3; m++) {
            if(m==n) continue;

            /* ro_mp1:
             * for the m dimension, we add 1 to the origin
             * this represents the side of the wall that is outside of the cell i,j,k 
             */

            ro_mp1 = dim(pdim[m][0],pdim[m][1],pdim[m][2]);
            ro_mm1 = dim(ndim[m][0],ndim[m][1],ndim[m][2]);

            /* Flux from wall centered to the east/north/top */
            H_vel = (vel[m][ro]*af[m][ro] +
                     vel[m][ro_p1]*af[m][ro_p1]) / 2;
            if (af[m][ro] < 0.01)
              H_vel = vel[m][ro_p1]*af[m][ro_p1];
            else if (af[m][ro_p1] < 0.01)
              H_vel = vel[m][ro] * af[m][ro];

            if(H_vel >= 0)
              upwind = vel[n][ro];
            else
              upwind = vel[n][ro_mp1];

            if ((af[m][ro]>0.01 || af[m][ro_p1]>0.01) && 
                 af[n][ro_mp1]>0.01)
              Q_W += 2/del[m] * H_vel * ((upwind - vel[n][ro]));
            /* ro_mm1:
             * for the m dimension, we subract 1 from the origin
             * this represents the side of the wall that is inside the cell i,j,k
             *
             * ro_nmm1:
             * For the n dimension, we add 1
             * for the m dimension, we subtract 1
             * this represents the side of the wall that is outside of the cell i,j,k
             */
            switch(m) {
            case 0:
              ro_nmm1 = dim(0, pdim[n][1], pdim[n][2]);
              break;
            case 1:
              ro_nmm1 = dim(pdim[n][0], 0, pdim[n][2]);
              break;
            case 2:
              ro_nmm1 = dim(pdim[n][0], pdim[n][1], 0);
              break;
            }

            /* Flux from wall centered to the west/south/bottom */
            H_vel = (vel[m][ro_mm1]*af[m][ro_mm1] +
                     vel[m][ro_nmm1]*af[m][ro_nmm1]) / 2;
            if (af[m][ro_mm1] < 0.01)
              H_vel = vel[m][ro_nmm1]*af[m][ro_nmm1];
            else if (af[m][ro_nmm1] < 0.01)
              H_vel = vel[m][ro_mm1] * af[m][ro_mm1];

            if(H_vel >= 0)
              upwind = vel[n][ro_mm1];
            else
              upwind = vel[n][ro];

            if ((af[m][ro_mm1]>0.01 || af[m][ro_nmm1]>0.01) &&
                 af[n][ro_mm1]>0.01)
              Q_W += 2/del[m] * H_vel * ((vel[n][ro] - upwind));

            /* Viscocity calc */
            vis[m] = 0;
            if(af[n][ro_mp1] > 0.01)
              vis[m] += ((af[m][ro]+af[m][ro_p1])/2) * (vel[n][ro_mp1] - vel[n][ro]);
            if(af[n][ro_mm1] > 0.01)
              vis[m] -= ((af[m][ro_mm1]+af[m][ro_nmm1])/2) * (vel[n][ro] - vel[n][ro_mm1]);

          }

          sum_fv = (FV(i,j,k) + FV(i+odim[n][0],j+odim[n][1],k+odim[n][2]));
          delp   = (P(i,j,k)  -  P(i+odim[n][0],j+odim[n][1],k+odim[n][2]));
          if(FV(i+odim[n][0],j+odim[n][1],k+odim[n][2]) < 0.000001) delp=0; /* ADDED 2/27/16 testing */

          Flux = (Q_C + Q_W) / sum_fv;
          
          if(solver->turbulence_nu != NULL) {
            nu = (solver->turbulence_nu(solver,i,j,k) + 
                  solver->turbulence_nu(solver,i+odim[n][0],j+odim[n][1],k+odim[n][2]))/2;
          }
          else
            nu = solver->nu;
          solver->nu_max = max(nu, solver->nu_max);
                 
          Viscocity = nu * (vis[0]/pow(del[0],2) + vis[1]/pow(del[1],2) + vis[2]/pow(del[2],2));

          /* deleted from this code 6/18
           * sum_fv/2 * delp: this created discontinuity at pressure boundaries */

          delv = solver->delt * ( /*(sum_fv/2) * */ (1/del[n]) * delp / solver->rho +
                 solver->gx * odim[n][0] + solver->gy * odim[n][1] + solver->gz * odim[n][2] -
                 Flux + Viscocity ); 

          /* if(fabs(delv) < solver->epsi * solver->dzro * solver->delt / (solver->rho * del[n]))
            delv = 0; uncomment to eliminate spurious velocity currents */

          switch(n) {
          case 0:
            if(i != IMAX-2)  U(i,j,k) = UN(i,j,k) + delv;
            break;
          case 1:
            if(j != JMAX-2)  V(i,j,k) = VN(i,j,k) + delv;
            break;
          case 2:
            if(k != KMAX-2)  W(i,j,k) = WN(i,j,k) + delv;
            break;
          }
          
        
        }
      
      }
    }
  }

  return 0;
#undef dim
 }
 