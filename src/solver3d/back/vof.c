/* 3dvof.c
 *
 * default solver based on VOF and fractional area/volume obstacles
 */

#include <stdio.h>
#include <math.h>

#include "vof.h"
#include "solver.h"
#include "mesh.h"

/* macros to reduce code needed for access to array elements */
#define U(i, j, k) solver->mesh->u[mesh_index(solver->mesh, i, j, k)]
#define V(i, j, k) solver->mesh->v[mesh_index(solver->mesh, i, j, k)]
#define W(i, j, k) solver->mesh->w[mesh_index(solver->mesh, i, j, k)]
#define UN(i, j, k) mesh_n->u[mesh_index(mesh_n, i, j, k)]
#define VN(i, j, k) mesh_n->v[mesh_index(mesh_n, i, j, k)]
#define WN(i, j, k) mesh_n->w[mesh_index(mesh_n, i, j, k)]
#define P(i, j, k) solver->mesh->P[mesh_index(solver->mesh, i, j, k)]

#define VOF(i, j, k) solver->mesh->vof[mesh_index(solver->mesh, i, j, k)]
#define VOF_N(i, j, k) mesh_n->vof[mesh_index(mesh_n, i, j, k)]
#define N_VOF(i, j, k) solver->mesh->n_vof[mesh_index(solver->mesh, i, j, k)]
#define FV(i, j, k) solver->mesh->fv[mesh_index(solver->mesh, i, j, k)]
#define AE(i, j, k) solver->mesh->ae[mesh_index(solver->mesh, i, j, k)]
#define AN(i, j, k) solver->mesh->an[mesh_index(solver->mesh, i, j, k)]
#define AT(i, j, k) solver->mesh->at[mesh_index(solver->mesh, i, j, k)]

#define PETA(i, j, k) solver->mesh->peta[mesh_index(solver->mesh, i, j, k)]
#define BETA(i, j, k) solver->mesh->beta[mesh_index(solver->mesh, i, j, k)]
#define TANTH(i, j, k) solver->mesh->tanth[mesh_index(solver->mesh, i, j, k)]

#define DELX solver->mesh->delx
#define DELY solver->mesh->dely
#define DELZ solver->mesh->delz
#define RDX solver->mesh->rdx
#define RDY solver->mesh->rdy
#define RDZ solver->mesh->rdz
#define IMAX solver->mesh->imax
#define JMAX solver->mesh->jmax
#define KMAX solver->mesh->kmax

int vof_init_solver(struct solver_data *solver) {

  solver->loop = vof_loop;
  solver->boundaries = vof_boundaries;
  solver->special_boundaries = vof_special_boundaries;
  solver->pressure = vof_pressure;
  solver->velocity = vof_velocity;
  solver->vfconv = vof_vfconv;
  solver->petacal = vof_petacal;
  solver->betacal = vof_betacal;

  return 0;
}

#define dim(i,j,k) i+3*(j+k*3)

int vof_velocity(struct solver_data *solver) {
  double vel[3][27];
  double af[3][27];
  double Flux[3];

  double Hu_E, u_E_CD, u_E_upwind, Q_E;
  int n;
  const int pdim[3][3] = { {  2,1,1 }, { 1,2,1 }, { 1,1,2 } };
  const int ndim[3][3] = { {  0,1,1 }, { 1,0,1 }, { 1,1,0 } };
  const int odim[3][3] = { {  1,0,0 }, { 0,1,0 }, { 0,0,1 } };
  int ro_p1, ro_m1; 
  const int ro=dim(1,1,1);

  const double del[3] = { DELX, DELY, DELZ };

  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
          
        if (FV(i,j,k) == 0.0) continue;

        for(n=0; n<3; n++) {
          for(m=0; m<3; m++) {
            for(o=0; o<3; o++) {
              vel[0][dim(n,m,o)] = UN(i-1+m,j-1+n,k-1+o);
              vel[1][dim(n,m,o)] = VN(i-1+m,j-1+n,k-1+o);
              vel[2][dim(n,m,o)] = WN(i-1+m,j-1+n,k-1+o);

              af[0][dim(n,m,o)] = AE(i-1+m,j-1+n,k-1+o);
              af[1][dim(n,m,o)] = AN(i-1+m,j-1+n,k-1+o);
              af[2][dim(n,m,o)] = AT(i-1+m,j-1+n,k-1+o);
            }
          }
        }

        for(n=0; n<3; n++) {

          if(VOF(i,j,k) + VOF(i+odim[n][0],j+odim[n][1],k+odim[n][2]) < solver->emf 
             && solver->nmat == 1) {
            U(i,j,k) = 0;
            continue;
          }

          ro_p1 = dim(pdim[n][0], pdim[n][1], pdim[n][2]);
          ro_m1 = dim(ndim[n][0], ndim[n][1], ndim[n][2]);

          /* cell centered fluxes */
          Q_C = 0;

          /* Flux from cell centered to the east/north/top */
          H_vel = (vel[n][ro] * af[n][ro] + 
                   vel[n][ro_p1] * af[n][ro_p1]) / 2;

          CD = (vel[n][ro] +
                vel[n][ro_p1]) / 2
          
          if (CD >= 0) 
            upwind = vel[n][ro]; 
          else
            upwind = vel[n][ro_p1];

          if(af[n][ro] > 0.01)
            Q_C += 2/del[n] * H_vel * ((1-solver->alpha)*(CD - vel[n][ro]) +
                                         solver->alpha*(upwind - vel[n][ro]));
          
          /* Flux from cell centered to the west/south/bottom */
          H_vel = (vel[n][ro] * af[n][ro] + 
                   vel[n][ro_m1] * af[n][ro_m1]) / 2;

          CD = (vel[n][ro] +
                vel[n][ro_m1]) / 2
          
          if (CD >= 0) 
            upwind = vel[n][ro_m1]; 
          else
            upwind = vel[n][ro];

          if(af[n][ro] > 0.01)
            Q_C += 2/del[n] * H_vel * ((1-solver->alpha)*(vel[n][ro] - CD) +
                                         solver->alpha*(vel[n][ro] - upwind));
 
          

          /* Viscocity Calculation */
          vis[n] = 0;
          if(af[n][ro_p1] > 0.01) 
            vis[n] += ((af[n][ro]+af[n][ro_p1]) / 2) * (vel[n][ro_p1]-vel[n][ro]);
          if(af[n][ro_m1] > 0.01) 
            vis[n] -= ((af[n][ro]+af[n][ro_m1]) / 2) * (vel[n][ro]-vel[n][ro_m1]);

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

            CD = (vel[n][ro] + vel[n][ro_mp1])/2

            if(H_vel >= 0)
              upwind = vel[n][ro];
            else
              upwind = vel[n][ro_mp1];

            if ((af[m][ro]>0.01 || af[m][ro_p1]>0.01) && 
                 af[n][ro_mp1]>0.01)
              Q_W += 2/del[m] * H_vel * ((1-alpha)*(CD - vel[n][ro]) + 
                                          alpha*(upwind - vel[n][ro]));
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
              ro_nmm1 = dim(pdim[n][0], 0, pdim[n][1]);
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

            CD = (vel[n][ro] + vel[n][ro_mp1])/2

            if(H_vel >= 0)
              upwind = vel[n][ro_mm1];
            else
              upwind = vel[n][ro];

            if ((af[m][ro]>0.01 || af[m][ro_p1]>0.01) &&
                 af[n][ro_mp1]>0.01)
              Q_W += 2/del[m] * H_vel * ((1-alpha)*(vel[n][ro] - CD) + 
                                          alpha*(vel[n][ro] - upwind));

            /* Viscocity calc */
            vis[m] = 0;
            if(af[n][ro_mp1] > 0.01)
              vis[m] += ((af[m][ro]+af[m][ro_p1])/2) * (vel[n][ro_mp1] - vel[n][ro]);
            if(af[n][ro_mm1] > 0.01)
              vis[m] -= ((af[m][ro_mm1]+af[m][ro_nmm1])/2) * (vel[n][ro] - vel[n][ro_mm1]);

          }

          sum_fv = (FV(i,j,k) + FV(i+odim[n][0],j+odim[n][1],k+odim[n][2]));
          delp   = (P(i,j,k)  -  P(i+odim[n][0],j+odim[n][1],k+odim[n][3]));

          Flux = (Q_E + Q_W) / sum_fv;
                            
          Viscocity = NU * (vis[0]/pow(del[0],2) + vis[1]/pow(del[1],2) + vis[2]/pow(del[2],2));

          delv = solver->delt * ( (sum_fv/2) * RDX * delp / solver->rho +
                 solver->gx * odim[n][0] + solver->gy * odim[n][1] + solver->gz * odim[n][2] -
                 Flux + Viscocity );

          switch(n) {
          case 0:
            U(i,j,k) = UN(i,j,k) + delv;
            break;
          case 1:
            V(i,j,k) = VN(i,j,k) + delv;
            break;
          case 2:
            W(i,j,k) = WN(i,j,k) + delv;
            break;
          }
        
        }
      
      }
    }
  }
}



int vof_velocity(struct solver_data *solver) {
  double Hu_E, u_E_CD, u_E_upwind, Q_E;
  
  for(i=0; i<IMAX; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
          
        if (FV(i,j,k) == 0.0) continue;

        /* Calculate U fluxes */

        /* Flux from cell centered to the east */
        Hu_E = (AE(i,j,k)*UN(i,j,k)+AE(i+1,j,k)*UN(i+1,j,k))/2;
        u_E_CD = (UN(i,j,k)+UN(i+1,j,k))/2;
        
        if (u_E_CD >= 0) 
          u_E_upwind = UN(i,j,k);
        else
          u_E_upwind = UN(i+1,j,k);

        if(AE(i,j,k) > 0.01)
          Q_E = 2/DELX * Hu_E * ((1-solver->alpha)*(u_E_CD - UN(i,j,k)) +
                                         alpha*(u_E_upwind - UN(i,j,k)));
        else Q_E=0;
        
        /* flux from cell centered to the west */
        Hu_C = (AE(i,j,k)*UN(i,j,k)+AE(i-1,j,k)*UN(i-1,j,k))/2;
        u_C_CD = (UN(i,j,k)+UN(i-1,j,k))/2;

        if(u_C_CD >=0)
          u_C_upwind = UN(i-1,j,k);
        else
          u_C_upwind = UN(i,j,k);


        if(AE(i-1,j,k) > 0.01)
          Q_C = 2/delx * Hu_C * ((1-alpha)*(UN(i,j,k) - u_C_CD) + 
                               alpha*(UN(i,j,k) - u_C_upwind));
        else
          Q_C=0;

        /* Flux from wall centered to the north */
        Hv_ne = (AN(i,j,k)*VN(i,j,k)+AN(i+1,j,k)*VN(i+1,j,k))/2;
        if (AN(i,j,k) < 0.01)
          Hv_ne = AN(i+1,j,k)*VN(i+1,j,k);
        else if (AN(i+1,j,k) < 0.01)
          Hv_ne = AN(i,j,k)*VN(i,j,k);

        u_ne_CD = (UN(i,j,k)+UN(i,j+1,k))/2;
        v_ne = (VN(i,j,k)+VN(i+1,j,k))/2;

        if(v_ne >= 0)
          u_ne_upwind = UN(i,j,k);
        else
          u_ne_upwind = UN(i,j+1,k);

        if ((AN(i,j,k)<0.01 and AN(i+1,j,k)<0.01) or AE(i,j+1,k)<0.01)
          Q_ne=0;
        else
          Q_ne = 2/dely * Hv_ne * ((1-alpha)*(u_ne_CD - UN(i,j,k)) + 
                                 alpha*(u_ne_upwind - UN(i,j,k)));
         
        /* Flux from wall centered to the south */
        Hv_se = (AN(i,j-1,k)*VN(i,j-1,k)+AN(i+1,j-1,k)*VN(i+1,j-1,k))/2;
        if (AN(i,j-1,k) < 0.01)
          Hv_se = AN(i+1,j-1,k)*VN(i+1,j-1,k);
        else if (AN(i+1,j-1,k) < 0.01)
          Hv_se = AN(i,j-1,k)*VN(i,j-1,k);

        u_se_CD = (UN(i,j,k)+UN(i,j-1,k))/2;
        v_se = (VN(i,j-1,k)+VN(i+1,j-1,k))/2;

        if(v_se >= 0)
          u_se_upwind = UN(i,j-1,k);
        else
          u_se_upwind = UN(i,j,k);

        if((AN(i,j-1,k)<0.01 && AN(i+1,j-1,k)<0.01) || AE(i,j-1,k)<0.01)
          Q_se=0;
        else
          Q_se = 2/dely * Hv_se * ((1-alpha)*(UN(i,j,k) - u_se_CD) + 
                                 alpha*(UN(i,j,k) - u_se_upwind));

        /* Flux from wall centered on top */
        Hw_te = (AT(i,j,k)*WN(i,j,k)+AT(i+1,j,k)*WN(i+1,j,k))/2;
        if (AT(i,j,k) < 0.01)
          Hw_te = AT(i+1,j,k)*WN(i+1,j,k);
        else if (AT(i+1,j,k) < 0.01)
          Hw_te = AT(i,j,k)*WN(i,j,k);

        u_te_CD = (UN(i,j,k)+UN(i,j,k+1))/2;
        w_te = (WN(i,j,k)+WN(i+1,j,k))/2;

        if(w_te >= 0)
          u_te_upwind = UN(i,j,k);
        else
          u_te_upwind = UN(i,j,k+1);

        if ((AT(i,j,k)<0.01 && AT(i+1,j,k)<0.01) || AE(i,j,k+1)<0.01)
          Q_te=0;
        else
          Q_te = 2/delz * Hw_te * ((1-alpha)*(u_te_CD - UN(i,j,k)) + 
                                 alpha*(u_te_upwind - UN(i,j,k)));
         
        /* Flux from wall centered on bottom */
        Hw_be = (AT(i,j,k-1)*WN(i,j,k-1)+AT(i+1,j,k-1)*WN(i+1,j,k-1))/2;
        if (AT(i,j,k-1) < 0.01)
          Hw_te = AT(i+1,j,k-1)*WN(i+1,j,k-1);
        else if (AT(i+1,j,k-1) < 0.01)
          Hw_te = AT(i,j,k-1)*WN(i,j,k-1);

        u_be_CD = (UN(i,j,k-1)+UN(i,j,k))/2;
        w_be = (WN(i,j,k-1)+WN(i+1,j,k-1))/2;

        if(w_be >= 0)
          u_be_upwind = UN(i,j,k-1);
        else
          u_be_upwind = UN(i,j,k);

        if ((AT(i,j,k-1)<0.01 && AT(i+1,j,k-1)<0.01) || AE(i,j,k-1)<0.01)
          Q_be=0;
        else
          Q_be = 2/delz * Hw_be * ((1-alpha)*(UN(i,j,k) - u_be_CD) + 
                                 alpha*(UN(i,j,k) - u_be_upwind));
          
        fUxyz = (Q_E + Q_C + Q_ne + Q_se + Q_te + Q_be) / 
                (2 * (FV(i,j,k) + FV(i+1,j,k))/2);

        /* Calculate V fluxes */

        /* Flux from cell centered to the north */
        Hv_N = (AN(i,j,k)*VN(i,j,k)+AN(i,j+1,k)*VN(i,j+1,k))/2;
        v_N_CD = (VN(i,j,k)+VN(i,j+1,k))/2;
        
        if (v_N_CD >= 0) 
          v_N_upwind = VN(i,j,k);
        else
          v_N_upwind = VN(i,j+1,k);

        if(AE(i,j,k) > 0.01)
          R_N = 2/DELX * Hv_N * ((1-solver->alpha)*(v_N_CD - VN(i,j,k)) +
                                         alpha*(v_N_upwind - VN(i,j,k)));
        else R_N=0;
        
        /* flux from cell centered to the south */
        Hv_C = (AN(i,j,k)*VN(i,j,k)+AN(i,j-1,k)*VN(i,j-1,k))/2;
        v_C_CD = (VN(i,j,k)+VN(i,j-1,k))/2;

        if(v_C_CD >=0)
          v_C_upwind = VN(i,j-1,k);
        else
          v_C_upwind = VN(i,j,k);


        if(AN(i,j-1,k) > 0.01)
          R_C = 2/DELX * Hv_C * ((1-alpha)*(VN(i,j,k) - v_C_CD) + 
                               alpha*(VN(i,j,k) - v_C_upwind));
        else
          R_C=0;

        /* Flux from wall centered to the east 
         * START HERE*/
        Hu_ne = (AE(i,j,k)*UN(i,j,k)+AE(i,j+1,k)*UN(i,j+1,k))/2;
        if (AE(i,j,k) < 0.01)
          Hu_ne = AE(i,j+1,k)*UN(i,j+1,k);
        else if (AN(i+1,j,k) < 0.01)
          Hv_ne = AN(i,j,k)*VN(i,j,k);

        u_ne_CD = (UN(i,j,k)+UN(i,j+1,k))/2;
        v_ne = (VN(i,j,k)+VN(i+1,j,k))/2;

        if(v_ne >= 0)
          u_ne_upwind = UN(i,j,k);
        else
          u_ne_upwind = UN(i,j+1,k);

        if ((AN(i,j,k)<0.01 and AN(i+1,j,k)<0.01) or AE(i,j+1,k)<0.01)
          Q_ne=0;
        else
          Q_ne = 2/dely * Hv_ne * ((1-alpha)*(u_ne_CD - UN(i,j,k)) + 
                                 alpha*(u_ne_upwind - UN(i,j,k)));
         
        /* Flux from wall centered to the south */
        Hv_se = (AN(i,j-1,k)*VN(i,j-1,k)+AN(i+1,j-1,k)*VN(i+1,j-1,k))/2;
        if (AN(i,j-1,k) < 0.01)
          Hv_se = AN(i+1,j-1,k)*VN(i+1,j-1,k);
        else if (AN(i+1,j-1,k) < 0.01)
          Hv_se = AN(i,j-1,k)*VN(i,j-1,k);

        u_se_CD = (UN(i,j,k)+UN(i,j-1,k))/2;
        v_se = (VN(i,j-1,k)+VN(i+1,j-1,k))/2;

        if(v_se >= 0)
          u_se_upwind = UN(i,j-1,k);
        else
          u_se_upwind = UN(i,j,k);

        if((AN(i,j-1,k)<0.01 && AN(i+1,j-1,k)<0.01) || AE(i,j-1,k)<0.01)
          Q_se=0;
        else
          Q_se = 2/dely * Hv_se * ((1-alpha)*(UN(i,j,k) - u_se_CD) + 
                                 alpha*(UN(i,j,k) - u_se_upwind));

        /* Flux from wall centered on top */
        Hw_te = (AT(i,j,k)*WN(i,j,k)+AT(i+1,j,k)*WN(i+1,j,k))/2;
        if (AT(i,j,k) < 0.01)
          Hw_te = AT(i+1,j,k)*WN(i+1,j,k);
        else if (AT(i+1,j,k) < 0.01)
          Hw_te = AT(i,j,k)*WN(i,j,k);

        u_te_CD = (UN(i,j,k)+UN(i,j,k+1))/2;
        w_te = (WN(i,j,k)+WN(i+1,j,k))/2;

        if(w_te >= 0)
          u_te_upwind = UN(i,j,k);
        else
          u_te_upwind = UN(i,j,k+1);

        if ((AT(i,j,k)<0.01 && AT(i+1,j,k)<0.01) || AE(i,j,k+1)<0.01)
          Q_te=0;
        else
          Q_te = 2/delz * Hw_te * ((1-alpha)*(u_te_CD - UN(i,j,k)) + 
                                 alpha*(u_te_upwind - UN(i,j,k)));
         
        /* Flux from wall centered on bottom */
        Hw_be = (AT(i,j,k-1)*WN(i,j,k-1)+AT(i+1,j,k-1)*WN(i+1,j,k-1))/2;
        if (AT(i,j,k-1) < 0.01)
          Hw_te = AT(i+1,j,k-1)*WN(i+1,j,k-1);
        else if (AT(i+1,j,k-1) < 0.01)
          Hw_te = AT(i,j,k-1)*WN(i,j,k-1);

        u_be_CD = (UN(i,j,k-1)+UN(i,j,k))/2;
        w_be = (WN(i,j,k-1)+WN(i+1,j,k-1))/2;

        if(w_be >= 0)
          u_be_upwind = UN(i,j,k-1);
        else
          u_be_upwind = UN(i,j,k);

        if ((AT(i,j,k-1)<0.01 && AT(i+1,j,k-1)<0.01) || AE(i,j,k-1)<0.01)
          Q_be=0;
        else
          Q_be = 2/delz * Hw_be * ((1-alpha)*(UN(i,j,k) - u_be_CD) + 
                                 alpha*(UN(i,j,k) - u_be_upwind));
          
        fUxyz = (Q_E + Q_C + Q_ne + Q_se + Q_te + Q_be) / 
                (2 * (FV(i,j,k) + FV(i+1,j,k))/2);



      }
    }
  }
}

int vof_loop(struct solver_data *solver) {

  struct mesh_data *mesh_n; /* describes the mesh at the previous 
                             * timestep.  used for any explicit
                             * calculations. */

  mesh_n = mesh_init_copy(solver->mesh);
  if(mesh_n == NULL)
    return 1;

  if(solver->petacal != NULL)
    solver->petacal(solver);
  
  solver->boundaries(solver);
  
  if(solver->betacal != NULL)
    solver->betacal(solver);

  vof_hydrostatic(solver);

  while(solver->t < solver->endt) {

    solver->velocity(solver);

    solver->boundaries(solver);
    
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);

    solver->iter = 0;

    while(solver->iter < solver->niter) {
    
      if(solver->pressure(solver) == 0) break;

      solver->boundaries(solver);
      if(solver->special_boundaries != NULL)
        solver->special_boundaries(solver);
    
      solver->iter++;
    }

    solver->boundaries(solver);
    solver->vfconv(solver);
    solver->boundaries(solver);

    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);
    if(solver->petacal != NULL)
      solver->petacal(solver);

    if(solver->iter > solver->niter) {
      printf("timestep: %lf pressure did not converge\n", solver->t);
    }
    else {
      printf("timestep: %lf convergence in %ld iterations\n", solver->t, solver->iter);
    }


    solver->t += solver->delt;
    mesh_copy_data(mesh_n, solver->mesh);
  }

  return 0;
}

int vof_boundaries(struct solver_data *solver);

int vof_special_boundaries(struct solver_data *solver);

int vof_pressure(struct solver_data *solver);

int vof_velocity(struct solver_data *solver);

int vof_vfconv(struct solver_data *solver);

int vof_petacal(struct solver_data *solver);

int vof_betacal(struct solver_data *solver);

int vof_hydrostatic(struct solver_data *solver);
