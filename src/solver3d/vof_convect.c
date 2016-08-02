
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
#include "vof.h"
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

  if(surf_flag) {
    if(VOF_N(acceptor[0], acceptor[1], acceptor[2]) > emf &&
       VOF_N(donor_m[0], donor_m[1], donor_m[2]) > emf) {
      acceptor_d[0] = donor[0];
      acceptor_d[1] = donor[1];
      acceptor_d[2] = donor[2];
    }
  }
  else {
    acceptor_d[0] = acceptor[0];
    acceptor_d[1] = acceptor[1];
    acceptor_d[2] = acceptor[2];
  }
  
  VOFdm = max(VOF_N(donor_m[0], donor_m[1], donor_m[2]),
              VOF_N(donor[0], donor[1], donor[2]));
  if(a < solver->emf) VOFdm = 1.0;

  if(rb > emf && ra > emf && rd > emf) {
    CF = max((VOFdm-VOF_N(acceptor_d[0], acceptor_d[1], acceptor_d[2]))*fabs(v) - 
             (VOFdm-VOF_N(donor[0],donor[1], donor[2]))*del, 0.0);
    dVOF  = VOF_N(acceptor_d[0], acceptor_d[1], acceptor_d[2])*fabs(v) + CF;                  
    
  /*  # check - you can't give more than you got
    # also this would mean the timestep is too long */
    dVOF = min(dVOF,VOF_N(donor[0],donor[1], donor[2],j,k)*del*rd/rb);

    if(v < 0) dVOF *= -1.0; /* adjust sign */
    
    /* test stability */
    if(fabs(v) / del * (rb/rd) * VOF_N(donor[0],donor[1], donor[2]) > 0.8 * ra) 
      solver->vof_flag = 1;
    
  }
  else dVOF = 0;

  return dVOF; 

}

int vof_mpi_vfconv(struct solver_data *solver) {
  long int i,j,k, ia, iad, id, idm, ja, jad, jd, jdm, ka, kad, kd, kdm;
  double vchg = 0.0;
  double rb, ra, rd, aedm, andm, atdm, vx, vy, vz;
  double dVOF, VOFdm, CF;

#define emf solver->emf

  solver->vof_flag = 0;
  
  if(solver->t > 0) {
  /* this code only executes after the first timestep */
    for(i=0; i<IRANGE-1; i++) {
      for(j=0; j<JMAX-1; j++) {
        for(k=0; k<KMAX-1; k++) { 

          if(FV(i,j,k)==0) continue;
          
          vx=U(i,j,k)*solver->delt;
          vy=V(i,j,k)*solver->delt;
          vz=W(i,j,k)*solver->delt;
          
          if(fabs(vx) > 0.5 * DELX 
          || fabs(vy) > 0.5 * DELY 
          || fabs(vz) > 0.5 * DELZ) {
            solver->vof_flag = 1;
          }

          if(vx >=0) {
            ia=i+1;
            id=i;
            idm=max(i-1,0);

            aedm = AE(idm,j,k);
          }
          else if(vx < 0) {
            ia=i;
            id=i+1;
            idm=min(i+2,IRANGE-1);

            aedm = AE(idm-1,j,k);
          }
          rb = AE(i,j,k);
          ra = FV(ia,j,k);
          rd = FV(id,j,k);

          iad = ia;

          /*
          # the variables are:
          # id = cell we are convecting from (i_donor)
          # ia = cell we are convecting to (i_acceptor)
          # idm = correction to move id inside the mesh if it is on a bdry
          # so idm becomes the cell we are convecting to

          # if the cell we are convecting from has a surface parallel to U */
          if(N_VOF(id,j,k) == south   || N_VOF(id,j,k) == north ||
             N_VOF(id,j,k) == bottom  || N_VOF(id,j,k) == top ) iad = id;
            /* # then we set iad to equal the cell we are convecting from
            # this will effectively make convection = 0 */

          /* # if the acceptor cell is emtpy in the previous timestep
          # or if the donor cell is empty in the previous timestep */
          if(VOF_N(ia,j,k) < emf || VOF_N(idm,j,k) < emf) iad = ia;

          /* # VOFdm is equal to the maximum donor VOF between idm and id
          # which is just id if we are interior to the mesh */
          VOFdm = max(VOF_N(idm,j,k), VOF_N(id,j,k));
          if(aedm < emf) VOFdm = 1.0;

        /*  # dVOF = VOF in the acceptor cell * velocity +
          #       [ (donor VOF - acceptor VOF) * velocity -
          #         (donor VOF - donor VOF) * velocity      ]
          # where the last term in [ ] must be at least 0
          # it is basically equal to the difference in VOF
          # between the donor and acceptor, multiplied by the velocity */

          if(rb > emf && ra > emf && rd > emf) {
            CF = max((VOFdm-VOF_N(iad,j,k))*fabs(vx)-(VOFdm-VOF_N(id,j,k))*DELX, 0.0);
            dVOF  = VOF_N(iad,j,k)*fabs(vx) + CF;                  
            
          /*  # check - you can't give more than you got
            # also this would mean the timestep is too long */
            dVOF = min(dVOF,VOF_N(id,j,k)*DELX*rd/rb);

          /*  # donor acceptor calc */
            VOF(id,j,k) = VOF(id,j,k) - dVOF*RDX*(rb/rd);
            VOF(ia,j,k) = VOF(ia,j,k) + dVOF*RDX*(rb/ra);
            
            /* test stability */
            if(fabs(vx) * RDX * (rb/rd) * VOF_N(id,j,k) > 0.8 * ra) {
              printf("Overconvecting f:   %e > %e.  Adjusting timestep\n", fabs(vx) * RDX * (rb/rd) , 0.8 * ra);
              solver->vof_flag = 1;
            }
            
          }

          /* repeat the calculation for the y axis */

          if(N_VOF(i,j,kd) == west   || N_VOF(i,j,kd) == east ||
             N_VOF(i,j,kd) == north  || N_VOF(i,j,kd) == south) kad=kd;
          if(vy >= 0) {
            ja = j+1;
            jd = j;
            jdm = max(j-1,0);
            andm = AN(i,jdm,k);
          }
          else {
            ja = j;
            jd = j+1;
            jdm = min(j+2,JMAX-1);
            andm = AN(i,jdm-1,k);
          }

          rb = AN(i,j,k);
          ra = FV(i,ja,k);
          rd = FV(i,jd,k);

          jad = ja;

          if(N_VOF(i,jd,k) == west || N_VOF(i,jd,k) == east ||
             N_VOF(i,jd,k) == top  || N_VOF(i,jd,k) == bottom) jad=jd;

          if(VOF_N(i,ja,k) < emf || VOF_N(i,jdm,k) < emf) jad=ja;

          VOFdm = max(VOF_N(i,jdm,k),VOF_N(i,jd,k));

          if(andm < emf) VOFdm = 1.0;

          if (rb > emf && ra > emf && rd > emf) {
            CF  = max((VOFdm-VOF_N(i,jad,k))*fabs(vy)-(VOFdm-VOF_N(i,jd,k))*DELY,0.0);
            dVOF  = VOF_N(i,jad,k)*fabs(vy) + CF;
                  
            dVOF  = min(dVOF,VOF_N(i,jd,k)*DELY*rd/rb);

            VOF(i,jd,k) = VOF(i,jd,k) - dVOF*RDY*(rb/rd);
            VOF(i,ja,k) = VOF(i,ja,k) + dVOF*RDY*(rb/ra);
            
            if(fabs(vy) * RDY * (rb/rd) * VOF_N(i, jd, k) > 0.8 * ra) {
              printf("Overconvecting f:   %e > %e.  Adjusting timestep\n", fabs(vy) * RDZ * (rb/rd) , 0.8 * ra);
              solver->vof_flag = 1;
            }

          } 
          /* repeat the calculation for the z axis */

          if(vz >= 0) {
            ka = k+1;
            kd = k;
            kdm = max(k-1,0);
            atdm = AT(i,j,kdm);
          }
          else {
            ka = k;
            kd = k+1;
            kdm = min(k+2,KMAX-1);
            atdm = AT(i,j,kdm-1);
          }

          rb = AT(i,j,k);
          ra = FV(i,j,ka);
          rd = FV(i,j,kd);

          kad = ka;

          if(N_VOF(i,j,kd) == west   || N_VOF(i,j,kd) == east ||
             N_VOF(i,j,kd) == north  || N_VOF(i,j,kd) == south) kad=kd;

          if(VOF_N(i,j,ka) < emf || VOF_N(i,j,kdm) < emf) kad=ka;

          VOFdm = max(VOF_N(i,j,kdm),VOF_N(i,j,kd));

          if(atdm < emf) VOFdm = 1.0;

          if (rb > emf && ra > emf && rd > emf) {
            CF = max((VOFdm-VOF_N(i,j,kad))*fabs(vz)-(VOFdm-VOF_N(i,j,kd))*DELZ,0.0);
            dVOF = VOF_N(i,j,kad)*fabs(vz) + CF;
                  
            dVOF  = min(dVOF,VOF_N(i,j,kd)*DELZ*rd/rb);

            VOF(i,j,kd) = VOF(i,j,kd) - dVOF*RDZ*(rb/rd);
            VOF(i,j,ka) = VOF(i,j,ka) + dVOF*RDZ*(rb/ra);
                        
            if(fabs(vz) * RDZ * (rb/rd) * VOF_N(i,j,kd) > 0.8 * ra) {
              printf("Overconvecting f:   %e > %e.  Adjusting timestep\n", fabs(vz) * RDZ * (rb/rd) , 0.8 * ra);
              solver->vof_flag = 1;
            }

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
          /* # this code executes if this cell is not a free surface
          # as a note, this is just the kind of code that could be inlined
          # the potential cases of not being a free surface: */
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