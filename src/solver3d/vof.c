/* vof.c
 *
 * default solver based on VOF and fractional area/volume obstacles
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
#include "track.h"
#include "vof_baffles.h"

#include "vof_macros.h"

struct mesh_data *mesh_n; /* describes mesh at previous timestep for explicit calcs */

#ifdef DEBUG
float u(struct solver_data *solver, long int i,long int j,long int k) { return U(i,j,k); }
float v(struct solver_data *solver, long int i,long int j,long int k) { return V(i,j,k); }
float w(struct solver_data *solver, long int i,long int j,long int k) { return W(i,j,k); }
float un(struct solver_data *solver, long int i,long int j,long int k) { return UN(i,j,k); }
float vn(struct solver_data *solver, long int i,long int j,long int k) { return VN(i,j,k); }
float wn(struct solver_data *solver, long int i,long int j,long int k) { return WN(i,j,k); }
float p(struct solver_data *solver, long int i,long int j,long int k) { return P(i,j,k); }
float vof(struct solver_data *solver, long int i,long int j,long int k) { return VOF(i,j,k); }
float n_vof(struct solver_data *solver, long int i,long int j,long int k) { return N_VOF(i,j,k); }
float ae(struct solver_data *solver, long int i,long int j,long int k) { return AE(i,j,k); }
float an(struct solver_data *solver, long int i,long int j,long int k) { return AN(i,j,k); }
float at(struct solver_data *solver, long int i,long int j,long int k) { return AT(i,j,k); }
float fv(struct solver_data *solver, long int i,long int j,long int k) { return FV(i,j,k); }
float peta(struct solver_data *solver, long int i,long int j,long int k) { return PETA(i,j,k); }
float debug_d(struct solver_data *solver, long int i,long int j,long int k) { return D(i,j,k); }

void track_cell(struct solver_data *solver, long int i,long int j,long int k) {
  printf("U(e/n/t) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i,j,k),V(i,j,k),W(i,j,k));
  printf("U(w/s/b) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i-1,j,k),V(i,j-1,k),W(i,j,k-1));  
  printf("P: %lf    VOF: %lf    N_VOF: %d\n",P(i,j,k),VOF(i,j,k),N_VOF(i,j,k));
}
#endif

int vof_setup_solver(struct solver_data *solver) {
  
  solver->init = vof_init_solver;
  solver->kill = vof_kill_solver;
  solver->loop = vof_loop;
  solver->boundaries = vof_boundaries;
  solver->special_boundaries = vof_special_boundaries;
  solver->pressure = vof_pressure_gmres;
  solver->velocity = vof_velocity;
  solver->vfconv = vof_vfconv;
  solver->petacal = vof_petacal;
  solver->betacal = vof_betacal;
  solver->deltcal = vof_deltcal;
  solver->write = vof_write;
  solver->output = vof_output;
  
  return 0;
}

int vof_init_solver(struct solver_data *solver) {
 
  mesh_n = mesh_init_copy(solver->mesh);
  if(mesh_n == NULL)
    return 1;
    
  mesh_set_array(solver->mesh, "D", 1.0, -1, 0, 0, 0, 0, 0);
  
  if(solver->csq > 0 )
    solver->rdtexp = 2.0 * sqrt(fabs(solver->csq))/min(min(DELX,DELY),DELZ);
  else
    solver->rdtexp = 10000000000;
  
  if(vof_pressure_init(solver) == 1) return 1;
  
  
  return 0;
}

int vof_kill_solver(struct solver_data *solver) {
  mesh_free(solver->mesh);

  return 0;
}

int vof_betacal(struct solver_data *solver) {
  long int i,j,k;
  double abe, abw, abn, abs, abt, abb, xx;

#define emf solver->emf

   for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        if(FV(i,j,k) < emf) continue;

        abe = AE(i,j,k);
        abw = AE(i-1,j,k);
        abn = AN(i,j,k);
        abs = AN(i,j-1,k);
        abt = AT(i,j,k);
        abb = AT(i,j,k-1);

        xx = 2.0 * solver->delt * ( RDX * ( abe * RDX * 0.5 + abw * RDX * 0.5) +
                                    RDY * ( abn * RDY * 0.5 + abs * RDY * 0.5) +
                                    RDZ * ( abt * RDZ * 0.5 + abb * RDZ * 0.5)) / solver->rho;

        /* xx = xx / FV(i,j,k); */
        BETA(i,j,k) = solver->omg / xx;

      }
    }
  } 

#undef emf

  return 0;
}

int vof_petacal(struct solver_data *solver) {
  long int i,j,k,l,m,n,ii,jj,kk;
  int mobs, inf, infcr, iobs;
  double vf, fxm, fxp, fym, fyp, fzm, fzp; 
  double vfxm, vfxp, vfym, vfyp, vfzm, vfzp;
  double dmx, dmin, amn, bpd, dd, sdis;
  enum cell_boundaries N_VOF_F;

#define emf solver->emf
#define emf_c solver->emf_c

  const double nemf = -1.0 * emf;
 
  for(i=0; i<IMAX; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
        PETA(i,j,k) = 1.0;
        N_VOF(i,j,k) = none;
        if(i==0 || j==0 || k==0 || i==IMAX-1 || j==JMAX-1 || k==KMAX-1 || FV(i,j,k) == 0)
          N_VOF(i,j,k) = 0;
        else if(VOF(i+1,j,k) >= emf && VOF(i,j+1,k) >= emf 
          && VOF(i-1,j,k) >= emf && VOF(i,j-1,k) >= emf
          && VOF(i,j,k+1) >= emf && VOF(i,j,k-1) >= emf) {
          N_VOF(i,j,k) = 0;          
        }
      }
    }
  }
 
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
/*
       if(i==241 && j==48 && k==42) {
          printf("break\n");
        }  */
        /*# first consider all the loop exit conditions*/
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

        /* MOVED TO ABOVE  if(VOF(i+1,j,k) >= emf && VOF(i,j+1,k) >= emf 
          && VOF(i-1,j,k) >= emf && VOF(i,j-1,k) >= emf
          && VOF(i,j,k+1) >= emf && VOF(i,j,k-1) >= emf) {
          N_VOF(i,j,k) = 0;
          continue; * # this is the case where the fluid is not part of a free surf *
        } */
        if(N_VOF(i,j,k) == 0) continue;

        /*# now calculate the partial derivatives of F
        # this code calculates how much F changes in each axis
        # the goal is to determine where the free surface lies
        # and the slope of the free surface */
               
               fxm=VOF(i,j,k);
               fxp=VOF(i,j,k);
               fym=VOF(i,j,k);
               fyp=VOF(i,j,k);
               fzm=VOF(i,j,k);
               fzp=VOF(i,j,k);

               if (AE(i-1,j,k)>emf) fxm=VOF(i-1,j,k);
               if (AE(i,j,k)>emf)   fxp=VOF(i+1,j,k);

               if (AN(i,j-1,k)>emf) fym=VOF(i,j-1,k);
               if (AN(i,j,k)>emf)   fyp=VOF(i,j+1,k);

               if (AT(i,j,k-1)>emf) fzm=VOF(i,j,k-1);
               if (AT(i,j,k)>emf)   fzp=VOF(i,j,k+1);

               mobs=1;
               inf=1;
               iobs=1;
               vf=0.0;

               vfxm=0.0;
               vfxp=0.0;

               for (kk=1;kk<4;kk++)
                {
                  n=k-2+kk;
                  for (jj=1;jj<4;jj++)
                    {
                       m=j-2+jj;
                       vfxm=vfxm+VOF(i-1,m,n);
                       vfxp=vfxp+VOF(i+1,m,n);
                    }
                }

               if (FV(i-1,j,k)==0.0 && solver->gx > nemf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fxm>=emf))
                 {
                    inf=inf+1;
                    if (vfxm>vf) N_VOF(i,j,k)=west;
                    if (N_VOF(i,j,k)==west) vf=vfxm;
                 }

               iobs=1;
               if (FV(i+1,j,k)==0.0 && solver->gx < emf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fxp>=emf))
                 {
                    inf=inf+1;
                    if (vfxp>vf) N_VOF(i,j,k)=east;
                    if (N_VOF(i,j,k)==east) vf=vfxp;
                 }
               iobs=1;
                                /* z-axis **/
               vfzm=0.0;
               vfzp=0.0;
               for (ii=1;ii<4;ii++)
                 {
                   l=i-2+ii;
                   for (jj=1;jj<4;jj++)
                     {
                       m=j-2+jj;
                       vfzm=vfzm+VOF(l,m,k-1);
                       vfzp=vfzp+VOF(l,m,k+1);
                     }
                 }
/*
       if(i==241 && j==48 && k==42) {
          printf("break\n");
        }  */
        
               if (FV(i,j,k-1)==0.0 && solver->gz > nemf) iobs=2; 
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fzm>=emf))
                 {
                   inf=inf+1;
                   if (vfzm>vf) N_VOF(i,j,k)=bottom;
                   if (N_VOF(i,j,k)==bottom) vf=vfzm;
                 }
               iobs=1;

               if (FV(i,j,k+1)==0.0 && solver->gz < emf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fzp>=emf))
                 {
                   inf=inf+1;
                   if (vfzp>vf) N_VOF(i,j,k)=top;
                   if (N_VOF(i,j,k)==top) vf=vfzp;
                 }
               iobs=1;
                                /*y-axis**/

               vfym=0.0;
               vfyp=0.0;

               for (kk=1;kk<4;kk++)
                 {
                   n=k-2+kk;
                   for (ii=1;ii<4;ii++)
                     {
                       l=i-2+ii;
                       vfym=vfym+VOF(l,j-1,n); 
                       vfyp=vfyp+VOF(l,j+1,n);
                     }
                 }

               if (FV(i,j-1,k)==0.0  && solver->gy > nemf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fym>=emf))
                 {
                   inf=inf+1;
                   if (vfym>vf) N_VOF(i,j,k)=south;
                   if (N_VOF(i,j,k)==south) vf=vfym;
                 }
               iobs=1;

               if (FV(i,j+1,k)==0.0  && solver->gy < emf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fyp>=emf))
                 {
                   inf=inf+1;
                   if (vfyp>vf) N_VOF(i,j,k)=north;
                   if (N_VOF(i,j,k)==north) vf=vfyp;
                 }
               iobs=1;

              
               /*if(N_VOF(i,j,k) != bottom) {
                printf("breakpoint\n");
               }*/

                                    /***/

               /* check if it is not a free surface, but is bounded by an obstacle */
               /* essentially we have fluid in each direction that isn't an obstacle */
               infcr=8-mobs;
               if ((inf==infcr)&&(infcr>1)) N_VOF(i,j,k)=0; 
               

      }
    }
  }

    for (i=1; i<IMAX-1; i++) {
     for (j=1; j<JMAX-1; j++) {
      for (k=1; k<KMAX-1; k++) {
          if ((N_VOF(i,j,k)>0.0)&&(N_VOF(i,j,k)<8.0)&&(FV(i,j,k)!=0)) /* changed from 
                                                                       N_VOF < 7 on 9/12 */
          {
           l=i; m=j; n=k;
           switch (N_VOF(i,j,k))
           {
            case west:
              dd=DELX; l=i-1;
            break;

            case east:
             dd=DELX;  l=i+1;
            break;

            case south:
             dd=DELY;  m=j-1;
            break;

            case north:
             dd=DELY;  m=j+1;
            break;

            case bottom:
             dd=DELZ;  n=k-1;
            break;

            case top:
             dd=DELZ;  n=k+1;
            break;

            case none:
              P(i,j,k)=0.1666667*(P(i+1,j,k)+P(i-1,j,k)+
                                    P(i,j+1,k)+P(i,j-1,k)+
                                    P(i,j,k+1)+P(i,j,k-1));
            break; 

          }

           if (N_VOF(i,j,k)!=none)
             {
               /* CHANGED 9/12 to allow N_VOF across boundary */
               sdis=VOF(i,j,k)*dd+VOF(l,m,n)*dd*0.5;
               if(FV(l,m,n) < emf) sdis = dd*0.5 + VOF(i,j,k) * dd;
               sdis=max(sdis,0.5*dd);
               PETA(i,j,k)=dd/sdis;

               if ((FV(l,m,n) >= emf)&&(N_VOF(l,m,n)!=0.0)) PETA(i,j,k)=1.0;
               if (PETA(i,j,k)>2.0) PETA(i,j,k)=2.0;
               if (PETA(i,j,k)<0.0) PETA(i,j,k)=0.0;
             }
          }
        }
      }
    }

  for(i=0; i<IMAX-1; i++) {
    for(j=0; j<JMAX-1; j++) {
      for(k=0; k<KMAX-1; k++) {

        N_VOF_F=N_VOF(i,j,k);

        if(N_VOF_F==0 || FV(i,j,k) < emf)
          continue;

        if(N_VOF_F>7) {
          P(i,j,k) = 0.0;
          continue;
        }

        l=i;
        m=j;
        n=k;

        switch(N_VOF_F) {
        case west:
          l=i-1;
          dmx=DELX;
          dmin=0.5*(dmx+DELX);
          amn=AE(l,j,k);
          break;

        case east:
          l=i+1;
          dmx=DELX;
          dmin=0.5*(dmx+DELX);
          amn=AE(i,j,k);
          break;

        case south:
          m=j-1;
          dmx=DELY;
          dmin=0.5*(dmx+DELY);
          amn=AN(i,m,k);
          break;

        case north:
          m=j+1;
          dmx=DELY;
          dmin=0.5*(dmx+DELY);
          amn=AN(i,j,k);
          break;

         case bottom:
          n=k-1;
          dmx=DELZ; 
          dmin=0.5*(dmx+DELZ);
          amn=AT(i,j,n);
          break;

        case top:
          n=k+1;
          dmx=DELZ;
          dmin=0.5*(dmx+DELZ);
          amn=AT(i,j,k);
          break;
          
        case none:
          continue;

       }

        if(N_VOF(l,m,n) > 0)
          continue;

        if(amn < emf)
          continue;

        /* calculate relaxation factor and store it in PETA */
        bpd = 1.0 / PETA(l,m,n) - BETA(l,m,n) * (1.0-PETA(i,j,k)) *
               amn * /* / FV(l,m,n) * */ solver->delt/(dmin*dmx) / solver->rho;
        
        PETA(l,m,n) = min(1.8/solver->omg, 1.0/bpd);

      }
    }
  }

  return 0;

#undef emf
#undef emf_c
}

int vof_vfconv(struct solver_data *solver) {
  long int i,j,k, ia, iad, id, idm, ja, jad, jd, jdm, ka, kad, kd, kdm;
  double vchg = 0.0;
  double rb, ra, rd, aedm, andm, atdm, vx, vy, vz;
  double fx, fx1, fy, fy1, fz, fz1, fdm;

#define emf solver->emf

  solver->vof_flag = 0;
  
  if(solver->t > 0 /*emf + solver->delt*/) {
  /* this code only executes after the first timestep */
    for(i=0; i<IMAX-1; i++) {
      for(j=0; j<JMAX-1; j++) {
        for(k=0; k<KMAX-1; k++) { /* corrected from MAX to MAX-1 6/28/14 */

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
            idm=min(i+2,IMAX-1);

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

          /* # fdm is equal to the maximum donor VOF between idm and id
          # which is just id if we are interior to the mesh */
          fdm = max(VOF_N(idm,j,k), VOF_N(id,j,k));
          if(aedm < emf) fdm = 1.0;

        /*  # fx1 = VOF in the acceptor cell * velocity +
          #       [ (donor VOF - acceptor VOF) * velocity -
          #         (donor VOF - donor VOF) * velocity      ]
          # where the last term in [ ] must be at least 0
          # it is basically equal to the difference in VOF
          # between the donor and acceptor, multiplied by the velocity */

          if(rb > emf && ra > emf && rd > emf) {
            fx1 = VOF_N(iad,j,k)*fabs(vx) + 
                  max((fdm-VOF_N(iad,j,k))*fabs(vx)-(fdm-VOF_N(id,j,k))*DELX, 0.0);
            
          /*  # check - you can't give more than you got
            # also this would mean the timestep is too long */
            fx = min(fx1,VOF_N(id,j,k)*DELX*rd/rb);

          /*  # donor acceptor calc */
            VOF(id,j,k) = VOF(id,j,k) - fx*RDX*(rb/rd);
            VOF(ia,j,k) = VOF(ia,j,k) + fx*RDX*(rb/ra);
            
            /* test stability */
            if(fabs(vx) * RDX * (rb/rd) * VOF_N(id,j,k) > 0.8 * ra) {
              printf("Overconvecting f:   %e > %e.  Adjusting timestep\n", fabs(vx) * RDX * (rb/rd) , 0.8 * ra);
              solver->vof_flag = 1;
            }
            
          }

          /* repeat the calculation for the y axis */

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

          fdm = max(VOF_N(i,jdm,k),VOF_N(i,jd,k));

          if(andm < emf) fdm = 1.0;

          if (rb > emf && ra > emf && rd > emf) {
            fy1 = VOF_N(i,jad,k)*fabs(vy) + 
                  max((fdm-VOF_N(i,jad,k))*fabs(vy)-(fdm-VOF_N(i,jd,k))*DELY,0.0);
            fy  = min(fy1,VOF_N(i,jd,k)*DELY*rd/rb);

            VOF(i,jd,k) = VOF(i,jd,k) - fy*RDY*(rb/rd);
            VOF(i,ja,k) = VOF(i,ja,k) + fy*RDY*(rb/ra);
            
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

          fdm = max(VOF_N(i,j,kdm),VOF_N(i,j,kd));

          if(atdm < emf) fdm = 1.0;

          if (rb > emf && ra > emf && rd > emf) {
            fz1 = VOF_N(i,j,kad)*fabs(vz) + 
                  max((fdm-VOF_N(i,j,kad))*fabs(vz)-(fdm-VOF_N(i,j,kd))*DELZ,0.0);
            fz  = min(fz1,VOF_N(i,j,kd)*DELZ*rd/rb);

            VOF(i,j,kd) = VOF(i,j,kd) - fz*RDZ*(rb/rd);
            VOF(i,j,ka) = VOF(i,j,ka) + fz*RDZ*(rb/ra);
                        
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
  for(i=1; i<IMAX-1; i++) {
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

int vof_velocity(struct solver_data *solver) {
  double vel[3][27];
  double af[3][27];
  double vis[3];
  double Flux, Viscocity, Q_C, Q_W, H_vel, CD, upwind, sum_fv, delp, delv, resi, nu;

  long int i,j,k;
  int n,m,o;

#define dim(i,j,k) i+3*(j+k*3)
  const int pdim[3][3] = { {  2,1,1 }, { 1,2,1 }, { 1,1,2 } };
  const int ndim[3][3] = { {  0,1,1 }, { 1,0,1 }, { 1,1,0 } };
  const int odim[3][3] = { {  1,0,0 }, { 0,1,0 }, { 0,0,1 } };
  int ro_p1, ro_m1, ro_mp1, ro_mm1, ro_nmm1; 
  
  const int ro=dim(1,1,1);

  const double del[3] = { DELX, DELY, DELZ };

#pragma omp parallel for shared (solver) \
            private(i,j,k,vel,af,vis,Flux,Viscocity,Q_C,Q_W,H_vel,CD,upwind,sum_fv,delp, \
                    delv,resi,nu,n,m,o,ro_p1,ro_m1,ro_mp1,ro_mm1,ro_nmm1) \
            collapse(3) schedule(static)
  for(i=1; i<IMAX-1; i++) {
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

          CD = (vel[n][ro] +
                vel[n][ro_p1]) / 2;
          
          if (CD >= 0) 
            upwind = vel[n][ro]; 
          else
            upwind = vel[n][ro_p1];

          if(af[n][ro_p1] > 0.01)
            Q_C += 2/del[n] * H_vel * ((1-solver->alpha)*(CD - vel[n][ro]) +
                                         solver->alpha*(upwind - vel[n][ro]));
          
          /* Flux from cell centered to the west/south/bottom */
          H_vel = (vel[n][ro] * af[n][ro] + 
                   vel[n][ro_m1] * af[n][ro_m1]) / 2;

          CD = (vel[n][ro] +
                vel[n][ro_m1]) / 2;
          
          if (CD >= 0) 
            upwind = vel[n][ro_m1]; 
          else
            upwind = vel[n][ro];

          if(af[n][ro_m1] > 0.01)
            Q_C += 2/del[n] * H_vel * ((1-solver->alpha)*(vel[n][ro] - CD) +
                                         solver->alpha*(vel[n][ro] - upwind));
 
          

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

            CD = (vel[n][ro] + vel[n][ro_mp1])/2;

            if(H_vel >= 0)
              upwind = vel[n][ro];
            else
              upwind = vel[n][ro_mp1];

            if ((af[m][ro]>0.01 || af[m][ro_p1]>0.01) && 
                 af[n][ro_mp1]>0.01)
              Q_W += 2/del[m] * H_vel * ((1-solver->alpha)*(CD - vel[n][ro]) + 
                                          solver->alpha*(upwind - vel[n][ro]));
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

            CD = (vel[n][ro] + vel[n][ro_mm1])/2;

            if(H_vel >= 0)
              upwind = vel[n][ro_mm1];
            else
              upwind = vel[n][ro];

            if ((af[m][ro_mm1]>0.01 || af[m][ro_nmm1]>0.01) &&
                 af[n][ro_mm1]>0.01)
              Q_W += 2/del[m] * H_vel * ((1-solver->alpha)*(vel[n][ro] - CD) + 
                                          solver->alpha*(vel[n][ro] - upwind));

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
          resi   = (DN(i,j,k)  +  DN(i+odim[n][0],j+odim[n][1],k+odim[n][2])) / 2;
          switch(n) {
          case 0:
            resi  *= UN(i,j,k);
            break;
          case 1:
            resi *= VN(i,j,k);
            break;
          case 2:
            resi *= WN(i,j,k);
            break;
          }
          if(solver->t < solver->emf || solver->iter > solver->niter || solver->p_flag==1 ) resi = 0.0;

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
                 Flux + Viscocity /* - resi /solver->rho */ );  /* uncomment to use residual as volume source */

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


int vof_loop(struct solver_data *solver) {
  double t_n;

  
  mesh_copy_data(mesh_n, solver->mesh);
    
  if(solver->petacal != NULL)
    solver->petacal(solver);
  
  solver->boundaries(solver);
  if(solver->special_boundaries != NULL)
    solver->special_boundaries(solver);
  
  if(solver->betacal != NULL)
    solver->betacal(solver);

  /*vof_hydrostatic(solver);*/
  
  solver->turbulence_loop(solver);

  while(solver->t < solver->endt) {
    
      
    solver->iter = 0;
    solver->resimax = 0;
    solver->nu_max = solver->nu;
    solver->delt_n = solver->delt;

    solver->velocity(solver);
    
    /* printf("velocity\n");
    track_cell(solver, 17,2,1); */
    
    solver->wall_shear(solver); 

    solver->boundaries(solver);
    
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);
    /* printf("boundaries\n");
    track_cell(solver, 17,2,1); */

    while(solver->iter < solver->niter) {
    
      if(solver->pressure(solver) == 0) break;

      solver->boundaries(solver);
      if(solver->special_boundaries != NULL)
        solver->special_boundaries(solver);
    
      solver->iter++;
    }
    /* vof_pressure_test(solver); */
    
    t_n = solver->t;
    solver->t += solver->delt;

    /* printf("pressure\n");
    track_cell(solver, 17,2,1); */

    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);
      
    solver->turbulence_loop(solver);
    solver->vfconv(solver);
    
    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver); 

    if(solver->petacal != NULL)
      solver->petacal(solver);
      
    if(solver->deltcal != NULL) {
      if(solver->deltcal(solver) == 0) 
        mesh_copy_data(mesh_n, solver->mesh);
      else 
        solver->t = t_n;
    }

    /* printf("boundaries\n");
    track_cell(solver, 17,2,1);  */
    solver->output(solver);
         
    solver->write(solver); 
  }

  solver->kill(solver);
  solver->turbulence_kill(solver);

  return 0;
}

int vof_output(struct solver_data *solver) {
	time_t current;
	double elapsed;

  if(solver->iter >= solver->niter) {
    printf("timestep: %lf | delt: %lf | pressure did not converge\n", solver->t, solver->delt);
  }
  else {
    printf("timestep: %lf | delt %lf | convergence in %ld iterations.\n", solver->t, solver->delt, solver->iter);
  }
  printf("Max residual %lf | epsi %lf", solver->resimax, solver->epsi);
  if(solver->pressure == vof_pressure)
  	printf(" | omega %lf", solver->omg_final);
  printf("\n");
  
  printf("max u, v, w: %lf, %lf, %lf\n",solver->umax,solver->vmax,solver->wmax);  
  
  printf("max nu: %lf\n",solver->nu_max);
  
  vof_baffles_output(solver);
  
  time(&current);
  elapsed = difftime(current, solver->start_time);
  printf("Elapsed time: %.0lf s | Timesteps per second: %lf\n",elapsed,solver->t / elapsed);
  
  
  printf("\n");
  
  return 0;
}

int vof_deltcal(struct solver_data *solver) {
  double delt, delt_conv, dt_U, alpha, dp, dv;
  int ret = 0;
  double dtvis;
  double mindx;
  long int i,j,k;
  int nan_flag = 0;
  
  #ifdef DEBUG
  long int umax_cell[3], vmax_cell[3], wmax_cell[3];
  #endif
  
  
  solver->umax = 0;
  solver->vmax = 0;
  solver->wmax = 0;
  /* dtvis = 0.5 * (pow(DELX,2) * pow(DELY,2) * pow(DELZ,2)) / 
                       (solver->nu_max * (pow(DELX,2) + pow(DELY,2) + pow(DELZ,2)) ); */
  
  dtvis = 0.25  / (solver->nu_max * (1 / pow(DELX,2) + 1 / pow(DELY,2) + 1 /pow(DELZ,2)) );
  
  if(solver->t < solver->emf) return 0;
  
  delt = solver->delt_n;
  
  
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {        
        if(isnan(U(i,j,k)) || isnan(V(i,j,k)) || isnan(W(i,j,k))) {
          nan_flag = 1;
          printf("divide by zero error in cell %ld %ld %ld | ",i,j,k);
          if(isnan(U(i,j,k))) printf("east");
          if(isnan(V(i,j,k))) printf("north");
          if(isnan(W(i,j,k))) printf("top");
          printf(" of cell\nwriting previous timestep and exiting\n");
        }
      
#ifdef DEBUG
        if(fabs(U(i,j,k)) > solver->umax) {
          umax_cell[0]=i; umax_cell[1]=j; umax_cell[2]=k;
        }
        if(fabs(V(i,j,k)) > solver->vmax) {
          vmax_cell[0]=i; vmax_cell[1]=j; vmax_cell[2]=k;
        }
        if(fabs(W(i,j,k)) > solver->wmax) {
          wmax_cell[0]=i; wmax_cell[1]=j; wmax_cell[2]=k;
        }        
#endif
        
        solver->umax = max(fabs(U(i,j,k)), solver->umax);
        solver->vmax = max(fabs(V(i,j,k)), solver->vmax);
        solver->wmax = max(fabs(W(i,j,k)), solver->wmax);
      }
    }
  }

#ifdef DEBUG
  printf("Max u: %lf in cell %ld %ld %ld\n", solver->umax, umax_cell[0], umax_cell[1], umax_cell[2]);
  printf("Max v: %lf in cell %ld %ld %ld\n", solver->vmax, vmax_cell[0], vmax_cell[1], vmax_cell[2]);
  printf("Max w: %lf in cell %ld %ld %ld\n", solver->wmax, wmax_cell[0], wmax_cell[1], wmax_cell[2]);  
#endif

  if(solver->vof_flag == 1) {
    delt = solver->delt * 0.67;
    ret = 1;
    solver->con *= 0.975;
 
#pragma omp parallel for shared(solver) private(i,j,k)
    for(i=1; i<IMAX-1; i++) {
      for(j=1; j<JMAX-1; j++) {
        for(k=1; k<KMAX-1; k++) {
          U(i,j,k) = 0;
          V(i,j,k) = 0;
          W(i,j,k) = 0;
          P(i,j,k) = PN(i,j,k);
          VOF(i,j,k) = VOF_N(i,j,k);
          N_VOF(i,j,k) = N_VOF_N(i,j,k);
          
        }
      }
    }
              
  }
  
  if(nan_flag == 1) { /* divide by zero error - exit */
    for(i=1; i<IMAX-1; i++) {
      for(j=1; j<JMAX-1; j++) {
        for(k=1; k<KMAX-1; k++) {
          U(i,j,k) = UN(i,j,k);
          V(i,j,k) = VN(i,j,k);
          W(i,j,k) = WN(i,j,k);
          P(i,j,k) = PN(i,j,k);
          VOF(i,j,k) = VOF_N(i,j,k);
          N_VOF(i,j,k) = N_VOF_N(i,j,k);
          
        }
      }
    }    
    vof_write_timestep(solver);
    exit(1);
  }
  
  if(solver->iter > 60) delt *= 0.99; 
  if(solver->iter < 20) delt *= 1.01; 

  dv = 0;
  delt_conv = delt * 100;
  for(i=0; i<IMAX; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {  
        /*
        if((N_VOF(i,j,k) > 7 && N_VOF(i+1,j,k) > 0) || (N_VOF(i,j,k) > 0 && N_VOF(i+1,j,k) > 7)) {
          if(U(i,j,k) > 0) { // && N_VOF(i+1,j,k) > 7) {
            dt_U = 0.5 * FV(i+1,j,k) * DELX * FV(i,j,k) / (fabs(U(i,j,k)) * AE(i,j,k) * VOF(i,j,k));
          } else if (U(i,j,k) < 0) { // && N_VOF(i,j,k) > 7) {
            dt_U = 0.5 * FV(i+1,j,k) * DELX * FV(i,j,k) / (fabs(U(i,j,k)) * AE(i,j,k) * VOF(i+1,j,k));
          } 
        }
        dp = P(i,j,k) - P(min(IMAX-1,i+1),j,k);
        if(i < IMAX-1) dp += solver->gx * solver->rho * DELX;
        if(dp * U(i,j,k) > solver->emf) dv = solver->delt * ( (1/DELX) * dp / solver->rho);
        else dv = 0;
        dv *= 1.25; */
      
        dt_U = solver->con * min(1, min(FV(i,j,k),FV(min(IMAX-1,i+1),j,k)) / AE(i,j,k)) * DELX/(fabs(dv + U(i,j,k)));
      //}
        if(AE(i,j,k) > solver->emf && !isnan(dt_U))
          delt_conv = min(delt_conv, dt_U);
/*
        if((N_VOF(i,j,k) > 7 && N_VOF(i,j+1,k) > 0) || (N_VOF(i,j,k) > 0 && N_VOF(i,j+1,k) > 7)) {
          if(V(i,j,k) > 0) { // && N_VOF(i,j+1,k) > 7) {
            dt_U = 0.5 * FV(i,j+1,k) * DELY * FV(i,j,k) / (fabs(V(i,j,k)) * AN(i,j,k) * VOF(i,j,k));
          } else if (V(i,j,k) < 0) { // && N_VOF(i,j,k) > 7) {
            dt_U = 0.5 * FV(i,j+1,k) * DELY * FV(i,j,k) / (fabs(V(i,j,k)) * AN(i,j,k) * VOF(i,j+1,k));
          } 
        } 
        dp = P(i,j,k) - P(i, min(JMAX-1,j+1), k);
        if(j < JMAX-1) dp += solver->gy * solver->rho * DELY;
        if(dp * V(i,j,k) > solver->emf) dv = solver->delt * ( (1/DELY) * dp / solver->rho); // dv is additive to calc stability limit based on pressure gradient effect on velocity 
        else dv = 0;
        dv *= 1.25; */
              
        dt_U = solver->con * min(1, min(FV(i,j,k),FV(i,min(JMAX-1,j+1),k)) / AN(i,j,k)) * DELY/(fabs(dv + V(i,j,k)));
      //}
        if(AN(i,j,k) > solver->emf && !isnan(dt_U))
          delt_conv = min(delt_conv, dt_U);
          /*
        if((N_VOF(i,j,k) > 7 && N_VOF(i,j,k+1) > 0) || (N_VOF(i,j,k) > 0 && N_VOF(i,j,k+1) > 7)) {					
          if(W(i,j,k) > 0) { // && N_VOF(i,j,k+1) > 7) {
            dt_U = 0.5 * FV(i,j,k+1) * DELZ * FV(i,j,k) / (fabs(W(i,j,k)) * AT(i,j,k) * VOF(i,j,k));
          } else if (W(i,j,k) < 0) { // && N_VOF(i,j,k) > 7) {
            dt_U = 0.5 * FV(i,j,k+1) * DELZ * FV(i,j,k) / (fabs(W(i,j,k)) * AT(i,j,k) * VOF(i,j,k+1));
          } 
        } 
        dp = P(i,j,k) - P(i,j, min(KMAX-1,k+1));
        if(k < KMAX-1) dp += solver->gz * solver->rho * DELZ;
        if(dp * W(i,j,k) > solver->emf) dv = solver->delt * ( (1/DELZ) * dp / solver->rho); // dv is additive to calc stability limit based on pressure gradient effect on velocity 
        else dv = 0;		
        dv *= 1.25;		*/
      
        dt_U = solver->con * min(1, min(FV(i,j,k),FV(i,j,min(KMAX-1,k+1))) / AT(i,j,k)) * DELZ/(fabs(dv + W(i,j,k)));
      //}
        if(AT(i,j,k) > solver->emf && !isnan(dt_U))
          delt_conv = min(delt_conv, dt_U);					
					
			}
		}
	}
  printf("maximum timestep for convective stability: %lf\n",delt_conv);
  
  delt = min(delt, delt_conv);
  
	delt = min(delt, 0.8 * dtvis);
   
  if(solver->delt_n != delt) {
    printf("timestep adjusted from %lf to %lf\n",solver->delt_n,delt);

    solver->delt = delt;

    if(delt < solver->delt_min) {
      printf("timestep too small to continue.  writing current timestep and exiting...");
      vof_write_timestep(solver);
      exit(1);
    }
  }

/*   Leaves a ridiculous amount of files on your drive
#ifdef DEBUG
  if(delt / solver->delt_n < 0.7 && solver->vof_flag != 1) {
    printf("large change in timestep - saving\n");
    vof_write_timestep(solver);
  }
#endif */
  
  /* now adjust alpha */
  alpha = solver->alpha;
  alpha = max(alpha, 1.5*solver->umax*solver->delt/DELX);
  alpha = max(alpha, 1.5*solver->vmax*solver->delt/DELY);
  alpha = max(alpha, 1.5*solver->wmax*solver->delt/DELZ);
  alpha = min(alpha, 1.0);
  /* printf("alpha adjusted from %lf to %lf\n",solver->alpha, alpha); */
  solver->alpha = alpha;
  
  /* adjust epsi */
  /* solver->epsi = solver->epsi * solver->delt / solver->delt_n; */
  /* solver->epsi = 0.0001 / solver->delt; */
  mindx = min(DELX,DELY);
  mindx = min(DELY,DELZ);
  solver->epsi = 1 * solver->rho * mindx / (solver->delt * solver->dzro);
                /* (solver->dzro * 0.01 * max(pow(solver->delt * 100, 1.2),0.001)); */
  
  delt_conv = delt_conv * 0.5 / solver->con;
  solver->omg = solver->omg_init * delt / delt_conv + (1 - delt / delt_conv);
  solver->betacal(solver);
  if(solver->petacal != NULL)
    solver->petacal(solver);
  
  return ret;
}


int vof_write(struct solver_data *solver) {
  static double write_flg = 0;
    
  if(solver->t >= write_flg) {
  
		vof_write_timestep(solver);
    
    write_flg = solver->t + solver->writet;
  }

  return 0;
}

int vof_write_timestep(struct solver_data * solver) {
	int write_step;
	
	write_step = track_add(solver->t);
	track_write();
	
	vtk_write_P(solver->mesh,write_step);
	vtk_write_U(solver->mesh,write_step);
	vtk_write_vof(solver->mesh,write_step);
	vtk_write_n_vof(solver->mesh,write_step);
	
	csv_write_P(solver->mesh,solver->t);
	csv_write_U(solver->mesh,solver->t);
	csv_write_vof(solver->mesh,solver->t);
	csv_write_n_vof(solver->mesh,solver->t);
	
	if(solver->mesh->turbulence_model != NULL) {
		vtk_write_k(solver->mesh,write_step);
		vtk_write_E(solver->mesh,write_step);
		csv_write_k(solver->mesh,solver->t);
		csv_write_E(solver->mesh,solver->t);
	}
	
	vof_baffles_write(solver);
	
	vof_vorticity(solver);
	
	vtk_write_vorticity(solver->mesh,write_step);
	csv_write_vorticity(solver->mesh,solver->t);
	
	return 1;
}

