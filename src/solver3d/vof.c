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
float pp(struct solver_data *solver, long int i,long int j,long int k) { return P(i,j,k); }
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
  solver->convect = vof_convect;
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
  static int omg = 0;
  
  if(omg == solver->omg) return 0;
  else omg = solver->omg;

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

int vof_loop(struct solver_data *solver) {
  double t_n;

  
  mesh_copy_data(mesh_n, solver->mesh);

  if(solver->betacal != NULL)
    solver->betacal(solver);    
  if(solver->petacal != NULL)
    solver->petacal(solver); 
  
  solver->boundaries(solver);
  if(solver->special_boundaries != NULL)
    solver->special_boundaries(solver);
  

  /*vof_hydrostatic(solver);*/
  
  solver->turbulence_loop(solver);

  while(solver->t < solver->endt) {
    
      
    solver->iter = 0;
    solver->resimax = 0;
    solver->nu_max = solver->nu;
    solver->delt_n = solver->delt;
    solver->vchgt  = 0;

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
    solver->convect(solver);
    
    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver); 

      
    if(solver->deltcal != NULL) {
      if(solver->deltcal(solver) == 0) 
        mesh_copy_data(mesh_n, solver->mesh);
      else 
        solver->t = t_n;
    }

    if(solver->betacal != NULL)
      solver->betacal(solver);
    if(solver->petacal != NULL)
      solver->petacal(solver);
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
  
  if(solver->vchgt > solver->emf) {
    printf("Fluid volume lost: %lf L | Flow change: %lf L/s\n",solver->vchgt*1000,solver->vchgt*1000/solver->delt);
  }
  
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
  
  if(solver->pressure == vof_pressure_gmres) {
    if(solver->iter > 150) delt *= 0.975; 
    if(solver->iter < 100) delt *= 1.025; 
  } else {
    if(solver->iter > 60) delt *= 0.99; 
    if(solver->iter < 20) delt *= 1.01; 
  }

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

