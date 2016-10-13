/* vof.c
 *
 * default solver based on VOF and fractional area/volume obstacles
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

void track_cell(struct solver_data *solver, long int i,long int j,long int k) {
  printf("U(e/n/t) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i,j,k),V(i,j,k),W(i,j,k));
  printf("U(w/s/b) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i-1,j,k),V(i,j-1,k),W(i,j,k-1));  
  printf("P: %lf    VOF: %lf    N_VOF: %d\n",P(i,j,k),VOF(i,j,k),N_VOF(i,j,k));
}
#endif

int vof_mpi_setup_solver(struct solver_data *solver) {
  
  solver->init = vof_mpi_init_solver;
  solver->kill = vof_mpi_kill_solver;
  solver->loop = vof_mpi_loop;
  solver->boundaries = vof_boundaries;
  solver->special_boundaries = vof_special_boundaries;
  solver->pressure = vof_pressure_gmres_mpi;
  solver->velocity = vof_mpi_velocity_upwind;
  solver->convect = vof_mpi_convect; //_depreciated;
  solver->nvof = vof_mpi_nvof;
  solver->deltcal = vof_mpi_deltcal;
  solver->write = vof_mpi_write;
  solver->output = vof_mpi_output;
  
  return 0;
}

int vof_mpi_init_solver(struct solver_data *solver) {
 
  mesh_n = mesh_mpi_init_copy(solver->mesh);
  if(mesh_n == NULL)
    return 1;
  
  return 0;
}

int vof_mpi_kill_solver(struct solver_data *solver) {
  mesh_free(solver->mesh);
  PetscEnd();

  exit(0);
}

int vof_mpi_loop(struct solver_data *solver) {
  double t_n;

  
  mesh_mpi_copy_data(mesh_n, solver->mesh);

  if(solver->nvof != NULL)
    solver->nvof(solver); 
  solver_sendrecv_edge_int(solver, solver->mesh->n_vof);
  
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
    
    solver->wall_shear(solver); 

    solver->boundaries(solver);
    
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);

    solver_sendrecv_edge(solver, solver->mesh->u);
    solver_sendrecv_edge(solver, solver->mesh->v);
    solver_sendrecv_edge(solver, solver->mesh->w);


    solver->pressure(solver);
    if(solver->iter < 1) solver->iter = 1;

    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);
  
    
    t_n = solver->t;
    solver->t += solver->delt;

    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);

    solver_sendrecv_edge(solver, solver->mesh->u);
    solver_sendrecv_edge(solver, solver->mesh->v);
    solver_sendrecv_edge(solver, solver->mesh->w); 
    solver_sendrecv_edge(solver, solver->mesh->P);

    solver->turbulence_loop(solver);

    solver->convect(solver);
    
    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver); 
    solver_sendrecv_edge(solver, solver->mesh->vof);
      
    if(solver->deltcal != NULL) {
      if(solver->deltcal(solver) == 0) 
        mesh_mpi_copy_data(mesh_n, solver->mesh);
      else 
        solver->t = t_n;
    }

    if(solver->nvof != NULL)
      solver->nvof(solver);
    solver_sendrecv_edge_int(solver, solver->mesh->n_vof);

    solver->output(solver);
         
    solver->write(solver); 
  }

  solver->turbulence_kill(solver);
  solver->kill(solver);

  return 0;
}

int vof_mpi_output(struct solver_data *solver) {
  time_t current;
  double elapsed;

  if(solver->rank > 0) return 0;

  if(solver->iter >= solver->niter) {
    printf("timestep: %lf | delt: %lf | pressure did not converge\n", solver->t, solver->delt);
  }
  else {
    printf("timestep: %lf | delt %lf | convergence in %ld iterations.\n", solver->t, solver->delt, solver->iter);
  }
  printf("Max residual %lf | Convergence reason: %s", solver->resimax, solver->conv_reason_str);
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

int vof_mpi_deltcal(struct solver_data *solver) {
  double delt, delt_conv, dt_U, dv;
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
  
  dtvis = 0.25  / (solver->nu_max * (1 / pow(DELX,2) + 1 / pow(DELY,2) + 1 /pow(DELZ,2)) );
  
  if(solver->t < solver->emf) return 0;
  
  delt = solver->delt_n;
  
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {        
        if(isnan(U(i,j,k)) || isnan(V(i,j,k)) || isnan(W(i,j,k))) {
          nan_flag = 1;
          if(!solver->rank) printf("divide by zero error in cell %ld %ld %ld | ",i,j,k);
          if(isnan(U(i,j,k))) printf("east");
          if(isnan(V(i,j,k))) printf("north");
          if(isnan(W(i,j,k))) printf("top");
          if(!solver->rank) printf(" of cell\nwriting previous timestep and exiting\n");
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

  solver->umax = solver_mpi_max(solver,solver->umax);
  solver->vmax = solver_mpi_max(solver,solver->vmax);
  solver->wmax = solver_mpi_max(solver,solver->wmax);

  solver->vof_flag = solver_mpi_max(solver, solver->vof_flag);
  if(solver->vof_flag == 1) {
    delt = solver->delt * 0.67;
    ret = 1;
    solver->con *= 0.975;
 
    for(i=1; i<IRANGE-1; i++) {
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
    for(i=1; i<IRANGE-1; i++) {
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
    vof_mpi_write_timestep(solver);
    PetscEnd();
    exit(1);
  }
  
  if(solver->iter > 150) delt *= 0.975; 
  if(solver->iter < 100) delt *= 1.025; 

  dv = 0;
  delt_conv = delt * 100;
  for(i=0; i<IRANGE; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {  
        dt_U = solver->con * min(1, min(FV(i,j,k),FV(min(IRANGE-1,i+1),j,k)) / AE(i,j,k)) * DELX/(fabs(dv + U(i,j,k)));

        if(AE(i,j,k) > solver->emf && !isnan(dt_U))
          delt_conv = min(delt_conv, dt_U);
              
        dt_U = solver->con * min(1, min(FV(i,j,k),FV(i,min(JMAX-1,j+1),k)) / AN(i,j,k)) * DELY/(fabs(dv + V(i,j,k)));

        if(AN(i,j,k) > solver->emf && !isnan(dt_U))
          delt_conv = min(delt_conv, dt_U);

        dt_U = solver->con * min(1, min(FV(i,j,k),FV(i,j,min(KMAX-1,k+1))) / AT(i,j,k)) * DELZ/(fabs(dv + W(i,j,k)));

        if(AT(i,j,k) > solver->emf && !isnan(dt_U))
          delt_conv = min(delt_conv, dt_U);					
 #ifdef DEBUG         
        if(delt_conv < 0.0001) {

          printf("");

        }
  #endif
      }
    }
  }

  if(!solver->rank) printf("maximum timestep for convective stability: %lf\n",delt_conv);
  
  delt = min(delt, delt_conv);
  delt = min(delt, 0.8 * dtvis);

  delt = solver_mpi_min(solver, delt);
   
  if(solver->delt_n != delt) {
    if(!solver->rank) printf("timestep adjusted from %lf to %lf\n",solver->delt_n,delt);

    solver->delt = delt;

    if(delt < solver->delt_min) {
      if(!solver->rank) printf("timestep too small to continue.  writing current timestep and exiting...");
      vof_mpi_write_timestep(solver);
      vof_mpi_kill_solver(solver);
    }
  }
  
  mindx = min(DELX,DELY);
  mindx = min(DELY,DELZ);
  
  delt_conv = delt_conv * 0.5 / solver->con;
  
  return ret;
}


int vof_mpi_write(struct solver_data *solver) {
  static double write_flg = 0;
    
  if(solver->t >= write_flg) {
  
    vof_mpi_write_timestep(solver);
    
    write_flg = solver->t + solver->writet;
  }

  return 0;
}

int vof_mpi_write_timestep(struct solver_data * solver) {
  int write_step;
  struct kE_data *kE;
  
  solver_mpi_gather(solver, solver->mesh->P);
  solver_mpi_gather(solver, solver->mesh->u);
  solver_mpi_gather(solver, solver->mesh->v);
  solver_mpi_gather(solver, solver->mesh->w);
  solver_mpi_gather(solver, solver->mesh->vof);
  solver_mpi_gather_int(solver, solver->mesh->n_vof);

  if(solver->mesh->turbulence_model != NULL) {
    kE = solver->mesh->turbulence_model;
    solver_mpi_gather(solver, kE->k);
    solver_mpi_gather(solver, kE->E);
  }

  if(!solver->rank) {
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
  }

  vof_baffles_write(solver);
  
  vof_vorticity(solver);
  
  if(!solver->rank) {
    vtk_write_vorticity(solver->mesh,write_step);
    csv_write_vorticity(solver->mesh,solver->t);
  }
  
  return 1;
}

