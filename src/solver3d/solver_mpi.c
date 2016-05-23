/*
 * solver_mpi.c
 *
 * extends solver.c to allow mpi paralellism
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <petscksp.h>
#include <mpi.h>

#include "mesh.h"
#include "readfile.h"
#include "solver.h"
#include "readsolver.h"
#include "laminar.h"
#include "kE.h"
#include "vof.h"
#include "mesh_mpi.h"
#include "vof_macros.h"

int solver_mpi_range(struct solver_data *solver) {
	long int range, start;
	int size, rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  range = (IMAX + (size - 1)) / size;
  start = range * rank;
  start -= 1;
  start = min(start, 0);
  range += 1;
	if(start + range > IMAX) range = IMAX - start;
  solver->mesh->i_range = range;
  solver->mesh->i_start = start;
}

int solver_mpi(struct solver_data *solver, double timestep, double delt)
{
	long int range, start;
	int size, rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if(rank > 0) {
    if(solver_mpi_high_rank(solver)) {
      return 1;
    }
    return 0;
  }
  
  solver = solver_init_empty();
  if(solver == NULL) return 1;

  vof_setup_solver(solver);
  
  if(solver_load(solver, "solverfile", "meshfile", "initials")==1)
    return 1;

  if(solver_init_complete(solver)==1)
    return 1;

  if(solver->turbulence_read != NULL) solver->turbulence_read("turbulencefile");
  
  solver->init(solver);
  solver->turbulence_init(solver);
  
  if(mesh_load_csv(solver->mesh, 0) == 1) return 1;
  solver_initial_values(solver);
  solver->petacal(solver);

  if(timestep > solver->emf) {
    solver->t = timestep;
    csv_read_U_p_vof(solver->mesh, timestep);
    solver->turbulence_load_values(solver);
  }
  else
    solver->write(solver); 
  
  if (delt > solver->emf) solver->delt = delt;
  
  track_read();
  
	solver_broadcast_all(solver);
	mesh_broadcast_all(solver->mesh);
  
  solver_mpi_range(solver);
  
  if(kE_check(solver)) kE_broadcast(solver);
  
  solver->send_all(solver);
  
  if(solver_run(solver)==1)
    return 1;

  return(0);

}

int solver_mpi_high_rank(struct solver_data *solver) {
	long int range, start;
	int size, rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	solver = solver_init_empty();
	if(solver == NULL) return 1;
	
	vof_setup_solver(solver);
	
	solver_broadcast_all(solver);
	mesh_broadcast_all(solver->mesh);
	
  if(solver_check(solver) == 1) {
    return(1);
  }
  solver->ready = 1;
  
  solver_mpi_range(solver);
  
  mesh_mpi_init_complete(solver->mesh, range, start);
  
  if(kE_check(solver)) kE_broadcast(solver);
  
  solver->init(solver);
  solver->turbulence_init(solver);
	
	solver_recv_all(solver);
	
  if(solver_run(solver)==1)
    return 1;

  return(0);
}


int solver_recv_all(struct solver_data *solver) {
	struct kE_data *kE;
	
	solver_mpi_recv(solver, solver->mesh->fv, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->ae, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->an, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->at, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->vof, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->P, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->u, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->v, 0, 0, IRANGE);
	solver_mpi_recv(solver, solver->mesh->w, 0, 0, IRANGE);
	
	if(kE_check(solver))  {
		kE = solver->mesh->turbulence_model;
		solver_mpi_recv(solver, kE->k, 0, 0, IRANGE);
	}
}

int solver_send_recv_edge(struct solver_data *solver) {
  /* each process will send it's eastern most edge to the next process */
	long int range, start, rank;
	int size, n;
	struct kE_data *kE;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(!rank) {
    solver_mpi_send(solver, solver->mesh->vof, 1, IRANGE-1, 1);
    solver_mpi_send(solver, solver->mesh->P, 1, IRANGE-1, 1);
    solver_mpi_send(solver, solver->mesh->u, 1, IRANGE-1, 1);
    solver_mpi_send(solver, solver->mesh->v, 1, IRANGE-1, 1);
    solver_mpi_send(solver, solver->mesh->w, 1, IRANGE-1, 1);		
		if(kE_check(solver))  {
			kE = solver->mesh->turbulence_model;
      solver_mpi_send(solver, kE->k, 1, IRANGE-1, 1);		
		}
  } else if(rank + 1 < size) {  
    /* internal cases */
    solver_mpi_sendrecv(solver, rank + 1, solver->mesh->vof, IRANGE-1, 1, 
                                rank - 1, solver->mesh->vof, 0, 1);
    solver_mpi_sendrecv(solver, rank + 1, solver->mesh->P, IRANGE-1, 1, 
                                rank - 1, solver->mesh->P, 0, 1);
    solver_mpi_sendrecv(solver, rank + 1, solver->mesh->u, IRANGE-1, 1, 
                                rank - 1, solver->mesh->u, 0, 1);
    solver_mpi_sendrecv(solver, rank + 1, solver->mesh->v, IRANGE-1, 1, 
                                rank - 1, solver->mesh->v, 0, 1);
    solver_mpi_sendrecv(solver, rank + 1, solver->mesh->w, IRANGE-1, 1, 
                                rank - 1, solver->mesh->w, 0, 1);
		if(kE_check(solver))  {
			kE = solver->mesh->turbulence_model;
      solver_mpi_sendrecv(solver, rank + 1, kE->k, IRANGE-1, 1, 
                                  rank - 1, kE->k, 0, 1);
		}
   } else {
    solver_mpi_recv(solver, solver->mesh->vof, rank - 1, 0, 1);	
    solver_mpi_recv(solver, solver->mesh->P, rank - 1, 0, 1);	
    solver_mpi_recv(solver, solver->mesh->u, rank - 1, 0, 1);	
    solver_mpi_recv(solver, solver->mesh->v, rank - 1, 0, 1);	
    solver_mpi_recv(solver, solver->mesh->w, rank - 1, 0, 1);	
		if(kE_check(solver))  {
			kE = solver->mesh->turbulence_model;
      solver_mpi_recv(solver, kE->k, rank - 1, 0, 1);	
		}
   }
  
}

int solver_send_all(struct solver_data *solver) {
	long int range, start, rank;
	int size, n;
	struct kE_data *kE;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	range = (IMAX + (size - 1)) / size;
	
	for(n=0; n<size; n++) {
		start = range * n;
		start -= 1;	
		
		if(start + range + 1 > IMAX) range = IMAX - start - 1;
		
		solver_mpi_send(solver, solver->mesh->fv, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->ae, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->an, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->at, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->vof, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->P, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->u, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->v, n, start, range + 1);
		solver_mpi_send(solver, solver->mesh->w, n, start, range + 1);
		
		if(kE_check(solver))  {
			kE = solver->mesh->turbulence_model;
			solver_mpi_send(solver, kE->k, n, start, range + 1);
		}
	}
}

int solver_mpi_send(struct solver_data *solver, double *data, int to, long int i_start, long int i_range) {
	MPI_SEND(data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_DOUBLE, to, 1, MPI_COMM_WORLD);
}

int solver_mpi_recv(struct solver_data *solver, double *data, int from, long int i_start, long int i_range) {
	MPI_RECV(data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int solver_mpi_sendrecv(struct solver_data *solver, int to, double *send, long int send_i_start, long int send_i_range, int from, double *recv, long int recv_i_start, long int recv_i_range) {
	MPI_Sendrecv(send[mesh_index(solver->mesh,send_i_start,0,0)], send_i_range * JMAX * KMAX, 
               MPI_DOUBLE, to, 1,
               recv[mesh_index(solver->mesh,recv_i_start,0,0)], recv_i_range * JMAX * KMAX,
               MPI_DOUBLE, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int solver_broadcast_all(struct solver_data *solver) {
	int rank, ic, turb, autot;
	
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	MPI_Bcast(&solver->gx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->gy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->gz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->rho, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->nu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->endt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->writet, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&solver->delt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if(!rank) {
		if(kE_check(solver)) turb = 1;
		else turb = 0;
		MPI_Bcast(&turb, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
		if(solver->deltcal != NULL) autot = 1;
		else autot = 0;
		MPI_Bcast(&autot, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Bcast(&turb, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&autot, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		if(turb) kE_setup(solver);
		else laminar_setup(solver);
		
		if(!autot) solver->deltcal = NULL;
	}
	
	return 0;
}