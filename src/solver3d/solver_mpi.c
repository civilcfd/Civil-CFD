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

int solver_mpi(struct solver_data *solver) {
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
  
  /* CODE FOR RANGE */
  range = (IMAX + (size - 1)) / size;
  start = range * rank;
  start -= 1;
  range += 1;
	if(start + range > IMAX) range = IMAX - start;
  
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
	
	solver_mpi_recv(solver, solver->mesh->fv, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->ae, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->an, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->at, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->vof, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->P, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->u, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->v, 0, ISTART, IRANGE);
	solver_mpi_recv(solver, solver->mesh->w, 0, ISTART, IRANGE);
	
	if(kE_check(solver))  {
		kE = solver->mesh->turbulence_model;
		solver_mpi_recv(solver, kE->k, 0, ISTART, IRANGE);
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
		
		solver_mpi_send(solver, solver->mesh->fv, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->ae, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->an, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->at, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->vof, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->P, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->u, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->v, 0, start, range + 1);
		solver_mpi_send(solver, solver->mesh->w, 0, start, range + 1);
		
		if(kE_check(solver))  {
			kE = solver->mesh->turbulence_model;
			solver_mpi_send(solver, kE->k, 0, start, range + 1);
		}
	}
}

int solver_mpi_send(struct solver_data *solver, double *data, int to, long int i_start, long int i_range) {
	MPI_SEND(data[i_start * JMAX * KMAX], i_range * JMAX * KMAX, MPI_DOUBLE, to, 1, MPI_COMM_WORLD);
}

int solver_mpi_recv(struct solver_data *solver, double *data, int from, long int i_start, long int i_range) {
	MPI_RECV(data[i_start * JMAX * KMAX], i_range * JMAX * KMAX, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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