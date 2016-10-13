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
#include "mesh_mpi.h"
#include "vof_macros.h"
#include "solver_mpi.h"
#include "vof_mpi.h"
#include "csv.h"
#include "track.h"

int solver_mpi_range(struct solver_data *solver) {
  long int range, start;
  int size, rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  range = (IMAX + (size - 1)) / size;
  start = range * rank;
  start -= 1;
  start = max(start, 0);
  range += 2;
  if(start + range > IMAX) range = IMAX - start;
  if(!rank) range--;
  solver->mesh->i_range = range;
  solver->mesh->i_start = start;

  return 0;
}

int solver_mpi_init_comm(struct solver_data *solver) {
  int colour;
  
  if(solver->rank % 2 == 0) { /* even rank */
    colour = solver->rank;

    if(solver->rank + 1 < solver->size)
      MPI_Comm_split(MPI_COMM_WORLD, colour, solver->rank, &solver->comm_downstream);
    else 
      MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, solver->rank, &solver->comm_downstream);

    colour--;
    if(colour > 0)
      MPI_Comm_split(MPI_COMM_WORLD, colour, solver->rank, &solver->comm_upstream);
    else
      MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, solver->rank, &solver->comm_upstream);

  } else {
    colour = solver->rank - 1;

    MPI_Comm_split(MPI_COMM_WORLD, colour, solver->rank, &solver->comm_upstream);
    
    colour++;
    if(solver->rank + 1 < solver->size)
      MPI_Comm_split(MPI_COMM_WORLD, colour, solver->rank, &solver->comm_downstream);
    else 
      MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, solver->rank, &solver->comm_downstream);

  }

  return 0;

}

int solver_mpi(struct solver_data *solver, double timestep, double delt)
{
  int size, rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if(rank > 0) {
    if(solver_mpi_high_rank(solver, timestep)) {
      return 1;
    }
    return 0;
  }
  
  solver = solver_init_empty();
  if(solver == NULL) return 1;
  solver->rank = rank;
  solver->size = size;
  solver_mpi_init_comm(solver);

  vof_mpi_setup_solver(solver);
  
  if(solver_load(solver, "solver.xml")==1)
    return 1;

  solver_mpi_range(solver);
  if(solver_mpi_init_complete(solver)==1)
    return 1;

  if(solver->turbulence_read != NULL) solver->turbulence_read("solver.xml");
  
  solver->init(solver);
  solver->turbulence_init(solver);
  
  if(mesh_load_csv(solver->mesh, 0) == 1) return 1;
  solver_initial_values(solver);
  solver->nvof(solver);

  if(timestep > solver->emf) {
    solver->t = timestep;
    csv_read_U_p_vof(solver->mesh, timestep);
    solver->turbulence_load_values(solver);
  }
  
  if (delt > solver->emf) solver->delt = delt;
  
  solver_broadcast_all(solver);
  mesh_broadcast_all(solver->mesh);
  
  if(kE_check(solver)) kE_broadcast(solver);
  
  solver_send_all(solver);
  if(timestep < solver->emf) solver->write(solver);
  track_read();
  
  if(solver_run(solver)==1)
    return 1;

  return(0);

}

int solver_mpi_high_rank(struct solver_data *solver, double timestep) {
  int size, rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  solver = solver_init_empty();
  if(solver == NULL) return 1;
  solver->rank = rank;
  solver->size = size;
  solver_mpi_init_comm(solver);
  
  vof_mpi_setup_solver(solver);
  
  solver_broadcast_all(solver);
  mesh_broadcast_all(solver->mesh);
  solver_mpi_range(solver);
  
  if(solver_check(solver) == 1) {
    return(1);
  }
  solver->ready = 1;
  
  if(solver_mpi_init_complete(solver)==1)
    return 1;
  
  if(kE_check(solver)) kE_broadcast(solver);
  
  solver->init(solver);
  solver->turbulence_init(solver);
  
  solver_recv_all(solver);
  if(timestep < solver->emf) solver->write(solver);
  
  if(solver_run(solver)==1)
    return 1;

  return(0);
}

int solver_mpi_init_complete(struct solver_data *solver) {

  if(solver_check(solver) == 1) {
    return(1);
  }

  if(mesh_mpi_init_complete(solver->mesh) == 1) {
    return(1);
  }

  solver->ready = 1;

  return 0;

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
    solver_mpi_recv(solver, kE->E, 0, 0, IRANGE);
    solver_mpi_recv(solver, kE->nu_t, 0, 0, IRANGE);
    kE_copy(solver);
  }

  return 0;
}

int solver_sendrecv_edge(struct solver_data *solver, double *data) {
  /* each process will communicate its eastern edge with the next process' western edge */
  int rank;
  int size;
  MPI_Request requests[4];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(!rank) {
    solver_mpi_sendrecv(solver, rank + 1, data, IRANGE-2, 1,
                                rank + 1, data, IRANGE-1, 1);
  } else if(rank + 1 < size) {  
    /* internal cases */
    solver_mpi_isend(solver, data, rank-1, 1, 1, &requests[0]);
    solver_mpi_isend(solver, data, rank+1, IRANGE-2, 1, &requests[1]);
    solver_mpi_irecv(solver, data, rank-1, 0, 1, &requests[2]);
    solver_mpi_irecv(solver, data, rank+1, IRANGE-1, 1, &requests[3]);
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

    /* solver_mpi_sendrecv(solver, rank - 1, data, 1, 1, 
                                rank - 1, data, 0, 1);
    solver_mpi_sendrecv(solver, rank + 1, data, IRANGE-2, 1,
                                rank + 1, data, IRANGE-1, 1); */
   } else {
    solver_mpi_sendrecv(solver, rank - 1, data, 1, 1, 
                                rank - 1, data, 0, 1);
   }
  
   return 0;
}

int solver_sendrecv_delu(struct solver_data *solver) {
  double *ds = solver->mesh->delu_downstream;
  double *us = solver->mesh->delu_upstream;

  if(!solver->rank) {
    solver_mpi_sendrecv_replace(solver, ds, 0, 1, 1, 1);
  } else if(solver->rank + 1 < solver->size) {
    if(solver->rank % 2) { /* odd valued cases */ 
      solver_mpi_sendrecv_replace(solver, us, 0, 1, solver->rank - 1, solver->rank - 1);
      solver_mpi_sendrecv_replace(solver, ds, 0, 1, solver->rank + 1, solver->rank + 1);
    } else
    {
      solver_mpi_sendrecv_replace(solver, ds, 0, 1, solver->rank + 1, solver->rank + 1);
      solver_mpi_sendrecv_replace(solver, us, 0, 1, solver->rank - 1, solver->rank - 1);
    }
  } else {
    solver_mpi_sendrecv_replace(solver, us, 0, 1, solver->rank - 1, solver->rank - 1);
  }

  return 0;
}

int solver_sendrecv_edge_int(struct solver_data *solver, int *data) {
  /* each process will communicate its eastern edge with the next process' western edge */
  /* TODO: instead of sendrecv, do series of non-blocking sends followed by receives for faster execution */
  int rank;
  int size;
  MPI_Request requests[4];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(!rank) {
    solver_mpi_sendrecv_int(solver, rank + 1, data, IRANGE-2, 1,
                                 rank + 1, data, IRANGE-1, 1);
  } else if(rank + 1 < size) {  
    /* internal cases */
    solver_mpi_isend_int(solver, data, rank-1, 1, 1, &requests[0]);
    solver_mpi_isend_int(solver, data, rank+1, IRANGE-2, 1, &requests[1]);
    solver_mpi_irecv_int(solver, data, rank-1, 0, 1, &requests[2]);
    solver_mpi_irecv_int(solver, data, rank+1, IRANGE-1, 1, &requests[3]);
    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
    /* internal cases *
    solver_mpi_sendrecv_int(solver, rank - 1, data, 1, 1, 
                                rank - 1, data, 0, 1);
    solver_mpi_sendrecv_int(solver, rank + 1, data, IRANGE-2, 1,
                                rank + 1, data, IRANGE-1, 1); */
   } else {
    solver_mpi_sendrecv_int(solver, rank - 1, data, 1, 1, 
                                    rank - 1, data, 0, 1);
   }

  return 0; 
}

int solver_mpi_gather_int(struct solver_data *solver, int *data) {
  long int range;
  int i;

  static int initialized = 0;
  static int *cnts, *displs;
  int disp = 0;

  range = (IMAX + (solver->size - 1)) / solver->size;
  range *= JMAX;
  range *= KMAX;

  if (!initialized) {
    cnts = malloc(solver->size * sizeof(int));
    displs = malloc(solver->size * sizeof(int));
    for (i = 0; i < solver->size; i++)
    {
      cnts[i] = (i != 0) ? range : 0;
      displs[i] = disp;
      disp += range;
    }
  }

  if (!solver->rank)
    MPI_Gatherv(MPI_IN_PLACE, cnts[solver->rank], MPI_INT,
    data, cnts, displs, MPI_INT,
    0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(data, cnts[solver->rank], MPI_INT,
    data, cnts, displs, MPI_INT,
    0, MPI_COMM_WORLD);

  return 0;
}

int solver_mpi_gather(struct solver_data *solver, double *data) {
  long int range;
  int i;

  static int initialized = 0;
  static int *cnts, *displs;
  int disp = 0;

  range = (IMAX + (solver->size - 1)) / solver->size;
  range *= JMAX; 
  range *= KMAX;
  
  if (!initialized) {
	  cnts = malloc(solver->size * sizeof(int));
	  displs = malloc(solver->size * sizeof(int));
	  for (i = 0; i < solver->size; i++)
	  {
		  cnts[i] = (i != 0) ? range : 0;
		  displs[i] = disp;
		  disp += range;
	  }
  }
  
  if (!solver->rank)
	  MPI_Gatherv(MPI_IN_PLACE, cnts[solver->rank], MPI_DOUBLE,
	  data, cnts, displs, MPI_DOUBLE,
	  0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(data, cnts[solver->rank], MPI_DOUBLE,
	  data, cnts, displs, MPI_DOUBLE,
	  0, MPI_COMM_WORLD);
  
  return 0;
}

int solver_send_all(struct solver_data *solver) {
  long int range, start;
  int size, rank, n;
  struct kE_data *kE;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  range = (IMAX + (size - 1)) / size;
  
  for(n=1; n<size; n++) {
    start = range * n;
    start -= 1;	
    start = max(start, 0);
    
    if(start + range + 2 > IMAX) range = IMAX - start - 2;
    
    solver_mpi_send(solver, solver->mesh->fv, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->ae, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->an, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->at, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->vof, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->P, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->u, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->v, n, start, range + 2);
    solver_mpi_send(solver, solver->mesh->w, n, start, range + 2);
    
    if(kE_check(solver))  {
      kE = solver->mesh->turbulence_model;
      solver_mpi_send(solver, kE->k, n, start, range + 2);
      solver_mpi_send(solver, kE->E, n, start, range + 2);
      solver_mpi_send(solver, kE->nu_t, n, start, range + 2);
    }
  }

  return 0;
}

int solver_mpi_isend_int(struct solver_data *solver, int *data, int to, long int i_start, long int i_range, MPI_Request *request) {
  MPI_Isend(&data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_INT, to, 1, MPI_COMM_WORLD, request);

  return 0;
}

int solver_mpi_irecv_int(struct solver_data *solver, int *data, int from, long int i_start, long int i_range, MPI_Request *request) {
  MPI_Irecv(&data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_INT, from, 1, MPI_COMM_WORLD, request);

  return 0;
}

int solver_mpi_isend(struct solver_data *solver, double *data, int to, long int i_start, long int i_range, MPI_Request *request) {
  MPI_Isend(&data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_DOUBLE, to, 1, MPI_COMM_WORLD, request);

  return 0;
}

int solver_mpi_irecv(struct solver_data *solver, double *data, int from, long int i_start, long int i_range, MPI_Request *request) {
  MPI_Irecv(&data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, request);

  return 0;
}

int solver_mpi_send(struct solver_data *solver, double *data, int to, long int i_start, long int i_range) {
  MPI_Send(&data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_DOUBLE, to, 1, MPI_COMM_WORLD);

  return 0;
}

int solver_mpi_recv(struct solver_data *solver, double *data, int from, long int i_start, long int i_range) {
  MPI_Recv(&data[mesh_index(solver->mesh,i_start,0,0)], i_range * JMAX * KMAX, MPI_DOUBLE, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return 0;
}

int solver_mpi_sendrecv(struct solver_data *solver, int to, double *send, long int send_i_start, long int send_i_range, int from, double *recv, long int recv_i_start, long int recv_i_range) {
  MPI_Sendrecv(&send[mesh_index(solver->mesh,send_i_start,0,0)], send_i_range * JMAX * KMAX, 
               MPI_DOUBLE, to, 1,
               &recv[mesh_index(solver->mesh,recv_i_start,0,0)], recv_i_range * JMAX * KMAX,
               MPI_DOUBLE, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
  return 0;
}

int solver_mpi_sendrecv_replace(struct solver_data *solver, double *data, long int start, long int range, int to, int from) {
  MPI_Sendrecv_replace(&data[mesh_index(solver->mesh,start,0,0)], 
                       range * JMAX * KMAX, MPI_DOUBLE, to, 1, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  
  return 0;
}

int solver_mpi_sendrecv_int(struct solver_data *solver, int to, int *send, long int send_i_start, long int send_i_range, int from, int *recv, long int recv_i_start, long int recv_i_range) {
  MPI_Sendrecv(&send[mesh_index(solver->mesh,send_i_start,0,0)], send_i_range * JMAX * KMAX, 
               MPI_INT, to, 1,
               &recv[mesh_index(solver->mesh,recv_i_start,0,0)], recv_i_range * JMAX * KMAX,
               MPI_INT, from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  return 0;
}

double solver_mpi_sum(struct solver_data *solver, double x) {
  double ret;

  if(solver->size == 1) return x;

  MPI_Allreduce(&x, &ret, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return ret;
}

double solver_mpi_max(struct solver_data *solver, double x) {
  double ret;

  if(solver->size == 1) return x;

  MPI_Allreduce(&x, &ret, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return ret;
}

double solver_mpi_min(struct solver_data *solver, double x) {
  double ret;

  if(solver->size == 1) return x;

  MPI_Allreduce(&x, &ret, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  return ret;
}

int solver_broadcast_all(struct solver_data *solver) {
  int rank, turb, autot;
  
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
