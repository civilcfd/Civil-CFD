#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petscksp.h>

#include "solver_data.h"
#include "solver.h"
#include "vof.h"
#include "intersections.h"
#include "mesh.h"
#include "readfile.h"
#include "volfract.h"
#include "vtk.h"
#include "csv.h"
#include "kE.h"
#include "track.h"
#include "vof_macros.h"
#include "solver_mpi.h"


int main(int argc, char *argv[])
{
  struct solver_data *solver = NULL;
  double timestep, delt;
  int aborted_run;
  int rank, size;
  int ret = 0;

	PetscInitialize(NULL, NULL, NULL, NULL);
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
  
		
  printf("solver3d: 3d solver to accompany the Civil CFD gui\n");


  #ifdef DEBUG
    printf("%s %s\n", argv[0], argv[1]);
  #endif
  if (argc > 1) timestep = atof(argv[1]);
  else timestep = 0;
  if (argc > 2) delt = atof(argv[2]);
  else delt = -1;
  
  if(size > 1) {
    ret = solver_mpi(solver, timestep, delt);
    PetscFinalize();
    return ret;
  }
  
  aborted_run = 0;
  if (argc > 2) {
    if(strcmp(argv[2], "abort") == 0) {
      aborted_run = 1;
      timestep = 0;
      printf("this run will be aborted after writing timestep 0\n");
    }
  }

  solver = solver_init_empty();
  if(solver == NULL) return 1;
  solver->size = 1;
  solver->rank = 1;

  vof_setup_solver(solver);
  
  if(solver_load(solver, "solverfile", "meshfile", "initials")==1)
    return 1;

  if(solver_init_complete(solver)==1)
    return 1;

  if(solver->turbulence_read != NULL) solver->turbulence_read("turbulencefile");
  
  solver->init(solver);
  solver->turbulence_init(solver);
  if (delt > solver->emf) solver->delt = delt;
  
  if(mesh_load_csv(solver->mesh, 0) == 1) return 1;
  solver_initial_values(solver);
  solver->petacal(solver);

  if(timestep != 0) {
    solver->t = timestep;
    csv_read_U_p_vof(solver->mesh, timestep);
    solver->turbulence_load_values(solver);
    solver->deltcal(solver);
  }
  else
    solver->write(solver); 
  
  if(aborted_run) return 0;
  
  
  track_read();
  solver_mpi_range(solver);
  
  if(solver_run(solver)==1)
    return 1;

  return(0);
}
