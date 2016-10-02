#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <petscksp.h>

#include "solver_data.h"
#include "solver.h"
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
  int rank, size;
  int ret = 0;

	PetscInitialize(NULL, NULL, NULL, NULL);
#ifdef _WIN32
	PetscOptionsSetValue(NULL, "-log_summary", "true");
#else
	PetscOptionsSetValue("-log_summary","true");
#endif
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
  else {
    printf("error: Run must be invoked using mpirun/mpiexec with at least 2 processes\n");
  }


  return(0);
}
