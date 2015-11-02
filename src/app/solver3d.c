#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int main(int argc, char *argv[])
{
  struct solver_data *solver;
  double timestep;
  int aborted_run;

  printf("solver3d: 3d solver to accompany the Civil CFD gui\n");

  /* if (argc<2) {
    printf("usage: vofsolver <timestep> <delt>\n");
    return(1);
  } */

  #ifdef DEBUG
    printf("%s %s\n", argv[0], argv[1]);
  #endif
  if (argc > 1) timestep = atof(argv[1]);
  else timestep = 0;
  
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

  vof_setup_solver(solver);
  
  if(solver_load(solver, "solverfile", "meshfile", "initials")==1)
    return 1;

  if(solver_init_complete(solver)==1)
    return 1;

  if(solver->turbulence_read != NULL) solver->turbulence_read("turbulencefile");
  
  solver->init(solver);
  solver->turbulence_init(solver);
  
  mesh_load_csv(solver->mesh, 0); 
  solver_initial_values(solver);

  if(timestep != 0) {
    csv_read_U(solver->mesh,timestep);
    csv_read_P(solver->mesh,timestep);
    csv_read_vof(solver->mesh,timestep);
    solver->t = timestep;
  }
  else
    solver->write(solver); 
  
  if(aborted_run) return 0;
  
  if (argc > 2) solver->delt = atof(argv[2]);
  
  track_read();
  
  if(solver_run(solver)==1)
    return 1;

  return(0);
}
