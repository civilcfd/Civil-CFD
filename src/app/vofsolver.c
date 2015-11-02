#include <stdio.h>
#include <stdlib.h>

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

int main(int argc, char *argv[])
{
  struct solver_data *solver;
  double timestep;

  printf("vofsolver: 3d VOF solver\n");

  if (argc<2) {
    printf("usage: vofsolver <source file> <timestep>\n");
    return(1);
  }

  #ifdef DEBUG
    printf("%s %s\n", argv[0], argv[1]);
  #endif

  solver = solver_init_empty();
  if(solver == NULL) return 1;

  if(solver_load(solver, NULL, argv[1], NULL)==1)
    return 1;

  vof_setup_solver(solver);
  kE_setup(solver);

  if(solver_init_complete(solver)==1)
    return 1;

  solver->endt=10.0;
  solver->gz = 0; 
  solver->writet=0.05;
  solver->delt=0.01;

  mesh_set_array(solver->mesh, "u", 1.0, -1, 0, 0, 0, 0, 0);
  /* mesh_set_array(solver->mesh, "v", 0.5,  0, 20, 0, 40, 0, 10); sloshing tank */
  mesh_set_array(solver->mesh, "w", 0.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(solver->mesh, "w", 0.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(solver->mesh, "P", 0.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(solver->mesh, "vof", 1.0, -1, 0, 0, 0, 0, 0);
  /* mesh_set_array(solver->mesh, "vof", 1.0, 0, 20, 0, 40, 0, 10); sloshing tank */
  /* mesh_set_array(solver->mesh, "vof", 1.0, 0, 10, 0, 20, 0, 20); dambreak */
  mesh_set_array(solver->mesh, "fv", 1.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(solver->mesh, "ae", 1.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(solver->mesh, "an", 1.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(solver->mesh, "at", 1.0, -1, 0, 0, 0, 0, 0);

  mesh_load_csv(solver->mesh, 0); 

  if(argc==3) {
    timestep = atof(argv[2]);
    if(timestep != 0) {
      mesh_set_array(solver->mesh, "u", 0.0, -1, 0, 0, 0, 0, 0);
      mesh_set_array(solver->mesh, "w", 0.0, -1, 0, 0, 0, 0, 0);
      mesh_set_array(solver->mesh, "w", 0.0, -1, 0, 0, 0, 0, 0);
      mesh_set_array(solver->mesh, "P", 0.0, -1, 0, 0, 0, 0, 0);
      mesh_set_array(solver->mesh, "vof", 0.0, -1, 0, 0, 0, 0, 0);

      csv_read_U(solver->mesh,timestep);
      csv_read_P(solver->mesh,timestep);
      csv_read_vof(solver->mesh,timestep);
    
      solver->t = timestep;
    }
  }

  solver->init(solver);
  solver->turbulence_init(solver);
  
  if(solver_run(solver)==1)
    return 1;

  return(0);
}
