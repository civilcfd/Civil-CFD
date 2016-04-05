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
  long int i,j,k;

  printf("inspect_cell: diagnostic tool that outputs relevant data for a cell\n");

  if (argc<4) {
    printf("usage: inspect_cell <timestep> <i> <j> <k>\n");
    return(1);
  } 

  #ifdef DEBUG
    printf("%s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4]);
  #endif
  timestep = atof(argv[1]);
  i = atol(argv[2]);
  j = atol(argv[3]);
  k = atol(argv[4]);

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

  csv_read_U(solver->mesh,timestep);
  csv_read_P(solver->mesh,timestep);
  csv_read_vof(solver->mesh,timestep);
  solver->petacal(solver);
  solver->t = timestep;

  printf("\ntimestep %lf     cell: %ld %ld %ld\n\ncell data:\n\n",timestep, i, j, k);
  
  printf("FV: %lf\n",FV(i,j,k));
  printf("P: %lf      VOF: %lf     N_VOF: %d\n\nedge data:\n\n",P(i,j,k),VOF(i,j,k),N_VOF(i,j,k));
  printf("A(e/n/t) %ld %ld %ld: %lf %lf %lf\n",i,j,k,AE(i,j,k),AN(i,j,k),AT(i,j,k));
  printf("A(w/s/b) %ld %ld %ld: %lf %lf %lf\n\n",i,j,k,AE(i-1,j,k),AN(i,j-1,k),AT(i,j,k-1));  
  printf("U(e/n/t) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i,j,k),V(i,j,k),W(i,j,k));
  printf("U(w/s/b) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i-1,j,k),V(i,j-1,k),W(i,j,k-1));  

  return(0);
}
