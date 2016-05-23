/*
 * laminar.c
 *
 * sets up the laminar model */

#include <stdio.h>

#include "laminar.h"
#include "solver.h"
 
int laminar_read(char *filename) { return 0; }
int laminar_load_values(char *filename) { return 0; }
int laminar_write(char *filename) { return 0; }

int laminar_init(struct solver_data *solver) { return 0; }

int laminar_loop(struct solver_data *solver) { return 0; }

int laminar_kill(struct solver_data *solver) { return 0; }

int laminar_setup(struct solver_data *solver) {

  solver->turbulence_init = laminar_init;
  solver->turbulence_loop = laminar_loop;
  solver->turbulence_kill = laminar_kill;
  solver->turbulence_read = laminar_read;
  solver->turbulence_write = laminar_write;
  solver->wall_shear = laminar_wall_shear;
  solver->turbulence_nu   = laminar_nu;
  solver->turbulence_load_values = laminar_load_values;
    
  solver->mesh->turbulence_model = NULL;
  
  return 0;
}

int laminar_wall_shear(struct solver_data *solver) {
/* write code for laminar wall shear */
  return 0;
}

double laminar_nu(struct solver_data *solver, long int i, long int j, long int k) {
  return solver->nu;
}
