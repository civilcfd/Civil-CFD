/*
 * solver_data.h
 *
 * defines data structures and functions for a CFD solver
 * actual functions for pressure / velocity coupling, etc appear
 * in other source files and are called by function pointers
 *
 */

#ifndef _SOLVER_DATA_H
#define _SOLVER_DATA_H

#include <time.h>
#include "mesh.h"

struct ic_data {
  /* Describes initial conditions
   * They are referenced in a list (not linked for simplicity)
   * based on a string describing the variable 
   * the list ends when param == "end" */
   
   char param[256];
   int dims;
   double value[3];
};

struct solver_data {

  struct ic_data ic[16]; /* describe up to 16 initial conditions */

	/* MPI DATA */
	int rank;

  double emf; 
  double emf_c;

  double csq;
  double rdtexp;
  
  double min_vof;
  double max_vof;

  double vchgt; /* change in VOF throughout the solution */

  double delt, delt_n, delt_min; /* timestep */
  double endt; /* solution completion criteria */
  double writet;

  long int iter;
  long int niter;
  int vof_flag; /* flag to indicate excessive VOF advection */
  int p_flag;
  
  double con;

  double nu; /* kinematic viscocity m2/s */
  double nu_max; /* maximum effective viscocity, including turbulence */
  double rho; /* density kg/m3 */
  double gx;   /* gravity m2/s */
  double gy;
  double gz;
  
  double resistart; /*initial residual*/
  double resimax; /*finalresidual */
  
  double umax, vmax, wmax;

  double epsi; /* pressure iteration convergence criteria */
  double dzro; /* scaling factor for pressure convergence test */
  double omg;  /* over-relaxation factor for pressure iteration */
  double omg_init; /* starting omg factor */
  double omg_final; /* omg factor at end of pressure iterations */
  double alpha; /* donor cell fluxing (0 means central differencing */
  
  double vof_height; /* initial fluid level */
  double vof_delay;  /* delay to reach this height */
  
	time_t start_time;
  /*
   * COMMENTED OUT - INCLUDED IN THE MESH DEFINITION
   * KEPT IN CASE I CHANGE MY MIND 
   *

  double *peta; 
  
  int nvrm;       * number of void regions *
  double *volvr;  * void region volume *
  double *nvr;    * void region labels *
  double *pvr;    * void region pressure *

  double *u;  
  double *v;
  double *w;

  double *P;

  double *vof;
  double *nvof;
  double *peta;
  double *tanth;
  */

  double t;

  /* functions to define the solution
   * these should be assigned by the application
   * default solvers will be built as part of this library */
  int (*init)(struct solver_data *solver);
  int (*kill)(struct solver_data *solver);
  int (*loop)(struct solver_data *solver);
  int (*boundaries)(struct solver_data *solver);
  int (*special_boundaries)(struct solver_data *solver);
  int (*pressure)(struct solver_data *solver);
  int (*velocity)(struct solver_data *solver);
  int (*vfconv)(struct solver_data *solver);
  int (*petacal)(struct solver_data *solver);
  int (*betacal)(struct solver_data *solver);
  int (*deltcal)(struct solver_data *solver);
  int (*write)(struct solver_data *solver);
  int (*output)(struct solver_data *solver);
  int (*wall_shear)(struct solver_data *solver);
  
  
  int (*turbulence_init)(struct solver_data *solver);
  int (*turbulence_loop)(struct solver_data *solver);
  int (*turbulence_kill)(struct solver_data *solver);
  double (*turbulence_nu)(struct solver_data *solver, long int i, long int j, long int k);
  int (*turbulence_read)(char *filename);
  int (*turbulence_write)(char *filename);
  int (*turbulence_load_values) (struct solver_data *solver);

  /* solution mesh */

  struct mesh_data *mesh;

  int ready;
};

#endif
