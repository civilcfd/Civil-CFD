/*
 * writesolver.c
 *
 * functions to write a solver into a file
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "solver.h"
#include "readsolver.h"
#include "laminar.h"
#include "kE.h"

int write_solver(struct solver_data *solver, char *filename) {
  FILE *fp;
  
  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to write_solver\n");
    return(1);
  }

  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: cannot open %s in write_solver\n",filename);
    return(1);
  }


  fprintf(fp,"gravity %e %e %e\n",solver->gx, solver->gy, solver->gz);
  fprintf(fp,"nu %e\n",solver->nu);
  fprintf(fp,"rho %e\n",solver->rho);
  fprintf(fp,"t %e\n",solver->t);
  fprintf(fp,"delt %e\n",solver->delt);
  fprintf(fp,"writet %e\n",solver->writet);
  fprintf(fp,"endt %e\n",solver->endt);
  if(solver->deltcal == NULL) fprintf(fp,"autot 0\n");
  
  if(kE_check(solver)) fprintf(fp,"turbulence 1\n");
  else fprintf(fp,"turbulence 0\n");

  fclose(fp);

  return 0;
}

int write_initial(struct solver_data *solver, char *filename) {
  FILE *fp;
  double velocity[3], inside[3];
  double value;
  
  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to write_initial\n");
    return(1);
  }

  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: cannot open %s in write_initial\n",filename);
    return(1);
  }

  solver_get_initial_vector(solver, "velocity", velocity);
  fprintf(fp,"velocity %e %e %e\n",velocity[0],velocity[1],velocity[2]);

  while(solver_get_initial_vector(solver, "inside", inside) != 1); // make sure we start at the first instance
  while(solver_get_initial_vector(solver, "inside", inside) != 1) {
    fprintf(fp,"inside %e %e %e\n",inside[0],inside[1],inside[2]);
  }
  
  value = solver_get_initial_scalar(solver, "vof_height");  
  if(value>0) fprintf(fp,"vof_height %e\n",value);
  
  value = solver_get_initial_scalar(solver, "hydrostatic");    
  if(value > 0) fprintf(fp,"hydrostatic 1\n");
  
  value = solver_get_initial_scalar(solver, "kE_k");     
  if(value > 0 && kE_check(solver)) fprintf(fp,"kE_k %e\n",value);

  fclose(fp);

  return 0;
}
