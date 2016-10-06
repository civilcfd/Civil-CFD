/*
 * solver.h
 *
 * header for functions in solver.c
 */

#ifndef _SOLVER_H
#define _SOLVER_H

#include "solver_data.h"

struct solver_data *solver_init_empty();

int solver_init_complete(struct solver_data *solver);

int solver_check(struct solver_data *solver);

int solver_load(struct solver_data *solver, char *solverfile);

int solver_run(struct solver_data *solver);

int solver_set_value(struct solver_data *solver, char *param, int dims, double *vector);

int solver_store_initial(struct solver_data *solver, char *param, int dims, 
                   double *vector);

int solver_set_initial(struct solver_data *solver, char *param, int dims, 
                   double *vector);
                   
double solver_get_initial_scalar(struct solver_data *solver, char *param);
int solver_get_initial_vector(struct solver_data *solver, char *param, double *vector);

int solver_initial_values(struct solver_data *solver);                                

#endif
