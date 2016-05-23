/*
 * laminar.h
 *
 * sets up the laminar model */
 
#ifndef LAMINAR_H
#define LAMINAR_H

#include "solver_data.h"

int laminar_init(struct solver_data *solver);
int laminar_loop(struct solver_data *solver);
int laminar_kill(struct solver_data *solver);
int laminar_setup(struct solver_data *solver);
int laminar_wall_shear(struct solver_data *solver);
double laminar_nu(struct solver_data *solver, long int i, long int j, long int k);
int laminar_read(char *filename);
int laminar_write(char *filename);
int laminar_load_values(char *filename);

#endif
