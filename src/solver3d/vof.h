/*
 * vof.h
 *
 * 3d solver with VOF functionality
 * Fractional Areas / Volume used to define obstacles
 * This is the default solver
 */

#ifndef _VOF_H
#define _VOF_H

#include "solver_data.h"

int vof_init_solver(struct solver_data *solver); 
int vof_loop(struct solver_data *solver);
int vof_boundaries(struct solver_data *solver);
int vof_special_boundaries(struct solver_data *solver);
int vof_pressure(struct solver_data *solver);
int vof_velocity(struct solver_data *solver);
int vof_convect(struct solver_data *solver);
int vof_petacal(struct solver_data *solver);
int vof_hydrostatic(struct solver_data *solver);
int vof_betacal(struct solver_data *solver);
int vof_deltcal(struct solver_data *solver);
int vof_write(struct solver_data *solver);
int vof_write_timestep(struct solver_data * solver);
int vof_output(struct solver_data *solver);
int vof_setup_solver(struct solver_data *solver);
int vof_kill_solver(struct solver_data *solver);
int vof_pressure_test(struct solver_data *solver);
int vof_pressure_init(struct solver_data *solver);
int vof_vorticity(struct solver_data *solver);
int vof_pressure_sor(struct solver_data *solver, long int i, long int j, long int k);
int vof_pressure_mp(struct solver_data *solver);
int vof_pressure_gmres(struct solver_data *solver);
int vof_pressure_gmres_mpi(struct solver_data *solver);
#endif
