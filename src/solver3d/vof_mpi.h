/*
 * vof.h
 *
 * 3d solver with VOF functionality
 * Fractional Areas / Volume used to define obstacles
 * This is the default solver
 */

#ifndef _VOF_MPI_H
#define _VOF_MPI_H

#include "solver_data.h"

int vof_mpi_init_solver(struct solver_data *solver); 
int vof_mpi_loop(struct solver_data *solver);
int vof_mpi_boundaries(struct solver_data *solver);
int vof_mpi_special_boundaries(struct solver_data *solver);
int vof_mpi_pressure(struct solver_data *solver);
int vof_mpi_velocity(struct solver_data *solver);
int vof_mpi_vfconv(struct solver_data *solver);
int vof_mpi_petacal(struct solver_data *solver);
int vof_mpi_hydrostatic(struct solver_data *solver);
int vof_mpi_betacal(struct solver_data *solver);
int vof_mpi_deltcal(struct solver_data *solver);
int vof_mpi_write(struct solver_data *solver);
int vof_mpi_write_timestep(struct solver_data * solver);
int vof_mpi_output(struct solver_data *solver);
int vof_mpi_setup_solver(struct solver_data *solver);
int vof_mpi_kill_solver(struct solver_data *solver);
int vof_mpi_pressure_test(struct solver_data *solver);
int vof_mpi_pressure_init(struct solver_data *solver);
int vof_mpi_vorticity(struct solver_data *solver);
int vof_mpi_pressure_sor(struct solver_data *solver, long int i, long int j, long int k);
int vof_mpi_pressure_mp(struct solver_data *solver);
int vof_mpi_pressure_gmres(struct solver_data *solver);
int vof_mpi_pressure_gmres_mpi(struct solver_data *solver);
#endif
