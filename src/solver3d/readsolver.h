/* readsolver.h
 *
 * header filef or file i/o for solver_data
 */

#ifndef _READSOLVER_H
#define _READSOLVER_H

#include "solver.h"

int read_solver(struct solver_data *solver, char *filename);
int read_solver_xml(struct solver_data *solver, char *filename);
int read_initial(struct solver_data *solver, char *filename);

#endif
