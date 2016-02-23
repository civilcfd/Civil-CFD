/* markcells.h
 *
 * 
 */

#ifndef _MARKCELLS_H
#define _MARKCELLS_H

#include "mesh_data.h"
#include "stl.h"

#ifndef _VOF_MACROS_H
#define FV(i, j, k) mesh->fv[mesh_index(mesh, i, j, k)]
#define AE(i, j, k) mesh->ae[mesh_index(mesh, i, j, k)]
#define AN(i, j, k) mesh->an[mesh_index(mesh, i, j, k)]
#define AT(i, j, k) mesh->at[mesh_index(mesh, i, j, k)]
#endif

int markcells_initialize(struct mesh_data *mesh, 
                         struct stl_data *stl);

int markcells_dist_tri_point(double *p, double *v1, double *v2, double *v3);

int markcells_check(long int n);

#endif