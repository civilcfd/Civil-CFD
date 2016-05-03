/* mesh.h
 *
 * Definitions for mesh struct types
 */

#ifndef _MESH_H
#define _MESH_H

#include "mesh_data.h"
#include "stl.h"

#ifndef _VOF_MACROS_H
#define FV(i, j, k) mesh->fv[mesh_index(mesh, i, j, k)]
#define AE(i, j, k) mesh->ae[mesh_index(mesh, i, j, k)]
#define AN(i, j, k) mesh->an[mesh_index(mesh, i, j, k)]
#define AT(i, j, k) mesh->at[mesh_index(mesh, i, j, k)]
#endif

int mesh_baffle_create(struct mesh_data *mesh, int axis, int type, double value, long int pos);

int mesh_baffle_extent_a(struct mesh_data *mesh, int axis, long int extent_a_1, long int extent_a_2);

int mesh_baffle_extent_b(struct mesh_data *mesh, int axis, long int extent_b_1, long int extent_b_2);

int mesh_load_csv(struct mesh_data *mesh, int timestep);

struct mesh_data *mesh_init_copy(struct mesh_data *mesh_source);

int mesh_copy_data(struct mesh_data *mesh, struct mesh_data *mesh_source);

struct mesh_data *mesh_init_empty();

int moller_trumbore(double *r_o, double *r_d, double *v1,
                   double *v2, double *v3, double *pt);

int mesh_init_complete(struct mesh_data *mesh);

int mesh_check(struct mesh_data *mesh);

int mesh_set_value(struct mesh_data *mesh, char *param, int dims, 
                   double *vector);

int mesh_allocate(struct mesh_data *mesh, long int size);
int mesh_free(struct mesh_data *mesh);

long int mesh_index(struct mesh_data *mesh,
                    long int i, long int j, long int k);

int mesh_fill(struct mesh_data *mesh, struct stl_data *stl); 

int mesh_fill_vof(struct mesh_data *mesh, double *vector); 

int mesh_normalize(struct mesh_data *mesh);

int mesh_set_hydrostatic(struct mesh_data *mesh, double g, double rho);

int mesh_set_array(struct mesh_data *mesh, char *param, double value, 
                   long int imin, long int imax,
                   long int jmin, long int jmax,
                   long int kmin, long int kmax);

int mesh_sb_create(struct mesh_data *mesh, int wall, 
                   int type, double value, double turbulence);                    

int mesh_sb_extent_a(struct mesh_data *mesh, int wall, long int extent_a_1, long int extent_a_2);
int mesh_sb_extent_b(struct mesh_data *mesh, int wall, long int extent_b_1, long int extent_b_2);

int mesh_avratio(struct mesh_data *mesh, double avr_max);
int mesh_area_correct(double *a1, double *a2, double an1, double an2, double r);

#endif
