/* stl.h
 *
 * STL filetype header file
 * STL data is described by a simple struct with 2d arrays representing
 * a list of vertices and normals
 * 
 * v_n = the nth vertex in a facet
 * facets is the actual number of facets
 * 
 * a hard limit is put on the STL size, since gigantic STLs likely will
 * not work with this mesher anyways
 */

#ifndef _STL_H
#define _STL_H

#include "mesh_data.h"

#define MAX_FACETS 1000000

struct stl_data {
  int ready;

  double normal[MAX_FACETS][3];

  double v_1[MAX_FACETS][3];
  double v_2[MAX_FACETS][3];
  double v_3[MAX_FACETS][3];

  long int facets;

  char solid[1024];
};

struct stl_data *stl_init_empty();

int stl_set_value(struct stl_data *stl, int dims, 
                   char (*args)[256], int state, double *limits);

int stl_check(struct stl_data *stl);

int stl_free(struct stl_data *stl);

int stl_check_normals(struct mesh_data *mesh, struct stl_data *stl, 
                      long int i, long int j, long int k);
											
int stl_check_normals_point(struct mesh_data *mesh, struct stl_data *stl, 
                      double *pt);

#endif
