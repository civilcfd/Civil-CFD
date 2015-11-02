/* volfract.h
 *
 * header file for volfract.c 
 * mostly function definitions */

#ifndef VOLFRACT_H
#define VOLFRACT_H

#include "mesh.h"
#include "stl.h"

int count_bits(int n);

double tet_volume(double *a, double *b, double *c, double *d);

int volume_fractions(struct mesh_data *mesh, struct stl_data *stl); 

int tet_fraction(struct mesh_data *mesh, struct stl_data *stl,
                 long int i, long int j, long int k);

int face_hex_fraction(struct mesh_data *mesh, struct stl_data *stl,
                       long int i, long int j, long int k);
 
int line_pent_fraction(struct mesh_data *mesh, struct stl_data *stl,
                       long int i, long int j, long int k); 

int vertex_hex_fraction(struct mesh_data *mesh, struct stl_data *stl,
                        long int i, long int j, long int k); 

int vertex_pent_fraction(struct mesh_data *mesh, struct stl_data *stl,
                        long int i, long int j, long int k);
												
int vertex_qhull_fraction(struct mesh_data *mesh, struct stl_data *stl, long int i, long int j, long int k);												

#endif

