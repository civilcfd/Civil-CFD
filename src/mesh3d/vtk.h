/* vtk.h */

#ifndef VTK_H
#define VTK_H

#include "mesh.h"

int vtk_write_n_vof(struct mesh_data *mesh, int timestep);
int vtk_write_k(struct mesh_data *mesh, int timestep);
int vtk_write_E(struct mesh_data *mesh, int timestep);
int vtk_write_fv(struct mesh_data *mesh, int timestep);
int vtk_write_vof(struct mesh_data *mesh, int timestep);
int vtk_write_P(struct mesh_data *mesh, int timestep);
int vtk_write_U(struct mesh_data *mesh, int timestep);
int vtk_write_vorticity(struct mesh_data *mesh, int timestep);
int vtk_write_vector_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *v1, double *v2, double *v3);
int vtk_write_scalar_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars); 
int vtk_write_scalar_magnitude_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars); 
int vtk_write_integer_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          int *scalars); 
double double_swap(double d);
uint16_t int_swap(uint16_t n);
#endif
