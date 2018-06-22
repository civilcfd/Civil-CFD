/* vtk.h */

#ifndef VTK_XML_H
#define VTK_XML_H

#include <stdint.h>
#include "mesh.h"

int vtk_xml_write_n_vof(struct mesh_data *mesh, int timestep);
int vtk_xml_write_k(struct mesh_data *mesh, int timestep);
int vtk_xml_write_E(struct mesh_data *mesh, int timestep);
int vtk_xml_write_fv(struct mesh_data *mesh, int timestep);
int vtk_xml_write_vof(struct mesh_data *mesh, int timestep);
int vtk_xml_write_P(struct mesh_data *mesh, int timestep);
int vtk_xml_write_U(struct mesh_data *mesh, int timestep);
int vtk_xml_write_vorticity(struct mesh_data *mesh, int timestep);
int vtk_xml_write_vector_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *v1, double *v2, double *v3);
int vtk_xml_compressed_write_vector_grid(char *filename_vtk, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *v1, double *v2, double *v3);
int vtk_xml_write_scalar_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars); 
int vtk_xml_compressed_write_scalar_grid(char *filename_vtk, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars); 
int vtk_xml_write_scalar_magnitude_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars); 
int vtk_xml_write_integer_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          int *scalars); 
int vtk_xml_compressed_write_integer_grid(char *filename_vtk, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          int *scalars); 

void vtk_xml_remove(char *filename);
int vtk_xml_decompress(const char *cstr);
#endif
