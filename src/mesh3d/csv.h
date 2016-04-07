/* vtk.h */

#ifndef CSV_H 
#define CSV_H

#include "mesh.h"
int csv_write_n_vof(struct mesh_data *mesh, double timestep);
int csv_write_k(struct mesh_data *mesh, double timestep);
int csv_write_E(struct mesh_data *mesh, double timestep);
int csv_write_vof(struct mesh_data *mesh, double timestep);
int csv_read_vof(struct mesh_data *mesh, double timestep);
int csv_read_n_vof(struct mesh_data *mesh, double timestep);
int csv_write_U(struct mesh_data *mesh, double timestep);
int csv_read_U(struct mesh_data *mesh, double timestep);
int csv_write_P(struct mesh_data *mesh, double timestep);
int csv_read_P(struct mesh_data *mesh, double timestep);
int csv_write_fv(struct mesh_data *mesh, double timestep);
int csv_read_fv(struct mesh_data *mesh, double timestep);
int csv_write_vorticity(struct mesh_data *mesh, double timestep);

int csv_write_scalar_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double *scalars); 
int csv_compressed_write_scalar_grid(char *filename_csv, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double *scalars); 
int csv_write_integer_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          int *scalars); 
int csv_compressed_write_integer_grid(char *filename_csv, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          int *scalars); 
int csv_write_vector_grid(char *filename, char *dataset_name, long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2); 
int csv_compressed_write_vector_grid(char *filename_csv, char *dataset_name, long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2); 

long int csv_read_scalar_grid(char *filename,  
                          long int ni, long int nj, long int nk,
                          double *scalars);
long int csv_compressed_read_scalar_grid(char *filename_csv,  
                          long int ni, long int nj, long int nk,
                          double *scalars);

long int csv_read_integer_grid(char *filename,  
                          long int ni, long int nj, long int nk,
                          int *scalars);
long int csv_compressed_read_integer_grid(char *filename_csv,  
                          long int ni, long int nj, long int nk,
                          int *scalars);
                          
long int csv_read_vector_grid(char *filename,  
                          long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2); 
long int csv_compressed_read_vector_grid(char *filename_csv,  
                          long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2); 
                          
int csv_read_af(struct mesh_data *mesh, double timestep);
int csv_write_af(struct mesh_data *mesh, double timestep);

int csv_write_scalar_grid_paraview(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars);

void csv_remove(char *filename);
#endif
