/*
 * mesh_mpi.h:  function definitions for mesh_mpi.c 
 *
 */

#ifndef MESH_MPI_H
#define MESH_MPI_H


long int mesh_mpi_space(struct mesh_data *solver);
int mesh_broadcast_all(struct mesh_data *mesh);
int mesh_broadcast_constants(struct mesh_data *mesh);
int mesh_broadcast_special_boundaries(struct mesh_data *mesh);
int mesh_broadcast_baffles(struct mesh_data *mesh);
struct mesh_data *mesh_mpi_init_copy(struct mesh_data *mesh_source);
int mesh_mpi_copy_data(struct mesh_data *mesh, struct mesh_data *mesh_source);
int mesh_mpi_init_complete(struct mesh_data *mesh);

#endif