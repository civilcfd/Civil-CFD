/* mesh_mpi.c
 *
 * extends mesh.c to allow mpi paralellism
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <petscksp.h>

#include "stack.h"
#include "mesh.h"
#include "mesh_mpi.h"
#include "csv.h"

int mesh_mpi_init_complete(struct mesh_data *mesh, long int range, long int start) {
  long int size;

  if(mesh == NULL) {
    printf("error: null mesh passed to mesh_mpi_init_complete\n");
    return (1);
  }
  
  mesh_check(mesh);
  if(mesh->ready == 0) {
    printf("error: uninitialized mesh passed to mesh_mpi_init_complete\n");
    return (1);
  }
  
  mesh->i_range = range;
  mesh->i_start = start;
  size = mesh->i_range * mesh->jmax * mesh->kmax;

  if(size<=0) {
    printf("error: improper size mesh in mesh_mpi_init_complete\n");
    return (1);
  }

	mesh_allocate(mesh, size);
}

int mesh_broadcast_all(struct mesh_data *mesh) {
	int rank;
	
	rank = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	mesh_broadcast_constants(mesh);
	mesh_broadcast_special_boundaries(mesh);
	mesh_broadcast_baffles(mesh);

	return 0;
}

int mesh_broadcast_constants(struct mesh_data *mesh) {
	int rank;
	
	rank = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	
	MPI_Bcast(&mesh->imax, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->jmax, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->kmax, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&mesh->delx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->dely, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->delz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->rdx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->rdy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->rdz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&mesh->origin, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->inside, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	MPI_Bcast(&mesh->wb, 6, MPI_INT, 0, MPI_COMM_WORLD);

	return 0;
}

int mesh_broadcast_special_boundaries(struct mesh_data *mesh) {
  struct sb_data *sb;
  double value, turbulence;
  long int extent_a[2], extent_b[2];
  int type;
	int x, rank, sb_create = -1;	
		
	rank = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(!rank) {
		for(x=0; x < 6; x++) {
			for(sb = mesh->sb[x]; sb != NULL; sb = sb->next) {    
				MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD); /* tell processes to create special boundary x */
				MPI_Bcast(&sb->type, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&sb->extent_a, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&sb->extent_b, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&sb->value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&sb->turbulence, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
		}
		
		MPI_Bcast(&sb_create, 1, MPI_INT, 0, MPI_COMM_WORLD); /* stop signal */
  }
  else {
  
  	do {
			MPI_Bcast(&sb_create, 1, MPI_INT, 0, MPI_COMM_WORLD); 
			
			if(sb_create != -1) {
				MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&extent_a, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&extent_b, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&turbulence, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
				mesh_sb_create(mesh, sb_create, type, value, turbulence);
				mesh_sb_extent_a(mesh, sb_create, extent_a[0], extent_a[1]);
				mesh_sb_extent_b(mesh, sb_create, extent_b[0], extent_b[1]);
				
			}
  	
  	} while(sb_create != -1);
  	
  }
  
  return 0;
}
  

int mesh_broadcast_baffles(struct mesh_data *mesh) {
  struct baffle_data *baffle;
  int x, rank;
  int type;
  long int extent_a[2]; 
  long int extent_b[2];                     
  double value;      
  long int pos;
  
	rank = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  if(!rank) {
		for(x=0; x < 3; x++) {
			for(baffle = mesh->baffles[x]; baffle != NULL; baffle = baffle->next) { 
					MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD); /* tell processes to create baffle on axis, x */
					MPI_Bcast(&baffle->type, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->extent_a, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->extent_b, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->pos, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				}
			}
		
		x = -1;
		MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD); /* stop signal */
  }
  else {
  
  	do {
			MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD); 
			
			if(x != -1) {
				MPI_Bcast(&type, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&extent_a, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&extent_b, 2, MPI_LONG_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				MPI_Bcast(&pos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
				mesh_baffle_create(mesh, x, type, value, pos);
				mesh_baffle_extent_a(mesh, x, extent_a[0], extent_a[1]);
				mesh_baffle_extent_b(mesh, x, extent_b[0], extent_b[1]);
				
			}
  	
  	} while(x != -1);
  	
  }
  
  return 0;
}