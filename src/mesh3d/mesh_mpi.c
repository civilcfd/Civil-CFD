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

int mesh_mpi_init_complete(struct mesh_data *mesh, int rank) {
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
  
  if(!rank) size = mesh->imax * mesh->jmax * mesh->kmax;
  else size = mesh->i_range * mesh->jmax * mesh->kmax;

  if(size<=0) {
    printf("error: improper size mesh in mesh_mpi_init_complete\n");
    return (1);
  }

  mesh->delu_upstream = malloc(sizeof(double) * mesh->jmax * mesh->kmax);

  if(mesh->delu_upstream == NULL) {
    printf("error: memory could not be allocated for del_u in mesh_mpi_init_complete\n");
    return(1);
  }

  mesh->delu_downstream = malloc(sizeof(double) * mesh->jmax * mesh->kmax);

  if(mesh->delu_downstream == NULL) {
    printf("error: memory could not be allocated for del_u in mesh_mpi_init_complete\n");
    return(1);
  }

	mesh_allocate(mesh, size);
  return 0;
}

int mesh_broadcast_all(struct mesh_data *mesh) {
	int rank;
	
	rank = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	mesh_broadcast_constants(mesh);
	mesh_broadcast_special_boundaries(mesh);
	mesh_broadcast_baffles(mesh);

	return 0;
}

struct mesh_data *mesh_mpi_init_copy(struct mesh_data *mesh_source) {

  struct mesh_data *mesh;
  int i;

  mesh = mesh_init_empty();
  if (mesh==NULL || mesh_source==NULL) {
    printf("error: mesh_init_copy failed\n");
    return NULL;
  }

  mesh->imax = mesh_source->imax;
  mesh->jmax = mesh_source->jmax;
  mesh->kmax = mesh_source->kmax;

  mesh->i_range = mesh_source->i_range;
  mesh->i_start = mesh_source->i_start;

  mesh->delx = mesh_source->delx;
  mesh->dely = mesh_source->dely;
  mesh->delz = mesh_source->delz;
  mesh->rdx = mesh_source->rdx;
  mesh->rdy = mesh_source->rdy;
  mesh->rdz = mesh_source->rdz;


  mesh->origin[0] = mesh_source->origin[0];
  mesh->origin[1] = mesh_source->origin[1];
  mesh->origin[2] = mesh_source->origin[2];

  mesh->inside[0] = mesh_source->inside[0];
  mesh->inside[1] = mesh_source->inside[1];
  mesh->inside[2] = mesh_source->inside[2];

  for(i = 0; i < 6; i++) {
    mesh->wb[i] = mesh_source->wb[i];
    mesh->sb[i] = mesh_source->sb[i];
    
    if(i<3) mesh->baffles[i] = mesh_source->baffles[i];
  }

  mesh->ready = mesh_source->ready;

  if(mesh_mpi_init_complete(mesh, 1) == 1) {
    free(mesh);
    return NULL;
  }

  mesh_mpi_copy_data(mesh, mesh_source);

  return mesh;
}


int mesh_mpi_copy_data(struct mesh_data *mesh, struct mesh_data *mesh_source) {
  long int size;

  size = mesh->i_range * mesh->jmax * mesh->kmax;

  if(size <= 0) {
    printf("error: improper size in mesh_copy_data\n");
    return 1;
  }
  
  memcpy(mesh->P, mesh_source->P, size * sizeof(double));
  memcpy(mesh->D, mesh_source->D, size * sizeof(double));
  memcpy(mesh->u, mesh_source->u, size * sizeof(double));
  memcpy(mesh->v, mesh_source->v, size * sizeof(double));
  memcpy(mesh->w, mesh_source->w, size * sizeof(double));
  memcpy(mesh->u_omega, mesh_source->u_omega, size * sizeof(double));
  memcpy(mesh->v_omega, mesh_source->v_omega, size * sizeof(double));
  memcpy(mesh->w_omega, mesh_source->w_omega, size * sizeof(double));
  
  memcpy(mesh->vof, mesh_source->vof, size * sizeof(double));
  memcpy(mesh->n_vof, mesh_source->n_vof, 
          size * sizeof(enum cell_boundaries));
  
  memcpy(mesh->fv, mesh_source->fv, size * sizeof(double));
  memcpy(mesh->ae, mesh_source->ae, size * sizeof(double));
  memcpy(mesh->an, mesh_source->an, size * sizeof(double));
  memcpy(mesh->at, mesh_source->at, size * sizeof(double));

  memcpy(mesh->peta, mesh_source->peta, size * sizeof(double));
  memcpy(mesh->tanth, mesh_source->tanth, size * sizeof(double));
  memcpy(mesh->beta, mesh_source->beta, size * sizeof(double));

  return 0;
}

int mesh_broadcast_constants(struct mesh_data *mesh) {

	MPI_Bcast(&mesh->imax, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->jmax, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&mesh->kmax, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	
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
		
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(!rank) {
		for(x=0; x < 6; x++) {
			for(sb = mesh->sb[x]; sb != NULL; sb = sb->next) {    
				MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD); /* tell processes to create special boundary x */
				MPI_Bcast(&sb->type, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Bcast(&sb->extent_a[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&sb->extent_b[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
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
				MPI_Bcast(&extent_a[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&extent_b[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
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
  
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  
  if(!rank) {
		for(x=0; x < 3; x++) {
			for(baffle = mesh->baffles[x]; baffle != NULL; baffle = baffle->next) { 
					MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD); /* tell processes to create baffle on axis, x */
					MPI_Bcast(&baffle->type, 1, MPI_INT, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->extent_a[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->extent_b[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					MPI_Bcast(&baffle->pos, 1, MPI_LONG, 0, MPI_COMM_WORLD);
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
				MPI_Bcast(&extent_a[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
				MPI_Bcast(&extent_b[0], 2, MPI_LONG, 0, MPI_COMM_WORLD);
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