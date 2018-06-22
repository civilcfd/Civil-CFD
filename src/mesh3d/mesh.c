/* mesh.c
 *
 * mesh initialization and value population
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "stack.h"
#include "mesh.h"
#include "csv.h"

/* loads a mesh af and fv and sets it up from a csv file */
int mesh_load_csv(struct mesh_data *mesh, int timestep) {
  
  mesh_set_array(mesh, "fv", 0.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(mesh, "ae", 0.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(mesh, "an", 0.0, -1, 0, 0, 0, 0, 0);
  mesh_set_array(mesh, "at", 0.0, -1, 0, 0, 0, 0, 0);

  if(csv_read_af(mesh, timestep) == -1) return 1;
  if(csv_read_fv(mesh, timestep) == -1) return 1;

  return 0;
}

/* initializes a mesh of the same size as the argument 
 * copies all the data */
struct mesh_data *mesh_init_copy(struct mesh_data *mesh_source) {

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

  if(mesh_init_complete(mesh) == 1) {
    free(mesh);
    return NULL;
  }

  mesh_copy_data(mesh, mesh_source);

  return mesh;
}

int mesh_copy_data(struct mesh_data *mesh, struct mesh_data *mesh_source) {
  long int size;

  size = mesh->imax * mesh->jmax * mesh->kmax;

  if(size <= 0) {
    printf("error: improper size in mesh_copy_data\n");
    return 1;
  }
  
  memcpy(mesh->P, mesh_source->P, size * sizeof(double));
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

  return 0;
}

struct mesh_data *mesh_init_empty() {
  int i;
  struct mesh_data *mesh;
 
  mesh = malloc(sizeof(struct mesh_data));
  if (mesh == NULL) {
    printf("error: could not allocate mesh_data in mesh_init_empty\n");
    return(NULL);
  }

  mesh->imax = 20;
  mesh->jmax = 10;
  mesh->kmax = 10;

  mesh->delx = 0.075;
  mesh->dely = 0.075;
  mesh->delz = 0.075;
  mesh->rdx = 1.0 / mesh->delx;
  mesh->rdy = 1.0 / mesh->dely;
  mesh->rdz = 1.0 / mesh->delz;

  mesh->origin[0] = 0;
  mesh->origin[1] = 0;
  mesh->origin[2] = 0;

  mesh->inside[0] = 0.1;
  mesh->inside[1] = 0.1;
  mesh->inside[2] = 0.1;
  
  mesh->compress = 1;

  for(i = 0; i < 6; i++) {
    mesh->wb[i] = slip;
    mesh->sb[i] = NULL;
    if(i<3) mesh->baffles[i] = NULL;
  }

  mesh->P = NULL; 
  mesh->u = NULL;
  mesh->v = NULL;
  mesh->w = NULL;
  mesh->u_omega = NULL;
  mesh->v_omega = NULL;
  mesh->w_omega = NULL;
  
  mesh->vof = NULL;
  mesh->n_vof = NULL;

  mesh->fv = NULL;
  mesh->ae = NULL;
  mesh->an = NULL;
  mesh->at = NULL;

  mesh->turbulence_model = NULL;
  mesh->nut = NULL;
  
  mesh->ready = 0;

  return mesh;
}

#define max(a,b) (a>b?a:b)
#define min(a,b) (a<b?a:b)

int mesh_set_hydrostatic(struct mesh_data *mesh, double g, double rho) {
  long int i,j,k;


  mesh_set_array(mesh, "P", 0, -1, 0, 0, 0, 0, 0);

  for(i = 0; i < mesh->imax; i++) {
    for(j = 0; j < mesh->jmax; j++) {
      for(k = mesh->kmax-2; k > -1; k--) {
        if(mesh->vof[mesh_index(mesh,i,j,k)] == 1.0) {
          mesh->P[mesh_index(mesh,i,j,k)] = 
            mesh->P[mesh_index(mesh,i,j,k+1)] + mesh->delz * rho * g;
        }
        else  if(mesh->vof[mesh_index(mesh,i,j,k)] > 0.0) {
          mesh->P[mesh_index(mesh,i,j,k)] = 
            mesh->delz * rho * g * max(-0.5,mesh->vof[mesh_index(mesh,i,j,k)] - 0.5);          
        }
      }  
    }
  }

  return 0;
}

int mesh_set_array(struct mesh_data *mesh, char *param, double value, 
                   long int imin, long int imax,
                   long int jmin, long int jmax,
                   long int kmin, long int kmax) {
  long int i,j,k;
  double *p;

  /* negative value means fill entire mesh */
  if(imin < 0 || jmin < 0 || kmin < 0 ||
     imax < 0 || jmax < 0 || kmax < 0) {
    imin = jmin = kmin = 0;
    imax = mesh->imax;
    jmax = mesh->jmax;
    kmax = mesh->kmax;  
  }

  if(imax > mesh->imax) imax = mesh->imax;
  if(jmax > mesh->jmax) jmax = mesh->jmax;
  if(kmax > mesh->kmax) kmax = mesh->kmax;

  if(strcmp(param, "vof") == 0) p = mesh->vof;
  else if(strcmp(param, "P") == 0) p = mesh->P;
  else if(strcmp(param, "u") == 0) p = mesh->u;
  else if(strcmp(param, "v") == 0) p = mesh->v;
  else if(strcmp(param, "w") == 0) p = mesh->w;
  else if(strcmp(param, "u_omega") == 0) p = mesh->u_omega;
  else if(strcmp(param, "v_omega") == 0) p = mesh->v_omega;
  else if(strcmp(param, "w_omega") == 0) p = mesh->w_omega;
  else if(strcmp(param, "fv") == 0) p = mesh->fv;
  else if(strcmp(param, "ae") == 0) p = mesh->ae;
  else if(strcmp(param, "an") == 0) p = mesh->an;
  else if(strcmp(param, "at") == 0) p = mesh->at;
  else if(strcmp(param, "nut") == 0) p = mesh->nut;
  else {
    printf("warning: mesh_set_array could not recognize param %s\n",param);
    return 0;
  }

  for(i = imin; i < imax; i++) {
    for(j = jmin; j < jmax; j++) {
      for(k = kmin; k < kmax; k++) {
        p[mesh_index(mesh,i,j,k)] = value;
      }  
    }
  }

  return 0;
}

int mesh_sb_create(struct mesh_data *mesh, int wall, int type, double value, double turbulence) {
  struct sb_data *item, *sb;
  
  item = malloc(sizeof(struct sb_data));
  if(item == NULL) {
    printf("error: could not malloc in mesh_sb_create\n");
    return(1);
  }
  
  item->extent_a[0] = 0;
  item->extent_a[1] = 0;
  item->extent_b[0] = 0;
  item->extent_b[1] = 0;  
  item->type        = type;
  item->value       = value;
  item->turbulence  = turbulence;
  item->next        = NULL;
  
  if(mesh->sb[wall] == NULL) {
    mesh->sb[wall] = item;
  }
  else {
    for(sb = mesh->sb[wall]; sb != NULL; sb = sb->next) {
      if(sb->next == NULL) {
        sb->next = item;
        return 0;
      }
    }
  }
  
  return 0;
}

int mesh_sb_extent_a(struct mesh_data *mesh, int wall, long int extent_a_1, long int extent_a_2) {
  struct sb_data *sb;
  
  if(mesh->sb[wall] == NULL) {
    printf("error: attempting to set value on an empty special boundary\n");
    return(1);
  }
  
  sb = mesh->sb[wall];
  while(sb->next != NULL) sb = sb->next;
  
  sb->extent_a[0] = extent_a_1;
  sb->extent_a[1] = extent_a_2;
  
  return 0;
}

int mesh_sb_extent_b(struct mesh_data *mesh, int wall, long int extent_b_1, long int extent_b_2) {
  struct sb_data *sb;
  
  if(mesh->sb[wall] == NULL) {
    printf("error: attempting to set value on an empty special boundary\n");
    return(1);
  }
  
  sb = mesh->sb[wall];
  while(sb->next != NULL) sb = sb->next;
  
  sb->extent_b[0] = extent_b_1;
  sb->extent_b[1] = extent_b_2;
  
  return 0;
}
int mesh_baffle_create(struct mesh_data *mesh, int axis, int type, double value, long int pos) {
  struct baffle_data *item, *baffle;
  
  item = malloc(sizeof(struct baffle_data));
  if(item == NULL) {
    printf("error: could not malloc in mesh_baffle_create\n");
    return(1);
  }
  
  item->extent_a[0] = 0;
  item->extent_a[1] = 0;
  item->extent_b[0] = 0;
  item->extent_b[1] = 0;  
  item->type        = type;
  item->value       = value;
  item->pos         =  pos;
  item->next        = NULL;
  
  if(mesh->baffles[axis] == NULL) {
    mesh->baffles[axis] = item;
  }
  else {
    for(baffle = mesh->baffles[axis]; baffle != NULL; baffle = baffle->next) {
      if(baffle->next == NULL) {
        baffle->next = item;
        return 0;
      }
    }
  }
  
  return 0;
}

int mesh_baffle_extent_a(struct mesh_data *mesh, int axis, long int extent_a_1, long int extent_a_2) {
  struct baffle_data *baffle;
  
  if(mesh->baffles[axis] == NULL) {
    printf("error: attempting to set value on an empty special boundary\n");
    return(1);
  }
  
  baffle = mesh->baffles[axis];
  while(baffle->next != NULL) baffle = baffle->next;
  
  baffle->extent_a[0] = extent_a_1;
  baffle->extent_a[1] = extent_a_2;
  
  return 0;
}

int mesh_baffle_extent_b(struct mesh_data *mesh, int axis, long int extent_b_1, long int extent_b_2) {
  struct baffle_data *baffle;
  
  if(mesh->baffles[axis] == NULL) {
    printf("error: attempting to set value on an empty special boundary\n");
    return(1);
  }
  
  baffle = mesh->baffles[axis];
  while(baffle->next != NULL) baffle = baffle->next;
  
  baffle->extent_b[0] = extent_b_1;
  baffle->extent_b[1] = extent_b_2;
  
  return 0;
}

int mesh_set_value(struct mesh_data *mesh, char *param, int dims, 
                   double *vector) {
  int wall, axis;

  #ifdef DEBUG
    printf("mesh_set_value: %s, %d\n",param,dims);
  #endif

  if(mesh == NULL || param == NULL || vector == NULL) {
    printf("error: null values passed to mesh_set_value\n");
    return(1);
  }

  if(strcmp(param, "cells")==0) {
    if(dims != 3) {
      printf("error in source file: cells requires 3 arguments\n");
      return(1);
    }

    mesh->imax = vector[0];
    mesh->jmax = vector[1];
    mesh->kmax = vector[2];
  }
  else if (strcmp(param, "del")==0) {
    if(dims != 3) {
      printf("error in source file: del requires 3 arguments\n");
      return(1);
    }

    mesh->delx = vector[0];
    mesh->dely = vector[1];
    mesh->delz = vector[2];
    mesh->rdx = 1.0 / mesh->delx;
    mesh->rdy = 1.0 / mesh->dely;
    mesh->rdz = 1.0 / mesh->delz;
  }
  else if (strcmp(param, "inside")==0) {
    if(dims != 3) {
      printf("error in source file: inside requires 3 arguments\n");
      return(1);
    }

    mesh->inside[0] = vector[0];
    mesh->inside[1] = vector[1];
    mesh->inside[2] = vector[2];
  }
  else if (strcmp(param, "origin")==0) {
    if(dims != 3) {
      printf("error in source file: origin requires 3 arguments\n");
      return(1);
    }

    mesh->origin[0] = vector[0];
    mesh->origin[1] = vector[1];
    mesh->origin[2] = vector[2];
  } 
  else if (strncmp(param, "wall", 4)==0) {
    if(dims != 1) {
      printf("error in source file: wall_* requires 1 argument\n");
      return(1);
    }

    if(strcmp(param, "wall_north")==0) mesh->wb[3] = vector[0];
    else if(strcmp(param, "wall_south")==0) mesh->wb[2] = vector[0];
    else if(strcmp(param, "wall_east")==0) mesh->wb[1] = vector[0];
    else if(strcmp(param, "wall_west")==0) mesh->wb[0] = vector[0];
    else if(strcmp(param, "wall_top")==0) mesh->wb[5] = vector[0];
    else if(strcmp(param, "wall_bottom")==0) mesh->wb[4] = vector[0];
    else {
      printf("error in source file: %s unrecognized wall_* command\n",
             param);
      return(1);
    }
  }
  else if (strncmp(param, "sb_", 3)==0) {
    if(dims != 3) {
      printf("error in source file: sb_* requires 3 arguments\n");
      return(1);
    }

    if(strcmp(param, "sb_north")==0) wall = 3;
    else if(strcmp(param, "sb_south")==0) wall = 2;
    else if(strcmp(param, "sb_east")==0) wall = 1;
    else if(strcmp(param, "sb_west")==0) wall = 0;
    else if(strcmp(param, "sb_top")==0) wall = 5;
    else if(strcmp(param, "sb_bottom")==0) wall = 4;
    else {
      printf("error in source file: %s unrecognized sb_* command\n",
             param);
      return(1);
    }
    
    mesh_sb_create(mesh, wall, vector[0], vector[1], vector[2]);
  }
  else if (strncmp(param, "sbextent_a_", 11)==0) {
    if(dims != 2) {
      printf("error in source file: sbextent_a_* requires 2 arguments\n");
      return(1);
    }

    if(strcmp(param, "sbextent_a_north")==0) wall = 3;
    else if(strcmp(param, "sbextent_a_south")==0) wall = 2;
    else if(strcmp(param, "sbextent_a_east")==0) wall = 1;
    else if(strcmp(param, "sbextent_a_west")==0) wall = 0;
    else if(strcmp(param, "sbextent_a_top")==0) wall = 5;
    else if(strcmp(param, "sbextent_a_bottom")==0) wall = 4;
    else {
      printf("error in source file: %s unrecognized sbextent_a_* command\n",
             param);
      return(1);
    }
    
    mesh_sb_extent_a(mesh, wall, vector[0], vector[1]);
  }
  else if (strncmp(param, "sbextent_b_", 11)==0) {
    if(dims != 2) {
      printf("error in source file: sbextent_b_* requires 2 arguments\n");
      return(1);
    }

    if(strcmp(param, "sbextent_b_north")==0) wall = 3;
    else if(strcmp(param, "sbextent_b_south")==0) wall = 2;
    else if(strcmp(param, "sbextent_b_east")==0) wall = 1;
    else if(strcmp(param, "sbextent_b_west")==0) wall = 0;
    else if(strcmp(param, "sbextent_b_top")==0) wall = 5;
    else if(strcmp(param, "sbextent_b_bottom")==0) wall = 4;
    else {
      printf("error in source file: %s unrecognized sbextent_b_* command\n",
             param);
      return(1);
    }
    
    mesh_sb_extent_b(mesh, wall, vector[0], vector[1]);
  }
  /* baffles */
  else if (strncmp(param, "baffle_extent_a_", 16)==0) {
    if(dims != 2) {
      printf("error in source file: baffle_extent_a_* requires 2 arguments\n");
      return(1);
    }

    if(strcmp(param, "baffle_extent_a_x")==0) axis = 0;
    else if(strcmp(param, "baffle_extent_a_y")==0) axis = 1;
    else if(strcmp(param, "baffle_extent_a_z")==0) axis = 2;
    else {
      printf("error in source file: %s unrecognized baffle_extent_a_* command\n",
             param);
      return(1);
    }
    
    mesh_baffle_extent_a(mesh, axis, vector[0], vector[1]);
  }
  else if (strncmp(param, "baffle_extent_b_", 16)==0) {
    if(dims != 2) {
      printf("error in source file: baffleextent_b_* requires 2 arguments\n");
      return(1);
    }

    if(strcmp(param, "baffle_extent_b_x")==0) axis = 0;
    else if(strcmp(param, "baffle_extent_b_y")==0) axis = 1;
    else if(strcmp(param, "baffle_extent_b_z")==0) axis = 2;
    else {
      printf("error in source file: %s unrecognized baffle_extent_b_* command\n",
             param);
      return(1);
    }
    
    mesh_baffle_extent_b(mesh, axis, vector[0], vector[1]);
  }
  else if (strncmp(param, "baffle_", 6)==0) {
    if(dims != 3) {
      printf("error in source file: baffle_* requires 3 arguments\n");
      return(1);
    }

    if(strcmp(param, "baffle_x")==0) axis = 0;
    else if(strcmp(param, "baffle_y")==0) axis = 1;
    else if(strcmp(param, "baffle_z")==0) axis = 2;
    else {
      printf("error in source file: %s unrecognized baffle_* command\n",
             param);
      return(1);
    }
    
    mesh_baffle_create(mesh, axis, vector[0], vector[1], vector[2]);
  }
  else if(strncmp(param, "end", 3)==0) {
    return(0);
  }
  else {
    printf("error in source file: unrecognized command: %s\n", param);
    return(1);
  }

  return(0);
}

int mesh_check(struct mesh_data *mesh) {
  /* check if mesh is ready for solution
   * set mesh->ready to 1 and return TRUE 
   *
   * note: does not check cells
   * this function should be run before allocating cells */
  int flg=1;
  int i;

  if(mesh->imax == -1) flg = 0;
  if(mesh->jmax == -1) flg = 0;
  if(mesh->kmax == -1) flg = 0;

  if(mesh->delx == -1) flg = 0;
  if(mesh->dely == -1) flg = 0;
  if(mesh->delz == -1) flg = 0;

  if(mesh->origin[0] == -1) flg = 0;
  if(mesh->origin[1] == -1) flg = 0;
  if(mesh->origin[2] == -1) flg = 0;

  if(mesh->inside[0] == -1) flg = 0;
  if(mesh->inside[1] == -1) flg = 0;
  if(mesh->inside[2] == -1) flg = 0;


  if(flg == 1)
    mesh->ready = 1;
  else
    mesh->ready = 0;

  #ifdef DEBUG
    printf("mesh_check: flg=%d\n",flg);
    printf("imax, jmax, kmax: %ld, %ld, %ld\n", mesh->imax, mesh->jmax, mesh->kmax);
    printf("delx, dely, delz: %lf, %lf, %lf\n", mesh->delx, mesh->dely, mesh->delz);
    printf("origin: %lf, %lf, %lf\n", mesh->origin[0], mesh->origin[1], mesh->origin[2]);
    printf("inside: %lf, %lf, %lf\n", mesh->inside[0], mesh->inside[1], mesh->inside[2]);
    
    for(i = 0; i < 6; i++) {
      printf("wb %d: %d\n", i, mesh->wb[i]);
    }

    printf("ready: %d\n", mesh->ready);
  #endif

  return flg;
}


int mesh_init_complete(struct mesh_data *mesh) {

  long int size;

  if(mesh == NULL) {
    printf("error: null mesh passed to mesh_init_complete\n");
    return (1);
  }
  
  mesh_check(mesh);
  if(mesh->ready == 0) {
    printf("error: uninitialized mesh passed to mesh_init_complete\n");
    return (1);
  }
  
  size = mesh->imax * mesh->jmax * mesh->kmax;

  mesh->i_range = mesh->imax;
  mesh->i_start = 0;

  if(size<=0) {
    printf("error: improper size mesh in mesh_init_complete\n");
    return (1);
  }

	return mesh_allocate(mesh, size);
}

int mesh_allocate(struct mesh_data *mesh, long int size) {
  /* allocate memory */
  mesh->P = malloc(sizeof(double) * size);
  
  if(mesh->P == NULL) {
    printf("error: memory could not be allocated for P in mesh_init_complete\n");
    return (1);
  }
  
  mesh->u = malloc(sizeof(double) * size);

  if(mesh->u == NULL) {
    printf("error: memory could not be allocated for u in mesh_init_complete\n");
    return(1);
  }
  
  mesh->v = malloc(sizeof(double) * size);

  if(mesh->v == NULL) {
    printf("error: memory could not be allocated for v in mesh_init_complete\n");
    return(1);
  }

  mesh->w = malloc(sizeof(double) * size);

  if(mesh->w == NULL) {
    printf("error: memory could not be allocated for w in mesh_init_complete\n");
    return(1);
  }
  
  /* vorticity */
  mesh->u_omega = malloc(sizeof(double) * size);
  
  if(mesh->u_omega == NULL) {
    printf("error: memory could not be allocated for u_omega in mesh_init_complete\n");
    return(1);
  }
  
  mesh->v_omega = malloc(sizeof(double) * size);

  if(mesh->v_omega == NULL) {
    printf("error: memory could not be allocated for v_omega in mesh_init_complete\n");
    return(1);
  }

  mesh->w_omega = malloc(sizeof(double) * size);

  if(mesh->w_omega == NULL) {
    printf("error: memory could not be allocated for w_omega in mesh_init_complete\n");
    return(1);
  }
 
  mesh->vof = malloc(sizeof(double) * size);

  if(mesh->vof == NULL) {
    printf("error: memory could not be allocated for vof in mesh_init_complete\n");
    return(1);
  }

  mesh->n_vof = malloc(sizeof(enum cell_boundaries) * size);

  if(mesh->n_vof == NULL) {
    printf("error: memory could not be allocated for n_vof in mesh_init_complete\n");
    return(1);
  } 

  mesh->fv = malloc(sizeof(double) * size);

  if(mesh->fv == NULL) {
    printf("error: memory could not be allocated for fv in mesh_init_complete\n");
    return(1);
  }

  mesh->ae = malloc(sizeof(double) * size);

  if(mesh->ae == NULL) {
    printf("error: memory could not be allocated for ae in mesh_init_complete\n");
    return(1);
  }

  mesh->an = malloc(sizeof(double) * size);

  if(mesh->an == NULL) {
    printf("error: memory could not be allocated for an in mesh_init_complete\n");
    return(1);
  }

  mesh->at = malloc(sizeof(double) * size);

  if(mesh->at == NULL) {
    printf("error: memory could not be allocated for at in mesh_init_complete\n");
    return(1);
  }

  mesh->nut  = malloc(sizeof(double) * size);
   if(mesh->nut == NULL) {
    printf("error: memory could not be allocated for nut in mesh_init_complete\n");
    return(1);
  } 
  
  return(0);
}

int mesh_free(struct mesh_data *mesh) {
  int i;
  struct sb_data *sb;
  struct sb_data *sb_free;

  free(mesh->P);
  free(mesh->u);
  free(mesh->v);
  free(mesh->w);
  free(mesh->u_omega);
  free(mesh->v_omega);
  free(mesh->w_omega);

  free(mesh->vof);
  free(mesh->n_vof);

  free(mesh->fv);
  free(mesh->ae);
  free(mesh->an);
  free(mesh->at);
  
  free(mesh->nut);

  for(i = 0; i < 6; i++) {
    sb = mesh->sb[i];
    
    while(sb != NULL) {
      sb_free = sb;
      sb = sb->next;
      
      free(sb_free);
    }
    
    mesh->sb[i] = NULL;
  }

  return(0);
}

#ifdef _WIN32
long int mesh_index(struct mesh_data *mesh,
	long int i, long int j, long int k)
#else
inline long int mesh_index(struct mesh_data *mesh, 
                    long int i, long int j, long int k) 
#endif
{
  /* return i+j*mesh->imax+k*mesh->imax*mesh->jmax; */
  /* return i + mesh->imax * (j + k * mesh->jmax); */
  return k + mesh->kmax * (j + i * mesh->jmax);
}

int mesh_fill_vof(struct mesh_data *mesh, double *vector) {
  /* fill the mesh from the inside vector point until either a VOF > 0 or FV < 1 is encountered */

  struct point_data *stack = NULL, *p;

  stack = stack_push(stack, vector[0]/mesh->delx, vector[1]/mesh->dely, vector[2]/mesh->delz);

  while((p = stack_pop(stack)) != NULL) {
    
    if(p->i >= mesh->imax || p->j >= mesh->jmax || p->k >= mesh->kmax || 
       p->i < 0           || p->j < 0           || p->k < 0 ) {
      free(p);
      continue;
    }

    if(mesh->fv[mesh_index(mesh, p->i, p->j, p->k)]  > 0.0 && 
       mesh->vof[mesh_index(mesh, p->i, p->j, p->k)] < 1.0) {
    

      mesh->vof[mesh_index(mesh, p->i, p->j, p->k)] = 1.0;
      
      stack = stack_push(stack, p->i + 1, p->j, p->k);
      stack = stack_push(stack, p->i - 1, p->j, p->k);

      stack = stack_push(stack, p->i, p->j + 1, p->k);
      stack = stack_push(stack, p->i, p->j - 1, p->k);

      stack = stack_push(stack, p->i, p->j, p->k + 1);
      stack = stack_push(stack, p->i, p->j, p->k - 1);


    }
    free(p);
  
  }

  return 0;
}

int mesh_fill(struct mesh_data *mesh, struct stl_data *stl) {

  struct point_data *stack = NULL, *p;

  stack = stack_push(stack, (mesh->inside[0]-mesh->origin[0])/mesh->delx, 
    (mesh->inside[1]-mesh->origin[1])/mesh->dely, (mesh->inside[2]-mesh->origin[2])/mesh->delz);

  while((p = stack_pop(stack)) != NULL) {
    
    if(p->i >= mesh->imax || p->j >= mesh->jmax || p->k >= mesh->kmax || 
       p->i < 0           || p->j < 0           || p->k < 0 ) {
      free(p);
      continue;
    }

    if(mesh->fv[mesh_index(mesh, p->i, p->j, p->k)]  == 0.0) {
    
      if(stl_check_normals(mesh, stl, p->i, p->j, p->k)) {

        mesh->fv[mesh_index(mesh, p->i, p->j, p->k)] = 1.0;
        
        stack = stack_push(stack, p->i + 1, p->j, p->k);
        stack = stack_push(stack, p->i - 1, p->j, p->k);

        stack = stack_push(stack, p->i, p->j + 1, p->k);
        stack = stack_push(stack, p->i, p->j - 1, p->k);

        stack = stack_push(stack, p->i, p->j, p->k + 1);
        stack = stack_push(stack, p->i, p->j, p->k - 1);

      }
      else {
        printf("warning: mesh_fill cell %ld %ld %ld has reverse normals.  skipped.\n", p->i, p->j, p->k);
      }
    }
    free(p);
  
  }

  return 0;
}

int mesh_avratio(struct mesh_data *mesh, double avr_max) {
	
	double avr, avr_max_obs, avr_fix, r;
	long int im1, jm1, km1;
	int corrected;
	const double emf = 0.001;	
  long int i, j, k, size;
  
  double *ae_n;
  double *an_n;
  double *at_n;  
  
  size = mesh->imax * mesh->jmax * mesh->kmax;
  
  ae_n = malloc(sizeof(double) * size);
  if(ae_n == NULL) {
    printf("error: memory could not be allocated for ae_n in mesh_avratio\n");
    return(1);
  }
  an_n = malloc(sizeof(double) * size);
  if(an_n == NULL) {
    printf("error: memory could not be allocated for an_n in mesh_avratio\n");
    return(1);
  }
  at_n = malloc(sizeof(double) * size);
  if(at_n == NULL) {
    printf("error: memory could not be allocated for at_n in mesh_avratio\n");
    return(1);
  }
  
  for(i = 0; i < mesh->imax; i++) {
    for(j = 0; j < mesh->jmax; j++) {
      for(k = 0; k < mesh->kmax; k++) {  
      	ae_n[mesh_index(mesh,i,j,k)] = AE(i,j,k);
      	an_n[mesh_index(mesh,i,j,k)] = AN(i,j,k);
      	at_n[mesh_index(mesh,i,j,k)] = AT(i,j,k);
  		}
  	}
  }

	avr_max_obs = 0;
	avr_fix = 0;
	
  for(i = 0; i < mesh->imax; i++) {
    for(j = 0; j < mesh->jmax; j++) {
      for(k = 0; k < mesh->kmax; k++) {
				if(FV(i,j,k) < emf) continue; 
				if(FV(i,j,k) > (1-emf)) continue;
				
				corrected = 0;
				im1 = max(i-1, 0);
				jm1 = max(j-1, 0);
				km1 = max(k-1, 0);
				/*
				avr = AE(i,j,k)/min(FV(i,j,k),FV(i+1,j,k));
				if(!isnan(avr) && avr > avr_max) {
					corrected = 1;
					avr_fix = avr;
					AE(i,j,k) = AE(i,j,k) * avr_max / avr;

				}				
				
				avr = AN(i,j,k)/min(FV(i,j,k),FV(i,j+1,k));
				if(!isnan(avr) && avr > avr_max) {
					corrected = 1;
					avr_fix = avr;
					AN(i,j,k) = AN(i,j,k) * avr_max / avr;
				}				
				
				avr = AT(i,j,k)/min(FV(i,j,k),FV(i,j,k+1));
				if(!isnan(avr) && avr > avr_max) {
					corrected = 1;
					avr_fix = avr;
					AT(i,j,k) = AT(i,j,k) * avr_max / avr;
				} */
				

				avr = max(AE(i,j,k)/FV(i,j,k), AE(max(0,i-1),j,k)/FV(i,j,k));
				if(avr > avr_max + emf && !isnan(avr)) {
					corrected = 1;
					avr_fix = avr;
					r = avr / avr_max;
					FV(i,j,k) = FV(i,j,k) * r;
					
					/* correct all related areas perpendicular */
					mesh_area_correct(&AN(i,j,k), &AN(i,jm1,k), an_n[mesh_index(mesh,i,j,k)], 
														an_n[mesh_index(mesh,i,jm1,k)], r);		
					mesh_area_correct(&AT(i,j,k), &AT(i,j,km1), at_n[mesh_index(mesh,i,j,k)], 
														at_n[mesh_index(mesh,i,j,km1)], r);	
				  avr_max_obs = max(avr_fix,avr_max_obs);							

				}
				
				avr = max(AN(i,j,k)/FV(i,j,k), AN(i,max(0,j-1),k)/FV(i,j,k));
				if(avr > avr_max + emf && !isnan(avr)) {
					corrected = 1;
					avr_fix = avr;
					r = avr / avr_max;
					FV(i,j,k) = FV(i,j,k) * r;	
					
					/* TESTING: correct all related areas perpendicular */
					mesh_area_correct(&AE(i,j,k), &AE(im1,j,k), ae_n[mesh_index(mesh,i,j,k)], 
														ae_n[mesh_index(mesh,im1,j,k)], r);		
					mesh_area_correct(&AT(i,j,k), &AT(i,j,km1), at_n[mesh_index(mesh,i,j,k)], 
														at_n[mesh_index(mesh,i,j,km1)], r);		
				  avr_max_obs = max(avr_fix,avr_max_obs);			
					
				}			

				avr = max(AT(i,j,k)/FV(i,j,k), AT(i,j,max(0,k-1))/FV(i,j,k));
				if(avr > avr_max + emf && !isnan(avr)) {
					corrected = 1;
					avr_fix = avr;
					r = avr / avr_max;
					FV(i,j,k) = FV(i,j,k) * r;	
					
					/* TESTING: correct all related areas perpendicular */	
					mesh_area_correct(&AE(i,j,k), &AE(im1,j,k), ae_n[mesh_index(mesh,i,j,k)], 
														ae_n[mesh_index(mesh,im1,j,k)], r);		
					mesh_area_correct(&AN(i,j,k), &AN(i,jm1,k), an_n[mesh_index(mesh,i,j,k)], 
														an_n[mesh_index(mesh,i,jm1,k)], r);	
				  avr_max_obs = max(avr_fix,avr_max_obs);
				} 

				if(corrected == 1) {
					printf("Corrected AV ratio in cell %ld %ld %ld from %e to %e\n", i, j, k, avr_fix, avr_max);

				}
				
			}
		}
	}		
	
	if(avr_max_obs > avr_max)
		printf("Maximum observed avr: %e\n", avr_max_obs);
	
	free(ae_n);
	free(an_n);
	free(at_n);
	
	return 0;

}

int mesh_area_correct(double *a1, double *a2, double an1, double an2, double r) {
	const double emf = 0.001;
  double ave_a, del_a;
  
	if(an1 < emf) {
		*a2 = max(*a2, an2 * r);						
	} 
	else if(an2 < emf) {
		*a1 = max(*a1, an1 * r);
	} 
	else {
		ave_a = (an1 + an2) / 2;
		del_a = ave_a * r - ave_a;
		*a1 = max(*a1, del_a + an1);
		*a2 = max(*a2, del_a + an2);
	}
	
	return 0;
}

int mesh_normalize(struct mesh_data *mesh) {

  const double emf = 0.01;

  long int i, j, k;

  for(i = 0; i < mesh->imax; i++) {
    for(j = 0; j < mesh->jmax; j++) {
      for(k = 0; k < mesh->kmax; k++) {

        if(mesh->fv[mesh_index(mesh, i, j, k)] > 1-emf)
        {
          mesh->fv[mesh_index(mesh, i, j, k)] = 1;
          
          mesh->ae[mesh_index(mesh, i, j, k)] = 1;
          mesh->an[mesh_index(mesh, i, j, k)] = 1;
          mesh->at[mesh_index(mesh, i, j, k)] = 1;

          if(i>0) mesh->ae[mesh_index(mesh, i-1, j, k)] = 1;
          if(j>0) mesh->an[mesh_index(mesh, i, j-1, k)] = 1;
          if(k>0) mesh->at[mesh_index(mesh, i, j, k-1)] = 1;
        }      

        if(mesh->ae[mesh_index(mesh, i, j, k)] < emf)
          mesh->ae[mesh_index(mesh, i, j, k)] = 0;
        if(mesh->ae[mesh_index(mesh, i, j, k)] > 1-emf)
          mesh->ae[mesh_index(mesh, i, j, k)] = 1;

        if(mesh->an[mesh_index(mesh, i, j, k)] < emf)
          mesh->an[mesh_index(mesh, i, j, k)] = 0;
        if(mesh->an[mesh_index(mesh, i, j, k)] > 1-emf)
          mesh->an[mesh_index(mesh, i, j, k)] = 1;

        if(mesh->at[mesh_index(mesh, i, j, k)] < emf)
          mesh->at[mesh_index(mesh, i, j, k)] = 0;
        if(mesh->at[mesh_index(mesh, i, j, k)] > 1-emf)
          mesh->at[mesh_index(mesh, i, j, k)] = 1;      
      
			}
		}
	}

  for(i = 0; i < mesh->imax; i++) {
    for(j = 0; j < mesh->jmax; j++) {
      for(k = 0; k < mesh->kmax; k++) {
        
        if(mesh->fv[mesh_index(mesh, i, j, k)] < emf)
        {
          mesh->fv[mesh_index(mesh, i, j, k)] = 0;
          
          mesh->ae[mesh_index(mesh, i, j, k)] = 0;
          mesh->an[mesh_index(mesh, i, j, k)] = 0;
          mesh->at[mesh_index(mesh, i, j, k)] = 0;

          if(i>0) mesh->ae[mesh_index(mesh, i-1, j, k)] = 0;
          if(j>0) mesh->an[mesh_index(mesh, i, j-1, k)] = 0;
          if(k>0) mesh->at[mesh_index(mesh, i, j, k-1)] = 0;
        }



      }
    }
  }

  return 0;
}
