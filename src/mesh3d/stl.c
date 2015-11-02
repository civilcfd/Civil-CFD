/* stl.c
 *
 * implementation of STL functionality 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "stl.h"
#include "mesh.h"

struct stl_data *stl_init_empty() {
  struct stl_data *stl;
  
  stl = malloc(sizeof(struct stl_data));
  if(stl == NULL) {
    printf("error: could not allocate stl_data in stl_init_empty\n");
    return(NULL);
  }

  stl->solid[0] = 0;
  stl->facets = 0;
  stl->ready = 0;

  return stl;
}

int stl_check(struct stl_data *stl) {
  int flg;
  long int i;

  flg = 1;

  if(stl->solid[0] == 0) flg=0;
  if(stl->facets == 0) flg=0;

  if(flg == 1) stl->ready=1;
  else stl->ready=0;

  #ifdef DEBUG
    printf("stl solid: %s\n", stl->solid);
    printf("facets: %ld\n", stl->facets);
    
    for(i=0; i < stl->facets; i++) {
      printf("facet: %ld\n",i);
      printf("vertex 1: %lf %lf %lf\n", stl->v_1[i][0], 
             stl->v_1[i][1], stl->v_1[i][2]);
      printf("vertex 2: %lf %lf %lf\n", stl->v_2[i][0], 
             stl->v_2[i][1], stl->v_2[i][2]);
      printf("vertex 3: %lf %lf %lf\n", stl->v_3[i][0], 
             stl->v_3[i][1], stl->v_3[i][2]);
    }

    printf("ready: %d\n", stl->ready);

  #endif

  return(flg);
}

int stl_set_value(struct stl_data *stl, int dims, 
                   char (*args)[256], int state) {
  
  /* states are as follows:
   * -1 = error
   * 0 = default state, not inside facet 
   * 1 = inside facet
   * 2 = inside outer loop / vertex 1
   * 3 = vertex 2
   * 4 = vertex 3
   * 5 = endloop
   * 6 = endfacet
   */
  char *p;
  long int facet;

  #ifdef DEBUG
    printf("stl_set_value: %s %s %s %s %s, %d\n",args[0],
           args[1], args[2], args[3], args[4], dims);
    printf("state: %d\n",state);
  #endif

  if(stl == NULL || args == NULL ) {
    printf("error: null values passed to mesh_set_value\n");
    return(1);
  }

  facet = stl->facets;

  if(strcmp(args[0], "solid")==0) {
    if(state != 0) {
      printf("error in stl file: solid in wrong part of file\n");
      return(-1);
    }

    strcpy(stl->solid, args[1]);
  }
  else if(strcmp(args[0], "facet")==0) {
    if(state != 0) {
      printf("error in stl file: facet in wrong part of file\n");
      return(-1);
    }

    if(strcmp(args[1], "normal")!=0) {
      printf("error in stl file: facet without normal\n");
      return(-1);
    }

    stl->normal[facet][0] = strtod(args[2], &p);
    stl->normal[facet][1] = strtod(args[3], &p);
    stl->normal[facet][2] = strtod(args[4], &p);
  
    state++;
  }
  else if(strcmp(args[0], "outer")==0) {
    if(state != 1) {
      printf("error in stl file: outer loop in wrong part of file\n");
      return(-1);
    }

    if(strcmp(args[1], "loop")!=0) {
      printf("error in stl file: outer without loop\n");
      return(-1);
    }
  
    state++;
  }
  else if(strcmp(args[0], "vertex")==0) {

    if(state == 2)
    {
      stl->v_1[facet][0] = strtod(args[1], &p);
      stl->v_1[facet][1] = strtod(args[2], &p);
      stl->v_1[facet][2] = strtod(args[3], &p);
    }
    else if(state == 3)
    {
      stl->v_2[facet][0] = strtod(args[1], &p);
      stl->v_2[facet][1] = strtod(args[2], &p);
      stl->v_2[facet][2] = strtod(args[3], &p);
    }
    else if(state == 4)
    {
      stl->v_3[facet][0] = strtod(args[1], &p);
      stl->v_3[facet][1] = strtod(args[2], &p);
      stl->v_3[facet][2] = strtod(args[3], &p);
    }
    else
    { 
      printf("error in stl file: vertex outside of outer loop\n");
      return(-1);
    }

    state++;
  } 
  else if(strcmp(args[0], "endloop")==0) {
    if(state != 5) {
      printf("error in stl file: endloop in wrong part of file\n");
      return(-1);
    }

    state++;
  } 
  else if(strcmp(args[0], "endfacet")==0) {
    if(state != 6) {
      printf("error in stl file: endfacet in wrong part of file\n");
      return(-1);
    }

    state=0;
    stl->facets++;
  } 
  
  return state;
}

int stl_free(struct stl_data *stl) {
  free(stl);
  return(0);
}

int stl_check_normals(struct mesh_data *mesh, struct stl_data *stl, 
                      long int i, long int j, long int k) {


  int facing, n, x;

  double r_o[3], pt_int[3], f;
  double normals[3][3] = { { 1, 0, 0 },
                                 { 0, 1, 0 },
                                 { 0, 0, 1} };

  int flg = 0;
  double f_min = 9999999;

  r_o[0] = mesh->origin[0] + mesh->delx * i + mesh->delx/2;
  r_o[1] = mesh->origin[1] + mesh->dely * j + mesh->dely/2;
  r_o[2] = mesh->origin[2] + mesh->delz * k + mesh->delz/2;

  for(n=0; n < stl->facets; n++) {
    
    for(x=0; x < 3; x++) { /* iterate through each normal direction */

      if((facing = moller_trumbore(r_o, normals[x], stl->v_1[n],
                                   stl->v_2[n], stl->v_3[n], pt_int)) != 0) {
      
        f = pt_int[x] - r_o[x];

        
        if( ( f > 0 && facing == -1) ||
            ( f < 0 && facing ==  1) ) {
          /* in this case, the normal is facing away from our point
           * check if this is the closest intersection
           * if so, set the flag to 1 */
          if(fabs(f) < f_min) {
            flg = 1;
            f_min = fabs(f);
          }
        }
        else {
          /* in this case, the normal is facing towards our point
           * check if this is the closest intersection
           * if so, set the flag to 0 */
          if(fabs(f) < f_min) {
            flg = 0;
            f_min = fabs(f);
          }
        }
      
      }
    }
  }

  if(flg == 0)
    return 1; /* success */
  else
    return 0;
}

int stl_check_normals_point(struct mesh_data *mesh, struct stl_data *stl, 
                      double *pt) {


  int facing, n, x;

  double r_o[3], pt_int[3], f;
  double normals[3][3] = { { 1, 0, 0 },
                                 { 0, 1, 0 },
                                 { 0, 0, 1} };

  int flg = 0;
  double f_min = 9999999;

  r_o[0] = pt[0];
  r_o[1] = pt[1];
  r_o[2] = pt[2];

  for(n=0; n < stl->facets; n++) {
    
    for(x=0; x < 3; x++) { /* iterate through each normal direction */

      if((facing = moller_trumbore(r_o, normals[x], stl->v_1[n],
                                   stl->v_2[n], stl->v_3[n], pt_int)) != 0) {
      
        f = pt_int[x] - r_o[x];

        
        if( ( f > 0 && facing == -1) ||
            ( f < 0 && facing ==  1) ) {
          /* in this case, the normal is facing away from our point
           * check if this is the closest intersection
           * if so, set the flag to 1 */
          if(fabs(f) < f_min) {
            flg = 1;
            f_min = fabs(f);
          }
        }
        else {
          /* in this case, the normal is facing towards our point
           * check if this is the closest intersection
           * if so, set the flag to 0 */
          if(fabs(f) < f_min) {
            flg = 0;
            f_min = fabs(f);
          }
        }
      
      }
    }
  }

  if(flg == 0)
    return 1; /* success */
  else
    return 0;
}

