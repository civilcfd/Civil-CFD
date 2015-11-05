/* track.c
 *
 * tracks the solution steps and saves it to a file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "track.h"
#include "vtk.h"
#include "vof.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "readfile.h"

#include "vof_macros.h"

static struct track_data *track=NULL;
static struct track_data *p_track=NULL;

int track_rewind() {
  p_track = track;
  return 0;
}

int track_next() {
  if(p_track != NULL) p_track = p_track->next;
  if(p_track != NULL) return p_track->n;
  else return -1;
}

double track_t() {
  if(p_track !=NULL) return p_track->t;
  else return -1;
}

int track_get_n() {
  if(p_track !=NULL) return p_track->n;
  else return -1;
}

int track_add(double t) {
  struct track_data *node, *x;
  
  node = malloc(sizeof(struct track_data));
  if(node==NULL) {
    printf("error: could not allocate memory in track_add()\n");
    return -1;
  }
  node->t = t;
  
  if(track == NULL) {
    node->n = 0;
    track = node;
    track->next = NULL;
    track->first = node;
    return 0;
  }
  
  
  node->next = NULL;
  node->first = track;
  for(x = track; x->next != NULL; x = x->next) {
  }
  
  node->n = x->n + 1;
  x->next = node;
  
  return node->n;
}

int track_add_n(double t, int n) {
  struct track_data *node, *x;
  
  node = malloc(sizeof(struct track_data));
  if(node==NULL) {
    printf("error: could not allocate memory in track_add()\n");
    return -1;
  }
  node->t = t;
  
  if(track == NULL) {
    node->n = 0;
    track = node;
    track->next = NULL;
    track->first = node;
    return 0;
  }
  
  
  node->next = NULL;
  node->first = track;
  for(x = track; x->next != NULL; x = x->next) {
  }
  
  node->n = n;
  x->next = node;
  
  return node->n;
}

int track_write() {
  FILE *fp;
  struct track_data *x;
  
  fp = fopen("tracking", "w");

  if(fp == NULL) {
    printf("error: cannot open tracking file\n");
    return(1);
  }

  for(x = track; x != NULL; x = x->next) {
    fprintf(fp,"%4.3lf %d\n",x->t, x->n);
  }

  fclose(fp);

  return 0;
}

int track_read() {
  FILE *fp;
  struct track_data *x, *del;
  int nargs;
  char text[256];
  char args[2][256];

  track_empty();

  fp = fopen("tracking", "r");

  if(fp == NULL) {
    printf("warning: cannot open tracking file to read\n");
    return 1;
  }

  while(!feof(fp)) {
    if(!fgets(text, sizeof(text), fp))
    { 
      if(!feof(fp)) {
        printf("error: fgets in read_mesh\n");
        return(1);
      }
      else break;
    }

    nargs = read_args(text, 2, args);

    if(nargs < 2) break;
    
    track_add_n(strtod(args[0],NULL),atoi(args[1]));
  }

  /* while(!feof(fp))
  {
    fscanf(fp, "%le %d", &t, &n);
    track_add_n(t,n);
  } */

  fclose(fp);

  return 0;
}

int track_delete(int n) {
  struct track_data *x, *del;
  
  if(track == NULL) {
    return 0;
  }
  
  if(track->n == n) {
    x = track;
    track = track->next;
    free(x);
    return 1;
  }
  
  for(x = track; x->next != NULL; x = x->next) {
    if(x->next->n == n) {
      del = x->next;
      x->next = x->next->next;
      free(del);
      return 1;
    }
  }

  return 0;
}

int track_empty() {
  struct track_data *x, *del;
  
  x = track;
  track = NULL;
  while(x != NULL) {
    del = x;
    x = x->next;
    free(del);
  }
  
  return 0;
}

