/* readfile.c
 *
 * Opens the mesh sourcefile.  Also writes to a sourcefile.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "readfile.h"
#include "mesh.h"

const char *wall_names[] = { "west", "east", "south", "north", "bottom", "top" };
const char *axis_names[] = { "x", "y", "z" };

int write_mesh(struct mesh_data *mesh, char *filename) {
  FILE *fp;
  int i;
  struct sb_data *sb;
  struct baffle_data *baffle;
  
  if(filename == NULL || mesh == NULL) {
    printf("error: passed null arguments to write_mesh\n");
    return(1);
  }

  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: cannot open %s in write_mesh\n",filename);
    return(1);
  }
  
  fprintf(fp,"cells %ld %ld %ld\n",mesh->imax,mesh->jmax,mesh->kmax);
  fprintf(fp,"del %lf %lf %lf\n",mesh->delx,mesh->dely,mesh->delz);
  fprintf(fp,"inside %lf %lf %lf\n",mesh->inside[0],mesh->inside[1],mesh->inside[2]);
  fprintf(fp,"origin %lf %lf %lf\n",mesh->origin[0],mesh->origin[1],mesh->origin[2]);
  
  for(i=0; i<6; i++) {
    fprintf(fp,"wall_%s %d\n",wall_names[i],mesh->wb[i]);
  }

  for(i=0; i<6; i++) {
    sb = mesh->sb[i];
    
    while(sb != NULL) {
      fprintf(fp,"sb_%s %d %lf %lf\n",wall_names[i],sb->type,sb->value,sb->turbulence);
      fprintf(fp,"sbextent_a_%s %ld %ld\n",wall_names[i],sb->extent_a[0],sb->extent_a[1]);
      fprintf(fp,"sbextent_b_%s %ld %ld\n",wall_names[i],sb->extent_b[0],sb->extent_b[1]);

      sb = sb->next;
    }
  }
  
  for(i=0; i<3; i++) {
    baffle = mesh->baffles[i];
    
    while(baffle != NULL) {
      fprintf(fp,"baffle_%s %d %lf %ld\n",axis_names[i],baffle->type,baffle->value,baffle->pos);
      fprintf(fp,"baffle_extent_a_%s %ld %ld\n",axis_names[i],baffle->extent_a[0],baffle->extent_a[1]);
      fprintf(fp,"baffle_extent_b_%s %ld %ld\n",axis_names[i],baffle->extent_b[0],baffle->extent_b[1]);

      baffle = baffle->next;
    }
  }
  
  fprintf(fp,"end\n");
  
  fclose(fp);
  
  return 0;
}


int read_mesh(struct mesh_data *mesh, char *filename)
{
  FILE *fp;

  char text[1024];
  char args[4][256];
  double vector[3];

  char *s;
  int i, nargs;

  if(filename == NULL || mesh == NULL) {
    printf("error: passed null arguments to read_mesh\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in read_mesh\n",filename);
    return(1);
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

    nargs = read_args(text, 4, args);

    if(nargs > 0) nargs--;
    for(i = nargs; i < 3; i++) {
      vector[i] = 0;
    }

    for(i = 0; i < nargs; i++) {
      /* sscanf(args[i], "%lf", &vector[i]); */
      vector[i] = strtod(args[i+1],&s); 
    }

    #ifdef DEBUG
      printf("read_mesh: read %d args: %s %lf %lf %lf\n", nargs, args[0], 
             vector[0], vector[1], vector[2]);
    #endif

    if (args[0][0] != '#' && args[0][0] != 0) 
      mesh_set_value(mesh, args[0], i, vector); 
    args[0][0] = 0;
  }

  fclose(fp);

  return 0;
}

int read_stl(struct stl_data *stl, char *filename, double *limits) {
  FILE *fp;
  
  char text[80];
  char args[5][256];
  int nargs;
  
  if(filename == NULL || stl == NULL) {
    printf("error: passed null arguments to read_stl\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in read_stl\n",filename);
    return(1);
  }
  
  if(!fgets(text, sizeof(text), fp)) { 
      if(!feof(fp)) {
        printf("error: fgets in read_stl\n");
        return(1);
      }
  }
  fclose(fp);
  
  nargs = read_args(text, 5, args);
  
  if(strcmp(args[0], "solid") == 0) {
    /* detected ASCII stl file */
    return read_stl_ascii(stl, filename, limits);
  } else {
    return read_stl_binary(stl, filename, limits);
  }
  
}

int read_stl_binary(struct stl_data *stl, char *filename, double *limits) {
  FILE *fp;
  char text[80];
  float buf[12];
  int n, count, i;
  long int res;
  
  if(filename == NULL || stl == NULL) {
    printf("error: passed null arguments to read_stl\n");
    return(1);
  }

  fp = fopen(filename, "rb");

  if(fp == NULL) {
    printf("error: cannot open %s in read_stl\n",filename);
    return(1);
  }
  
  res = fread(&text, 80, 1, fp);
  if(res != 1) {
    printf("error: misformed binary stl file %s\n",filename);
    return(1);
  }  
  
  res = fread(&count, 4, 1, fp);
  if(res != 1) {
    printf("error: misformed binary stl file %s\n",filename);
    return(1);
  }
  
  i=0;
  stl->facets = 0;
  while(!feof(fp) && i < count) {
    i++;
  
    if(fread(&buf, 4, 12, fp) != 12) {
      printf("error: misformed binary stl file %s\n",filename);
      return(1);
    }
    if(fread(&n, 2, 1, fp) != 1) {
      printf("error: misformed binary stl file %s\n",filename);
      return(1);
    }
    
    stl->normal[stl->facets][0] = buf[0];
    stl->normal[stl->facets][1] = buf[1];
    stl->normal[stl->facets][2] = buf[2];
    stl->v_1[stl->facets][0]    = buf[3];
    stl->v_1[stl->facets][1]    = buf[4];
    stl->v_1[stl->facets][2]    = buf[5];
    stl->v_2[stl->facets][0]    = buf[6];
    stl->v_2[stl->facets][1]    = buf[7];
    stl->v_2[stl->facets][2]    = buf[8];
    stl->v_3[stl->facets][0]    = buf[9];
    stl->v_3[stl->facets][1]    = buf[10];
    stl->v_3[stl->facets][2]    = buf[11];  

    if((stl->v_1[stl->facets][0] > limits[0] && stl->v_2[stl->facets][0] > limits[0] && stl->v_3[stl->facets][0] > limits[0]) ||
       (stl->v_1[stl->facets][1] > limits[1] && stl->v_2[stl->facets][1] > limits[1] && stl->v_3[stl->facets][1] > limits[1]) ||
       (stl->v_1[stl->facets][2] > limits[2] && stl->v_2[stl->facets][2] > limits[2] && stl->v_3[stl->facets][2] > limits[2])) {
    } else {
      stl->facets++;
    }
  
  }  
  fclose(fp);
  
  strncpy(stl->solid, "binary", 6);
  
  return 0;
}

int read_stl_ascii(struct stl_data *stl, char *filename, double *limits) {

  FILE *fp;

  char text[1024];
  char args[5][256];

  int nargs;

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
  int state=0; 
  
  if(filename == NULL || stl == NULL) {
    printf("error: passed null arguments to read_stl\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in read_stl\n",filename);
    return(1);
  }

  while(!feof(fp)) {
    if(!fgets(text, sizeof(text), fp))
    { 
      if(!feof(fp)) {
        printf("error: fgets in read_stl\n");
        return(1);
      }
      else break;
    }

    nargs = read_args(text, 5, args);
   
    #ifdef DEBUG
      printf("read_stl: read %d args: %s %s %s %s %s\n", nargs,  
             args[0], args[1], args[2], args[3], args[4]);
    #endif

    state = stl_set_value(stl, nargs, args, state, limits); 
  }

  fclose(fp);

  return 0;
}

#ifndef min
#define min(a,b)  (a>b?b:a)
#endif

int read_args(char *text, int nargs, char (*args)[256])
{
  char *s, *p;
  int n, i;

  s = trimwhitespace(text);
  n = 0;
  
  while((p = strtok_r(s," ", &s))!=NULL) {

    /* strncpy(args[n], p, min(strlen(p)+1, sizeof(args[n])-1));

    p = trimwhitespace(args[n]); */
    p = trimwhitespace(p);
    strncpy(args[n], p, min(strlen(p)+1, sizeof(args[n])-1));

    n++;
    if(n>=nargs) break;
  }

  for(i=n; i<nargs; i++) {
    args[i][0] = 0;
  }

  return n;
}


char *trimwhitespace(char *str)
{
  char *end;

  /* Trim leading space */
  while(isspace(*str)) str++;

  if(*str == 0)  /* All spaces? */
    return str;

  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  
  /* Write new null terminator */
  *(end+1) = 0;
  
  return str;
}

#ifdef _WIN32
char* strtok_r(
    char *str, 
    const char *delim, 
    char **nextp)
{
    char *ret;

    if (str == NULL)
    {
        str = *nextp;
    }

    str += strspn(str, delim);

    if (*str == '\0')
    {
        return NULL;
    }

    ret = str;

    str += strcspn(str, delim);

    if (*str)
    {
        *str++ = '\0';
    }

    *nextp = str;

    return ret;
}
#endif
