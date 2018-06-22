/*
 * readsolver.c
 *
 * functions to load a solver into memory
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <libxml/xpath.h>

#include "solver.h"
#include "readsolver.h"
#include "laminar.h"
#include "kE.h"
#include "readfile.h"

/* list of one dimensional solver properties */
const char *solver_properties_double[] = { "nu", "rho", "t", "delt", "writet", "endt", 
                                           "autot", "abstol", "reltol", "end" };

int read_solver_xml(struct solver_data *solver, char *filename) {
  xmlXPathContext *xpathCtx;
  xmlDoc *doc;
  double vector[3];
  char buf[256];
  int i;

  xmlInitParser();
  LIBXML_TEST_VERSION
 
  doc = xmlParseFile(filename);
  if(doc == NULL) {
    printf("Could not open %s\n",filename);
    return 1;
  }

  xpathCtx = xmlXPathNewContext(doc);
  if(doc == NULL) {
    printf("Could not create xpath context\n");
    return 1;
  }

  i = 0;
  while(strcmp("end", solver_properties_double[i]) != 0) {
    sprintf(buf, "/Case/Solver/Methods/%s", solver_properties_double[i]);

    if(!read_xmlpath_double(&vector[0], buf, xpathCtx)) {
      printf("Could not evaluate %s\n", solver_properties_double[i]);
    }
    else {
      printf("Read %s: %lf\n", buf, vector[0]);
      solver_set_value(solver, solver_properties_double[i], 1, vector);
    }
 
    i++;
  }

  if(read_xmlpath_double_vector(vector, "/Case/Solver/Methods/Gravity", "x", "y", "z", xpathCtx)) {
    solver_set_value(solver, "gravity", 3, vector);
  }

  if(!read_xmlpath_str(buf, "/Case/Solver/Methods/turbulence", xpathCtx)) {
    printf("Could not evaluate turbulence");
  }
  else {
    printf("Read /Case/Solver/Methods/turbulence: %s\n", buf);

    if(!strcmp(buf, "kE")) {
      vector[0] = 1;
    }
    else {
      vector[0] = 0;
    }

    solver_set_value(solver, "turbulence", 1, vector);
  }


  return 0;
}


int read_initial_xml(struct solver_data *solver, char *filename) {
  xmlXPathContext *xpathCtx;
  xmlDoc *doc;
  double vector[3];
  char buf[256];
  int i;

  xmlInitParser();
  LIBXML_TEST_VERSION
 
  doc = xmlParseFile(filename);
  if(doc == NULL) {
    printf("Could not open %s\n",filename);
    return 1;
  }

  xpathCtx = xmlXPathNewContext(doc);
  if(doc == NULL) {
    printf("Could not create xpath context\n");
    return 1;
  }

  if(read_xmlpath_double_vector(vector, "/Case/Solver/Initial/Velocity", "u", "v", "w", xpathCtx)) {
    solver_store_initial(solver, "velocity", 3, vector);
  }

  if(!read_xmlpath_double(&vector[0], "/Case/Solver/Initial/vof_height", xpathCtx)) {
    printf("Could not evaluate /Case/Solver/Initial/vof_height\n");
    vector[0] = 0;
  }
  else {
    printf("Read /Case/Solver/Initial/vof_height: %lf\n", vector[0]);
    solver_store_initial(solver, "vof_height", 1, vector);
  }

  if(!read_xmlpath_double(&vector[0], "/Case/Solver/Initial/hydrostatic", xpathCtx)) {
    printf("Could not evaluate /Case/Solver/Initial/hydrostatic\n");
    vector[0] = 0;
  }
  else {
    printf("Read /Case/Solver/Initial/hydrostatic: %lf\n", vector[0]);
    solver_store_initial(solver, "hydrostatic", 1, vector);
  }

  /* read inside points */
  i = 1;
  while(1) {
    sprintf(buf, "/Case/Solver/Initial/Inside/Point%d",i);
    
    if(read_xmlpath_double_vector(vector, buf, "x", "y", "z", xpathCtx)) {
      solver_store_initial(solver, "inside", 3, vector);
    }
    else break;

    i++;
  }

  if(!read_xmlpath_double(&vector[0], "/Case/Solver/Initial/kE_k", xpathCtx)) {
    printf("Could not evaluate /Case/Solver/Initial/kE_k\n");
    vector[0] = 0;
  }
  else {
    printf("Read /Case/Solver/Initial/kE_k: %lf\n", vector[0]);
    solver_store_initial(solver, "kE_k", 1, vector);
  }

  solver_store_initial(solver, "end", 0, NULL);

  return 0;
}

int read_solver(struct solver_data *solver, char *filename)
{
  /* code  to read solver variables
   * this includes delt, endt, gx, etc.
   *
   */
  FILE *fp;

  char text[1024];
  char args[4][256];
  double vector[3];

  char *s;
  int i, nargs;

  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to read_solver\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in read_solver\n",filename);
    return(1);
  }

  while(!feof(fp)) {
    if(!fgets(text, sizeof(text), fp))
    { 
      if(!feof(fp)) {
        printf("error: fgets in read_solver\n");
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
      printf("read_solver: read %d args: %s %lf %lf %lf\n", nargs, args[0], 
             vector[0], vector[1], vector[2]);
    #endif

    if (args[0][0] != '#' && args[0][0] != 0) 
      solver_set_value(solver, args[0], i, vector); 
    args[0][0] = 0;
  }

  fclose(fp);

  return 0;
}

int read_initial(struct solver_data *solver, char *filename)
{
  /* code  to read solver variables
   * this includes delt, endt, gx, etc.
   *
   */
  FILE *fp;

  char text[1024];
  char args[4][256];
  double vector[3];

  char *s;
  int i, nargs;

  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to read_initial\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in read_initial\n",filename);
    return(1);
  }

  while(!feof(fp)) {
    if(!fgets(text, sizeof(text), fp))
    { 
      if(!feof(fp)) {
        printf("error: fgets in read_initial\n");
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
      printf("read_initial: read %d args: %s %lf %lf %lf\n", nargs, args[0], 
             vector[0], vector[1], vector[2]);
    #endif

    if (args[0][0] != '#' && args[0][0] != 0) 
      solver_store_initial(solver, args[0], i, vector); 
    args[0][0] = 0;
  }

  solver_store_initial(solver, "end", 0, NULL);

  fclose(fp);

  return 0;
}


