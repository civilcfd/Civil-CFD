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

int read_solver_xml(struct solver_data *solver, char *filename) {
  xmlXPathContext *xpathCtx;
  xmlXPathObject  *xpathObj;
  xmlNode *node;
  xmlDoc *doc;
  double vector[3];

  xmlInitParser();
  LIBXML_TEST_VERSION
 
  doc = xmlParseFile(filename);
  if(doc === NULL) {
    printf("Could not open %s\n",filename);
    return 1;
  }

  xpathCtx = xmlXPathNewContext(doc);
  if(doc == NULL) {
    printf("Could not create xpath context\n");
    return 1;
  }

  xpathObj = xmlXPathEvalExpression( (xmlChar*) "/Case/Solver/Methods/nu", xpathCtx);
  if(xpathObj == NULL) {
    printf("Could not evaluate nu\n");
  }
  else {
    node = xpathObj->nodesetval->nodeTab[0];
    vector[0] = atof(xmlNodeGetContent(node));

    printf("nu: %lf\n",vector[0]);

  }


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


