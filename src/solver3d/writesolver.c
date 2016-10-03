/*
 * writesolver.c
 *
 * functions to write a solver into a file
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>

#include "vof_mpi.h"
#include "writesolver.h"
#include "solver.h"
#include "readsolver.h"
#include "laminar.h"
#include "kE.h"
#include "../mesh3d/readfile.h"


#define MY_ENCODING "UTF-8"
int write_solver_xml(struct solver_data *solver, char *filename) {
  
  int rc;
  xmlTextWriterPtr writer;
  xmlDocPtr doc;
  xmlNodePtr node;
  char buf[256];
  
  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to write_solver\n");
    return(1);
  }


  /* Create a new XML DOM tree, to which the XML document will be
    * written */
  doc = xmlNewDoc(BAD_CAST XML_DEFAULT_VERSION);
  if (doc == NULL) {
    printf
        ("testXmlwriterTree: Error creating the xml document tree\n");
    return 1;
  }

  /* Create a new XML node, to which the XML document will be
    * appended */
  node = xmlNewDocNode(doc, NULL, BAD_CAST "Case", NULL);
  if (node == NULL) {
    printf("testXmlwriterTree: Error creating the xml node\n");
    return 1;
  }

  /* Make ELEMENT the root node of the tree */
  xmlDocSetRootElement(doc, node);

  /* Create a new XmlWriter for DOM tree, with no compression. */
  writer = xmlNewTextWriterTree(doc, node, 0);
  if (writer == NULL) {
    printf("testXmlwriterTree: Error creating the xml writer\n");
    return 1;
  }
  xmlTextWriterSetIndent(writer, 1);

  /* Start the document with the xml default for the version,
    * encoding ISO 8859-1 and the default for the standalone
    * declaration. */
  rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
  if (rc < 0) {
      printf("testXmlwriterTree: Error at xmlTextWriterStartDocument\n");
      return 1;
  }

  rc = xmlTextWriterStartElement(writer, BAD_CAST "Solver");
  if (rc < 0) {
      printf("testXmlwriterTree: Error at xmlTextWriterStartElement\n");
      return 1;
  }


  xmlTextWriterStartElement(writer, BAD_CAST "Methods");

	if(solver->pressure == vof_pressure_gmres) sprintf(buf,"gmres");
	else if(solver->pressure == vof_pressure_gmres_mpi) sprintf(buf,"gmres_mpi");

  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "implicit", "%s", buf);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "nu", "%e", solver->nu);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "rho", "%e", solver->rho);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "t", "%e", solver->t);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "delt", "%e", solver->delt);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "writet", "%e", solver->writet);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "endt", "%e", solver->endt);

  rc = xmlTextWriterStartElement(writer, BAD_CAST "Gravity");
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "x", "%e", solver->gx);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "y", "%e", solver->gy);
  rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "z", "%e", solver->gz);
  rc = xmlTextWriterEndElement(writer);

  if(solver->deltcal == NULL)
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "autot", "%d", 0);
  else
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "autot", "%d", 1);
  
  if(kE_check(solver)) 
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "turbulence", "%s", "kE");
  else 
    rc = xmlTextWriterWriteFormatElement(writer, BAD_CAST "turbulence", "%s", "laminar");

  xmlTextWriterEndElement(writer);

  write_mesh_xml(solver->mesh, writer);
  write_initial_xml(solver, writer);
  if(kE_check(solver)) kE_write_xml(writer);

  xmlTextWriterEndElement(writer);
  xmlTextWriterEndDocument(writer);
  xmlFreeTextWriter(writer);
  xmlSaveFileEnc(filename, doc, MY_ENCODING);
  xmlFreeDoc(doc);

  return 0;
}


int write_initial_xml(struct solver_data *solver, xmlTextWriterPtr writer) {
  double velocity[3], inside[3];
  double value;
  int n = 1;
  char buf[256]; 
 
  xmlTextWriterStartElement(writer, BAD_CAST "Initial");

  solver_get_initial_vector(solver, "velocity", velocity);
  xmlTextWriterStartElement(writer, BAD_CAST "Velocity");
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "u", "%e", velocity[0]);
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "v", "%e", velocity[1]);
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "w", "%e", velocity[2]);
  xmlTextWriterEndElement(writer);

  value = solver_get_initial_scalar(solver, "vof_height");
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "vof_height", "%e", value);  
  
  value = solver_get_initial_scalar(solver, "hydrostatic");    
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "hydrostatic", "%e", value); 

  while(solver_get_initial_vector(solver, "inside", inside) != 1); // make sure we start at the first instance
  xmlTextWriterStartElement(writer, BAD_CAST "Inside");
  while(solver_get_initial_vector(solver, "inside", inside) != 1) {
    sprintf(buf, "Point%d", n); n++;
    xmlTextWriterStartElement(writer, BAD_CAST buf);
    xmlTextWriterWriteFormatElement(writer, BAD_CAST "x", "%e", inside[0]);
    xmlTextWriterWriteFormatElement(writer, BAD_CAST "y", "%e", inside[1]);
    xmlTextWriterWriteFormatElement(writer, BAD_CAST "z", "%e", inside[2]);
    xmlTextWriterEndElement(writer);
  }
  xmlTextWriterEndElement(writer);


  value = solver_get_initial_scalar(solver, "kE_k");     
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "kE_k", "%e", value);

  xmlTextWriterEndElement(writer);

  return 0;
}

int write_solver(struct solver_data *solver, char *filename) {
  FILE *fp;
  
  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to write_solver\n");
    return(1);
  }

  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: cannot open %s in write_solver\n",filename);
    return(1);
  }

	if(solver->pressure == vof_pressure_gmres) fprintf(fp,"gmres 1\n");
	else if(solver->pressure == vof_pressure_gmres_mpi) fprintf(fp,"gmres_mpi 1\n");

  fprintf(fp,"gravity %e %e %e\n",solver->gx, solver->gy, solver->gz);
  fprintf(fp,"nu %e\n",solver->nu);
  fprintf(fp,"rho %e\n",solver->rho);
  fprintf(fp,"t %e\n",solver->t);
  fprintf(fp,"delt %e\n",solver->delt);
  fprintf(fp,"writet %e\n",solver->writet);
  fprintf(fp,"endt %e\n",solver->endt);
  if(solver->deltcal == NULL) fprintf(fp,"autot 0\n");
  
  if(kE_check(solver)) fprintf(fp,"turbulence 1\n");
  else fprintf(fp,"turbulence 0\n");

  fclose(fp);

  return 0;
}

int write_initial(struct solver_data *solver, char *filename) {
  FILE *fp;
  double velocity[3], inside[3];
  double value;
  
  if(filename == NULL || solver == NULL) {
    printf("error: passed null arguments to write_initial\n");
    return(1);
  }

  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: cannot open %s in write_initial\n",filename);
    return(1);
  }

  solver_get_initial_vector(solver, "velocity", velocity);
  fprintf(fp,"velocity %e %e %e\n",velocity[0],velocity[1],velocity[2]);
  
  value = solver_get_initial_scalar(solver, "vof_height");  
  if(value>0) fprintf(fp,"vof_height %e\n",value);
  
  value = solver_get_initial_scalar(solver, "hydrostatic");    
  if(value > 0) fprintf(fp,"hydrostatic 1\n");

  while(solver_get_initial_vector(solver, "inside", inside) != 1); // make sure we start at the first instance
  while(solver_get_initial_vector(solver, "inside", inside) != 1) {
    fprintf(fp,"inside %e %e %e\n",inside[0],inside[1],inside[2]);
  }
  
  value = solver_get_initial_scalar(solver, "kE_k");     
  if(value > 0 && kE_check(solver)) fprintf(fp,"kE_k %e\n",value);

  fclose(fp);

  return 0;
}
