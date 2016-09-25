/* readfile.h
 * 
 * header file for file i/o 
 */

#ifndef _READFILE_H
#define _READFILE_H

#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include <libxml/xpath.h>

#include "stl.h"
#include "mesh.h"

int string_index(char **str, char *str2);

int read_xmlpath_double_vector(double *vector, char *path, char *x, char *y, char *z, xmlXPathContext *xpathCtx);
int read_xmlpath_double(double *x, char *path, xmlXPathContext *xpathCtx);
int read_xmlpath_int(int *x, char *path, xmlXPathContext *xpathCtx);
int read_xmlpath_str(char *x, char *path, xmlXPathContext *xpathCtx);

int read_mesh(struct mesh_data *mesh, char *filename);
int read_mesh_xml(struct mesh_data *mesh, char *filename);

int write_mesh(struct mesh_data *mesh, char *filename);
int write_mesh_xml(struct mesh_data *mesh, xmlTextWriterPtr writer);

int read_stl(struct stl_data *stl, char *filename, double *limits);
int read_stl_ascii(struct stl_data *stl, char *filename, double *limits);
int read_stl_binary(struct stl_data *stl, char *filename, double *limits);

char *trimwhitespace(char *str);

int read_args(char *text, int nargs, char (*args)[256]);

#ifdef _WIN32
char* strtok_r(
    char *str, 
    const char *delim, 
    char **nextp);
#endif

#endif
