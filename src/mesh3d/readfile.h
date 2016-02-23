/* readfile.h
 * 
 * header file for file i/o 
 */

#ifndef _READFILE_H
#define _READFILE_H

#include "stl.h"
#include "mesh.h"

int read_mesh(struct mesh_data *mesh, char *filename);

int write_mesh(struct mesh_data *mesh, char *filename);

int read_stl(struct stl_data *stl, char *filename, double *limits);

char *trimwhitespace(char *str);

int read_args(char *text, int nargs, char (*args)[256]);

#ifdef __MINGW32__
char* strtok_r(
    char *str, 
    const char *delim, 
    char **nextp);
#endif

#endif
