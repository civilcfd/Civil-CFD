/* intersections.h
 *
 * header file with function definitions */

#ifndef _INTERSECTIONS_H
#define _INTERSECTIONS_H

#include "mesh.h"
#include "stl.h"

int intersect_area_fractions(struct mesh_data *mesh, 
                             struct stl_data *stl);

int check_intersect_tri(double *pt1, double *pt2, double *pt3, 
                        double *linept, double *vect,
                        double *pt_int);

int check_same_clock_dir(double *pt1, double *pt2, double *pt3, 
                         double *norm);

int check_coplanar(double *v0, double *v1, double *v2, double *p);

#endif
