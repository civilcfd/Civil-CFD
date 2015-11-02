#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <libqhull/libqhull.h>
#include <libqhull/mem.h>
#include <libqhull/qset.h>
#include "qh_interface.h"

double qhull_volume(double (*v_list)[3], int v_max) {
	
	double *points;
	int i, x;
	
	points = malloc(sizeof(double) * v_max * 3);
	
	for(i=0; i<v_max; i++) {
		for(x=0; x<3; x++) {
			points[i * 3 + x] = v_list[i][x];
		}
	}
	
	return qhull_get_volume(points, v_max);
}

double qhull_get_pent_volume(double *p0, double *p1, double *p2, double *p3,
                             double *p4, double *p5, double *p6, double *p7)
{
  double points[24];

  points[0] = p0[0];
  points[1] = p0[1];
  points[2] = p0[2];

  points[3] = p1[0];
  points[4] = p1[1];
  points[5] = p1[2];
  
  points[6] = p2[0];
  points[7] = p2[1];
  points[8] = p2[2];
  
  points[9] = p3[0];
  points[10] = p3[1];
  points[11] = p3[2];
  
  points[12] = p4[0];
  points[13] = p4[1];
  points[14] = p4[2];

  points[15] = p5[0];
  points[16] = p5[1];
  points[17] = p5[2];

  points[18] = p6[0];
  points[19] = p6[1];
  points[20] = p6[2];

  points[21] = p7[0];
  points[22] = p7[1];
  points[23] = p7[2];

  return qhull_get_volume(points, 8);
}

double qhull_get_volume(double *points, int numpoints) {
  double volume;
  int exitcode, dim=3;
  boolT ismalloc = 0;
  char flags[256] = "qhull s FA";

  /* qh_init_A(stdin, stdout, stderr, argc, argv);  * sets qh qhull_command */
  exitcode= qh_new_qhull(dim, numpoints, points, ismalloc, flags, NULL, stderr);
  /*setjmp(qh errexit); * simple statement for CRAY J916 */
  if (!exitcode)  
    volume = qh totvol;
  else
    volume = -1.0;

  qh NOerrexit = True;
  qh_freeqhull(!qh_ALL);

  return volume;
}

#ifdef TEST_INTERFACE

int main(int argc, char *argv[]) {

   coordT points[8*3] = { 0,   0.2,    0,
                         0,   0.042,  0,
                         0,   0.2,    0.2,
                         0,   0.147,  0.2,
                         0.1, 0.2,    0.2,
                         0.2, 0.2,    0,
                         0.2, 0.147,  0,
                         0.2, 0.2,    0.1 };
  
  printf("Total Volume: %lf\n", qhull_get_volume(points, 8));
  printf("If above number ~= 0.002202 then the interface is working\n");

  return 0;
}

#endif
