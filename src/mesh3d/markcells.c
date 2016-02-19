/* markcells.c
 *
 * creates a static structure to mark which mesh cells should be checked
 * for intersections with the stl file
 */
#include <stdio.h>
#include <math.h>
#include <omp.h>


#include "intersections.h"
#include "stl.h"
#include "mesh.h"
#include "vector_macros.h"

#include "markcells.h"

static int *marked_cells;

int markcells_initialize(struct mesh_data *mesh, 
                         struct stl_data *stl) {

	long int size, i, j, k, n, x, count;
  double v_1[3], v_2[3], v_3[3], norm[3], p[3], dist, v[3], s[3], dotp, min_dist;
  
	#ifndef __MINGW32__
  const double emf = 0.000001;
	#else
	const double emf = 0.0001;
	#endif
		
	size = mesh->imax * mesh->jmax * mesh->kmax;
	min_dist = sqrt(pow(mesh->delx,2) + pow(mesh->dely,2) + pow(mesh->delz,2)) + emf;
	count = 0;
	
	if(size < 1) return 1;
	
	marked_cells = malloc(sizeof(int) * size);
	
	if(marked_cells == NULL) {
    printf("error: memory could not be allocated for P in mesh_init_complete\n");
    return (1);	
	}
	
	for(i=0; i<size; i++) {
		marked_cells[i] = 0;
	}
	
#pragma omp parallel for shared (marked_cells, count) private(i, j, k, n, v_1, v_2, v_3, \
							 dotp, norm, p, dist, v, s, x) collapse(3) 
  for(i=0; i < mesh->imax; i++) {
    for(j=0; j< mesh->jmax; j++) {
      for(k=0; k < mesh->kmax; k++) {
      
					for(n=0; n < stl->facets; n++) {
					
						for(x=0; x<3; x++) {
							v_1[x] = stl->v_1[n][x];
							v_2[x] = stl->v_2[n][x];
							v_3[x] = stl->v_3[n][x];
							norm[x]= stl->normal[n][x];
						}
						
						p[0] = mesh->origin[0] + mesh->delx * i + mesh->delx/2;		
						p[1] = mesh->origin[1] + mesh->dely * j + mesh->dely/2;		
						p[2] = mesh->origin[2] + mesh->delz * k + mesh->delz/2;	
						
						vector_subtract(v, p, v_1);					
						dotp = inner_product(v, norm);
						dotp = dotp * -1.0;
						vector_multiply(s, norm, dotp);
						dist = vector_magnitude(s);
						/* vector_add(v, s, p); */
						
						/* ok s is the vector from the point to the plane described by the triangle */
						/* and dist is the minimum distance from the point to that plane 
							 need to fix such that dist is min distance to triangle v_1, v_2, v_3 */
						
						if(dist < min_dist) {
							count++;
							marked_cells[mesh_index(mesh,i,j,k)] = 1;
							break;
						}
					}      	
      
      }
    }
  }
	
	printf("Total number of marked cells: %ld / %ld \n", count, size);
	
	return 0;
}


