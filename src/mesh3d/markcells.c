/* markcells.c
 *
 * creates a static structure to mark which mesh cells should be checked
 * for intersections with the stl file
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
 

#include "intersections.h"
#include "stl.h"
#include "mesh.h"
#include "vector_macros.h"

#include "markcells.h"

static int *marked_cells;

int markcells_check(long int n) {
	return(marked_cells[n]);
}

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
    printf("error: memory could not be allocated for marked_cells in markcells_intialize\n");
    return (1);	
	}
	
	for(i=0; i<size; i++) {
		marked_cells[i] = 0;
	}
	
#pragma omp parallel for shared (marked_cells, count) private(i, j, k, n, v_1, v_2, v_3, \
							 dotp, norm, p, dist, v, s, x) schedule(static)
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

						dist = markcells_dist_tri_point(p, v_1, v_2, v_3);
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

int markcells_dist_tri_point(double *p, double *v1, double *v2, double *v3) {
	double a, b, c, d, e, f, det, s, t, sqrDistance, tmp0, tmp1, numer, denom, invDet;
	double BB[3], E0[3], E1[3], DD[3];

	/* rewrite triangle in normal form */
	/* B = TRI(1,:); */
	vector_copy(BB, v1);
	/* E0 = TRI(2,:)-B; */
	vector_subtract(E0, v2, BB);
	/* %E0 = E0/sqrt(sum(E0.^2)); %normalize vector
	E1 = TRI(3,:)-B;*/
	vector_subtract(E1, v3, BB);
	/* %E1 = E1/sqrt(sum(E1.^2)); %normalize vector*/


	/* D = B - P;
	a = dot(E0,E0);
	b = dot(E0,E1);
	c = dot(E1,E1);
	d = dot(E0,D);
	e = dot(E1,D);
	f = dot(D,D); */
	vector_subtract(DD, BB, p);
	a = inner_product(E0, E0);
	b = inner_product(E0, E1);
	c = inner_product(E1, E1);
	d = inner_product(E0, DD);
	e = inner_product(E1, DD);
	f = inner_product(DD, DD);

	/* det = a*c - b*b; % do we have to use abs here?
	s   = b*e - c*d;
	t   = b*d - a*e; */
	det = a * c - b * b;
	s = b * e - c * d;
	t = b * d - a * e;

	/* % Terible tree of conditionals to determine in which region of the diagram
	% shown above the projection of the point into the triangle-plane lies. */

	if(( s + t ) <= det) {
		if(s<0) {
			if(t<0) { /*region 4*/
				if(d<0) {
					t=0;
					if((d * -1.0) >= a) {
						s = 1;
						sqrDistance = a + 2*d + f;
					} else {
						s = -1.0 * d / a;
						sqrDistance = d * s + f;
					}
				} else {
					s = 0;
					if(e >= 0) {
						t = 0;
						sqrDistance = f;
					} else {
						if( -1.0 * e >= c) {
							t = 1;
							sqrDistance = c + 2*e + f;
						} else {
							t = -1.0 * e / c;
							sqrDistance = e * t + f;
						}
					}
				} /* end of region 4 */
		 	} else { /* region 3 */
				s = 0;
				if(e >= 0) {
					t = 0;
					sqrDistance =f;
				} else {
					if( -1.0 * e >= c ) {
						t = 1;
						sqrDistance = c + 2*e + f;
					} else {
						t = -1.0 * e / c;
						sqrDistance = e * t + f;
					}
				}
			} /* end of region 3 */	 
		} else {
			if(t < 0) { /* region 5 */
				t = 0;
				if(d >= 0) {
					s = 0;
					sqrDistance = f;
				} else {
					if(-1.0 * d >= a) {
						s = 1;
						sqrDistance = a + 2*d + f; 
					} else {
						s = -1.0 * d / a;
						sqrDistance = d * s + f;
					}
				} /* end region 5 */
			} else { /* region 0 */
				invDet = 1/det;
				s = s*invDet;
				t = t*invDet;
				sqrDistance = s*(a*s + b*t + 2*d) \
										+ t*(b*s + c*t + 2*e) + f;				
			}
		}
	} else {
		if(s < 0) {  /* region 2 */ 
			tmp0 = b + d;
			tmp1 = c + d;
			if(tmp1 > tmp0) {
				numer = tmp1 - tmp0;
				denom = a - 2*b + c;
				if(numer >= denom) {
					s = 1;
					t = 0;
					sqrDistance = a + 2*d + f;
				} else {
					s = numer / denom;
					t = 1 - s;
					sqrDistance = s * (a*s + b*t + 2*d) \
											+ t * (b*s + c*t + 2*e) + f;
				}
			} else {
				s = 0;
				if(tmp1 <= 0) {
					t = 1;
					sqrDistance = c + 2*e + f;
				} else {
					if(e >= 0) {
						t=0;
						sqrDistance = f;
					} else {
						t = -1.0 * e / c;
						sqrDistance = e * t + f;
					}
				}
			} /* end of region 2*/
		} else {
			if(t < 0) {
				/* region 6 */
				tmp0 = b + e;
				tmp1 = a + d;
				if(tmp1 > tmp0) {
					numer = tmp1 - tmp0;
					denom = a - 2*b + c;
					if(numer >= denom) {
						t = 1;
						s = 0;
						sqrDistance = c + 2*e + f;
					} else {
						t = numer / denom;
						s = 1 - t;
						sqrDistance = s * (a*s + b*t + 2*d) \
											  + t * (b*s + c*t + 2*e) + f;
					}
				} else {
					t = 0;
					if(tmp1 <= 0) {
						s = 1;
						sqrDistance = a + 2*d + f;
					} else {
						if(d >= 0) {
							s = 0;
							sqrDistance = f;
						} else {
							s = -1.0 * d/a;
							sqrDistance = d*s + f;
						}
					}
				}/*end of region 6 */
			} else { /* region 1 */
				numer = c + e - b - d;
				if(numer <= 0) {
					s = 0;
					t = 1;
					sqrDistance = c + 2*e + f;
				} else {
					denom = a - 2*b + c;
					if(numer >= denom) {
						s = 1;
						t = 0;
						sqrDistance = a + 2*d + f;
					} else {
						s = numer / denom;
						t = 1-s;
						sqrDistance = s * (a*s + b*t + 2*d) \
												+ t * (b*s + c*t + 2*e) + f;
					}
				} /*end of region 1 */
			}
		}
	}
/*			
			end
		end
	end

	% account for numerical round-off error
	if (sqrDistance < 0)
		sqrDistance = 0;
	end */
	
	if(sqrDistance < 0) {
		sqrDistance = 0;
	}
	
	return(sqrt(sqrDistance)); 
	
	/*

	dist = sqrt(sqrDistance);

	if nargout>1
		PP0 = B + s*E0 + t*E1;
	end*/
}




