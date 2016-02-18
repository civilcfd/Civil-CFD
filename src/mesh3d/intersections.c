/* intersections.c
 *
 * this code iterates through the mesh and calculates area fractions
 * then volume fractions are found
 */
#include <stdio.h>
#include <math.h>
#include <omp.h>


#include "intersections.h"
#include "stl.h"
#include "mesh.h"

#include "vector_macros.h"

double round_pres(double f,int pres);

int intersect_area_fractions(struct mesh_data *mesh, 
                             struct stl_data *stl) {

  long int i, j, k, n, x;
  int a, s, facing, flg;

  double pt_int[3], origin[3], r_o[3];
  double v_1[3], v_2[3], v_3[3];

  double x_af[3][4] = { { 0, 0, 0, 0},
                        { 0, 0, 0, 0},
                        { 0, 0, 0, 0} };
  
  int intersect[4] = { 0, 0, 0, 0 };
  
  double o_af[3][5][3] = { { { 1, 0, 0}, /* first this is ae */
                             { 1, 1, 0},
                             { 1, 1, 1},
                             { 1, 0, 1},
                             { 1, 0, 0} },
                           { { 0, 1, 0}, /* second, this is an */
                             { 0, 1, 1},
                             { 1, 1, 1},
                             { 1, 1, 0},
                             { 0, 1, 0} },
                           { { 0, 0, 1}, /* last, this is at */
                             { 1, 0, 1},
                             { 1, 1, 1},
                             { 0, 1, 1},
                             { 0, 0, 1} } };

  double del[3] = { mesh->delx, mesh->dely, mesh->delz };
  double o_n[3]; /* this is the normal vector of the line segment
                  * being checked for collisions */
  double sgn_n;  /* this is the sign of the normal vector*/
  int i_n;       /* this represents the axis of the normal vector */ 


  double f;
  int f0, f1;
  
  int abort = 0;

	#ifndef __MINGW32__
  const double emf = 0.000001;
	#else
	const double emf = 0.0001;
	#endif
  
  
#pragma omp parallel for shared (mesh) private(i, j, k, n, x, a, s, facing, flg, pt_int, \
								 origin, r_o, v_1, v_2, v_3, x_af, intersect, o_n, \
								 sgn_n, i_n, f, f0, f1) collapse(3)
  /* iterate through the mesh
   * we double calculate each line segment.  this could be optimized out
   * in the future, but would require more storage. */
  for(i=0; i < mesh->imax; i++) {
    for(j=0; j< mesh->jmax; j++) {
      for(k=0; k < mesh->kmax; k++) {

      for(a=0; a<4; a++) {
      	for(s=0; s<3; s++) {
							x_af[s][a] = 0; 
							intersect[a] = 0;
				}
			}
      
#pragma omp flush(abort)
      	if(!abort) {
				
					origin[0] = mesh->origin[0] + mesh->delx * i;
					origin[1] = mesh->origin[1] + mesh->dely * j;
					origin[2] = mesh->origin[2] + mesh->delz * k;

					/* need to populate ae, an and at, in that order */
					for(n=0; n < stl->facets; n++) {
					
					

						/* iterate over each side in this order: ae, an, at */
						for(s=0; s<3; s++) {

							for(a=0; a<4; a++) {
								/* r_o vector is the relative origin for this
								 * line segment. */ 
								r_o[0] = origin[0] + mesh->delx * o_af[s][a][0];
								r_o[1] = origin[1] + mesh->dely * o_af[s][a][1];
								r_o[2] = origin[2] + mesh->delz * o_af[s][a][2];

								/* o_n is the normal vector for this line segment */
								o_n[0] = o_af[s][a+1][0] - o_af[s][a][0];
								o_n[1] = o_af[s][a+1][1] - o_af[s][a][1];
								o_n[2] = o_af[s][a+1][2] - o_af[s][a][2];

								/* sgn_n is the sign of the vector.  this is important
								 * as the calculations change for a different sign */
								sgn_n = o_n[0] + o_n[1] + o_n[2];
								/* i_n is the axis of the normal vector */
								i_n = fabs(o_n[0]) * 0 + fabs(o_n[1]) * 1 + fabs(o_n[2]) * 2;

								if(fabs(sgn_n)>1 || i_n > 2) { 
									printf("error: normal vector not orthoganol\n");
									abort = 1;
									#pragma omp flush (abort)
									continue;
								}
							
								for(x=0; x<3; x++) {
									v_1[x] = stl->v_1[n][x];
									v_2[x] = stl->v_2[n][x];
									v_3[x] = stl->v_3[n][x];
								}
								/*
								if(i==7 && j==1 && k==3 && a==2 && s==0) {
									printf("7 1 3: a==2, r_o %lf %lf %lf o_n %lf %lf %lf\n",
												 r_o[0], r_o[1], r_o[2], o_n[0], o_n[1], o_n[2]);

									printf("7 1 3: v_1 %lf %lf %lf\n",v_1[0], v_1[1], v_1[2]);
									printf("7 1 3: v_1 %lf %lf %lf\n",v_2[0], v_2[1], v_2[2]);
									printf("7 1 3: v_1 %lf %lf %lf\n",v_3[0], v_3[1], v_3[2]);

								}*/

								/* check for intersection using the moller_trubore
								 * algorithm (google for information)
								 * standard graphics algorithm to find intersection of triangle
								 * and line segment
								 *
								 * inputs are  { r_o, o_n: origin and vector of line segment
								 *             { stl->v_x: 3 points, in order, for the triangle
								 * outputs are { facing  : true if triangle normal dot product with
								 *                         o_n is negative
								 *               pt_int  : cartesian coordinate of intersecting point */
								if((facing = moller_trumbore(r_o, o_n, v_1, 
												 v_2, v_3, pt_int)) != 0)  {
								
									/*if(i==7 && j ==1 && (k ==3 || k == 0)) 

										printf("Raw intersection %ld %ld %ld: %.20lf r_o %.10lf i_n %d s %d a %d\n", i,j,k,
										pt_int[i_n], r_o[i_n], i_n, s, a);*/

									f = (pt_int[i_n] - r_o[i_n])
											/ del[i_n];
									f *= sgn_n;
									if(f > -1.0 * emf && f < (1+emf)) {
										if(f<emf) f = emf;
										if(f>1-emf) f = 1-emf;
										x_af[s][a] = f * facing * -1;
									}
							
								}
							}

						} 
				
					}

					for(s=0; s<3; s++) {

	#ifdef DEBUG
						for(a=0; a<4; a++) {
							if(x_af[s][a] !=0) 
								printf("Intersection at %ld %ld %ld, face %d, side %d, fractional position %lf\n",
											 i, j, k, s, a, x_af[s][a]);
						}
	#endif
						/* now calculate ae/an/at */
						/* first check for multiple intersections */

						f = 0;
						flg = 0;
						for(a=0; a<4; a++) {
							if(x_af[s][a] != 0) {
							
								intersect[flg] = a;
								flg++;
							}

						}
						if(flg == 1) {
							printf("warning: single intersection in cell: %ld %ld %ld\n", i,j,k);
							printf("single intersections are ignored\n");
							f = 1;
						}
						if(flg >2 ) {
							/* first check if these intersections are at a shared vertex */
							/* make correction to this special case if necessary 
							 * the correction is to use the next vertex in the series. 
							 * the area is then calculated using the adjacent intersection case
							 * consider changing to the opposite intersection case! */
							x = 0;
							while (x<flg)
							{

								if(fabs(x_af[s][intersect[x]]) <= emf && 
									 fabs(x_af[s][intersect[(x-1 < 0 ? flg-1 : x-1)]]) >= 1-emf) {
									if(x+1 < flg) {
										/* intersect[0] = intersect[x];
										intersect[1] = intersect[x+1]; */
										f0 = intersect[x];
										f1 = intersect[x+1];

										/* in some cases the above code can create mismatched
										 * sides.  we need a check and correction */
										if(x_af[s][f0] * x_af[s][f1] > 0) {
											switch(x) {
											case 0:
												intersect[0] = intersect[1];
												intersect[1] = intersect[2];
												break;
											case 1:
												intersect[0] = intersect[0];
												intersect[1] = intersect[2];
												break;
											case 2:
												intersect[0] = intersect[1];
												intersect[1] = intersect[3];
												break; 
											}
										}
										/* no correction required */
										else {
											intersect[0] = f0;
											intersect[1] = f1;
										}

										flg=2;
										x=5;
										break;
									}
									else if(x-2 > -1) {
										f0 = intersect[x-2];
										f1 = intersect[x];
									
										/* again, correct for mismatching if necessary */
										if(x_af[s][f0] * x_af[s][f1] > 0) {
											switch(x) {
											case 2:
												intersect[0] = intersect[0]; /* for readability */
												intersect[1] = intersect[1];
												break;
											case 3:
												printf("error: cannot correct mismatch in cell %ld %ld %ld.  exiting.\n", i, j, k);
												abort = 1;
												#pragma omp flush (abort)
												break;
											}
										}
										else {
											intersect[0] = f0;
											intersect[1] = f1;
										}
										flg=2;
										x=5;
										break;
									}
								}
								x++;
							}
							/* second, another special case exists where one intersection is 
							 * at the end of the side (i.e. fractional position is 1-emf or
							 * emf ).  in this case, we should ignore this edge */
							if(x!=5) {
								x = 0;
								while (x<flg)
								{
									if(fabs(x_af[s][intersect[x]]) >= (1 - emf) ||
										 fabs(x_af[s][intersect[x]]) <= emf) {
										switch(x) {
										case 0:
											intersect[x]   = intersect[x+1];
											intersect[x+1] = intersect[x+2];
											break;
										case 1:
											intersect[x]   = intersect[x+1];
											break;
										}
										flg=2;
										x=5;
										break;
									}
									/* else
									**	printf("fabs: %e %d\n",fabs(x_af[s][intersect[x]]),  fabs(x_af[s][intersect[x]]) >= (1 - emf) );
									**/
									x++;
								}
							}
							/*if(flg >2  &&
								 (x_af[s][intersect[1]]-x_af[s][intersect[0]] <= -1 + 2.0 * emf || 
									x_af[s][intersect[1]]-x_af[s][intersect[0]] >= 1 - 2.0 * emf)) {
								 * special case of line crossing two opposing vertices *
								 intersect[1] = intersect[2];
								 flg = 2;
							}*/
							
							
							/* there is another special case where three adjacent sides have a total of 4 intersections, 
							and the side in the middle of the three sides has 2 intersections.  in this case, we can just
							delete the 2 intersections in the middle and make it a 2 intersection case */

							if(x!=5 && flg==4) {
								for(x = 0; x<3; x++) {
									if(intersect[x] == intersect[x+1] &&
											x_af[s][intersect[x]] > 0 &&
											x_af[s][intersect[x+1]] < 0) {
										switch(x) {
										case 0:
											intersect[0] = intersect[2];
											intersect[1] = intersect[3];
											flg=2;
											x=5;
											break;
										case 1:
											intersect[1] = intersect[3];
											flg=2;
											x=5;
											break;
										case 2:
											flg=2;
											x=5;
											break;
										}
										break;
									}
								}
							}
							
							if(x!=5) {
								printf("error: multiple intersections in cell: %ld %ld %ld\n", i,j,k);
								abort = 1;
								#pragma omp flush (abort)

							}
						}

					
						/* this is the adjacent intersection case */
						if((intersect[1] - intersect[0]) == 1 && flg==2) {
							f = ((1- fabs(x_af[s][intersect[0]])) * fabs(x_af[s][intersect[1]]))/2;
							printf("%ld %ld %ld: f %lf intersect[0] %d intersect[1] %d\n", i, j, k, 
											f, intersect[0], intersect[1]);
							if(x_af[s][intersect[0]] <= emf && x_af[s][intersect[1]] > 0)
								f = 1-f;
							if(f > emf && f < 1-emf &&
								 x_af[s][intersect[0]] * x_af[s][intersect[1]] > 0) {
								printf("error: triangle orientation mismatch in cell %ld %ld %ld\n", i,j,k);
								abort = 1;
								#pragma omp flush (abort)
							}

						}
						/* in this case, the intersections are on two opposing sides */
						else if((intersect[1] - intersect[0])  == 2 && flg==2) {
							f = (fabs(x_af[s][intersect[0]]) + (1-fabs(x_af[s][intersect[1]])))/2;
							if(x_af[s][intersect[0]] > 0 && x_af[s][intersect[1]] < 0)
								f = 1 - f; 
						
							else if(!(x_af[s][intersect[0]] < 0 && x_af[s][intersect[1]] > 0)) {
								printf("error: direction mismatch in cel: %ld %ld %ld\n", i, j, k);
								abort = 1;
								#pragma omp flush (abort)
							}

							if(!((intersect[0] == 0 && intersect[1] == 2) ||
									 (intersect[0] == 1 && intersect[1] == 3))) { 
								printf("error: intersection point mismatch in cell: %ld %ld %ld\n", i, j, k);
								printf("intersection points are: %d %d\n", intersect[0], intersect[1]);
								abort = 1;
								#pragma omp flush (abort)
							}

						}
						/* again, adjacent intersection, but the intersections are side
						 * 3 and side 0 */
						else if((intersect[1] - intersect[0]) == 3 && flg==2) {
							f = (fabs(x_af[s][intersect[0]]) * (1-fabs(x_af[s][intersect[1]])))/2;
							printf("%ld %ld %ld: f %lf intersect[0] %d intersect[1] %d\n", i, j, k, 
											f, intersect[0], intersect[1]);
							if(x_af[s][intersect[0]] > 0 && x_af[s][intersect[1]] <= emf) 
								f = 1-f;
							if(f > emf && f < 1-emf && 
								 x_af[s][intersect[0]] * x_af[s][intersect[1]] > 0) {
								printf("error: triangle orientation mismatch in cell %ld %ld %ld\n", i,j,k);
								abort = 1;
								#pragma omp flush (abort)

							}

						}
				
						switch (s) {
						case 0:
							mesh->ae[mesh_index(mesh, i, j, k)] = f;
							break;
						case 1:
							mesh->an[mesh_index(mesh, i, j, k)] = f;
							break;
						case 2:
							mesh->at[mesh_index(mesh, i, j, k)] = f;
							break;
						}

		#ifdef DEBUG
						if(f > 0 && f < 1) {
							printf("Area fraction in cell: %ld %ld %ld equals %lf | side: %d\n",
										 i,j,k,f,s);
							printf("Intersect[0] Intersect[1] flg: %d %d %d\n",intersect[0],intersect[1],flg);
						}
		#endif

						/* now clear the intersections once they have been used to
						 * calculate ae/an/at */
						for(a=0; a<4; a++) {
							x_af[s][a] = 0; 
							intersect[a] = 0;
						}
					}
				}
			}
		}
	}
	
	return(abort);

}

#define SAME_CLOCKNESS  1
#define DIFF_CLOCKNESS  0

int check_same_clock_dir(double *pt1, double *pt2, double *pt3, 
                         double *norm)
{  
   double testi, testj, testk;
   double dotprod;

   /* normal of triangle 
    * calculated based on cross product of two sides */
   testi = (((pt2[1] - pt1[1])*(pt3[2] - pt1[2])) - 
            ((pt3[1] - pt1[1])*(pt2[2] - pt1[2])));
   testj = (((pt2[2] - pt1[2])*(pt3[0] - pt1[0])) - 
            ((pt3[2] - pt1[2])*(pt2[0] - pt1[0])));
   testk = (((pt2[0] - pt1[0])*(pt3[1] - pt1[1])) - 
            ((pt3[0] - pt1[0])*(pt2[1] - pt1[1])));

   // Dot product with triangle normal
   dotprod = testi*norm[0] + testj*norm[1] + testk*norm[2];

   //answer
   if(dotprod < 0) return DIFF_CLOCKNESS;
   else return SAME_CLOCKNESS;
}

int check_intersect_tri(double *pt1, double *pt2, double *pt3, 
                        double *linept, double *vect,
                        double *pt_int)
{
   double V1x, V1y, V1z;
   double V2x, V2y, V2z;
   double norm[3];
   double dotprod;
   double t;

   // vector form triangle pt1 to pt2
   V1x = pt2[0] - pt1[0];
   V1y = pt2[1] - pt1[1];
   V1z = pt2[2] - pt1[2];

   /* vector form triangle pt2 to pt3 */
   V2x = pt3[0] - pt2[0];
   V2y = pt3[1] - pt2[1];
   V2z = pt3[2] - pt2[2];

   /* vector normal of triangle */
   norm[0] = V1y*V2z-V1z*V2y;
   norm[1] = V1z*V2x-V1x*V2z;
   norm[2] = V1x*V2y-V1y*V2x;

   /* dot product of normal and line's 
    * vector if zero line is parallel to triangle */
   dotprod = norm[0]*vect[0] + norm[1]*vect[1] + norm[2]*vect[2];

   if(dotprod > 0.000001)
   {
      /* Find point of intersect to triangle plane.
         find t to intersect point */
      t = -(norm[0]*(linept[0]-pt1[0]) + 
            norm[1]*(linept[1]-pt1[1]) + 
            norm[2]*(linept[2]-pt1[2])) /
           (norm[0]*vect[0] + norm[1]*vect[1] + norm[2]*vect[2]);

      /* if ds is neg line started past triangle so can't hit triangle. */
      if(t < 0) return 0;

      pt_int[0] = linept[0] + vect[0]*t;
      pt_int[1] = linept[1] + vect[1]*t;
      pt_int[2] = linept[2] + vect[2]*t;

      if(check_same_clock_dir(pt1, pt2, pt_int, norm) == SAME_CLOCKNESS)
      {
         if(check_same_clock_dir(pt2, pt3, pt_int, norm) == SAME_CLOCKNESS)
         {
            if(check_same_clock_dir(pt3, pt1, pt_int, norm) == SAME_CLOCKNESS)
            {
               /* answer in pt_int is insde triangle */
               return 1;
            }
         }
      }
   }
   return 0;
}

int moller_trumbore(double *r_o, double *r_d, double *v1,
                   double *v2, double *v3, double *pt) {
  double e1[3], e2[3], r[3], s[3], q[3];
  double f, u, v, t, a;
  int facing;      /* 1 for front facing
                    * -1 for rear facing
                    * 0 for no intersection */
  const double emf = 0.000001;
  const double em10 =    0.0000000001; /* accounts for fp errors more tightly */
  const double em10_n = -0.0000000001;

  vector_subtract(e2, v3, v1);
  vector_subtract(e1, v2, v1);
  cross_product(r, r_d, e2);

  vector_subtract(s, r_o, v1);

  a = inner_product(e1, r);
  f = 1 / a;

  cross_product(q, s, e1);

  u = inner_product(s, r);
  facing = 1;

  if(a > emf) { 
    if (u<em10_n || u>(a+em10)) return 0;

    v = inner_product(r_d, q);
    if (v<em10_n || (u+v)>(a+em10)) return 0;
  }
  else if(a < (emf * -1.0)) {
    facing = -1;

    if(u>em10 || u<(a-em10)) return 0;

    v = inner_product(r_d, q);

    if(v>em10 || (u+v)<(a-em10)) return 0;
  }
  else return 0; /* parallel case */

  t=f*inner_product(e2, q);

  u = u * f; v = v * f;

  pt[0] = r_o[0] + r_d[0] * t;
  pt[1] = r_o[1] + r_d[1] * t;
  pt[2] = r_o[2] + r_d[2] * t;

  return facing;
}

double round_pres(double x, int prec)
{
double power = 1.0;
int i;

if (prec > 0)
for (i = 0; i < prec; i++)
power *= 10.0;
else if (prec < 0)
for (i = 0; i < prec; i++)
power /= 10.0;

if (x > 0)
x = floor(x * power + 0.5) /power;
else if (x < 0)
x = ceil(x * power - 0.5) /power;

x = (double) ((float) x);

if (x == -0)
x = 0;

return x;
}

int check_coplanar(double *v0, double *v1, double *v2, double *p) {
  double a[3], b[3], n[3], mag_n, d, res;
  const double emf = 0.001;

  vector_subtract(a, v1, v0);
  vector_subtract(b, v2, v0);
  cross_product(n, a, b);

  mag_n = pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2);
  mag_n = sqrt(mag_n);
  n[0] /= mag_n;
  n[1] /= mag_n;
  n[2] /= mag_n;

  d = inner_product(n, v0);

  res = inner_product(n, p) - d;

  if(fabs(res) > emf) return 0;
  else return 1;
}
