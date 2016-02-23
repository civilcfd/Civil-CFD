/* volfract.c
 *
 * calculate volume fractions for each cell
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "stl.h"
#include "mesh.h"
#include "intersections.h"
#include "volfract.h"
#include "qh_interface.h"

#define vector_add(a,b,c) \
  (a)[0] = (b)[0] + (c)[0]; \
  (a)[1] = (b)[1] + (c)[1]; \
  (a)[2] = (b)[2] + (c)[2];

#define vector_subtract(a,b,c) \
  (a)[0] = (b)[0] - (c)[0]; \
  (a)[1] = (b)[1] - (c)[1]; \
  (a)[2] = (b)[2] - (c)[2];

#define cross_product(a,b,c) \
  (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
  (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
  (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define inner_product(v,q) \
    ((v)[0] * (q)[0] + \
    (v)[1] * (q)[1] + \
    (v)[2] * (q)[2])


int volume_fractions(struct mesh_data *mesh,
                     struct stl_data *stl) {

  long int i, j, k;

  int flg;

  /* the methodology here is to look for primitives that match the 
   * shape of the plane intersected hex, and calculate volume based on
   * simple geometric formulae
   *
   * for more complicated, non-planar intersections, the qhull library
   * is used to calculate volume
   *
   * it may be better to use qhull for all volumes, but there could be
   * a significant performance hit
   *
   * i find a simple mesh with 1 million cells takes approx 10 mins with this
   * code on 2012 gamer hardware */
#pragma omp parallel for shared (mesh) private(i,j,k,flg) schedule(dynamic, 300)
	for(i=0; i < mesh->imax; i++) {
		for(j=0; j < mesh->jmax; j++) {
			for(k=0; k < mesh->kmax; k++) {
			
				if(!markcells_check(mesh_index(mesh,i,j,k))) continue;
			
				if((flg = tet_fraction(mesh, stl, i, j, k))==2) {
					if(!line_pent_fraction(mesh, stl, i, j, k))
						if(!face_hex_fraction(mesh, stl, i, j, k))
							if(!vertex_pent_fraction(mesh, stl, i, j, k))
								if(!vertex_hex_fraction(mesh, stl, i, j, k))
									if(!vertex_qhull_fraction(mesh, stl, i, j, k)) 
									{
										printf("warning: volume_fractions found intersections but no suitable primitives to calculate volume in cell %ld %ld %ld\n",i,j,k);
										flg=0;
									}

				/* TODO:
				 * Code should include area fraction corrections for some of the
				 * primitives.  For example, line_pent_fraction can correct
				 * the area fractions, which will be incorrectly approximated
				 * in "intersections.c" 
				 *
				 * It may make sense to calculate area fractions for all primitives
				 * and use it as a debug error check against intersections.c
				 * */

				}
				if(flg == 0) {
					mesh->fv[mesh_index(mesh, i, j, k)] = 0;
				}

			}
		}
	}

  return 0;
}

int vertex_qhull_fraction(struct mesh_data *mesh, struct stl_data *stl, long int i, long int j, long int k)
{
	double v_list[32][3];
	double pt[3], origin[3];
	int n, s, v_index, v, v_vertices;
	int flg, x, a, sgn, facing;
	double r_o[3], o_n[3], pt_int[3];
	double f;
	
	const double vertex_list[8][3] = { { 0, 0, 0 },
                                     { 0, 1, 0 },
                                     { 1, 1, 0 },
                                     { 1, 0, 0 },
                                     { 0, 0, 1 },
                                     { 1, 0, 1 },
                                     { 1, 1, 1 },
                                     { 0, 1, 1 } }; 

  const double normals[3][3] =     { { 1, 0, 0 },
                                     { 0, 1, 0 },
                                     { 0, 0, 1 } };

  const double del[3] = { mesh->delx, mesh->dely, mesh->delz };
	

	#ifndef __MINGW32__
  const double emf = 0.000001;
	#else
	const double emf = 0.0001;
	#endif

  origin[0] = mesh->origin[0] + mesh->delx * i;
  origin[1] = mesh->origin[1] + mesh->dely * j;
  origin[2] = mesh->origin[2] + mesh->delz * k;

  flg=0;
	v_index = 0;
	
	/* Step 1: Iterate through each vertex and check if it is inside the mesh */

	for(s=0; s<8; s++) { /* iterate through each vertex */

		pt[0] = origin[0] + vertex_list[s][0] * mesh->delx;
		pt[1] = origin[1] + vertex_list[s][1] * mesh->dely;
		pt[2] = origin[2] + vertex_list[s][2] * mesh->delz;

		if(stl_check_normals_point(mesh, stl, pt)) { /* this point is inside the mesh */
			v_list[v_index][0] = pt[0] - origin[0];
			v_list[v_index][1] = pt[1] - origin[1];
			v_list[v_index][2] = pt[2] - origin[2];
			v_index++;
		}
 
  }
	
	if(v_index == 8) return 0; /* all points are inside the mesh */
	v_vertices = v_index;
	
	/* Step 2: Now add any intersections to the vertex list */
	for(s=0; s<8; s++) { /* iterate through each vertex */

    for(x=0; x<3; x++) {
      r_o[x] = origin[x] + del[x] * vertex_list[s][x];
    }
		
    for(n=0; n < stl->facets; n++) {
  
      for(a=0; a<3; a++) { /* iterate through each axis */
        
        for(x=0; x<3; x++) {
          o_n[x] = normals[a][x] * 
               (vertex_list[s][x] > 0 ? -1.0 : 1.0);
        }

       
        sgn = o_n[a];
        if((facing = moller_trumbore(r_o, o_n, stl->v_1[n],
                                     stl->v_2[n], stl->v_3[n],
                                     pt_int)) != 0) {
          f = (pt_int[a] - r_o[a]) / del[a];
          f *= sgn;

          if(f > -1.0 * emf && f < (1+emf)) {
            
						flg = 0;
						for(v = 0; v < v_index; v++) {
							/* check for duplicates */
							for(x=0; x<3; x++) {
								if((pt_int[x] - origin[x]) > (v_list[v][x] - emf) && (pt_int[x]-origin[x]) < (v_list[v][x] + emf)) {
									flg++;

								}
							}
							if (flg<3) flg=0;
							else break;
							

						}
						
						if(flg < 3) { /* non-duplicate entry.  Add to list */
							v_list[v_index][0] = pt_int[0] - origin[0];
							v_list[v_index][1] = pt_int[1] - origin[1];
							v_list[v_index][2] = pt_int[2] - origin[2];
							v_index++;
							
							if(v_index >= 32) {
								printf("error: impossibly high number of intersections at cell %ld %ld %ld.  Skipping", i, j, k);
								return 0;

							}
						}

          }

        }

      }
    }
  }
	
	if(v_vertices == v_index) {
		printf("warning: vertices outside the mesh without intersections at cell %ld %ld %ld.", i, j, k);

	}
	
	f = qhull_volume(v_list, v_index);
	
	if(f == -1.0) {
		printf("error in qhull_get_volume for cell %ld %ld %ld. ignoring and will set volume to 0.0\n", i, j, k);
		f = 0;
	}
	
	f /= mesh->delx * mesh->dely * mesh->delz;
	mesh->fv[mesh_index(mesh,i,j,k)] = f;
	
#ifdef DEBUG
    printf("%ld %ld %ld: Found qhull volume, volume of %lf\n", i, j, k, 
      mesh->fv[mesh_index(mesh, i, j ,k)]);
#endif	

	return 1;
}

int vertex_pent_fraction(struct mesh_data *mesh, struct stl_data *stl,
                        long int i, long int j, long int k) {

  long int n, a, b, s;
  int x, y, sgn;
  int bits = 0;
  int facing, f_facing, marker;

  double pt_int[3], origin[3], r_o[3], o_n[3], dv[3], d_o[3];
  double r_v[3], vert[4][3], tet[3], v_x[3][3], v_list[6][3];
  double f, dx_dz, dy_dz, dx_dy, dz_dx, dy_dx, dz_dy, dx, dy, dz;
 
  /* use -2 as a marker so we can see which normal direction
   * does not have 2 intersections on the secondary vertex */
  double hex[3][2][3] = { { { -2, -2, -2 },
                            { -2, -2, -2 } },
                          { { -2, -2, -2 },
                            { -2, -2, -2 } },
                          { { -2, -2, -2 },
                            { -2, -2, -2 } } };

  const double vertex_list[8][3] = { { 0, 0, 0 },
                                     { 0, 1, 0 },
                                     { 1, 1, 0 },
                                     { 1, 0, 0 },
                                     { 0, 0, 1 },
                                     { 1, 0, 1 },
                                     { 1, 1, 1 },
                                     { 0, 1, 1 } }; 

  const double normals[3][3] =     { { 1, 0, 0 },
                                     { 0, 1, 0 },
                                     { 0, 0, 1 } };

  const double del[3] = { mesh->delx, mesh->dely, mesh->delz };

	#ifndef __MINGW32__
  const double emf = 0.000001;
	#else
	const double emf = 0.0001;
	#endif

  origin[0] = mesh->origin[0] + mesh->delx * i;
  origin[1] = mesh->origin[1] + mesh->dely * j;
  origin[2] = mesh->origin[2] + mesh->delz * k;

  marker = -1;

  for(s=0; s<8; s++) { /* iterate through each vertex */

    for(x=0; x<3; x++) {
      r_o[x] = origin[x] + del[x] * vertex_list[s][x];
    
      for(a=0; a<2; a++) {
        hex[x][a][0] = -2;
        hex[x][a][1] = -2;
        hex[x][a][2] = -2;
      }
    }

    bits=0;
      
    for(a=0; a<3; a++) { /* iterate through each secondary vertex */

      vert[a][0] = vertex_list[s][0];
      vert[a][1] = vertex_list[s][1];
      vert[a][2] = vertex_list[s][2];

      if(vertex_list[s][0] < emf) {
        r_v[0] = r_o[0] + normals[a][0] * del[0];
        vert[a][0] = vert[a][0] + normals[a][0];
      }
      else { 
        r_v[0] = r_o[0] - normals[a][0] * del[0];
        vert[a][0] = vert[a][1] - normals[a][0];
      }

      if(vertex_list[s][1] < emf) {
        r_v[1] = r_o[1] + normals[a][1] * del[1];
        vert[a][1] = vert[a][1] + normals[a][1];
      }
      else {
        r_v[1] = r_o[1] - normals[a][1] * del[1];
        vert[a][1] = vert[a][1] - normals[a][1];
      }

      if(vertex_list[s][2] < emf) {
        r_v[2] = r_o[2] + normals[a][2] * del[2];
        vert[a][2] = vert[a][2] + normals[a][2];
      }
      else {
        r_v[2] = r_o[2] - normals[a][2] * del[2];
        vert[a][2] = vert[a][2] - normals[a][2];
      }

      vector_subtract(d_o, r_v, origin);

      dv[0]  = fabs(r_o[0] - r_v[0]);
      dv[1]  = fabs(r_o[1] - r_v[1]);
      dv[2]  = fabs(r_o[2] - r_v[2]);

      /* iterate through each axis on each secondary vertex */
      x=0;
      for(b=0; b<3; b++) {  
           
        if(inner_product(dv, normals[b]) > emf)
          continue;

        sgn = inner_product(d_o, normals[b]) > emf ? -1 : 1; 

        o_n[0] = normals[b][0]; /* sgn;*/
        o_n[1] = normals[b][1]; /* sgn;*/
        o_n[2] = normals[b][2]; /* sgn;*/
        
        for(n=0; n < stl->facets; n++) {

          if((facing = moller_trumbore(r_v, o_n, stl->v_1[n],
                                       stl->v_2[n], stl->v_3[n],
                                       pt_int)) != 0) {
 
            f = (pt_int[0] - r_v[0])/del[0] +
                (pt_int[1] - r_v[1])/del[1] +
                (pt_int[2] - r_v[2])/del[2];
            f *= sgn;

            if(f > -1.0 * emf && f < (1+emf)) {
              hex[a][x][0] = /* pt_int[0];*/(pt_int[0] - r_v[0])/del[0]; 
              hex[a][x][1] = /* pt_int[1];*/(pt_int[1] - r_v[1])/del[1]; 
              hex[a][x][2] = /* pt_int[2];*/(pt_int[2] - r_v[2])/del[2];
              bits |= 1 << (a*2+x);
              
              f_facing = facing * sgn;
            }
       
          }

        }

        x++;
      }
    }
    if(count_bits(bits) == 4) {
    /* we know that there are 2 secondary verticies with intersections
     * in opposite axis
     * check to see if the axis direction without the secondary vertices
     * has an intersection.  if so, this is a pentagonal intersection
     */
      for(b=0; b<3; b++) {
      /* we left a marker in hex for this purpose */
        if(hex[b][0][0] == -2 && hex[b][0][1] == -2 && hex[b][0][2] == -2 &&
           hex[b][1][0] == -2 && hex[b][1][1] == -2 && hex[b][1][2] == -2) {
        /* set up normals and r_v to look for intersection */
          marker = b; /* store the marker value to prevent future look-ups */

          vector_subtract(d_o, r_o, origin);
 
          sgn = inner_product(d_o, normals[b]) > emf ? -1 : 1; 

          o_n[0] = normals[b][0]; /* sgn;*/
          o_n[1] = normals[b][1]; /* sgn;*/
          o_n[2] = normals[b][2]; /* sgn;*/
          
          for(n=0; n < stl->facets; n++) {

            if((facing = moller_trumbore(r_o, o_n, stl->v_1[n],
                                         stl->v_2[n], stl->v_3[n],
                                         pt_int)) != 0) {
   
              f = (pt_int[0] - r_o[0])/del[0] +
                  (pt_int[1] - r_o[1])/del[1] +
                  (pt_int[2] - r_o[2])/del[2];
              f *= sgn;

              if(f > -1.0 * emf && f < (1+emf)) {
                hex[b][0][0] = /* pt_int[0];*/(pt_int[0] - r_o[0])/del[0]; 
                hex[b][0][1] = /* pt_int[1];*/(pt_int[1] - r_o[1])/del[1]; 
                hex[b][0][2] = /* pt_int[2];*/(pt_int[2] - r_o[2])/del[2];
                
                f_facing = facing * sgn;
              }
              else marker = -1;
               
         
            }

          }


        }


      }
    }

    if(marker != -1) {
      vert[3][0] = vertex_list[s][0];
      vert[3][1] = vertex_list[s][1];
      vert[3][2] = vertex_list[s][2];
 
      break;
    }

    

  }

  if(marker != -1) {  
    /* we have a vertex centered pentagonal plane intersection
     * now we need to calc volume
     * we split the problem into 3 cases, based on which axis
     * does not contain secondary verticies with intersections on
     * opposing axis 
     * this is stored in the marker variable */

    switch(marker) {
    case 0: /* x axis */
      v_x[0][0] = vert[3][0] + hex[0][0][0];
      v_x[0][1] = vert[3][1];
      v_x[0][2] = vert[3][2];

      dy = 1;
      dx = fabs(hex[1][0][0] - hex[0][0][0]);
      dy_dx = dy/dx;

      if(vert[1][1] > emf) 
        v_x[1][1] = dy_dx * fabs(hex[1][0][0]) + vert[1][1];
      else
        v_x[1][1] = dy_dx * fabs(hex[1][0][0]) * -1.0;
      v_x[1][0] = vert[1][0];
      v_x[1][2] = vert[1][2];

      dz = 1;
      dx = fabs(hex[2][0][0] - hex[0][0][0]);
      dz_dx = dz/dx;

      if(vert[2][2] > emf)
        v_x[2][2] = dz_dx * fabs(hex[2][0][0]) + vert[2][2];
      else
        v_x[2][2] = dz_dx * fabs(hex[2][0][0]) * -1.0;
      v_x[2][0] = vert[2][0];
      v_x[2][1] = vert[2][1];
      
      break;
    case 1: /* y axis */
      dx = 1;
      dy = fabs(hex[0][0][1] - hex[1][0][1]);
      dx_dy = dx/dy;

      if(vert[0][0] > emf)
        v_x[0][0] = dx_dy * fabs(hex[0][0][1]) + vert[0][0];
      else
        v_x[0][0] = dx_dy * fabs(hex[0][0][1]) * -1.0;
      v_x[0][1] = vert[0][1];
      v_x[0][2] = vert[0][2];

      v_x[1][0] = vert[3][0];
      v_x[1][1] = vert[3][1] + hex[1][0][1];
      v_x[1][2] = vert[3][2];

      dz = 1;
      dy = fabs(hex[2][1][1] - hex[1][0][1]);
      dz_dy = dz/dy;

      if(vert[2][2] > emf)
        v_x[2][2] = dz_dy * fabs(hex[2][1][1]) + vert[2][2];
      else
        v_x[2][2] = dz_dy * fabs(hex[2][1][1]) * -1.0;
      v_x[2][0] = vert[2][0];
      v_x[2][1] = vert[2][1];

      break;
    case 2: /* z axis */
      dx = 1;
      dz = fabs(hex[0][1][2] - hex[2][0][2]);
      dx_dz = dx/dz;

      if(vert[0][0] > emf) 
        v_x[0][0] = dx_dz * fabs(hex[0][1][2]) + vert[0][0];
      else
        v_x[0][0] = dx_dz * fabs(hex[0][1][2]) * -1.0;
      v_x[0][1] = vert[0][1];
      v_x[0][2] = vert[0][2];

      dy = 1;
      dz = fabs(hex[1][1][2] - hex[2][0][2]);
      dy_dz = dy/dz;

      if(vert[1][1] > emf)
        v_x[1][1] = dy_dz * fabs(hex[1][1][2]) + vert[1][1];
      else
        v_x[1][1] = dy_dz * fabs(hex[1][1][2]) * -1.0;
      v_x[1][0] = vert[1][0];
      v_x[1][2] = vert[1][2];

      v_x[2][0] = vert[3][0];
      v_x[2][1] = vert[3][1];
      v_x[2][2] = vert[3][2] + hex[2][0][2];

      break;
    }
    
    y = 0;
    for(x = 0; x<3; x++) {
      if(x==marker) continue;
      vector_add(hex[x][0], hex[x][0], vert[x]);
      vector_add(hex[x][1], hex[x][1], vert[x]);

      tet[y] = fabs(tet_volume(v_x[x], hex[x][0], hex[x][1], vert[x]));

      /* build a list of vector to check if they are coplanar */
      v_list[2*y][0] = hex[x][0][0];
      v_list[2*y][1] = hex[x][0][1];
      v_list[2*y][2] = hex[x][0][2];
      v_list[2*y+1][0] = hex[x][1][0];
      v_list[2*y+1][1] = hex[x][1][1];
      v_list[2*y+1][2] = hex[x][1][2];

      y++;
    }

    tet[y] = fabs(tet_volume(v_x[0], v_x[1], v_x[2], vert[3]));

    if(check_coplanar(v_x[marker], v_list[0], v_list[1], v_list[2]) &&
       check_coplanar(v_x[marker], v_list[0], v_list[1], v_list[3]))
    /* since the intersections are co-planar, this is the simpler case where volume
     * can be calculated by breaking the geometry into tetrahedra */
      mesh->fv[mesh_index(mesh, i, j, k)] = 
        tet[2] - (tet[0] + tet[1]);

    else {
    /* calculate volume of a cut-plane intersection
     * this is the case where a pentagonal intersection is formed by multiple planes
     * meeting inside the solid, and needs to be handled uniquely */
    /* qhull will be used to calculate the volume based on a convex polyhedra
     * defined by the intersections and vertices */
      y=0;
      for(x=0; x<3; x++) {
        if(x==marker) continue;

        v_list[4+y][0] = vert[x][0];
        v_list[4+y][1] = vert[x][1];
        v_list[4+y][2] = vert[x][2];
        y++;
      }

      /* vertex list is as follows:
       * v_list[0 - 3] are the double intersections 90 degree to each other, centered around vertices
       * v_list[4 - 5] are the vertices with the double intersections above
       * vert[3] is the vertex with only a single intersection
       * v_x[2] is the intersection for the above vertex */
      f = qhull_get_pent_volume(v_list[0], v_list[1], v_list[2], v_list[3], v_list[4], v_list[5],
                              vert[3], v_x[marker]);

      if(f == -1.0) {
        printf("error in qhull_get_pent_volume for cell %ld %ld %ld. ignoring and will set volume to 0.0\n", i, j, k);
        f = 0.0;
      }
      mesh->fv[mesh_index(mesh, i, j, k)] = f;
    }

    if(f_facing < 0)
      mesh->fv[mesh_index(mesh, i, j, k)] = 1 -
        mesh->fv[mesh_index(mesh, i, j, k)];

#ifdef DEBUG
    printf("%ld %ld %ld: Found vertex_pent, volume of %lf\n", i, j, k, 
      mesh->fv[mesh_index(mesh, i, j ,k)]);
#endif

    return 1;
    
  }

  return 0;
}


int vertex_hex_fraction(struct mesh_data *mesh, struct stl_data *stl,
                        long int i, long int j, long int k) {

  long int n, a, b, s;
  int x, sgn;
  int bits = 0;
  int facing, f_facing;

  double pt_int[3], origin[3], r_o[3], o_n[3], dv[3], d_o[3];
  double r_v[3], vert[4][3], tet[4], v_x[3][3];
  double hex[3][2][3];
  double f, dx_dy, dz_dx, dy_dx, dx, dy, dz;
  
  const double vertex_list[4][3] = { { 0, 0, 0 },
                                     { 0, 1, 0 },
                                     { 1, 0, 0 },
                                     { 0, 0, 1 } };

  /* const double vertex_list[8][3] = { { 0, 0, 0 },
                                     { 0, 1, 0 },
                                     { 1, 1, 0 },
                                     { 1, 0, 0 },
                                     { 0, 0, 1 },
                                     { 1, 0, 1 },
                                     { 1, 1, 1 },
                                     { 0, 1, 1 } }; */

  const double normals[3][3] =     { { 1, 0, 0 },
                                     { 0, 1, 0 },
                                     { 0, 0, 1 } };

  const double del[3] = { mesh->delx, mesh->dely, mesh->delz };

  const double emf = 0.000001;

  origin[0] = mesh->origin[0] + mesh->delx * i;
  origin[1] = mesh->origin[1] + mesh->dely * j;
  origin[2] = mesh->origin[2] + mesh->delz * k;

  for(s=0; s<4; s++) { /* iterate through each vertex */

    for(x=0; x<3; x++) {
      r_o[x] = origin[x] + del[x] * vertex_list[s][x];
    }

    bits=0;
      
    for(a=0; a<3; a++) { /* iterate through each secondary vertex */

      vert[a][0] = vertex_list[s][0];
      vert[a][1] = vertex_list[s][1];
      vert[a][2] = vertex_list[s][2];

      if(vertex_list[s][0] < emf) {
        r_v[0] = r_o[0] + normals[a][0] * del[0];
        vert[a][0] = vert[a][0] + normals[a][0];
      }
      else { 
        r_v[0] = r_o[0] - normals[a][0] * del[0];
        vert[a][0] = vert[a][1] - normals[a][0];
      }

      if(vertex_list[s][1] < emf) {
        r_v[1] = r_o[1] + normals[a][1] * del[1];
        vert[a][1] = vert[a][1] + normals[a][1];
      }
      else {
        r_v[1] = r_o[1] - normals[a][1] * del[1];
        vert[a][1] = vert[a][1] - normals[a][1];
      }

      if(vertex_list[s][2] < emf) {
        r_v[2] = r_o[2] + normals[a][2] * del[2];
        vert[a][2] = vert[a][2] + normals[a][2];
      }
      else {
        r_v[2] = r_o[2] - normals[a][2] * del[2];
        vert[a][2] = vert[a][2] - normals[a][2];
      }

      vector_subtract(d_o, r_v, origin);

      dv[0]  = fabs(r_o[0] - r_v[0]);
      dv[1]  = fabs(r_o[1] - r_v[1]);
      dv[2]  = fabs(r_o[2] - r_v[2]);

      /* iterate through each axis on each secondary vertex */
      x=0;
      for(b=0; b<3; b++) {  
           
        if(inner_product(dv, normals[b]) > emf)
          continue;

        sgn = inner_product(d_o, normals[b]) > emf ? -1 : 1; 

        o_n[0] = normals[b][0]; /* sgn;*/
        o_n[1] = normals[b][1]; /* sgn;*/
        o_n[2] = normals[b][2]; /* sgn;*/
        
        for(n=0; n < stl->facets; n++) {

          if((facing = moller_trumbore(r_v, o_n, stl->v_1[n],
                                       stl->v_2[n], stl->v_3[n],
                                       pt_int)) != 0) {
 
            f = (pt_int[0] - r_v[0])/del[0] +
                (pt_int[1] - r_v[1])/del[1] +
                (pt_int[2] - r_v[2])/del[2];
            f *= sgn;

            if(f > -1.0 * emf && f < (1+emf)) {
              hex[a][x][0] = /* pt_int[0];*/(pt_int[0] - r_v[0])/del[0]; 
              hex[a][x][1] = /* pt_int[1];*/(pt_int[1] - r_v[1])/del[1]; 
              hex[a][x][2] = /* pt_int[2];*/(pt_int[2] - r_v[2])/del[2];
              bits |= 1 << (a*2+x);
              
              f_facing = facing * sgn;
            }
       
          }

        }

        x++;
      }
      if(bits >= 63) break;
  
    }

    
    if(bits >= 63) {
      vert[3][0] = vertex_list[s][0];
      vert[3][1] = vertex_list[s][1];
      vert[3][2] = vertex_list[s][2];
      break;
    }

  }

  if(bits >= 63) {  
    /* we have intersections on all three connecting axis */

    /* need to extend the hexagon such that it forms a large
     * tetrahedron, and find the vertices of each point on that
     * tetrahedron. the vertices are stored in v_x.
     *
     * the code also rotates the hex so that the primary vertex 
     * is located at 0, 0, 0 */
    
    dx = (vert[0][0] + hex[0][0][0]) - (vert[1][0] + hex[1][0][0]);
    dy = (vert[0][1] + hex[0][0][1]) - (vert[1][1] + hex[1][0][1]);
    dx_dy = fabs(dx / dy);

    if(vert[0][0] > emf)
      v_x[0][0] = dx_dy * fabs(hex[0][0][1]) + vert[0][0];
    else
      v_x[0][0] = dx_dy * fabs(hex[0][0][1]) * -1.0;
    v_x[0][1] = vert[0][1];
    v_x[0][2] = vert[0][2];

    dy_dx = 1.0 / dx_dy;
    if(vert[1][1] > emf)
      v_x[1][1] = dy_dx * fabs(hex[1][0][0]) + vert[1][1];
    else
      v_x[1][1] = dy_dx * fabs(hex[1][0][0]) * -1.0;
    v_x[1][0] = vert[1][0];
    v_x[1][2] = vert[1][2];

    dx = (vert[2][0] + hex[2][0][0]) - (vert[0][0] + hex[0][1][0]);
    dz = (vert[2][2] + hex[2][0][2]) - (vert[0][2] + hex[0][1][2]);
    dz_dx = fabs(dz / dx);

    if(vert[2][2] > emf) 
      v_x[2][2] = dz_dx * fabs(hex[2][0][0]) + vert[2][2];
    else
      v_x[2][2] = dz_dx * fabs(hex[2][0][0]) * -1.0;
    v_x[2][0] = vert[2][0];
    v_x[2][1] = vert[2][1];

    /* REWRITE
     *
     * dx_dy = fabs((hex[1][0][0]-hex[0][0][0])/(hex[1][0][1]-hex[0][0][1])); 
    v_x[0][0] = dx_dy * (hex[0][0][1] - vert[0][1]) + 1; *vert[0][0];*
    v_x[0][1] = 0; * vert[0][1]; *
    v_x[0][2] = 0; * vert[0][2]; *

    dy_dx = 1.0 / dx_dy;
    v_x[1][0] = 0; * vert[1][0]; *
    v_x[1][1] = dy_dx * (hex[1][0][0] - vert[1][0]) + 1;
    v_x[1][2] = 0; * vert[1][2]; *

    dz_dx = fabs((hex[2][0][2] - hex[0][1][2])/(hex[2][0][0] - hex[0][1][0]));
    v_x[2][0] = 0; * vert[2][0]; *
    v_x[2][1] = 0; * vert[2][1]; *
    v_x[2][2] = dz_dx * (hex[2][0][0] - vert[2][0]) + 1; */

    for(x = 0; x<3; x++) {
      vector_add(hex[x][0], hex[x][0], vert[x]);
      vector_add(hex[x][1], hex[x][1], vert[x]);
    }

    tet[0] = fabs(tet_volume(v_x[0], hex[0][0], hex[0][1], vert[0]));
    tet[1] = fabs(tet_volume(v_x[1], hex[1][0], hex[1][1], vert[1]));
    tet[2] = fabs(tet_volume(v_x[2], hex[2][0], hex[2][1], vert[2]));
    tet[3] = fabs(tet_volume(v_x[0], v_x[1], v_x[2], vert[3]));

    mesh->fv[mesh_index(mesh, i, j, k)] = 
      tet[3] - (tet[0] + tet[1] + tet[2]);

    if(f_facing < 0)
      mesh->fv[mesh_index(mesh, i, j, k)] = 1 -
        mesh->fv[mesh_index(mesh, i, j, k)];

#ifdef DEBUG
    printf("%ld %ld %ld: Found vertex_hex, volume of %lf\n", i, j, k, 
      mesh->fv[mesh_index(mesh, i, j ,k)]);
#endif

    return 1;
    
  }

  return 0;
}

int line_pent_fraction(struct mesh_data *mesh, struct stl_data *stl,
                       long int i, long int j, long int k) {
  long int n, a, s, v, x;
  int sgn;
  int bits = 0;
  int facing, f_facing;

  double pt_int[3], origin[3], o_n[3], pent[4][3];
  double v1[3], v2[3], dv[3];
  double f;

  const double vertex_list[4][3] = { { 0, 0, 0 },
                                     { 1, 1, 0 },
                                     { 1, 0, 1 },
                                     { 0, 1, 1 } };


  /* normals representing each face */
  const double normal_list[3][3] = { {  1,  0,  0 },
                                     {  0,  1,  0 },
                                     {  0,  0,  1 } };
                 

  const double identity[3] = { 1, 1, 1 };

  const double del[3] = { mesh->delx, mesh->dely, mesh->delz };

  const double emf = 0.000001;

  origin[0] = mesh->origin[0] + mesh->delx * i;
  origin[1] = mesh->origin[1] + mesh->dely * j;
  origin[2] = mesh->origin[2] + mesh->delz * k;

  for(v=0; v<4; v++) { /* iterate through each vertex */

   for(s=0; s<3; s++) { /* iterate through each secondary vertex */
    
      /* bits stores intersection locations as follows:
       * bit 0: first intersection on primary vertex
       * bit 1: second intersection on primary vertex
       * bit 2: first intersection on secondary vertex
       * bit 3: second intersection on secondary vertex
       *
       * first and second intersection must be ordered
       * for example, x and then z, never z and then x */
      bits=0;
     
      /* set the primary vertex 
       * this cannot be done in the parent loop because we later
       * store the cartesian coordinate in this variable */ 
      v1[0] = vertex_list[v][0];
      v1[1] = vertex_list[v][1];
      v1[2] = vertex_list[v][2];
 
      /* set the secondary vertex */
      v2[0]  = v1[0] + normal_list[s][0] * (v1[0]>0? -1 : 1);
      v2[1]  = v1[1] + normal_list[s][1] * (v1[1]>0? -1 : 1);
      v2[2]  = v1[2] + normal_list[s][2] * (v1[2]>0? -1 : 1);

      /* determine v1 - v2 and store it.  we know that when testing
       * for line vertices, we need to check 2 normals per vertex.  
       * considering 3 orthogonal normals, the normal that is not 
       * included is dv */
      dv[0]  = fabs(v1[0] - v2[0]);
      dv[1]  = fabs(v1[1] - v2[1]);
      dv[2]  = fabs(v1[2] - v2[2]);

      /* adjust v1 and v2 for the cartesian coordinates */
      v1[0]  = del[0] * v1[0] + origin[0];
      v1[1]  = del[1] * v1[1] + origin[1];
      v1[2]  = del[2] * v1[2] + origin[2];
      
      v2[0]  = del[0] * v2[0] + origin[0];
      v2[1]  = del[1] * v2[1] + origin[1];
      v2[2]  = del[2] * v2[2] + origin[2];

      for(n=0; n < stl->facets; n++) {

        x = 0;
        for(a=0; a<3; a++) { /* iterate through each axis */
          
          if(inner_product(dv, normal_list[a]) == 1)
            continue;

          if(vertex_list[v][a] > 0) sgn = -1;
          else sgn = 1;
          
          o_n[0] = normal_list[a][0] * sgn;
          o_n[1] = normal_list[a][1] * sgn;
          o_n[2] = normal_list[a][2] * sgn;
   
          /* check for intersections on primary axis */
          if((facing = moller_trumbore(v1, o_n, stl->v_1[n],
                                       stl->v_2[n], stl->v_3[n],
                                       pt_int)) != 0) {
            f = (pt_int[0] - v1[0])/del[0] +
                (pt_int[1] - v1[1])/del[1] +
                (pt_int[2] - v1[2])/del[2];
            f *= sgn;

            if(f > -1.0 * emf && f < (1+emf)) {
              pent[x][0] = fabs(pt_int[0] - v1[0]); 
              pent[x][1] = fabs(pt_int[1] - v1[1]);
              pent[x][2] = fabs(pt_int[2] - v1[2]);
              bits |= 1 << x;
              f_facing = facing;
            }
          }

          /* check for intersections on secondary axis */
          if((facing = moller_trumbore(v2, o_n, stl->v_1[n],
                                       stl->v_2[n], stl->v_3[n],
                                       pt_int)) != 0) {
            f = (pt_int[0] - v2[0])/del[0] +
                (pt_int[1] - v2[1])/del[1] +
                (pt_int[2] - v2[2])/del[2];
            f *= sgn;

            if(f > -1.0 * emf && f < (1+emf)) {
              pent[x+2][0] = fabs(pt_int[0] - v2[0]); 
              pent[x+2][1] = fabs(pt_int[1] - v2[1]);
              pent[x+2][2] = fabs(pt_int[2] - v2[2]);
              bits |= 1 << (x+2);
              f_facing = facing;
            }

 
          }

          x++;

        }
      }
     
      if( (bits & (1 << 0)) && (bits & (1 << 1)) && 
          (bits & (1 << 2)) && (bits & (1 << 3))) break;
    }
    
    if( (bits & (1 << 0)) && (bits & (1 << 1)) && 
        (bits & (1 << 2)) && (bits & (1 << 3))) break;
  }

  if( (bits & (1 << 0)) && (bits & (1 << 1)) && 
      (bits & (1 << 2)) && (bits & (1 << 3))) {
    /* we have intersections on both primary and secondary vertices */

    f = ((inner_product(pent[0], identity) * 
          inner_product(pent[1], identity))) /2 +
        ((inner_product(pent[2], identity) * 
          inner_product(pent[3], identity))) /2;

    f /= 2.0;
    f *= inner_product(dv, del);


    mesh->fv[mesh_index(mesh, i, j, k)] = f;

    if(f_facing < 0)
      mesh->fv[mesh_index(mesh, i, j, k)] = (del[0] * del[1] * del[2]) -
        mesh->fv[mesh_index(mesh, i, j, k)];

    mesh->fv[mesh_index(mesh, i, j, k)] /= del[0] * del[1] * del[2];

#ifdef DEBUG
    printf("%ld %ld %ld: Found pent, volume of %lf\n", i, j, k, 
      mesh->fv[mesh_index(mesh, i, j ,k)]);
#endif

    return 1;
    
  }

  return 0;

 
}

int face_hex_fraction(struct mesh_data *mesh, struct stl_data *stl,
                       long int i, long int j, long int k) {
  long int n, a, s;
  int sgn;
  int bits = 0;
  int facing, f_facing;

  double pt_int[3], origin[3], r_o[3], o_n[3], hex[4][3];
  double f, face_area;

  /* normals representing each face */
  const double normal_list[3][3] = { {  1,  0,  0 },
                                     {  0,  1,  0 },
                                     {  0,  0,  1 } };
                 

  const double vertex_list[3][4][3] = 
                                   { { {  0,  0,  0},
                                       {  0,  1,  0},
                                       {  0,  1,  1},
                                       {  0,  0,  1} },
                                     { {  0,  0,  0},
                                       {  1,  0,  0},
                                       {  1,  0,  1},
                                       {  0,  0,  1} },
                                     { {  0,  0,  0},
                                       {  1,  0,  0},
                                       {  1,  1,  0},
                                       {  0,  1,  0} } };

  const double identity[3]   =       { 1, 1, 1};

  const double del[3] = { mesh->delx, mesh->dely, mesh->delz };

  const double emf = 0.000001;

  origin[0] = mesh->origin[0] + mesh->delx * i;
  origin[1] = mesh->origin[1] + mesh->dely * j;
  origin[2] = mesh->origin[2] + mesh->delz * k;

  for(s=0; s<3; s++) { /* iterate through each face */

    o_n[0] = normal_list[s][0];
    o_n[1] = normal_list[s][1];
    o_n[2] = normal_list[s][2];

    bits=0;

    for(n=0; n < stl->facets; n++) {

      for(a=0; a<4; a++) { /* iterate through each face vertex */
        
        r_o[0] = vertex_list[s][a][0] * del[0] + origin[0];
        r_o[1] = vertex_list[s][a][1] * del[1] + origin[1];
        r_o[2] = vertex_list[s][a][2] * del[2] + origin[2];
  
        sgn = o_n[0] + o_n[1] + o_n[2];
        if((facing = moller_trumbore(r_o, o_n, stl->v_1[n],
                                     stl->v_2[n], stl->v_3[n],
                                     pt_int)) != 0) {
          f = (pt_int[0] - r_o[0])/del[0] +
              (pt_int[1] - r_o[1])/del[1] +
              (pt_int[2] - r_o[2])/del[2];
          f *= sgn;

          if(f > -1.0 * emf && f < (1+emf)) {
            hex[a][0] = fabs(pt_int[0] - r_o[0]); 
            hex[a][1] = fabs(pt_int[1] - r_o[1]);
            hex[a][2] = fabs(pt_int[2] - r_o[2]);
            bits |= 1 << a;
            f_facing = facing;
          }

        }

      }
    }

    
    if( (bits & (1 << 0)) && (bits & (1 << 1)) && 
        (bits & (1 << 2)) && (bits & (1 << 3))) break;
  }

  if( (bits & (1 << 0)) && (bits & (1 << 1)) && 
      (bits & (1 << 2)) && (bits & (1 << 3))) {
    /* we have intersections on all three connecting axis */

    f = inner_product(identity, hex[0]) +
        inner_product(identity, hex[1]) +
        inner_product(identity, hex[2]) +
        inner_product(identity, hex[3]);

    f /= 4;

    switch(s) {
      case 0:
        face_area = del[1] * del[2];
        break;

      case 1:
        face_area = del[0] * del[2];
        break;

      case 2:
        face_area = del[0] * del[1];
        break;
    }

    mesh->fv[mesh_index(mesh, i, j, k)] = face_area * f;

    if(f_facing < 0)
      mesh->fv[mesh_index(mesh, i, j, k)] = (del[0] * del[1] * del[2]) -
        mesh->fv[mesh_index(mesh, i, j, k)];

    mesh->fv[mesh_index(mesh, i, j, k)] /= del[0] * del[1] * del[2];

#ifdef DEBUG
    printf("%ld %ld %ld: Found hex, volume of %lf\n", i, j, k, 
      mesh->fv[mesh_index(mesh, i, j ,k)]);
#endif

    return 1;
    
  }

  return 0;

}

int tet_fraction(struct mesh_data *mesh, struct stl_data *stl,
                 long int i, long int j, long int k) {

  long int n, a, s;
  int x, sgn;
  int bits = 0;
  int facing, f_facing;
  int flg;

  double pt_int[3], origin[3], r_o[3], o_n[3];
  double tet[4][3];
  double f;

  const double vertex_list[8][3] = { { 0, 0, 0 },
                                     { 0, 1, 0 },
                                     { 1, 1, 0 },
                                     { 1, 0, 0 },
                                     { 0, 0, 1 },
                                     { 1, 0, 1 },
                                     { 1, 1, 1 },
                                     { 0, 1, 1 } };

  const double normals[3][3] =     { { 1, 0, 0 },
                                     { 0, 1, 0 },
                                     { 0, 0, 1 } };

  const double del[3] = { mesh->delx, mesh->dely, mesh->delz };

  const double emf = 0.000001;

  origin[0] = mesh->origin[0] + mesh->delx * i;
  origin[1] = mesh->origin[1] + mesh->dely * j;
  origin[2] = mesh->origin[2] + mesh->delz * k;

  flg=0;

  for(s=0; s<8; s++) { /* iterate through each vertex */

    for(x=0; x<3; x++) {
      r_o[x] = origin[x] + del[x] * vertex_list[s][x];
    }

    bits=0;

    for(n=0; n < stl->facets; n++) {
  
      for(a=0; a<3; a++) { /* iterate through each axis */
        
        for(x=0; x<3; x++) {
          o_n[x] = normals[a][x] * 
               (vertex_list[s][x] > 0 ? -1.0 : 1.0);
        }

       
        sgn = o_n[a];
        if((facing = moller_trumbore(r_o, o_n, stl->v_1[n],
                                     stl->v_2[n], stl->v_3[n],
                                     pt_int)) != 0) {
          f = (pt_int[a] - r_o[a]) / del[a];
          f *= sgn;

          if(f > -1.0 * emf && f < (1+emf)) {
            flg = 2; /* flag to signal that there are intersections
                      * in this cell, even if it isn't for a tet */

            tet[a][0] = fabs(pt_int[0] - r_o[0]); 
            tet[a][1] = fabs(pt_int[1] - r_o[1]);
            tet[a][2] = fabs(pt_int[2] - r_o[2]);
            bits |= 1 << a;
            f_facing = facing;
          }

        }

      }
    }

    
    if( (bits & (1 << 0)) && (bits & (1 << 1)) && (bits & (1 << 2))) {
      tet[3][0] = 0.0; /*del[0] * vertex_list[s][0];*/
      tet[3][1] = 0.0; /*del[1] * vertex_list[s][1];*/
      tet[3][2] = 0.0; /*del[2] * vertex_list[s][2];*/
      break;
    }


  }

  if( (bits & (1 << 0)) && (bits & (1 << 1)) && (bits & (1 << 2))) {
    /* we have intersections on all three connecting axis */

    mesh->fv[mesh_index(mesh, i, j, k)] = 
             fabs(tet_volume(tet[0], tet[1], tet[2], tet[3] ));

    if(f_facing < 0)
      mesh->fv[mesh_index(mesh, i, j, k)] = (del[0] * del[1] * del[2]) -
        mesh->fv[mesh_index(mesh, i, j, k)];

    mesh->fv[mesh_index(mesh, i, j, k)] /= del[0] * del[1] * del[2];

#ifdef DEBUG
    printf("%ld %ld %ld: Found tet, volume of %lf\n", i, j, k, 
      mesh->fv[mesh_index(mesh, i, j ,k)]);
#endif

    return 1;
    
  }

  return flg;
}

double tet_volume(double *a, double *b, double *c, double *d) {

  double b_m_d[3], c_m_d[3], cross[3], a_m_d[3];
  double inner;

  vector_subtract(b_m_d, b, d);
  vector_subtract(c_m_d, c, d);

  cross_product(cross, b_m_d, c_m_d);

  vector_subtract(a_m_d, a, d);

  inner = inner_product(a_m_d, cross);

  return inner/6.0;
}

int count_bits(int n) {     
  unsigned int c; /* c accumulates the total bits set in v */
  for (c = 0; n; c++) 
    n &= n - 1; /* clear the least significant bit set */
  return c;
}
