/*
 * mesh_data.h
 *
 * defines data structures to prevent circular headers */

#ifndef _MESH_DATA_H
#define _MESH_DATA_H

/* define how 3d arrays are indexed */
enum cell_boundaries {
  east=1,
  west=2,
  north=3,
  south=4,
  top=5,
  bottom=6,
  none=7
};

enum wall_boundaries {
  slip=0,
  no_slip=1,
  zero_gradient=2
};

enum special_boundaries {
  fixed_velocity = 0,
  mass_outflow = 1,
  hgl = 2,
  weir = 3,
  wall = 4
};

enum baffle_type {
  flow = 0,
  barrier = 1,
  k = 2,
  swirl_angle = 3,
  v_deviation = 4
};

struct sb_data {
  /* Describes special boundaries
   * These are fixed domain boundaries that exist on wall boundaries
   * used to set features such as inflow and outflow conditions
   */

  enum special_boundaries type;

  long int extent_a[2]; /* describes the extent in the two dimensions parallel to the
                       * boundary plane */
  long int extent_b[2]; /* describes the extent in the two dimensions parallel to the
                       * boundary plane */                       
  double value;       /* sets the value at the boundary */
  double turbulence;  /* stores turbulence data at the boundary, if relevent */
  
  struct sb_data *next;
};


struct baffle_data {
  /* Describes baffles
   * These are fixed domain boundaries that exist inside the mesh
   * used to measure flow, induce backpressure, or block flow
   */

  enum baffle_type type;

  long int extent_a[2]; /* describes the extent in the two dimensions parallel to the
                       * boundary plane */
  long int extent_b[2]; /* describes the extent in the two dimensions parallel to the
                       * boundary plane */                       
  double value;       /* value dependent on type */
  long int pos;
  
  struct baffle_data *next;
};

struct mesh_data {
  /* Describes a single mesh
   * Code will eventually allow multiblock meshes
   * Whereby there will be multiple meshes, stiched together
   */

  long int imax, jmax, kmax;
  double delx, dely, delz;
  double rdx, rdy, rdz;
  double origin[3]; /* vector describing the coordinate of the 0,0,0
                     * point in the mesh */
  double inside[3]; /* a point inside the solution space */

  /* mesh edge boundaries */
  /* number   boundary
      0         west
      1         east
      2         south
      3         north
      4         bottom
      5         top */
  enum wall_boundaries wb[6];

  struct sb_data *sb[6]; /* linked lists of special boundaries for each wall */
  struct baffle_data *baffles[3]; /* linked lists of baffles for each axis */
  
  /* flag to indicate that the above properties are assigned
   * and mesh is ready for dynamic array allocation */
  int ready;
  
  /* flag to indicate whether output should be compressed */
  int compress;
  
  /* The next set of data describes 1D arrays of data that describe
   * the fluid properties.  The 1D arrays describe flattened
   * 3D arrays. */

  /* Calculated outputs describing fluid properties */
  double *P, *D;
  double *u, *v, *w;
  double *u_omega, *v_omega, *w_omega; /* vorticity */

  /* Volume of Fluid in each cell */
  double *vof;
  enum cell_boundaries *n_vof;

  /* Fractional area / volumes
   * fv = Fractional volume
   * ae = Fractional area to the east (x axis)
   * an = Fractional area to the north (y axis)
   * at = Fractional area to the top (z axis)
   */
  double *fv;
  double *ae, *an, *at;
  
  /* Coefficients used in calculations */
  double *peta, *tanth, *beta;
  
  /* Turbulence model declared void to allow flexibility
   * since different models use different variables 
   * Note the turbulence model routines are responsible for their 
   * own memory management */
   void *turbulence_model;
   double *nut; /* turbulent viscocity */
};

#endif
