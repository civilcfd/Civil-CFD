/* vof_boundary.c
 *
 * implement boundaries
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vtk.h"
#include "vof_mpi.h"
#include "solver.h"
#include "solver_mpi.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "vof_boundary.h"
#include "vof_baffles.h"

#include "vof_macros.h"

extern struct mesh_data *mesh_n; /* describes mesh at previous timestep for explicit calcs
                                    defined in vof.c */

  /* mesh edge boundaries */
  /* number   boundary
      0         west
      1         east
      2         south
      3         north
      4         bottom
      5         top */

enum special_boundaries vof_boundaries_check_inside_sb(struct solver_data *solver, long int a, long int b,
                                 int x) {
  struct sb_data *sb;
  
  switch(x) {
  case 0: /* west */
  case 1: /* east */
    break;
  case 2: /* south */
  case 3:
  case 4:
  case 5:
    a += ISTART;
    break;
  }
  
  for(sb = solver->mesh->sb[x]; sb != NULL; sb = sb->next) {   
      
        if(a >= sb->extent_a[0] && a <= sb->extent_b[0] &&
           b >= sb->extent_a[1] && b <= sb->extent_b[1])
          return sb->type;
  }    
   
  
  return wall;
}

int sboundary_setup(struct solver_data *solver, int x, long int *imin, long int *jmin, long int *kmin, 
                           long int *imax, long int *jmax, long int *kmax,
                           double min_1, double min_2, double max_1, double max_2) {
  switch(x) {
  case 0:
    *imin = 0;
    *imax = 0;
    break;
  case 1:
    *imin = IMAX-1;
    *imax = IMAX-1;
    break;
  case 2:
    *jmin = 0;
    *jmax = 0;
    break;
  case 3:
    *jmin = JMAX-1;
    *jmax = JMAX-1;
    break;    
  case 4:
    *kmin = 0;
    *kmax = 0;
    break;
  case 5:
    *kmin = KMAX-1;
    *kmax = KMAX-1;
    break;   
  }

  switch(x) {
  case 0:
  case 1:
    *jmin = min_1;
    *jmax = min(JMAX-1,max_1);
    *kmin = min_2;
    *kmax = min(KMAX-1,max_2);    
    break;
  case 2:
  case 3:
    *imin = min_1;
    *imax = min(IMAX-1,max_1);
    *kmin = min_2;
    *kmax = min(KMAX-1,max_2);   
    break;    
  case 4:
  case 5:
    *imin = min_1;
    *imax = min(IMAX-1,max_1);
    *jmin = min_2;
    *jmax = min(JMAX-1,max_2);   
    break;   
  }

  *imin = *imin - ISTART;
  *imin = max(*imin, 0);
  *imax = *imax - ISTART;
  *imax = min(*imax, IRANGE-1);

  return 0;
}

int boundary_hgl(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence) {
  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  double height, partial;
  long int coplanar[3] = {0, 0, 0};
  
  struct mesh_data *mesh = solver->mesh;
  
  sboundary_setup(solver, x, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
    
  /* value += mesh->delz; uncomment and everything is off-set by 1 cell vertically */
  value -= mesh->origin[2]; /* ADDED 7/27/2018 */

  switch(x) {
  case 0: /* west */
    /*imax++*/;
    coplanar[0] = 1;
    break;
  case 1: /* east */
    /*imin--*/;
    coplanar[0] = -1;
    break;
  case 2: /* south */
    /*jmax++*/;
    coplanar[1] = 1;
    break;
  case 3: /* north */
    /*jmin--*/;
    coplanar[1] = -1;
    break;
  case 4: /* bottom */
    /*kmax++*/;
    coplanar[2] = 1;
    break;
  case 5: /* top */
    /*kmin--*/;
    coplanar[2] = -1;
    break;
  }
  
  if(fmod(value,DELZ) < solver->emf) value -= solver->min_vof * 100;
 
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k = kmax-1; k > kmin-1; k--) {
        
        //if(FV(i,j,k) < 0.000001) continue;
      
        /* must set velocity to a Neumann boundary */
        
        /* first set U(i,j,k) equal to the next neighboring velocity inside the mesh */
        U(i,j,k) = U(i+coplanar[0], j+coplanar[1], k+coplanar[2]);
        V(i,j,k) = V(i+coplanar[0], j+coplanar[1], k+coplanar[2]);
        W(i,j,k) = W(i+coplanar[0], j+coplanar[1], k+coplanar[2]);
        
        /* now correct for east / north / top boundary */
        switch(x) {
        case 1:
          U(i,j,k) = U(i-2,j,k);
          U(i-1,j,k) = U(i-2,j,k); 
          break;
        case 3:
          V(i,j,k) = V(i,j-2,k);
          V(i,j-1,k) = V(i,j-2,k);
          break;
        case 5:
          W(i,j,k) = W(i,j,k-2);
          W(i,j,k-1) = W(i,j,k-2);
          break;
        }

        
        height = k * mesh->delz;
        if(value - height >= mesh->delz - solver->emf) {
          VOF(i,j,k) = 1.0;
          /*VOF(i+1,j,k) = 1.0;*/
        }
        else {
          partial = value - height;
          if(partial > solver->emf) {
            partial /= mesh->delz;
            VOF(i,j,k) = partial;  
            /*VOF(i+1,j,k) = partial;          */
          } 
          /* ADDED 01/03/2014 */
          else if(partial <= 0)
            VOF(i,j,k) = 0.0;
        }
        /* VOF(i+coplanar[0],j+coplanar[1],k+coplanar[2]) = (VOF(i,j,k) + VOF(i+2*coplanar[0],j+2*coplanar[1],k+2*coplanar[2]))/2; */
        
             
        /* ADDED 01/10/2014 */
        if(mesh->vof[mesh_index(mesh,i,j,k)] == 1.0) {
          mesh->P[mesh_index(mesh,i,j,k)] = 
            mesh->P[mesh_index(mesh,i,j,k+1)] + 
            (mesh->delz * solver->rho * fabs(solver->gz)) * 
            (0.5 + min(mesh->vof[mesh_index(mesh,i,j,k+1)],0.5));

        }
        else  if(mesh->vof[mesh_index(mesh,i,j,k)] > 0.0) {
          mesh->P[mesh_index(mesh,i,j,k)] = 
            mesh->delz * solver->rho * fabs(solver->gz) * max(-0.5,mesh->vof[mesh_index(mesh,i,j,k)] - 0.5);          
        }
        else if(mesh->vof[mesh_index(mesh,i,j,k)] <= 0.0) 
          mesh->P[mesh_index(mesh,i,j,k)] = 0.0;
          
        if(FV(i,j,k) < 0.000001) continue;
        P(i+coplanar[0],j+coplanar[1],k+coplanar[2]) = P(i,j,k); 
        VOF(i+coplanar[0],j+coplanar[1],k+coplanar[2]) = VOF(i,j,k); 


     }
    }
  }
  return 0;
}

double calc_flow(struct solver_data *solver, int x, long int imin, long int imax, 
                 long int jmin, long int jmax, long int kmin, long int kmax, double *area_ref) {
  long int i,j,k;
  
  double area_0, area, flow;
  area_0 = 0;
  area = 0;
  flow = 0;

  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {  
        switch(x) {
        case 0:
          area_0 = DELY * DELZ * AE(i+1,j,k);
          area_0 *= (VOF(i+1,j,k) + VOF(i+2,j,k))/2;
          
          flow += UN(i+1,j,k) * area_0;
          area += area_0;
          break;          
        case 1:
          area_0 = DELY * DELZ * AE(i-2,j,k);
          area_0 *= (VOF(i-1,j,k) + VOF(i-2,j,k))/2;
          
          flow += UN(i-2,j,k) * area_0;
          area += area_0;
          break;
        case 2:
          area_0 = DELX * DELZ * AN(i,j+1,k);
          area_0 *= (VOF(i,j+1,k) + VOF(i,j+2,k))/2;
          
          flow += VN(i,j+1,k) * area_0;
          area += area_0;
          break;          
        case 3:
          area_0 = DELX * DELZ * AN(i,j-2,k);
          area_0 *= (VOF(i,j-1,k) + VOF(i,j-2,k))/2;
          
          flow += VN(i,j-2,k) * area_0;
          area += area_0;
          break;          
        case 4:
          area_0 = DELX * DELY * AT(i,j,k+1);
          area_0 *= (VOF(i,j,k+1) + VOF(i,j,k+2))/2;
          
          flow += WN(i,j,k+1) * area_0;
          area += area_0;
          break;          
        case 5:
          area_0 = DELX * DELY * AT(i,j,k-2);
          area_0 *= (VOF(i,j,k-1) + VOF(i,j,k-2))/2;
          
          flow += WN(i,j,k-2) * area_0;
          area += area_0;
          break;     
        }        
      }
    }
  }
  *area_ref = area;

  flow = solver_mpi_sum(solver, flow);

  return flow;
}

int boundary_weir(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence) {
#define emf solver->emf  

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax, l, m, n;
  double flow, area, height, ave_height, count, head, sgn; 
  double coplanar[3] = {0, 0, 0};
  
  sboundary_setup(solver, x, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);

  flow = 0;
  area = 0;
  sgn  = 1.0;
  
  if(x>3) return 0;
  
  
  switch(x) {
  case 0: /* west */
    /*imax++*/;
    coplanar[0] = 1;
    break;
  case 1: /* east */
    /*imin--*/;
    coplanar[0] = -1;
    break;
  case 2: /* south */
    /*jmax++*/;
    coplanar[1] = 1;
    break;
  case 3: /* north */
    /*jmin--*/;
    coplanar[1] = -1;
    break;
  }
  
  
  switch(x) {
  case 0: /* west */
  case 2: /* south */
  case 4: /* bottom */
    sgn *= -1.0;
    break;
  }
  
  flow = calc_flow(solver, x, imin, imax, jmin, jmax, kmin, kmax, &area);
  flow *= sgn;
  
  count = 0;
  ave_height = 0;
  
  /* now calculate height of liquid */
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k = kmax-1; k > kmin-1; k--) {
        
        if(FV(i,j,k) < 0.000001) continue;

        l = i+coplanar[0];
        m = j+coplanar[1];
        n = k+coplanar[2];
        
        if(VOF(l,m,n) > (1-emf)) {
          height = k * DELZ;
          break;
        } else if(VOF(l,m,n) > emf) {
          if(VOF(l,m,max(0,n-1)) > emf) {
            height = max(k-1,0) * DELZ + VOF(l,m,n) * DELZ;
            break;
          }
        }
      }
      
      ave_height += height;
      count += 1.0;
    }
  } 
  
  ave_height = ave_height / count;
  
  /* first check water level is greater than weir crest, if not, v=0 */
  if(ave_height <= turbulence) {
    boundary_fixed_velocity(solver, x, min_1, min_2, max_1, max_2, 0, 0.00001);
  }
  else {

    if(flow < emf) flow = emf;
    head = pow(flow / (1.6 * value),0.667);
    boundary_hgl(solver, x, min_1, min_2, max_1, max_2, head + turbulence, 0.00001); 
    
    /*
    head = ave_height - solver->vof_height;
    flow = 1.6 * value * pow(head, 1.5);
    boundary_mass_outflow(solver, x, min_1, min_2, max_1, max_2, flow, 0);    */
  }

  return 0;
}

int boundary_mass_outflow(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence) {
#define emf solver->emf  

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  
  sboundary_setup(solver, x, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  double flow, area, area_0, flow_factor;
  flow = 0;
  area = 0;
  
  
  /* this is outflow so set value so that it is positive on an east/north/top boundary, and
  negative otherwise */
  value = fabs(value);
  
  switch(x) {
  case 0: /* west */
  case 2: /* south */
  case 4: /* bottom */
    value *= -1.0;
    break;
  }
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {  
        switch(x) {
        case 0:
          area_0 = DELY * DELZ * AE(i+1,j,k);
          area_0 *= (VOF(i+1,j,k) + VOF(i+2,j,k))/2;
          
          if(UN(i+1,j,k) * value > 0) flow += UN(i+1,j,k) * area_0;
          area += area_0;
          break;          
        case 1:
          area_0 = DELY * DELZ * AE(i-2,j,k);
          area_0 *= (VOF(i-1,j,k) + VOF(i-2,j,k))/2;
          
          if(UN(i-2,j,k) * value > 0) flow += UN(i-2,j,k) * area_0;
          area += area_0;
          break;
        case 2:
          area_0 = DELX * DELZ * AN(i,j+1,k);
          area_0 *= (VOF(i,j+1,k) + VOF(i,j+2,k))/2;
          
          if(VN(i,j+1,k) * value > 0) flow += VN(i,j+1,k) * area_0;
          area += area_0;
          break;          
        case 3:
          area_0 = DELX * DELZ * AN(i,j-2,k);
          area_0 *= (VOF(i,j-1,k) + VOF(i,j-2,k))/2;
          
          if(VN(i,j-2,k) * value > 0 ) flow += VN(i,j-2,k) * area_0;
          area += area_0;
          break;          
        case 4:
          area_0 = DELX * DELY * AT(i,j,k+1);
          area_0 *= (VOF(i,j,k+1) + VOF(i,j,k+2))/2;
          
          if(WN(i,j,k+1) * value > 0) flow += WN(i,j,k+1) * area_0;
          area += area_0;
          break;          
        case 5:
          area_0 = DELX * DELY * AT(i,j,k-2);
          area_0 *= (VOF(i,j,k-1) + VOF(i,j,k-2))/2;
          
          if(WN(i,j,k-2) * value > 0) flow += WN(i,j,k-2) * area_0;
          area += area_0;
          break;     
        }        
      }
    }
  }

  /* ADDED 7/27/18 */
  flow = solver_mpi_sum(solver, flow);
  area = solver_mpi_sum(solver, area);

  if(fabs(flow) < 0.1 * fabs(value) || flow * value < 0) { 
    /* in the case of no flow or reverse flow, we set a fixed velocity to start the solution */
    boundary_fixed_velocity(solver, x, min_1, min_2, max_1, max_2, value/area, turbulence);
    return 0;
  }
  
  flow_factor = flow / value;
  flow_factor = min(flow_factor, 1.25); /* maximum 25% deviation per timestep */
  flow_factor = max(flow_factor, 0.8); 

  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {        
        switch(x) {
        case 0:
        	if(VOF(i+1,j,k) > emf && FV(i,j,k) > emf) {
          	if(UN(i+1,j,k) * value > 0) U(i,j,k) = UN(i+1,j,k) / flow_factor;
          	else U(i,j,k) = 0.1 * value/area;
          	VOF(i,j,k) = 1.0;
          } 
          else {
          	U(i,j,k) = U(i+1,j,k);
          }
          V(i,j,k) = V(i+1,j,k);
          W(i,j,k) = W(i+1,j,k);
          /* V(i,j,k) = 0;
          W(i,j,k) = 0; */
          break;
        case 1:
        	if(VOF(i-1,j,k) > emf && FV(i,j,k) > emf) {
          	if(UN(i-2,j,k) * value > 0) U(i,j,k) = UN(i-2,j,k) / flow_factor;
          	else U(i,j,k) = 0.1 * value/area;
          	VOF(i,j,k) = 1.0;
          }
          else {
          	U(i,j,k) = U(i-2,j,k);
          }
          
          U(i-1,j,k) = U(i,j,k);
          
          V(i,j,k) = V(i-1,j,k);
          W(i,j,k) = W(i-1,j,k);
          /* V(i,j,k) = 0;
          V(i-1,j,k) = 0;
          W(i,j,k) = 0;
          W(i-1,j,k) = 0; */
          break;
        case 2:
        	if(VOF(i,j+1,k) > emf && FV(i,j,k) > emf) {
          	if(VN(i,j+1,k) * value > 0) V(i,j,k) = VN(i,j+1,k) / flow_factor;
          	else V(i,j,k) = 0.1 * value/area;
          	VOF(i,j,k) = 1.0;
          }
          else {
          	V(i,j,k) = V(i,j+1,k);
          }
          
          U(i,j,k) = U(i,j+1,k);
          W(i,j,k) = W(i,j+1,k);
          /* U(i,j,k) = 0;
          W(i,j,k) = 0; */
          break;
        case 3:
        	if(VOF(i,j-1,k) > emf && FV(i,j,k) > emf) {
          	if(VN(i,j-2,k) * value > 0) V(i,j,k) = VN(i,j-2,k) / flow_factor;
          	else V(i,j,k) = 0.1 * value/area;
          	VOF(i,j,k) = 1.0;
          }
          else {
          	V(i,j,k) = V(i,j-2,k);
          }
          V(i,j-1,k) = V(i,j,k);
          
          U(i,j,k) = U(i,j-1,k);
          W(i,j,k) = W(i,j-1,k);
          /* U(i,j,k) = 0;
          U(i,j-1,k) = 0;
          W(i,j,k) = 0;
          W(i,j-1,k) = 0; */
          break;
        case 4:
        	if(VOF(i,j,k+1) > emf && FV(i,j,k) > emf) {
          	if(WN(i,j,k+1) * value > 0) W(i,j,k) = WN(i,j,k+1) / flow_factor;
          	else W(i,j,k) = 0.1 * value/area;
          	VOF(i,j,k) = 1.0;
          }
          else {
          	W(i,j,k) = W(i,j,k+1);
          }
          
          U(i,j,k) = U(i,j,k+1);
          V(i,j,k) = V(i,j,k+1);
          /* U(i,j,k) = 0;
          V(i,j,k) = 0; */
          break;
        case 5:
        	if(VOF(i,j,k-1) > emf && FV(i,j,k) > emf) {
          	if(WN(i,j,k-2) * value > 0) W(i,j,k) = WN(i,j,k-2) / flow_factor;
          	else W(i,j,k) = 0.1 * value/area;
          	VOF(i,j,k) = 1.0;
          }
          else {
          	W(i,j,k) = W(i,j,k-2);
          }
          W(i,j,k-1) = W(i,j,k);
          
          U(i,j,k) = U(i,j,k-1);
          V(i,j,k) = V(i,j,k-1);
          /* U(i,j,k) = 0;
          U(i,j,k-1) = 0;
          V(i,j,k) = 0;
          V(i,j,k-1) = 0; */
          break;
        }        
      }
    }
  }
  
  return 0;
}
#undef emf
      
int boundary_fixed_velocity(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence) {

  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  
  sboundary_setup(solver, x, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) {   
        switch(x) {
        case 0:
          if(value > 0)
            VOF(i,j,k) = mesh_n->vof[mesh_index(solver->mesh,i,j,k)];
          U(i,j,k) = value;
          V(i,j,k) = 0;
          W(i,j,k) = 0;
          break;
        case 1:
          if(value < 0)
            VOF(i,j,k) = mesh_n->vof[mesh_index(solver->mesh,i,j,k)];
          U(i,j,k) = value;
          U(i-1,j,k) = value;
          V(i,j,k) = 0;
          W(i,j,k) = 0;
          break;
        case 2:
          if(value > 0)
            VOF(i,j,k) = mesh_n->vof[mesh_index(solver->mesh,i,j,k)];
          U(i,j,k) = 0;
          V(i,j,k) = value;
          W(i,j,k) = 0;
          break; 
        case 3:
          if(value < 0)
            VOF(i,j,k) = mesh_n->vof[mesh_index(solver->mesh,i,j,k)];
          U(i,j,k) = 0;
          V(i,j,k) = value;
          V(i,j-1,k) = value;
          W(i,j,k) = 0;         
          break;    
        case 4:
          if(value > 0)
            VOF(i,j,k) = mesh_n->vof[mesh_index(solver->mesh,i,j,k)];
          U(i,j,k) = 0;
          V(i,j,k) = 0;
          W(i,j,k) = value;
          break;   
        case 5:
          if(value < 0)
            VOF(i,j,k) = mesh_n->vof[mesh_index(solver->mesh,i,j,k)];
          U(i,j,k) = 0;
          V(i,j,k) = 0;
          W(i,j,k) = value;
          W(i,j,k-1) = value;
          break;   
        }   
    
      }
    }
  }

  return 0;
}

int vof_special_boundaries(struct solver_data *solver) {
  int x;
  struct sb_data *sb;
#define emf solver->emf  

  for(x=0; x < 6; x++) {
    for(sb = solver->mesh->sb[x]; sb != NULL; sb = sb->next) {    
        switch(sb->type) {
        case fixed_velocity:
          boundary_fixed_velocity(solver, x, sb->extent_a[0], sb->extent_a[1], 
                                     sb->extent_b[0], sb->extent_b[1], 
                                     sb->value, sb->turbulence);
          break;
        case mass_outflow:
          boundary_mass_outflow(solver, x, sb->extent_a[0], sb->extent_a[1], 
                                     sb->extent_b[0], sb->extent_b[1], 
                                     sb->value, sb->turbulence);
          break;
        case hgl:
         	boundary_hgl(solver, x, sb->extent_a[0], sb->extent_a[1], 
                                     sb->extent_b[0], sb->extent_b[1], 
                                     sb->value, sb->turbulence);
          break;        
      	case weir:
          boundary_weir(solver, x, sb->extent_a[0], sb->extent_a[1], 
                                     sb->extent_b[0], sb->extent_b[1], 
                                     sb->value, sb->turbulence);
          break;
          
        }
    }
  }

  return 0;
#undef emf
}

int vof_boundaries(struct solver_data *solver) {

  long int i,j,k;
  int bm[6], bmtot;
  enum cell_boundaries nff;
#define dim(i,j,k) i+3*(j+k*3)
  double dv, dA, delp;
  
  
  /* first boundaries on x-axis */
  for(j=0; j<JMAX; j++) {
    for(k=0; k<KMAX; k++) {

      /* west boundary */
      if(ISTART == 0) {
        if(vof_boundaries_check_inside_sb(solver, j, k, 0) == wall) {    
          switch(solver->mesh->wb[0]) {
          case slip:
                /* slip case is zero_gradient for parallel axis
                * and zero fixed value for perpendicular axis */
            U(0,j,k) = 0;         
            V(0,j,k) = V(1,j,k);  
            W(0,j,k) = W(1,j,k);  
            break;
          case no_slip:
                /* no slip case is velocity = -1.0 * interior velocity at boundary */
            U(0,j,k) = -1.0 * U(1,j,k); 
            V(0,j,k) = -1.0 * V(1,j,k);  
            W(0,j,k) = -1.0 * W(1,j,k);  
            break;
          case zero_gradient:
                /* zero gradient means that velocity at boundary is equal to interior velocity */
            U(0,j,k) = U(1,j,k);  
            V(0,j,k) = V(1,j,k);  
            W(0,j,k) = W(1,j,k);
            break;       
          }
        }
      
        P(0,j,k) = P(1,j,k);
        VOF(0,j,k) = VOF(1,j,k);
      }
      
      /* east boundary */
      if(IRANGE + ISTART == IMAX) {
        if(vof_boundaries_check_inside_sb(solver, j, k, 1) == wall) {    
          switch(solver->mesh->wb[1]) {
          case slip:
                /* slip case is zero_gradient for parallel axis
                * and zero fixed value for perpendicular axis */
            U(IRANGE-1,j,k) = 0;
            U(IRANGE-2,j,k) = U(IRANGE-1,j,k);
            V(IRANGE-1,j,k) = V(IRANGE-2,j,k);
            W(IRANGE-1,j,k) = W(IRANGE-2,j,k);
            break;
          case no_slip:
                /* no slip case is velocity = -1.0 * interior velocity at boundary */
            U(IRANGE-1,j,k) = -1.0 * U(IRANGE-3,j,k);
            U(IRANGE-2,j,k) = U(IRANGE-1,j,k);
            V(IRANGE-1,j,k) = -1.0 * V(IRANGE-2,j,k);
            W(IRANGE-1,j,k) = -1.0 * W(IRANGE-2,j,k);
            break;
          case zero_gradient:
                /* zero gradient means that velocity at boundary is equal to interior velocity */
            U(IRANGE-1,j,k) = U(IRANGE-3,j,k);
            U(IRANGE-2,j,k) = U(IRANGE-1,j,k);
            V(IRANGE-1,j,k) = V(IRANGE-2,j,k);
            W(IRANGE-1,j,k) = W(IRANGE-2,j,k);
            break;       
          }
        }    
      
        P(IRANGE-1,j,k) = P(IRANGE-2,j,k);
        VOF(IRANGE-1,j,k) = VOF(IRANGE-2,j,k);
      }

    } 
  }
  
  /* boundaries on y-axis */
  for(i=0; i<IRANGE; i++) {
    for(k=0; k<KMAX; k++) {
         
      /* south boundary */
      if(vof_boundaries_check_inside_sb(solver, i, k, 2) == wall) { 
        switch(solver->mesh->wb[2]) {
        case slip:
              /* slip case is zero_gradient for parallel axis
               * and zero fixed value for perpendicular axis */
          U(i,0,k) = U(i,1,k);  
          V(i,0,k) = 0;         
          W(i,0,k) = W(i,1,k); 
          break;
        case no_slip:
              /* no slip case is velocity = -1.0 * interior velocity at boundary */
          U(i,0,k) = -1.0 * U(i,1,k);  
          V(i,0,k) = -1.0 * V(i,1,k);  
          W(i,0,k) = -1.0 * W(i,1,k);  
          break;
        case zero_gradient:
              /* zero gradient means that velocity at boundary is equal to interior velocity */
          U(i,0,k) = U(i,1,k);  
          V(i,0,k) = V(i,1,k);  
          W(i,0,k) = W(i,1,k);  
          break;       
        }
      }

      P(i,0,k) = P(i,1,k);
      VOF(i,0,k) = VOF(i,1,k);
         
      /* north boundary */
      if(vof_boundaries_check_inside_sb(solver, i, k, 3) == wall) { 
        switch(solver->mesh->wb[3]) {
        case slip:
              /* slip case is zero_gradient for parallel axis
               * and zero fixed value for perpendicular axis */
          U(i,JMAX-1,k) = U(i,JMAX-2,k);
          V(i,JMAX-1,k) = 0;
          V(i,JMAX-2,k) = V(i,JMAX-1,k);
          W(i,JMAX-1,k) = W(i,JMAX-2,k);
          break;
        case no_slip:
              /* no slip case is velocity = -1.0 * interior velocity at boundary */
          U(i,JMAX-1,k) = -1.0 * U(i,JMAX-2,k);
          V(i,JMAX-1,k) = -1.0 * V(i,JMAX-3,k);
          V(i,JMAX-2,k) = V(i,JMAX-1,k);
          W(i,JMAX-1,k) = -1.0 * W(i,JMAX-2,k);
          break;
        case zero_gradient:
              /* zero gradient means that velocity at boundary is equal to interior velocity */
          U(i,JMAX-1,k) = U(i,JMAX-2,k);
          V(i,JMAX-1,k) = V(i,JMAX-3,k);
          V(i,JMAX-2,k) = V(i,JMAX-1,k);
          W(i,JMAX-1,k) = W(i,JMAX-2,k);
          break;       
        }
      }   
      
      P(i,JMAX-1,k) = P(i,JMAX-2,k);
      VOF(i,JMAX-1,k) = VOF(i,JMAX-2,k);
      
    } 
  }
  
    
  /* boundaries on z-axis */
  for(i=0; i<IRANGE; i++) {
    for(j=0; j<JMAX; j++) {
         
      /* bottom boundary */
      if(vof_boundaries_check_inside_sb(solver, i, j, 4) == wall) { 
        switch(solver->mesh->wb[4]) {
        case slip:
              /* slip case is zero_gradient for parallel axis
               * and zero fixed value for perpendicular axis */
          U(i,j,0) = U(i,j,1); 
          V(i,j,0) = V(i,j,1);  
          W(i,j,0) = 0;         
          break;
        case no_slip:
              /* no slip case is velocity = -1.0 * interior velocity at boundary */
          U(i,j,0) = -1.0 * U(i,j,1);  
          V(i,j,0) = -1.0 * V(i,j,1);  
          W(i,j,0) = -1.0 * W(i,j,1);  
          break;
        case zero_gradient:
              /* zero gradient means that velocity at boundary is equal to interior velocity */
          U(i,j,0) = U(i,j,1);  
          V(i,j,0) = V(i,j,1);  
          W(i,j,0) = W(i,j,1);  
          break;       
        }
      }

      P(i,j,0) = P(i,j,1);
      VOF(i,j,0) = VOF(i,j,1);
      
      /* top boundary */
      if(vof_boundaries_check_inside_sb(solver, i, j, 5) == wall) { 
        switch(solver->mesh->wb[5]) {
        case slip:
              /* slip case is zero_gradient for parallel axis
               * and zero fixed value for perpendicular axis */
          U(i,j,KMAX-1) = U(i,j,KMAX-2);
          V(i,j,KMAX-1) = V(i,j,KMAX-2);
          W(i,j,KMAX-1) = 0;
          W(i,j,KMAX-2) = W(i,j,KMAX-1);
          break;
        case no_slip:
              /* no slip case is velocity = -1.0 * interior velocity at boundary */
          U(i,j,KMAX-1) = -1.0 * U(i,j,KMAX-2);
          V(i,j,KMAX-1) = -1.0 * V(i,j,KMAX-2);
          W(i,j,KMAX-1) = -1.0 * W(i,j,KMAX-3);
          W(i,j,KMAX-2) = W(i,j,KMAX-1);
          break;
        case zero_gradient:
              /* zero gradient means that velocity at boundary is equal to interior velocity */
          U(i,j,KMAX-1) = U(i,j,KMAX-2);
          V(i,j,KMAX-1) = V(i,j,KMAX-2);
          W(i,j,KMAX-1) = W(i,j,KMAX-3);
          W(i,j,KMAX-2) = W(i,j,KMAX-1);
          break;       
        }
      }     
      
      P(i,j,KMAX-1) = P(i,j,KMAX-2);
      VOF(i,j,KMAX-1) = VOF(i,j,KMAX-2);
      
    } 
  }
  
  vof_baffles(solver);
        
  /* Free surface and sloped boundary conditions */

  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        if(FV(i,j,k) < solver->emf) {
          
          VOF(i,j,k) = 0;
          P(i,j,k) = 0;
          
          bm[0]=1.0; if (FV(i+1,j,k) < solver->emf) bm[0]=0.0;
          bm[1]=1.0; if (FV(i-1,j,k) < solver->emf) bm[1]=0.0;
          bm[2]=1.0; if (FV(i,j+1,k) < solver->emf) bm[2]=0.0;
          bm[3]=1.0; if (FV(i,j-1,k) < solver->emf) bm[3]=0.0;
          bm[4]=1.0; if (FV(i,j,k+1) < solver->emf) bm[4]=0.0;
          bm[5]=1.0; if (FV(i,j,k-1) < solver->emf) bm[5]=0.0;

          bmtot=bm[0]+bm[1]+bm[2]+bm[3]+bm[4]+bm[5];
          if (bmtot > solver->emf) {
            VOF(i,j,k)=(bm[0]*VOF(i+1,j,k) + bm[2]*VOF(i,j+1,k) + bm[4]*VOF(i,j,k+1) + 
                        bm[1]*VOF(i-1,j,k) + bm[3]*VOF(i,j-1,k) + bm[5]*VOF(i,j,k-1))/bmtot;
            P(i,j,k)  =(bm[0]*P(i+1,j,k) + bm[2]*P(i,j+1,k) + bm[4]*P(i,j,k+1) + 
                        bm[1]*P(i-1,j,k) + bm[3]*P(i,j-1,k) + bm[5]*P(i,j,k-1))/bmtot; 
          }
        
          continue;
        }

        nff = N_VOF(i,j,k);
        
        
        if(nff > 0 && nff < 8) { /* code applies to free surface */
      
          switch(nff) {
          case west:
            if(AE(i,j,k)   > 0 && N_VOF(i+1,j,k) != 0) U(i,  j,k) = U(i-1,j,  k);
            if(AN(i,j,k)   > 0 && N_VOF(i,j+1,k) != 0) V(i,  j,k) = V(i-1,j,  k);
            if(AN(i,j-1,k) > 0 && N_VOF(i,j-1,k) != 0) V(i,j-1,k) = V(i-1,j-1,k);
            if(AT(i,j,k)   > 0 && N_VOF(i,j,k+1) != 0) W(i,j,k  ) = W(i-1,j,  k);
            if(AT(i,j,k-1) > 0 && N_VOF(i,j,k-1) != 0) W(i,j,k-1) = W(i-1,j,k-1);
            break;
          case east:
            if(AE(i-1,j,k) > 0 && N_VOF(i-1,j,k) != 0) U(i-1,j,k) = U(i  ,j,  k);
            if(AN(i,j,k)   > 0 && N_VOF(i,j+1,k) != 0) V(i,  j,k) = V(i+1,j,  k);
            if(AN(i,j-1,k) > 0 && N_VOF(i,j-1,k) != 0) V(i,j-1,k) = V(i+1,j-1,k);
            if(AT(i,j,k)   > 0 && N_VOF(i,j,k+1) != 0) W(i,j,k  ) = W(i+1,j,  k);
            if(AT(i,j,k-1) > 0 && N_VOF(i,j,k-1) != 0) W(i,j,k-1) = W(i+1,j,k-1);
            break;
          case south:
            if(AE(i,j,k)   > 0 && N_VOF(i+1,j,k) != 0) U(i,  j,k) = U(i  ,j-1,k);
            if(AE(i-1,j,k) > 0 && N_VOF(i-1,j,k) != 0) U(i-1,j,k) = U(i-1,j-1,k);
            if(AN(i,j,k)   > 0 && N_VOF(i,j+1,k) != 0) V(i,  j,k) = V(i  ,j-1,k);
            if(AT(i,j,k)   > 0 && N_VOF(i,j,k+1) != 0) W(i,j,k  ) = W(i,  j-1,k);
            if(AT(i,j,k-1) > 0 && N_VOF(i,j,k-1) != 0) W(i,j,k-1) = W(i,  j-1,k-1);
            break;          
          case north:
            if(AE(i,j,k)   > 0 && N_VOF(i+1,j,k) != 0) U(i,  j,k) = U(i  ,j+1,k);
            if(AE(i-1,j,k) > 0 && N_VOF(i-1,j,k) != 0) U(i-1,j,k) = U(i-1,j+1,k);
            if(AN(i,j-1,k) > 0 && N_VOF(i,j-1,k) != 0) V(i,j-1,k) = V(i  ,j  ,k);
            if(AT(i,j,k)   > 0 && N_VOF(i,j,k+1) != 0) W(i,j,k  ) = W(i,  j+1,k);
            if(AT(i,j,k-1) > 0 && N_VOF(i,j,k-1) != 0) W(i,j,k-1) = W(i,  j+1,k-1);
            break;  
          case bottom:
            if(AE(i,j,k)   > 0 && N_VOF(i+1,j,k) != 0) U(i,  j,k) = U(i  ,j,  k-1);
            if(AE(i-1,j,k) > 0 && N_VOF(i-1,j,k) != 0) U(i-1,j,k) = U(i-1,j,  k-1);
            if(AN(i,j,k)   > 0 && N_VOF(i,j+1,k) != 0) V(i,  j,k) = V(i  ,j,  k-1);
            if(AN(i,j-1,k) > 0 && N_VOF(i,j-1,k) != 0) V(i,j-1,k) = V(i  ,j-1,k-1);
            if(AT(i,j,k)   > 0 && N_VOF(i,j,k+1) != 0) W(i,j,k  ) = W(i  ,j,  k-1);
            break;
          case top:
            if(AE(i,j,k)   > 0 && N_VOF(i+1,j,k) != 0) U(i,  j,k) = U(i  ,j,  k+1);
            if(AE(i-1,j,k) > 0 && N_VOF(i-1,j,k) != 0) U(i-1,j,k) = U(i-1,j,  k+1);
            if(AN(i,j,k)   > 0 && N_VOF(i,j+1,k) != 0) V(i,  j,k) = V(i  ,j,  k+1);
            if(AN(i,j-1,k) > 0 && N_VOF(i,j-1,k) != 0) V(i,j-1,k) = V(i  ,j-1,k+1);
            if(AT(i,j,k-1) > 0 && N_VOF(i,j,k-1) != 0) W(i,j,k-1) = W(i  ,j,  k);
            break;
          case none:
            break;
          }

        	dA = 0; 
        	
        	if(N_VOF(i+1,j,k) > 7 && AE(i,j,k) > solver->emf) {
        		if(N_VOF(i-1,j,k) > 7) {		
        			U(i,j,k) = (UN(i,j,k) + UN(i-1,j,k))/2;
        			U(i-1,j,k) = U(i,j,k);
        		}
        		else {
        			U(i,j,k) = U(i-1,j,k);
        		}
        		dA += AE(i,j,k) * pow(RDX,2);
        	}
        	
        	if(N_VOF(i-1,j,k) > 7 && AE(i-1,j,k) > solver->emf) {
        		if(N_VOF(i+1,j,k) > 7) {
        			U(i-1,j,k) = (UN(i,j,k) + UN(i-1,j,k))/2;
        			U(i,j,k) = U(i-1,j,k);
        		}
        		else {
        			U(i-1,j,k) = U(i,j,k);
        		}
        		dA += AE(i-1,j,k) * pow(RDX,2);
        	}
        
        	if(N_VOF(i,j+1,k) > 7 && AN(i,j,k) > solver->emf) {
        		if(N_VOF(i,j-1,k) > 7) {		
        			V(i,j,k) = (VN(i,j,k) + VN(i,j-1,k))/2;
        			V(i,j-1,k) = V(i,j,k);
        		}
        		else {
        			V(i,j,k) = V(i,j-1,k);
        		}
        		dA += AN(i,j,k) * pow(RDY,2);
        	}
        	
        	if(N_VOF(i,j-1,k) > 7 && AN(i,j-1,k) > solver->emf) {
        		if(N_VOF(i,j+1,k) > 7) {
        			V(i,j-1,k) = (VN(i,j,k) + VN(i,j-1,k))/2;
        			V(i,j,k) = V(i,j-1,k);
        		}
        		else {
        			V(i,j-1,k) = V(i,j,k);
        		}
        		dA += AN(i,j-1,k) * pow(RDY,2);
        	}
        	
        
        	if(N_VOF(i,j,k+1) > 7 && AT(i,j,k) > solver->emf) {
        		if(N_VOF(i,j,k-1) > 7) {		
        			W(i,j,k) = (WN(i,j,k) + WN(i,j,k-1))/2;
        			W(i,j,k-1) = W(i,j,k);
        		}
        		else {
        			W(i,j,k) = W(i,j,k-1);
        		}
        		dA += AT(i,j,k) * pow(RDZ,2);
        	}
        	
        	if(N_VOF(i,j,k-1) > 7 && AT(i,j,k-1) > solver->emf) {
        		if(N_VOF(i,j,k+1) > 7) {
        			W(i,j,k-1) = (WN(i,j,k) + WN(i,j,k-1))/2;
        			W(i,j,k) = W(i,j,k-1);
        		}
        		else {
        			W(i,j,k-1) = W(i,j,k);
        		}
        		dA += AT(i,j,k-1) * pow(RDZ,2);
        	}
        	
        	if(dA > 0) {
					
						dv  = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k)) +
								RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k)) +
								RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1));

						delp = -1.0 * dv * solver->rho / (dA * solver->delt);
						delp *= VOF(i,j,k);

						//if(solver->iter > 0) P(i,j,k) += delp;

						if(N_VOF(i+1,j,k) > 7 && AE(i,j,k) > solver->emf) 
							U(i,j,k)=U(i,j,k) + solver->delt* RDX * delp / (solver->rho);
						
						if(N_VOF(i-1,j,k) > 7 && AE(i-1,j,k) > solver->emf) 
							U(i-1,j,k)=U(i-1,j,k) - solver->delt* RDX * delp / (solver->rho);

						if(N_VOF(i,j+1,k) > 7 && AN(i,j,k) > solver->emf) 
							V(i,j,k)=V(i,j,k) + solver->delt* RDY * delp / (solver->rho);
						
						if(N_VOF(i,j-1,k) > 7 && AN(i,j-1,k) > solver->emf) 
							V(i,j-1,k)=V(i,j-1,k) - solver->delt* RDY * delp / (solver->rho);

						if(N_VOF(i,j,k+1) > 7 && AT(i,j,k) > solver->emf) 
							W(i,j,k)=W(i,j,k) + solver->delt* RDZ * delp / (solver->rho);
						
						if(N_VOF(i,j,k-1) > 7 && AT(i,j,k-1) > solver->emf) 
							W(i,j,k-1)=W(i,j,k-1) - solver->delt* RDZ * delp / (solver->rho);
					} 


#define emf solver->emf
   /* # set velocities in empty cells adjacent to partial fluid cells */
          if(solver->iter==0) {
        
            if(VOF(i+1,j,k) < emf) {
              if(VOF(i+1,j+1,k) < emf && AN(i+1,j,k) > emf)
                V(i+1,j,k) = VOF(i,j,k) * V(i,j,k);
              if(VOF(i+1,j-1,k) < emf && AN(i+1,j-1,k) > emf)
                V(i+1,j-1,k) = VOF(i,j,k) * V(i,j-1,k);
              
              if(VOF(i+1,j,k+1) < emf && AT(i+1,j,k) > emf)
                W(i+1,j,k) = VOF(i,j,k) * W(i,j,k);
              if(VOF(i+1,j,k-1) < emf && AT(i+1,j,k-1) > emf)
                W(i+1,j,k-1) = VOF(i,j,k) * W(i,j,k-1);
            }
           
            if(VOF(i-1,j,k) < emf) {
              if(VOF(i-1,j+1,k) < emf && AN(i-1,j,k) > emf)
                V(i-1,j,k) = VOF(i,j,k) * V(i,j,k);
              if(VOF(i-1,j-1,k) < emf && AN(i-1,j-1,k) > emf)
                V(i-1,j-1,k) = VOF(i,j,k) * V(i,j-1,k);
              
              if(VOF(i-1,j,k+1) < emf && AT(i-1,j,k) > emf)
                W(i-1,j,k) = VOF(i,j,k) * W(i,j,k);
              if(VOF(i-1,j,k-1) < emf && AT(i-1,j,k-1) > emf)
                W(i-1,j,k-1) = VOF(i,j,k) * W(i,j,k-1);
            }
                           
            if(VOF(i,j+1,k) < emf) {
              if(VOF(i+1,j+1,k) < emf && AE(i,j+1,k) > emf)
                U(i,j+1,k) = VOF(i,j,k) * U(i,j,k);
              if(VOF(i-1,j+1,k) < emf && AE(i-1,j+1,k) > emf)
                U(i-1,j+1,k) = VOF(i,j,k) * U(i-1,j,k);

              if(VOF(i,j+1,k+1) < emf && AT(i,j+1,k) > emf)
                W(i,j+1,k) = VOF(i,j,k) * W(i,j,k);
              if(VOF(i,j+1,k-1) < emf && AT(i-1,j+1,k) > emf)
                W(i,j+1,k-1) = VOF(i,j,k) * W(i-1,j,k);
            }

            if(VOF(i,j-1,k) < emf)
            {
              if(VOF(i+1,j-1,k) < emf  &&  AE(i,j-1,k) > emf)
                U(i,j-1,k) = VOF(i,j,k) * U(i,j,k);
              if(VOF(i-1,j-1,k) < emf  &&  AE(i-1,j-1,k) > emf)
                U(i-1,j-1,k) = VOF(i,j,k) * U(i-1,j,k);
             
              if(VOF(i,j-1,k+1) < emf  &&  AT(i,j-1,k) > emf)
                W(i,j-1,k) = VOF(i,j,k) * W(i,j,k);
              if(VOF(i,j-1,k-1) < emf  &&  AT(i,j-1,k-1) > emf)
                W(i,j-1,k-1) = VOF(i,j,k) * W(i,j,k-1);    
            
            }
            if(VOF(i,j,k+1) < emf)
            {
              if(VOF(i+1,j,k+1) < emf && AE(i,j,k+1) > emf)
                U(i,j,k+1) = VOF(i,j,k) * U(i,j,k);
              if(VOF(i-1,j,k+1) < emf && AE(i-1,j,k+1) > emf)
                U(i-1,j,k+1) = VOF(i,j,k) * U(i-1,j,k);
              
              if(VOF(i,j+1,k+1) < emf && AN(i,j,k+1) > emf)
                V(i,j,k+1) = VOF(i,j,k) * V(i,j,k);
              if(VOF(i,j-1,k+1) < emf && AN(i,j-1,k+1) > emf)
                V(i,j-1,k+1) = VOF(i,j,k) * V(i,j-1,k);
            }
            if(VOF(i,j,k-1) < emf)
            {
              if(VOF(i+1,j,k-1) < emf && AE(i,j,k-1) > emf)
                U(i,j,k-1) = VOF(i,j,k) * U(i,j,k);
              if(VOF(i-1,j,k-1) < emf && AE(i-1,j,k-1) > emf)
                U(i-1,j,k-1) = VOF(i,j,k) * U(i-1,j,k);
              
              if(VOF(i,j+1,k-1) < emf && AN(i,j,k-1) > emf)
                V(i,j,k-1) = VOF(i,j,k) * V(i,j,k);
              if(VOF(i,j-1,k-1) < emf && AN(i,j-1,k-1) > emf)
                V(i,j-1,k-1) = VOF(i,j,k) * V(i,j-1,k);
            } 
          } 
        }
      }
    }
  }

  vof_baffles(solver);
  
  return 0;
#undef emf
}
