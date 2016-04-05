/* vof_boundary.c
 *
 * implement boundaries
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "vtk.h"
#include "vof.h"
#include "solver.h"
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
  int i;
  struct sb_data *sb;
  
  
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
        
        if(FV(i,j,k) < 0.000001) continue;
        
        if(solver->p_flag == 0) {
          /* must set velocity to a Von Neumann boundary */
          
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
        VOF(i+coplanar[0],j+coplanar[1],k+coplanar[2]) = VOF(i,j,k); 
     
        if(solver->p_flag != 0) continue; 
        
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
        /* ADDED 01/03/2014 */
        else if(mesh->vof[mesh_index(mesh,i,j,k)] <= 0.0) 
          mesh->P[mesh_index(mesh,i,j,k)] = 0.0;
          
        P(i+coplanar[0],j+coplanar[1],k+coplanar[2]) = P(i,j,k); 
        /* now set value 1 cell interior to the mesh to be the same */
        /* P(i+coplanar[0],j+coplanar[1],k+coplanar[2]) = (P(i,j,k) + P(i+2*coplanar[0],j+2*coplanar[1],k+2*coplanar[2]))/2; */

        
        /* ADDED 02/27/2016 *
        if(FV(i,j,k) < 1.0 && FV(i,j,k) > 0.0) {
          if(i + coplanar[0] < IMAX-1 && j + coplanar[0] < JMAX - 1) {
            if(FV(i+coplanar[0],j+coplanar[0],k)
            P(i,j,k) = (P(i,j,k) + P( 
        }*/

     }
    }
  }

/*  if(solver->p_flag != 0) return;*/

  /* UNCOMMENT BELOW TO INCLUDE VELOCITY HEAD IN THE CALC, 
     CHANGING THIS BOUNDARY FROM HGL TO EGL */

  /* now add the velocity pressure
     the boundary assumes stagnation if upstream
     and static pressure if downstream  *
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k = kmin; k <= kmax; k++) {
        if(mesh->vof[mesh_index(mesh,i,j,k)] <= 0.0) continue;
        if(FV(i,j,k) < 0.000001) continue;
        
        if((U(i,j,k) * normal[0] * coplanar[0] + V(i,j,k) * normal[1] * coplanar[1]
           + W(i,j,k) * normal[2] * coplanar[2]) > 0.0) continue; * assume stagnation if upstream *
 
        l = i+coplanar[0];
        m = j+coplanar[1];
        n = k+coplanar[2];
 
        if(N_VOF(l,m,n) == 0) {
          rho_v2 = solver->rho *     (pow(U(i,j,k),2) * normal[0] + 
                                      pow(V(i,j,k),2) * normal[1] + 
                                      pow(W(i,j,k),2) * normal[2]) / 2;
          mesh->P[mesh_index(mesh,i,j,k)] += rho_v2;
        } else {
          if(k > 0) {
            l = i;
            m = j;
            n = k-1;
          
            sdis=VOF(i,j,k)*dd+VOF(l,m,n)*dd*0.5;
            if(FV(l,m,n) < emf) sdis = dd*0.5 + VOF(i,j,k) * dd;
            sdis=max(sdis,0.5*dd);
            peta=dd/sdis;

            if ((FV(l,m,n) >= emf)&&(N_VOF(l,m,n)!=0.0)) peta=1.0;
            if (peta>2.0) peta=2.0;
            if (peta<0.0) peta=0.0;
          
            P(i,j,k) = (1 - peta) * P(i,j,k-1);
          }
        }
      }
    }
  } */ /* For simplicity this is being removed.  The boundary is "hgl" which denotes stagnation pressure anyways */

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
  
  #pragma omp parallel for shared (solver, imin, jmin, kmin, imax, jmax, kmax, value, turbulence) \
              private(i,j,k) collapse(3) schedule(static)
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
        	if(vof_vof_height_boundary(solver))
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

  vof_vof_height_boundary(solver);
  return 0;
#undef emf
}

int vof_vof_height_boundary(struct solver_data *solver) {
  static int done = 0;
  double z_height;
  
  if(done == 1) return 1;
  
  if(solver->vof_delay < solver->emf) {
    done = 1;
    return 1;
  }
  
  if(solver->t < solver->vof_delay) { /* prior to delay, go single phase */
  
    mesh_set_array(solver->mesh, "vof", 1.0, -1, 0, 0, 0, 0, 0);
    
  }

  else { /* after, go multiphase */
  
    mesh_set_array(solver->mesh, "vof", 0.0, -1, 0, 0, 0, 0, 0);
    z_height = floor(solver->vof_height / solver->mesh->delz);

    mesh_set_array(solver->mesh, "vof", 1.0, 0, solver->mesh->imax, 
                                             0, solver->mesh->jmax, 
                                             0, z_height);  
                                             
    if(solver->vof_height / solver->mesh->delz - z_height > 0) {
      mesh_set_array(solver->mesh, "vof", solver->vof_height / solver->mesh->delz - z_height, 
                                             0, solver->mesh->imax, 
                                             0, solver->mesh->jmax, 
                                             z_height, z_height+1);        
    }

		/* Now set this for the previous timestep */
    mesh_set_array(mesh_n, "vof", 0.0, -1, 0, 0, 0, 0, 0);

    mesh_set_array(mesh_n, "vof", 1.0, 0, solver->mesh->imax, 
                                             0, solver->mesh->jmax, 
                                             0, z_height);  
                                             
    if(solver->vof_height / solver->mesh->delz - z_height > 0) {
      mesh_set_array(mesh_n, "vof", solver->vof_height / solver->mesh->delz - z_height, 
                                             0, solver->mesh->imax, 
                                             0, solver->mesh->jmax, 
                                             z_height, z_height+1);        
    }
  
  
    mesh_set_hydrostatic(solver->mesh, fabs(solver->gz), solver->rho);
		solver->delt = 0.001;

#ifdef DEBUG
    printf("VOF Boundary Delay timeout\n");
#endif
  
      done = 1;
    solver->boundaries(solver);
    if(solver->special_boundaries != NULL)
      solver->special_boundaries(solver);
      
    return 1;
  }
  
  return 0;
}

int vof_boundaries(struct solver_data *solver) {

  long int i,j,k,l,m,n;
  int bm[6], bmtot, nindex, flg;
  enum cell_boundaries nff;
#define dim(i,j,k) i+3*(j+k*3)
  const int ndim[3][3] = { {  0,1,1 }, { 1,0,1 }, { 1,1,0 } };
  const int odim[3][3] = { {  1,0,0 }, { 0,1,0 }, { 0,0,1 } };
  double denom, dv;
  
  
  
  /* first boundaries on x-axis */
  for(j=0; j<JMAX; j++) {
    for(k=0; k<KMAX; k++) {

      /* west boundary */
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
          if(solver->p_flag == 0) {
              /* zero gradient means that velocity at boundary is equal to interior velocity */
            U(0,j,k) = U(1,j,k);  
            V(0,j,k) = V(1,j,k);  
            W(0,j,k) = W(1,j,k);
          }
          break;       
        }
      }
      
      if(solver->p_flag == 0) {
        P(0,j,k) = P(1,j,k);
      }
      VOF(0,j,k) = VOF(1,j,k);
      
      /* east boundary */
      if(vof_boundaries_check_inside_sb(solver, j, k, 1) == wall) {    
        switch(solver->mesh->wb[1]) {
        case slip:
              /* slip case is zero_gradient for parallel axis
               * and zero fixed value for perpendicular axis */
          U(IMAX-1,j,k) = 0;
          U(IMAX-2,j,k) = U(IMAX-1,j,k);
          V(IMAX-1,j,k) = V(IMAX-2,j,k);
          W(IMAX-1,j,k) = W(IMAX-2,j,k);
          break;
        case no_slip:
              /* no slip case is velocity = -1.0 * interior velocity at boundary */
          U(IMAX-1,j,k) = -1.0 * U(IMAX-3,j,k);
          U(IMAX-2,j,k) = U(IMAX-1,j,k);
          V(IMAX-1,j,k) = -1.0 * V(IMAX-2,j,k);
          W(IMAX-1,j,k) = -1.0 * W(IMAX-2,j,k);
          break;
        case zero_gradient:
          if(solver->p_flag == 0) {
              /* zero gradient means that velocity at boundary is equal to interior velocity */
            U(IMAX-1,j,k) = U(IMAX-3,j,k);
            U(IMAX-2,j,k) = U(IMAX-1,j,k);
            V(IMAX-1,j,k) = V(IMAX-2,j,k);
            W(IMAX-1,j,k) = W(IMAX-2,j,k);
          }
          break;       
        }
      }    
      
      if(solver->p_flag == 0) {
        P(IMAX-1,j,k) = P(IMAX-2,j,k);
      }
      VOF(IMAX-1,j,k) = VOF(IMAX-2,j,k);
      
    } 
  }
  
  /* boundaries on y-axis */
  for(i=0; i<IMAX; i++) {
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
          if(solver->p_flag == 0) {
            U(i,0,k) = U(i,1,k);  
            V(i,0,k) = V(i,1,k);  
            W(i,0,k) = W(i,1,k);  
          }
          break;       
        }
      }
      if(solver->p_flag == 0) {
        P(i,0,k) = P(i,1,k);
      }
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
          if(solver->p_flag == 0) {
            U(i,JMAX-1,k) = U(i,JMAX-2,k);
            V(i,JMAX-1,k) = V(i,JMAX-3,k);
            V(i,JMAX-2,k) = V(i,JMAX-1,k);
            W(i,JMAX-1,k) = W(i,JMAX-2,k);
          }
          break;       
        }
      }   
      
      if(solver->p_flag == 0) {
        P(i,JMAX-1,k) = P(i,JMAX-2,k);
      }
      VOF(i,JMAX-1,k) = VOF(i,JMAX-2,k);
      
    } 
  }
  
    
  /* boundaries on z-axis */
  for(i=0; i<IMAX; i++) {
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
          if(solver->p_flag == 0) {
            U(i,j,0) = U(i,j,1);  
            V(i,j,0) = V(i,j,1);  
            W(i,j,0) = W(i,j,1);  
          }
          break;       
        }
      }
      if(solver->p_flag == 0) {
        P(i,j,0) = P(i,j,1);
      }
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
          if(solver->p_flag == 0) {
            U(i,j,KMAX-1) = U(i,j,KMAX-2);
            V(i,j,KMAX-1) = V(i,j,KMAX-2);
            W(i,j,KMAX-1) = W(i,j,KMAX-3);
            W(i,j,KMAX-2) = W(i,j,KMAX-1);
          }
          break;       
        }
      }     
      
      if(solver->p_flag == 0) {
        P(i,j,KMAX-1) = P(i,j,KMAX-2);
      }
      VOF(i,j,KMAX-1) = VOF(i,j,KMAX-2);
      
    } 
  }
  
  vof_vof_height_boundary(solver);
  vof_baffles(solver);
  if(solver->p_flag != 0) return 0;      
        
  /* Free surface and sloped boundary conditions */
  #pragma omp parallel for shared (solver, ndim, odim) private(i,j,k,l,m,n,bm,bmtot,nindex,nff,denom) \
            collapse(3) schedule(dynamic, 100)
  for(i=1; i<IMAX-1; i++) {
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
        
          /* REMOVED 08/13/2014
           * NOT SURE WHY THIS WAS EVER INCLUDED *
          flg=0;
          bmtot = 0.0;
          N_VOF(i,j,k)=0;
          P(i,j,k)=0;
          
          for(x=0; x<3; x++) {
          
            bm[x*2] = 0;
            l = i + odim[x][0];
            m = j + odim[x][1];
            n = k + odim[x][2];
            
            if(FV(l,m,n) > solver->emf) {
              /* check for free surface, if so, use neighbor cell as basis *
              if(VOF(l,m,n) > solver->emf && VOF(l,m,n) < solver->emf_c) {

                if(N_VOF(l,m,n) == (x+1)*2) /*west/south/bottom* {
                  VOF(i,j,k) = 1.0;
                  flg=1;
                }
                else if(N_VOF(l,m,n) == (x*2)+1) /*east/north/top* {
                  VOF(i,j,k) = 0.0;
                  flg=1;
                }
                else if(N_VOF(l,m,n) > 0 || N_VOF(l,m,n) < 7) {
                  VOF(i,j,k) = VOF(l,m,n);
                  flg=1;
                }
              }
              bm[x*2]=1;
            }
            if(flg==1) break;
            
            bm[x*2 + 1] = 0;
            l = i - odim[x][0];
            m = j - odim[x][1];
            n = k - odim[x][2];

            if(FV(l,m,n) > 0) {
              /* check for free surface, if so, use neighbor cell as basis *
              if(VOF(l,m,n) > solver->emf && VOF(l,m,n) < solver->emf_c) {

                if(N_VOF(l,m,n) == (x+1)*2) /*west/south/bottom* {
                  VOF(i,j,k) = 0.0;
                  flg=1;
                }
                else if(N_VOF(l,m,n) == (x*2)+1) /*east/north/top* {
                  VOF(i,j,k) = 1.0;
                  flg=1;
                }
                else if(N_VOF(l,m,n) > 0 || N_VOF(l,m,n) < 7) {
                  VOF(i,j,k) = VOF(l,m,n);
                  flg=1;
                }
              }
              bm[x*2 + 1]=1;
            }                      
            if(flg==1) break;
          
            bmtot += bm[x*2] + bm[x*2 + 1];
          }
          
          if(flg==1) continue;          

          if(bmtot <= 0)
            continue;
          
          VOF(i,j,k) = (bm[0]*VOF(i+1,j,k) + bm[2]*VOF(i,j+1,k) + bm[4]*VOF(i,j,k+1) + 
                        bm[1]*VOF(i-1,j,k) + bm[3]*VOF(i,j-1,k) + bm[5]*VOF(i,j,k-1))/bmtot;
          P(i,j,k)   = (bm[0]*P(i+1,j,k) + bm[2]*P(i,j+1,k) + bm[4]*P(i,j,k+1) + 
                        bm[1]*P(i-1,j,k) + bm[3]*P(i,j-1,k) + bm[5]*P(i,j,k-1))/bmtot;          
          */
          continue;
        }

        /* this code creates a von neumann boundary at the fluid surface
        # ( sort of )
        # it solves the continuity equation for each surface cell
        # and uses that to set the value for U and V
        #
        # if two axis are unknown
        # then the velocity is arbitrarily copied from one axis
        # and continuity is used to solve the other axis
        #
        # this leads me to question: if we are forcing continuity in
        # the free surface cell, should we skip the pressure correction? */

        nff = N_VOF(i,j,k);
        
        
        if(nff > 0 && nff < 8) { /* code applies to free surface */
 #pragma omp critical(free_surf_bdry) 
 {
 
          if(AE(i,j,k)   > emf && N_VOF(i+1,j,k) > 7) U(i,j,k)   = U(i-1,j,k) * AE(i-1,j,k) / AE(i,j,k);
          if(AE(i-1,j,k) > emf && N_VOF(i-1,j,k) > 7) U(i-1,j,k) = U(i,j,k)   * AE(i,j,k)   / AE(i-1,j,k); 
          if(AN(i,j,k)   > emf && N_VOF(i,j+1,k) > 7) V(i,j,k)   = V(i,j-1,k) * AN(i,j-1,k) / AN(i,j,k);
          if(AN(i,j-1,k) > emf && N_VOF(i,j-1,k) > 7) V(i,j-1,k) = V(i,j,k)   * AN(i,j,k)   / AN(i,j-1,k);
          if(AT(i,j,k)   > emf && N_VOF(i,j,k+1) > 7) W(i,j,k)   = W(i,j,k-1) * AT(i,j,k-1) / AT(i,j,k);
          if(AT(i,j,k-1) > emf && N_VOF(i,j,k-1) > 7) W(i,j,k-1) = W(i,j,k)   * AT(i,j,k)   / AT(i,j,k-1);

          flg = 0;
          nindex = 0;
          nff = N_VOF(i,j,k);
          dv  = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k)) +
              RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k)) +
              RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1));
          
          while(flg == 0) {
          
            switch(nff) {
            case west:
              if(N_VOF(i+1,j,k) > 7 && AE(i,j,k) > emf) {
                U(i,j,k) = U(i,j,k) - DELX * dv / AE(i,j,k);
                flg = 1;
                dv  = 0;
              }
              break;
            case east:
              if(N_VOF(i-1,j,k) > 7 && AE(i-1,j,k) > emf) {
                U(i-1,j,k) = U(i-1,j,k) - DELX * dv / AE(i-1,j,k);
                flg = 1;
                dv  = 0;
              }
              break;
            case south:
              if(N_VOF(i,j+1,k) > 7 && AN(i,j,k) > emf) {
                V(i,j,k) = V(i,j,k) - DELY * dv / AN(i,j,k);
                flg = 1;
                dv  = 0;
              }
              break;
            case north:
              if(N_VOF(i,j-1,k) > 7 && AN(i,j-1,k) > emf) {
                V(i,j-1,k) = V(i,j-1,k) - DELY * dv / AN(i,j-1,k);
                flg = 1;
                dv  = 0;
              }
              break;
            case bottom:
              if(N_VOF(i,j,k+1) > 7 && AT(i,j,k) > emf) {
                W(i,j,k) = W(i,j,k) - DELZ * dv / AT(i,j,k);
                flg = 1;
                dv  = 0;
              }
              break;
            case top:
              if(N_VOF(i,j,k-1) > 7 && AT(i,j,k-1) > emf) {
                W(i,j,k-1) = W(i,j,k-1) - DELZ * dv / AT(i,j,k-1);
                flg = 1;
                dv  = 0;
              }
              break;
            }
            
            if(flg == 0) {
              nff++;
              if(nff > 6) nff = 1;
            }
            
            nindex++;
            if(nindex > 6) flg = 1;
            
          }
 
/* OLD METHOD FROM 3dVOF - make sure to add that p_flag == 0 to run this
          switch(nff) {
          case west:
            if(AE(i,j,k)   > 0 && N_VOF(i+1,j,k) != 0) U(i,  j,k) = U(i-1,j,  k);
            if(AN(i,j,k)   > 0) V(i,  j,k) = V(i-1,j,  k);
            if(AN(i,j-1,k) > 0) V(i,j-1,k) = V(i-1,j-1,k);
            if(AT(i,j,k)   > 0) W(i,j,k  ) = W(i-1,j,  k);
            if(AT(i,j,k-1) > 0) W(i,j,k-1) = W(i-1,j,k-1);
            break;
          case east:
            if(AE(i-1,j,k) > 0 && N_VOF(i-1,j,k) != 0) U(i-1,j,k) = U(i  ,j,  k);
            if(AN(i,j,k)   > 0) V(i,  j,k) = V(i+1,j,  k);
            if(AN(i,j-1,k) > 0) V(i,j-1,k) = V(i+1,j-1,k);
            if(AT(i,j,k)   > 0) W(i,j,k  ) = W(i+1,j,  k);
            if(AT(i,j,k-1) > 0) W(i,j,k-1) = W(i+1,j,k-1);
            break;
          case south:
            if(AE(i,j,k)   > 0) U(i,  j,k) = U(i  ,j-1,k);
            if(AE(i-1,j,k) > 0) U(i-1,j,k) = U(i-1,j-1,k);
            if(AN(i,j,k)   > 0 && N_VOF(i,j+1,k) != 0) V(i,  j,k) = V(i  ,j-1,k);
            if(AT(i,j,k)   > 0) W(i,j,k  ) = W(i,  j-1,k);
            if(AT(i,j,k-1) > 0) W(i,j,k-1) = W(i,  j-1,k-1);
            break;          
          case north:
            if(AE(i,j,k)   > 0) U(i,  j,k) = U(i  ,j+1,k);
            if(AE(i-1,j,k) > 0) U(i-1,j,k) = U(i-1,j+1,k);
            if(AN(i,j-1,k) > 0 && N_VOF(i,j-1,k) != 0) V(i,j-1,k) = V(i  ,j  ,k);
            if(AT(i,j,k)   > 0) W(i,j,k  ) = W(i,  j+1,k);
            if(AT(i,j,k-1) > 0) W(i,j,k-1) = W(i,  j+1,k-1);
            break;  
          case bottom:
            if(AE(i,j,k)   > 0) U(i,  j,k) = U(i  ,j,  k-1);
            if(AE(i-1,j,k) > 0) U(i-1,j,k) = U(i-1,j,  k-1);
            if(AN(i,j,k)   > 0) V(i,  j,k) = V(i  ,j,  k-1);
            if(AN(i,j-1,k) > 0) V(i,j-1,k) = V(i  ,j-1,k-1);
            if(AT(i,j,k)   > 0 && N_VOF(i,j,k+1) != 0) W(i,j,k  ) = W(i  ,j,  k-1);
            break;
          case top:
            if(AE(i,j,k)   > 0) U(i,  j,k) = U(i  ,j,  k+1);
            if(AE(i-1,j,k) > 0) U(i-1,j,k) = U(i-1,j,  k+1);
            if(AN(i,j,k)   > 0) V(i,  j,k) = V(i  ,j,  k+1);
            if(AN(i,j-1,k) > 0) V(i,j-1,k) = V(i  ,j-1,k+1);
            if(AT(i,j,k-1) > 0 && N_VOF(i,j,k-1) != 0) W(i,j,k-1) = W(i  ,j,  k);
            break;
          case none:
            break;
          } 
           
          nindex = 0;
          while(nindex<8) {
            nindex++;

            switch(nff) {
            case west:
              if(N_VOF(i+1,j,k) > 7 && AE(i,j,k) > 0) {
                denom = -1.0 * RDX * AE(i,j,k);
                U(i,j,k) = ( RDX * (-1.0 * U(i-1,j,k) * AE(i-1,j,k))  +
                             RDY * (V(i,j,k) * AN(i,j,k) - V(i,j-1,k) * AN(i,j-1,k)) +
                             RDZ * (W(i,j,k) * AT(i,j,k) - W(i,j,k-1) * AT(i,j,k-1)) ) / denom;
              }
              break;
            case east:
              if(N_VOF(i-1,j,k) > 7 && AE(i-1,j,k) > 0) {
                denom = RDX * AE(i-1,j,k);
                U(i-1,j,k) = ( RDX * (U(i,j,k) * AE(i,j,k)) +
                               RDY * (V(i,j,k) * AN(i,j,k) - V(i,j-1,k) * AN(i,j-1,k)) +
                               RDZ * (W(i,j,k) * AT(i,j,k) - W(i,j,k-1) * AT(i,j,k-1)) ) / denom;
              }
              break;
            case south:
              if(N_VOF(i,j+1,k) > 7 && AN(i,j,k) > 0) {
                denom = -1.0 * RDY * AN(i,j,k);
                V(i,j,k) = ( RDX * (U(i,j,k) * AE(i,j,k) - U(i-1,j,k) * AE(i-1,j,k))  +
                             RDY * (-1.0 * V(i,j-1,k) * AN(i,j-1,k)) +
                             RDZ * (W(i,j,k) * AT(i,j,k) - W(i,j,k-1) * AT(i,j,k-1)) ) / denom;
              }
              break;
            case north:
              if(N_VOF(i,j-1,k) > 7 && AN(i,j-1,k) > 0) {
                denom = RDY * AN(i,j-1,k);
                V(i,j-1,k) = ( RDX * (U(i,j,k) * AE(i,j,k) - U(i-1,j,k) * AE(i-1,j,k)) +
                               RDY * (V(i,j,k) * AN(i,j,k)) +
                               RDZ * (W(i,j,k) * AT(i,j,k) - W(i,j,k-1) * AT(i,j,k-1)) ) / denom;
              }
              break;                
            case bottom:
              if(N_VOF(i,j,k+1) > 7 && AT(i,j,k) > 0) {
                denom = -1.0 * RDZ * AT(i,j,k);
                W(i,j,k) = ( RDX * (U(i,j,k) * AE(i,j,k) - U(i-1,j,k) * AE(i-1,j,k)) +
                             RDY * (V(i,j,k) * AN(i,j,k) - V(i,j-1,k) * AN(i,j-1,k)) +
                             RDZ * (-1.0 * W(i,j,k-1) * AT(i,j,k-1)) ) / denom;
              }
              break;
            case top:
              if(N_VOF(i,j,k-1) > 7 && AT(i,j,k-1) > 0) {
                denom = RDZ * AT(i,j,k-1);
                W(i,j,k-1) = ( RDX * (U(i,j,k) * AE(i,j,k) - U(i-1,j,k) * AE(i-1,j,k)) +
                               RDY * (V(i,j,k) * AN(i,j,k) - V(i,j-1,k) * AN(i,j-1,k)) +
                               RDZ * (W(i,j,k) * AT(i,j,k)) ) / denom;
              }
              break;  
            case none:
              break;               
            }
              
            nff++;
            if(nff>7) nff = 1; 
          } */
          
#define emf solver->emf

   /* # set velocities in empty cells adjacent to partial fluid cells */
          if(solver->p_flag==0 || solver->iter==0) {
        
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
  }

  vof_vof_height_boundary(solver);
  vof_baffles(solver);
  
  return 0;
#undef emf
}

