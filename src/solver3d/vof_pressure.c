/* vof_pressure.c
 *
 * implement pressure iteration
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
#include "track.h"

#include "vof_macros.h"

struct mesh_data *mesh_n; /* describes mesh at previous timestep for explicit calcs */

static double *up, *um, *vp, *vm, *wp, *wm;

int vof_pressure_init(struct solver_data *solver) {
  long int size;

  size = solver->mesh->imax * solver->mesh->jmax * solver->mesh->kmax;
  
  up = malloc(sizeof(double) * size);
  um = malloc(sizeof(double) * size);
  if(up == NULL || um == NULL) {
    printf("error: memory could not be allocated in vof_pressure_init\n");
    return(1);
  }
  
  vp = malloc(sizeof(double) * size);
  vm = malloc(sizeof(double) * size);
  if(vp == NULL || vm == NULL) {
    printf("error: memory could not be allocated in vof_pressure_init\n");
    return(1);
  }  
  
  wp = malloc(sizeof(double) * size);
  wm = malloc(sizeof(double) * size);
  if(wp == NULL || wm == NULL) {
    printf("error: memory could not be allocated in vof_pressure_init\n");
    return(1);
  }
  
  return 0;
  
}

/* iterate through the mesh and check for overpressurization.  relax it and report it to the user */
int vof_pressure_test(struct solver_data *solver) {
  long int i,j,k,l,m,n,x;
  double dp, dv, stabil_limit, ax, del, ux, gx, hydrostatic;
  
  
#define emf solver->emf

#pragma omp parallel for shared (solver) private(i,j,k,l,m,n,x,dp,dv,stabil_limit,ax,del,ux,gx,hydrostatic) \
            collapse(3) schedule(dynamic, 100)
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
      
      	for(x=0; x<6; x++) {
      		
      		l = i; m = j; n = k;
      		
      		switch(x) {
      		case 0:
      			l = i+1;
      			ux = U(i,j,k);
      			ax = AE(i,j,k);
      			del = DELX;
      			gx = solver->gx;
      			break;
      		case 1:
      			l = i-1;
      			ux = U(l,m,n);
      			ax = AE(l,m,n);
      			del = DELX;
      			gx = solver->gx;
      			break;
      		case 2: 
      			m = j+1;
      			ux = V(i,j,k);
      			ax = AN(i,j,k);
      			del = DELY;
      			gx = solver->gy;
      			break;
      		case 3:
      			m = j-1;
      			ux = V(l,m,n);
      			ax = AN(l,m,n);
      			del = DELY;
      			gx = solver->gy;
      			break;
      		case 4:
      			n = k+1;
      			ux = W(i,j,k);
      			ax = AT(i,j,k);
      			del = DELZ;
      			gx = solver->gz;
      			break;
      		case 5:
      			n = k-1;
      			ux = W(l,m,n);
      			ax = AT(l,m,n);
      			del = DELZ;
      			gx = solver->gz;
      			break;
      		default:
      			continue;      		
      		}
      		
      		if(VOF(i,j,k) < emf || VOF(l,m,n) < emf) continue;
      		/* if(N_VOF(i,j,k) != 0 || N_VOF(l,m,n) != 0) continue; */
      		if(FV(i,j,k) < emf || FV(l,m,n) < emf) continue;	
						
					dp = P(i,j,k) - P(l,m,n);
					
					hydrostatic = gx * del * solver->rho * (l-i + m-j + n-k);
					dp += hydrostatic;
					
					if(dp < emf) continue;
					
          dv = solver->delt * ( (1/del) * dp / solver->rho);
					
					stabil_limit = solver->con * (min(FV(i,j,k),FV(l,m,n)) / ax)  * del/solver->delt;
			
					if((l-i + m-j + n-k) * ux > emf)
						stabil_limit = min(max(stabil_limit - fabs(ux), emf), stabil_limit);
      				
      		if(fabs(dv) > stabil_limit) {
      			printf("warning: stability limit violation in cell %ld %ld %ld to %ld %ld %ld\ncorrecting dp from %lf to ", i, j, k, l, m, n, dp);
      		
            dp = stabil_limit / (solver->delt * (1/del) / solver->rho);
            dp -= hydrostatic;
            
            printf("%lf\n", dp);
            
            P(i,j,k) = P(l,m,n) + dp - emf;
      		}
      			

      	}
      
      }
    }
  }
  
  return 0;
}


int vof_pressure_mp(struct solver_data *solver) {
  
  long int i,j,k, imax, jmax, kmax, imin, jmin, kmin;
  
  
  static int direction = -1;
  
  if(direction == -1) {
    imin = 1;
    jmin = 1;
    kmin = 1;
    imax = IMAX-1;
    jmax = JMAX-1;
    kmax = KMAX-1;
    direction = 1;
  } else {
    imax = 0;
    jmax = 0;
    kmax = 0;
    imin = IMAX-2;
    jmin = JMAX-2;
    kmin = KMAX-2;
    direction = -1;
  }  
  
  
  
  solver->p_flag = 0;
  
#define emf solver->emf

#pragma omp parallel for shared (solver) private(i,j,k) \
            collapse(3) schedule(dynamic, 100)
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {  
      
        vof_pressure_sor(solver,i,j,k);
      }
    }
  }
  
 /*
#pragma omp parallel for shared (solver) private(i,j,k) \
            collapse(3) schedule(dynamic, 100) 
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
        
        if(N_VOF(i,j,k) == none) continue;
        
        if(AE(i,j,k) > emf && i < IMAX-1)
          U(i,j,k) += up[mesh_index(solver->mesh,i,j,k)];
          
        if(AN(i,j,k) > emf && j < JMAX-1)
          V(i,j,k) += vp[mesh_index(solver->mesh,i,j,k)];
          
        if(AT(i,j,k) > emf && k < KMAX-1)
          W(i,j,k) += wp[mesh_index(solver->mesh,i,j,k)];

      }
    } 
  }

#pragma omp parallel for shared (solver) private(i,j,k) \
            collapse(3) schedule(dynamic, 100)
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
        
        if(N_VOF(i,j,k) == none) continue;
        
        if(AE(i-1,j,k) > emf && i > 0)
          U(i-1,j,k) -= um[mesh_index(solver->mesh,i-1,j,k)];

        if(AN(i,j-1,k) > emf && j > 0)
          V(i,j-1,k) -= vm[mesh_index(solver->mesh,i,j-1,k)];

        if(AT(i,j,k-1) > emf && k > 0)
          W(i,j,k-1) -= wm[mesh_index(solver->mesh,i,j,k-1)];

      }
    } 
  }    */

  return solver->p_flag;
}

int vof_pressure_sor(struct solver_data *solver, long int i, long int j, long int k) {
  long int l,m,n;
  double plmn, delp;
  
  double g, del;
  double dp, dv;
  double ax, ux, stabil_limit;
  
  if(FV(i,j,k)<emf) return 0;

  if(VOF(i,j,k) < emf) return 0;
       
  if(N_VOF(i,j,k) != 0) {
    l = i;
    m = j;
    n = k;
    switch(N_VOF(i,j,k)) {
    case east:
      l=i+1;
      g = solver->gx;
      del = DELX;
      ax = AE(i,j,k); ux = U(i,j,k);
      break;
    case west:
      l=i-1;
      g = solver->gx;
      del = DELX;
      ax = AE(i-1,j,k); ux = U(i-1,j,k);
      break;
    case north:
      m=j+1;
      g = solver->gy;
      del = DELY;
      ax = AN(i,j,k); ux = V(i,j,k);
      break;
    case south:
      m=j-1;
      g = solver->gy;
      del = DELY;
      ax = AN(i,j-1,k); ux = V(i,j-1,k);
      break;
    case top:
      n=k+1;
      g = solver->gz;
      del = DELZ;
      ax = AT(i,j,k); ux = W(i,j,k);
      break;
    case bottom:
      n=k-1;
      g = solver->gz;
      del = DELZ;
      ax = AT(i,j,k-1); ux =  W(i,j,k-1);
      break;
    case none:
    default: /* ADDED 2/29/16 */
      return 0;
    }
 
    plmn = P(l,m,n);
  
    /* if the neighbor is also a surface cell, then set plm to psurf
     * this is the case for a very small floating droplet */
    if(N_VOF(l,m,n) != 0 && FV(i,j,k) != 0)
      plmn=0; /* psurf */

    if(FV(l,m,n) < emf) { 
      g = sqrt(pow(solver->gx,2) + pow(solver->gy,2) + pow(solver->gz,2));
      plmn = fabs(g) * solver->rho * (VOF(i,j,k) + 0.5) * del; 
    } 

    /* using the interpolation factor
     * we are setting delp to a mix of surface pressure and
     * internal pressure according to the formula below */
    delp = (1.0-PETA(i,j,k)) * plmn; 
  
    /* check if excessive pressure is building near the free surface and artificially release it */
    dp = plmn - delp;
    dv = solver->delt * ( (1/del) * dp / solver->rho);
    stabil_limit = solver->con * (min(FV(i,j,k),FV(l,m,n)) / ax)  * del/solver->delt;
  
    if((i-l + j-m + k-n) * ux > 0)
      stabil_limit = min(max(stabil_limit - fabs(ux), emf), stabil_limit / 2);
    else
      stabil_limit /= 2;      
  
    if(fabs(dv) > stabil_limit) {

      dp = stabil_limit / (solver->delt * (1/del) / solver->rho);
    
      /* printf("\n*** Excessive pressure gradient from %lf", delp - P(i,j,k)); */
    
      plmn = dp / PETA(i,j,k);
      P(l,m,n) = plmn;
      delp = (1-PETA(i,j,k)) * plmn;
    
    
      /* printf(" to %lf\n", delp - P(i,j,k));
      printf("dv: %lf    ux: %lf\n", dv, ux);
      printf("i,j,k: %ld %ld %ld\n", i, j, k);
      printf("l,m,n: %ld %ld %ld\n", l, m, n);
      printf("P: %lf   plmn: %lf\n\n", delp, plmn); */
    }
  
    delp = delp - P(i,j,k); 
  
    D(i,j,k) = 0;
  

  }
  else {

  /* if(i==229 && j==48 && k==42) {
    printf("break\n");
  } */
    D(i,j,k) = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k)) +
        RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k)) +
        RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1)); 
    D(i,j,k) = D(i,j,k) / FV(i,j,k); 
  
    /* de-foaming */
    if(VOF(i,j,k) < (1-emf)) {
      /* don't de-foam next to boundaries */
      if(!(FV(i+1,j,k) < emf || FV(i-1,j,k) < emf ||
           FV(i,j+1,k) < emf || FV(i,j-1,k) < emf ||
           FV(i,j,k+1) < emf || FV(i,j,k-1) < emf)) {
        D(i,j,k) = D(i,j,k) + (min(1000 * solver->epsi, 
                              0.1 * (1.0 - VOF(i,j,k)) / solver->delt)) / solver->rho;
      }
    }
  
    delp=-BETA(i,j,k)*D(i,j,k)*PETA(i,j,k);
  
    if(solver->iter > 50 && solver->iter < 75) delp /= (solver->omg - 0.1);
    if(solver->iter > 100 && solver->iter < 125) delp /= (solver->omg + 0.1);
  
    if(fabs(D(i,j,k)) > solver->epsi / FV(i,j,k)) 
    {
      solver->p_flag=1;
    }

  }

  P(i,j,k)=P(i,j,k)+delp;
  
  #pragma omp critical(pressure)
  {
  if(AE(i,j,k) > emf && i < IMAX-1)
    U(i,j,k)=U(i,j,k) + solver->delt* RDX * delp / (solver->rho );
  if(AE(i-1,j,k) > emf && i > 0)
    U(i-1,j,k)=U(i-1,j,k) - solver->delt* RDX * delp / (solver->rho );
  if(AN(i,j,k) > emf && j < JMAX-1)
    V(i,j,k)=V(i,j,k) + solver->delt * RDY * delp / (solver->rho );
  if(AN(i,j-1,k) > emf && j > 0)
    V(i,j-1,k)=V(i,j-1,k) - solver->delt * RDY * delp / (solver->rho );
  if(AT(i,j,k) > emf && k < KMAX-1)
    W(i,j,k)=W(i,j,k) + solver->delt * RDZ * delp / (solver->rho );
  if(AT(i,j,k-1) > emf && k > 0)
    W(i,j,k-1)=W(i,j,k-1) - solver->delt * RDZ * delp / (solver->rho );
  }
         
  /*           
  if(AE(i,j,k) > emf && i < IMAX-1)
    up[mesh_index(solver->mesh,i,j,k)]  =solver->delt * RDX * delp / (solver->rho );
  if(AE(i-1,j,k) > emf && i > 0)
    um[mesh_index(solver->mesh,i-1,j,k)]=solver->delt * RDX * delp / (solver->rho );
  if(AN(i,j,k) > emf && j < JMAX-1)
    vp[mesh_index(solver->mesh,i,j,k)]  =solver->delt * RDY * delp / (solver->rho );
  if(AN(i,j-1,k) > emf && j > 0)
    vm[mesh_index(solver->mesh,i,j-1,k)]=solver->delt * RDY * delp / (solver->rho );
  if(AT(i,j,k) > emf && k < KMAX-1)
    wp[mesh_index(solver->mesh,i,j,k)]  =solver->delt * RDZ * delp / (solver->rho );
  if(AT(i,j,k-1) > emf && k > 0)
    wm[mesh_index(solver->mesh,i,j,k-1)]=solver->delt * RDZ * delp / (solver->rho ); */

  D(i,j,k) = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k))+
        RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k))+
        RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1));
  D(i,j,k) = D(i,j,k) * solver->rho / FV(i,j,k);

  return 0;
}
#undef emf


int vof_pressure(struct solver_data *solver) {
  
  long int i,j,k,l,m,n, imax, jmax, kmax, imin, jmin, kmin;
  double plmn, delp;
  
  double g, del, dp, dv;
  double ax, ux, stabil_limit;
  
  
  solver->p_flag = 0;
  
#define emf solver->emf

  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {      
      
        
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
               
        if(N_VOF(i,j,k) != 0) {
          l = i;
          m = j;
          n = k;
          switch(N_VOF(i,j,k)) {
          case east:
            l=i+1;
            g = solver->gx;
            del = DELX;
            ax = AE(i,j,k); ux = U(i,j,k);
            break;
          case west:
            l=i-1;
            g = solver->gx;
            del = DELX;
            ax = AE(i-1,j,k); ux = U(i-1,j,k);
            break;
          case north:
            m=j+1;
            g = solver->gy;
            del = DELY;
            ax = AN(i,j,k); ux = V(i,j,k);
            break;
          case south:
            m=j-1;
            g = solver->gy;
            del = DELY;
            ax = AN(i,j-1,k); ux = V(i,j-1,k);
            break;
          case top:
            n=k+1;
            g = solver->gz;
            del = DELZ;
            ax = AT(i,j,k); ux = W(i,j,k);
            break;
          case bottom:
            n=k-1;
            g = solver->gz;
            del = DELZ;
            ax = AT(i,j,k-1); ux =  W(i,j,k-1);
            break;
          case none:
          default: /* ADDED 2/29/16 */
            continue;
          }
         
          plmn = P(l,m,n);
          
          /* if the neighbor is also a surface cell, then set plm to psurf
           * this is the case for a very small floating droplet */
          if(N_VOF(l,m,n) != 0 && FV(i,j,k) != 0)
            plmn=0; /* psurf */

          if(FV(l,m,n) < emf) { 
            g = sqrt(pow(solver->gx,2) + pow(solver->gy,2) + pow(solver->gz,2));
            plmn = fabs(g) * solver->rho * (VOF(i,j,k) + 0.5) * del; 
          } 
  
          /* using the interpolation factor
           * we are setting delp to a mix of surface pressure and
           * internal pressure according to the formula below */
          delp = (1.0-PETA(i,j,k)) * plmn; 
          
          /* check if excessive pressure is building near the free surface and artificially release it */
          dp = plmn - delp;
          dv = solver->delt * ( (1/del) * dp / solver->rho);
          stabil_limit = solver->con * (min(FV(i,j,k),FV(l,m,n)) / ax)  * del/solver->delt;
          
          if((i-l + j-m + k-n) * ux > 0)
          	stabil_limit = min(max(stabil_limit - fabs(ux), emf), stabil_limit / 2);
          else
          	stabil_limit /= 2;      
          
          if(fabs(dv) > stabil_limit) {

            dp = stabil_limit / (solver->delt * (1/del) / solver->rho);
            
            /* printf("\n*** Excessive pressure gradient from %lf", delp - P(i,j,k)); */
            
            plmn = dp / PETA(i,j,k);
            P(l,m,n) = plmn;
            delp = (1-PETA(i,j,k)) * plmn;
            
            
            /* printf(" to %lf\n", delp - P(i,j,k));
            printf("dv: %lf    ux: %lf\n", dv, ux);
            printf("i,j,k: %ld %ld %ld\n", i, j, k);
            printf("l,m,n: %ld %ld %ld\n", l, m, n);
            printf("P: %lf   plmn: %lf\n\n", delp, plmn); */
          }
          
          delp = delp - P(i,j,k); 
          
          D(i,j,k) = 0;
          

        }
        else {


          D(i,j,k) = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k)) +
              RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k)) +
              RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1)); 
          D(i,j,k) = D(i,j,k) / FV(i,j,k); 
          
          /* de-foaming */
          if(VOF(i,j,k) < (1-emf)) {
            /* don't de-foam next to boundaries */
            if(!(FV(i+1,j,k) < emf || FV(i-1,j,k) < emf ||
                 FV(i,j+1,k) < emf || FV(i,j-1,k) < emf ||
                 FV(i,j,k+1) < emf || FV(i,j,k-1) < emf)) {
              D(i,j,k) = D(i,j,k) + (min(1000 * solver->epsi, 
                                    (1.0 - VOF(i,j,k)) / solver->delt)) / solver->rho;
            }
          }
          
          delp=-BETA(i,j,k)*D(i,j,k)*PETA(i,j,k);
          
          if(solver->iter > 50 && solver->iter < 75) delp /= (solver->omg - 0.1);
          if(solver->iter > 100 && solver->iter < 125) delp /= (solver->omg + 0.1);
          
          if(fabs(D(i,j,k)) > solver->epsi / FV(i,j,k)) 
          {
            solver->p_flag=1;
          }

        }
        
        P(i,j,k)=P(i,j,k)+delp;

        if(AE(i,j,k) > emf && i < IMAX-1)
          U(i,j,k)=U(i,j,k) + solver->delt* RDX * delp / (solver->rho /* AE(i,j,k) */);
        if(AE(i-1,j,k) > emf && i > 0)
          U(i-1,j,k)=U(i-1,j,k) - solver->delt* RDX * delp / (solver->rho /* AE(i-1,j,k) */);
        if(AN(i,j,k) > emf && j < JMAX-1)
          V(i,j,k)=V(i,j,k) + solver->delt * RDY * delp / (solver->rho /* AN(i,j,k) */);
        if(AN(i,j-1,k) > emf && j > 0)
          V(i,j-1,k)=V(i,j-1,k) - solver->delt * RDY * delp / (solver->rho /* AN(i,j-1,k) */);
        if(AT(i,j,k) > emf && k < KMAX-1)
          W(i,j,k)=W(i,j,k) + solver->delt * RDZ * delp / (solver->rho /* AT(i,j,k) */);
        if(AT(i,j,k-1) > emf && k > 0)
          W(i,j,k-1)=W(i,j,k-1) - solver->delt * RDZ * delp / (solver->rho /* AT(i,j,k-1) */);
    
 

        D(i,j,k) = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k))+
              RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k))+
              RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1));
        D(i,j,k) = D(i,j,k) * solver->rho / FV(i,j,k);

      }
    }
  }


  return solver->p_flag;
#undef emf
}


int vof_pressure_old(struct solver_data *solver) {
  /* double D=0.0; */
  enum cell_boundaries ignore;
  
  long int i,j,k,l,m,n;
  double plmn, delp;
  
  double g, del, dp, dv;
  double rcsq;
  double sum_a, ax, ux, stabil_limit;
  
  /*enum cell_boundaries nfew, nfee, nfen, nfes, nfet, nfeb, nfe;*/

  rcsq = 1 / pow(solver->csq, 2);
  
  solver->p_flag = 0;
  
#define emf solver->emf

  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        ignore = none; /* ********* D = 0.0; */
        
        if(FV(i,j,k)<emf) continue;

        if(VOF(i,j,k) < emf) continue;
               
        if(N_VOF(i,j,k) != 0) {
          l = i;
          m = j;
          n = k;
          switch(N_VOF(i,j,k)) {
          case east:
            l=i+1;
            if(N_VOF(i-1,j,k) == 0) ignore = west;
            g = solver->gx;
            del = DELX;
            ax = AE(i,j,k); ux = U(i,j,k);
            break;
          case west:
            l=i-1;
            if(N_VOF(i+1,j,k) == 0) ignore = east;
            g = solver->gx;
            del = DELX;
            ax = AE(i-1,j,k); ux = U(i-1,j,k);
            break;
          case north:
            m=j+1;
            if(N_VOF(i,j-1,k) == 0) ignore = south;
            g = solver->gy;
            del = DELY;
            ax = AN(i,j,k); ux = V(i,j,k);
            break;
          case south:
            m=j-1;
            if(N_VOF(i,j+1,k) == 0) ignore = north;
            g = solver->gy;
            del = DELY;
            ax = AN(i,j-1,k); ux = V(i,j-1,k);
            break;
          case top:
            n=k+1;
            if(N_VOF(i,j,k-1) == 0) ignore = bottom;
            g = solver->gz;
            del = DELZ;
            ax = AT(i,j,k); ux = W(i,j,k);
            break;
          case bottom:
            n=k-1;
            if(N_VOF(i,j,k+1) == 0) ignore = top;
            g = solver->gz;
            del = DELZ;
            ax = AT(i,j,k-1); ux =  W(i,j,k-1);
            break;
          case none:
            continue;
          }
         
          /*
          nfew = N_VOF(i-1,j,k);
          nfee = N_VOF(i+1,j,k);
          nfes = N_VOF(i,j-1,k);
          nfen = N_VOF(i,j+1,k);
          nfeb = N_VOF(i,j,k-1);
          nfet = N_VOF(i,j,k+1);

          nfe = fmax(fmax(fmax(nfee,nfew),fmax(nfes,nfen)),fmax(nfeb,nfet));

          * psurf = PS(i,j,k)+PR(nfe); SURFACE PRESSURE CODE NOT WRITTEN */

          /* plm is the fluid pressure normal from the interface */
          plmn = P(l,m,n);
          
          /* if the neighbor is also a surface cell, then set plm to psurf
           * this is the case for a very small floating droplet */
          if(N_VOF(l,m,n) != 0 && FV(i,j,k) != 0)
            plmn=0; /* psurf */

          /* ADDED 9/12 */
          if(FV(l,m,n) < emf) { 
            g = sqrt(pow(solver->gx,2) + pow(solver->gy,2) + pow(solver->gz,2));
            plmn = fabs(g) * solver->rho * (VOF(i,j,k) + 0.5) * del; 
          } 
  
          /* using the interpolation factor
           * we are setting delp to a mix of surface pressure and
           * internal pressure according to the formula below */
          delp = (1.0-PETA(i,j,k)) * plmn; /*+ PETA(i,j,k) * psurf  - P(i,j,k); */
          /* if(delp < 0) {
            if(delp < (1.0-PETA(i,j,k)) *  fabs(g) * solver->rho * del * 0.5) {
              P(i,j,k) = (1.0-PETA(i,j,k)) *  fabs(g) * solver->rho * del * 0.5;
              continue;
            }
          } */
          dp = plmn - delp;
          dv = solver->delt * ( /*(sum_fv/2) * */ (1/del) * dp / solver->rho);
          /* dtmax = CON * min(1, min(FV(i,j,k),FV(i,min(JMAX-1,j+1),k)) / AN(i,j,k)) * DELY/fabs(V(i,j,k)); */
          /* TODO: limit dv based on stability limit */
          stabil_limit = solver->con * (min(FV(i,j,k),FV(l,m,n)) / ax)  * del/solver->delt;
          if(fabs(dv * (i-l + j-m + k-n) 
                  + ux) > stabil_limit ||
             fabs(dv) > stabil_limit / 4) {
            /* dv = 0.4 * (min(FV(i,j,k),FV(l,m,n)) / ax) * del / solver->delt; 
            dp = dv / (solver->delt * (1/del) / solver->rho); */
            dp /= 2;
            printf("\n*** Excessive pressure adjusted from %lf to ", delp);
            delp = plmn - dp;
            printf(" to %lf\n", delp);
            printf("dv: %lf    ux: %lf\n", dv, ux);
            printf("i,j,k: %ld %ld %ld\n", i, j, k);
            printf("l,m,n: %ld %ld %ld\n", l, m, n);
            printf("plmn: %lf\n\n", plmn);
          }
          
          delp = delp - P(i,j,k); 
          /* P(i,j,k) = delp; */
          
          D(i,j,k) = 0;
          /* continue; */
          
          /* CONFIRM IF THIS SHOULD BE INCLUDED *
          if(FV(l,m,n) < emf) { 
            P(i,j,k)=P(i,j,k)+delp;
            continue; 
          } */

        }
        else {

          D(i,j,k) = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k)) +
              RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k)) +
              RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1)); /* +
              (1 / solver->rho) * rcsq * FV(i,j,k) * pow(VOF(i,j,k),2) * 
              (P(i,j,k) - PN(i,j,k)) / solver->delt; Uncomment for compressibility */
          D(i,j,k) = D(i,j,k) /* * solver->rho*/ / FV(i,j,k); 
          
          /* de-foaming */
          if(VOF(i,j,k) < (1-emf))
            D(i,j,k) = D(i,j,k) + min(100 * solver->epsi, 
                                    0.1 * (1.0 - VOF(i,j,k)) / solver->delt);
          
          delp=-BETA(i,j,k)*D(i,j,k)*PETA(i,j,k);
          
          if(solver->iter > 50 && solver->iter < 75) delp /= (solver->omg - 0.1);
          if(solver->iter > 100 && solver->iter < 125) delp /= (solver->omg + 0.1);
          
          if(fabs(D(i,j,k)) > solver->epsi / FV(i,j,k)) 
          {
            solver->p_flag=1;
#ifdef DEBUG
            if(solver->iter > 200 && /* j>9 && j<11 && */ solver->iter % 100 == 0) {
              printf("\nCell %ld %ld %ld having trouble with convergence\n",i,j,k);
              printf("P: %lf  PN: %lf  VOF: %lf  delp: %lf D: %lf\n",P(i,j,k),PN(i,j,k),VOF(i,j,k),delp,D(i,j,k));
              printf("N_VOF: %d  PETA: %lf  BETA: %lf\n",N_VOF(i,j,k),PETA(i,j,k),BETA(i,j,k));
              printf("N_N_VOF: east %d west %d north %d south %d top %d bottom %d\n",
                      N_VOF(i+1,j,k), N_VOF(i-1,j,k), N_VOF(i,j+1,k), N_VOF(i,j-1,k),
                      N_VOF(i,j,k+1), N_VOF(i,j,k-1)); 
              printf("U: %lf %lf %lf\n",U(i,j,k),V(i,j,k),W(i,j,k));
              printf("U-1: %lf %lf %lf\n",U(i-1,j,k),V(i,j-1,k),W(i,j,k-1));
              printf("UN: %lf %lf %lf\n\n",UN(i,j,k),VN(i,j,k),WN(i,j,k));     
                if(solver->iter > 1001) {
                  vtk_write_P(solver->mesh,0);
                  vtk_write_U(solver->mesh,0);
                  vtk_write_vof(solver->mesh,0);
                  
              }         
            }
#endif                     
          }
          else;

        }
        ignore = none;
        
        sum_a = AE(i,j,k) + AE(i-1,j,k) + AN(i,j,k) + AN(i,j-1,k) + AT(i,j,k) + AT(i,j,k-1);
        
        P(i,j,k)=P(i,j,k)+delp;
            
        if(AE(i,j,k) > emf && i < IMAX-1)
          U(i,j,k)=U(i,j,k) + solver->delt* RDX * delp / (solver->rho /* AE(i,j,k) */);
        if(AE(i-1,j,k) > emf && i > 0)
          U(i-1,j,k)=U(i-1,j,k) - solver->delt* RDX * delp / (solver->rho /* AE(i-1,j,k) */);
        if(AN(i,j,k) > emf && j < JMAX-1)
          V(i,j,k)=V(i,j,k) + solver->delt * RDY * delp / (solver->rho /* AN(i,j,k) */);
        if(AN(i,j-1,k) > emf && j > 0)
          V(i,j-1,k)=V(i,j-1,k) - solver->delt * RDY * delp / (solver->rho /* AN(i,j-1,k) */);
        if(AT(i,j,k) > emf && k < KMAX-1)
          W(i,j,k)=W(i,j,k) + solver->delt * RDZ * delp / (solver->rho /* AT(i,j,k) */);
        if(AT(i,j,k-1) > emf && k > 0)
          W(i,j,k-1)=W(i,j,k-1) - solver->delt * RDZ * delp / (solver->rho /* AT(i,j,k-1) */);

            if(solver->iter > 200 && j>9 && j<11 && fabs(D(i,j,k)/solver->dzro) > solver->epsi && solver->iter % 100 == 0) {
              printf("After correction: P: %lf\n",P(i,j,k));
              printf("U: %lf %lf %lf\n",U(i,j,k),V(i,j,k),W(i,j,k));
              printf("U-1: %lf %lf %lf\n",U(i-1,j,k),V(i,j-1,k),W(i,j,k-1));            
            
            }
        D(i,j,k) = RDX*(AE(i,j,k)*U(i,j,k)-AE(i-1,j,k)*U(i-1,j,k))+
              RDY*(AN(i,j,k)*V(i,j,k)-AN(i,j-1,k)*V(i,j-1,k))+
              RDZ*(AT(i,j,k)*W(i,j,k)-AT(i,j,k-1)*W(i,j,k-1));
        D(i,j,k) = D(i,j,k) * solver->rho / FV(i,j,k);
        /* if(fabs(D(i,j,k) /solver->dzro) > solver->epsi) D(i,j,k) = 0; */

      }
    }
  }

  return solver->p_flag;
#undef emf
}
