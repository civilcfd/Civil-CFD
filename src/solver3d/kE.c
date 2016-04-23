/*
 * kE.c 
 *
 * Implements a k-epsilon turbulence model */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "vtk.h"
#include "vof.h"
#include "solver.h"
#include "readfile.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "vof_boundary.h"
 
#include "vof_macros.h" 
#include "vector_macros.h"
  
static struct kE_data kE, kE_n;

#define k(l,m,n) kE.k[mesh_index(solver->mesh, l, m, n)]
#define E(l,m,n) kE.E[mesh_index(solver->mesh, l, m, n)]
#define nu_t(l,m,n) kE.nu_t[mesh_index(solver->mesh, l, m, n)]
#define k_N(l,m,n) kE_n.k[mesh_index(solver->mesh, l, m, n)]
#define E_N(l,m,n) kE_n.E[mesh_index(solver->mesh, l, m, n)]
#define nu_t_N(l,m,n) kE_n.nu_t[mesh_index(solver->mesh, l, m, n)]
#define tau_x(l,m,n) kE.tau_x[mesh_index(solver->mesh, l, m, n)]
#define tau_y(l,m,n) kE.tau_y[mesh_index(solver->mesh, l, m, n)]
#define tau_z(l,m,n) kE.tau_z[mesh_index(solver->mesh, l, m, n)]

#ifdef DEBUG
float kk(struct solver_data *solver, long int i,long int j,long int k);
float kk(struct solver_data *solver, long int i,long int j,long int k) { return k(i,j,k); } 
float kk_n(struct solver_data *solver, long int i,long int j,long int k);
float kk_n(struct solver_data *solver, long int i,long int j,long int k) { return k_N(i,j,k); } 
float EE(struct solver_data *solver, long int i,long int j,long int k);
float EE(struct solver_data *solver, long int i,long int j,long int k) { return E(i,j,k); } 
float EE_n(struct solver_data *solver, long int i,long int j,long int k);
float EE_n(struct solver_data *solver, long int i,long int j,long int k) { return E_N(i,j,k); } 
#endif

int kE_set_value(char *param, int dims, 
                   double *vector) {


  if(param == NULL || vector == NULL) {
    printf("error: null values passed to kE_set_value\n");
    return(1);
  }

  if(strcmp(param, "length_scale")==0) {
    if(dims != 1) {
      printf("error in source file: length_scale requires 3 arguments\n");
      return(1);
    }

    kE.length_scale = vector[0];
    kE_n.length_scale = kE.length_scale;
  }
  else if (strcmp(param, "length")==0) {
    if(dims != 1) {
      printf("error in source file: length requires 1 arguments\n");
      return(1);
    }

    kE.length = vector[0];
    kE_n.length = kE.length;

  }
  else if (strcmp(param, "rough")==0) {
    if(dims != 1) {
      printf("error in source file: rough requires 1 arguments\n");
      return(1);
    }

    kE.rough = vector[0];
    kE_n.rough = kE.rough;

  }

	return 0;
}

int kE_read(char *filename)
{
  /* code  to read kE variables
   * this includes length, length_scale, rough
   *
   */
  FILE *fp;

  char text[1024];
  char args[4][256];
  double vector[3];

  char *s;
  int i, nargs;

  if(filename == NULL) {
    printf("error: passed null arguments to kE_read\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in kE_read\n",filename);
    return(1);
  }

  while(!feof(fp)) {
    if(!fgets(text, sizeof(text), fp))
    { 
      if(!feof(fp)) {
        printf("error: fgets in kE_read\n");
        return(1);
      }
      else break;
    }

    nargs = read_args(text, 4, args);

    if(nargs > 0) nargs--;
    for(i = nargs; i < 3; i++) {
      vector[i] = 0;
    }

    for(i = 0; i < nargs; i++) {
      /* sscanf(args[i], "%lf", &vector[i]); */
      vector[i] = strtod(args[i+1],&s); 
    }

    #ifdef DEBUG
      printf("kE_read: read %d args: %s %lf %lf %lf\n", nargs, args[0], 
             vector[0], vector[1], vector[2]);
    #endif

    if (args[0][0] != '#' && args[0][0] != 0) 
      kE_set_value(args[0], i, vector); 
    args[0][0] = 0;
  }

  fclose(fp);

  return 0;
}

/*
int kE_read(char *filename) {
  FILE *fp;
  char dummy[256];

  if(filename == NULL) {
    printf("error: passed null arguments to kE_read\n");
    return(1);
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    printf("error: cannot open %s in kE_read\n",filename);
    return(1);
  }

  fscanf(fp,"length_scale %le",&kE.length_scale);
  
  fscanf(fp,"%s %le",dummy,&kE.rough);
  kE_n.rough = kE.rough;

  fclose(fp);
  return 0;
}*/

int kE_write(char *filename) {
  FILE *fp;

  if(filename == NULL) {
    printf("error: passed null arguments to kE_write\n");
    return(1);
  }

  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: cannot open %s in kE_write\n",filename);
    return(1);
  }

  fprintf(fp,"length_scale %le\n",kE.length_scale);
  fprintf(fp,"length %le\n",kE.length);
  fprintf(fp,"rough %le\n",kE.rough);

  fclose(fp);
  return 0;
}

int kE_setup(struct solver_data *solver) {

  solver->turbulence_init = kE_init;
  solver->turbulence_loop = kE_loop_explicit;
  solver->turbulence_kill = kE_kill;
  solver->turbulence_read = kE_read;
  solver->turbulence_write = kE_write;
  solver->wall_shear = kE_wall_shear;
  solver->turbulence_nu   = kE_nu;
    
  solver->mesh->turbulence_model = &kE;
  
  kE.C1E = kE_n.C1E = 1.44;
  kE.C2E = kE_n.C2E = 1.92;
  kE.C3E = kE_n.C3E = -0.33;
  kE.C_mu = kE_n.C_mu = 0.09;
  kE.sigma_k = kE_n.sigma_k = 1.0;
  kE.sigma_E = kE_n.sigma_E = 1.3;
  kE.vonKarman = kE_n.vonKarman = 0.41;
  kE.length_scale = 0.038;
  kE.length = min(IMAX * DELX, JMAX * DELY);
  kE.length = min(kE.length, KMAX * DELZ);

  kE.rough = kE_n.rough = 0.00161; /* concrete */
  
  return 0;
}

int kE_init(struct solver_data *solver) {
  long int size = IMAX * JMAX * KMAX;

  kE.k = malloc(sizeof(double) * size);
  kE.E = malloc(sizeof(double) * size);
  kE.nu_t = malloc(sizeof(double) * size);
  kE_n.k = malloc(sizeof(double) * size);
  kE_n.E = malloc(sizeof(double) * size);
  kE_n.nu_t = malloc(sizeof(double) * size);
  
  if(kE.k == NULL || kE.E == NULL || kE.nu_t == NULL || 
     kE_n.k == NULL || kE_n.E == NULL || kE_n.nu_t == NULL) {
    printf("error: kE_init could not allocate k, E or nu_t\n");
    return 1;
  }
  
  /* allocate tau.  note we don't need to allocate for kE_n since tau is an intermediate calculated value */
  /* i.e. it is never explicitly calculated */
  kE.tau_x = malloc(sizeof(double) * size);
  kE.tau_y = malloc(sizeof(double) * size);
  kE.tau_z = malloc(sizeof(double) * size);
  if(kE.tau_x == NULL || kE.tau_y == NULL || kE.tau_z == NULL) {
    printf("error: kE_init could not allocate tau\n");
    return 1;
  }
  
   
  kE_set_internal(solver, 0.0001, 0);
  kE.length *= kE.length_scale;
 
  return 0;
}

int kE_check(struct solver_data *solver) {
  return (solver->mesh->turbulence_model == &kE);
}

int kE_boundary_fixed_k(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence) {
  long int i, j, k, imin, jmin, kmin, imax, jmax, kmax;
  
  sboundary_setup(solver, x, &imin, &jmin, &kmin, &imax, &jmax, &kmax, min_1, min_2, max_1, max_2);
  
  for(i=imin; i <= imax; i++) {
    for(j=jmin; j <= jmax; j++) {
      for(k=kmin; k <= kmax; k++) { 
        /* first check if this is an outflow boundary.  if so, skip it */
        switch(x) {
        case 0: /* west */
          if(U(i,j,k) < 0) continue;
          break;
        case 1: /* east */
          if(U(i,j,k) > 0) continue;
          break;
        case 2: /* south */
          if(V(i,j,k) < 0) continue;
          break;
        case 3: /* north */
          if(V(i,j,k) > 0) continue;
          break;
        case 4: /* bottom */
          if(W(i,j,k) < 0) continue;
          break;
        case 5: /* top */
          if(W(i,j,k) > 0) continue;
          break;
        }
      
        k(i,j,k) = turbulence;  
      }
    }
  }

  return 0;                            
}

int kE_special_boundaries(struct solver_data *solver) {
  int x;
  struct sb_data *sb;

  for(x=0; x < 6; x++) {
    for(sb = solver->mesh->sb[x]; sb != NULL; sb = sb->next) {    
        switch(sb->type) {
        case hgl:
        case mass_outflow:
        case fixed_velocity:
          kE_boundary_fixed_k(solver, x, sb->extent_a[0], sb->extent_a[1], 
                                     sb->extent_b[0], sb->extent_b[1], 
                                     sb->value, sb->turbulence);
          break;
        case weir:
          kE_boundary_fixed_k(solver, x, sb->extent_a[0], sb->extent_a[1], 
                                     sb->extent_b[0], sb->extent_b[1], 
                                     sb->value, 0.0001);
          break;
        case wall:
        	break;
        }
    }
  }

  return 0;
}

int kE_set_internal(struct solver_data *solver, double k, double E) {
  long int l,m,n;
  double E_limit;
  
  E_limit = kE.C_mu * pow(k, 1.5) / kE.length;

#pragma omp parallel for shared(solver,k,E,E_limit) private(l,m,n)  
  for(l=0; l<IMAX; l++) {
    for(m=0; m<JMAX; m++) {
      for(n=0; n<KMAX; n++) {
        k(l,m,n) = k;
        E(l,m,n) = max(E,E_limit);
        
        /* set shear stresses to zero upon initialization */
        tau_x(l,m,n) = 0;
        tau_y(l,m,n) = 0;
        tau_z(l,m,n) = 0;
        
        nu_t(l,m,n) = 0;
      }
    }
  }
      
  kE_boundaries(solver);
  kE_copy(solver);
  
  return 0;
}

int kE_copy(struct solver_data *solver) {
  long int l,m,n;

#pragma omp parallel for shared(solver) private(l,m,n)   
  for(l=0; l<IMAX; l++) {
    for(m=0; m<JMAX; m++) {
      for(n=0; n<KMAX; n++) {
        k_N(l,m,n) = k(l,m,n);
        E_N(l,m,n) = E(l,m,n);
        nu_t_N(l,m,n) = nu_t(l,m,n);
      }
    }
  }
      
  return 0;
}

double kE_nu(struct solver_data *solver, long int i, long int j, long int k) {
  const double nu = solver->nu / kE.sigma_k; 

  return nu_t(i,j,k) + nu;
}

int kE_boundaries(struct solver_data *solver) {

  /* sets boundaries on k and epsilon
   * does not set boundaries on U, see wall_shear */
  long int i,j,k, l, m, n;


#pragma omp parallel for shared(solver) private(i,j,k)  
  for(i=0; i<IMAX; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
      
      	tau_x(i,j,k) = 0;
      	tau_y(i,j,k) = 0;
      	tau_z(i,j,k) = 0;
      
        if(i>0 && i<IMAX-1 && j>0 && j<JMAX-1 && k>0 && k<KMAX-1) continue;
        /* set mesh boundaries */
        k(0,j,k) = k(1,j,k);
        k(IMAX-1,j,k) = k(IMAX-2,j,k);
        E(0,j,k) = E(1,j,k);
        E(IMAX-1,j,k) = E(IMAX-2,j,k);
        nu_t(0,j,k) = nu_t(1,j,k);
        nu_t(IMAX-1,j,k) = nu_t(IMAX-2,j,k);

        k(i,0,k) = k(i,1,k);
        k(i,JMAX-1,k) = k(i,JMAX-2,k);
        E(i,0,k) = E(i,1,k);
        E(i,JMAX-1,k) = E(i,JMAX-2,k);
        nu_t(i,0,k) = nu_t(i,1,k);
        nu_t(i,JMAX-1,k) = nu_t(i,JMAX-2,k);
        

        k(i,j,0) = k(i,j,1);
        k(i,j,KMAX-1) = k(i,j,KMAX-2);
        E(i,j,0) = E(i,j,1);
        E(i,j,KMAX-1) = E(i,j,KMAX-2);
        nu_t(i,j,0) = nu_t(i,j,1);
        nu_t(i,j,KMAX-1) = nu_t(i,j,KMAX-2);

      }
    }
  }

  kE_special_boundaries(solver);
  
  kE_tau(solver);

	/* free surface boundaries - the free surface acts as a von Neumann boundary */
	/* also takes care of cells outside of the solution space */
	/* ADDED 10/26 */
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
      
      	if(FV(i,j,k) < solver->emf || VOF(i,j,k) < solver->emf || N_VOF(i,j,k) > 6) { 
      		/* at an interior obstacle or empty cell, set k = 0 */
      		k(i,j,k) = 0;
      		E(i,j,k) = 0;
      		nu_t(i,j,k) = 0;
        	tau_x(i,j,k) = 0;
        	tau_y(i,j,k) = 0;
        	tau_z(i,j,k) = 0;
      		
      		continue;
      	}      	
      	
      	if(N_VOF(i,j,k) == 0 || FV(i,j,k) < (1 - solver->emf)) continue;

				l = i;
				m = j;
				n = k;
				switch(N_VOF(i,j,k)) {
				case east:
					l=i+1;
					break;
				case west:
					l=i-1;
					break;
				case north:
					m=j+1;
					break;
				case south:
					m=j-1;
					break;
				case top:
					n=k+1;
					break;
				case bottom:
					n=k-1;
					break;
				case none:
				case empty:
					continue;
				}
      
      		/*
          if(N_VOF(l,m,n) != 0 && FV(i,j,k) != 0) {
            k(i,j,k) = 0;
            E(i,j,k) = 0;
        		tau_x(i,j,k) = 0;
        		tau_y(i,j,k) = 0;
        		tau_z(i,j,k) = 0;
        		nu_t(i,j,k) = 0;
            continue;   		
          }*/
				k(i,j,k) = k(l,m,n);
				E(i,j,k) = E(l,m,n);
				nu_t(i,j,k) = nu_t(l,m,n);
      
      }
    }
  }

  return 0;
}

#ifdef _WIN32
double log_law(double u, double d, double mu, double rho, double rough)
#else
inline double log_law(double u, double d, double mu, double rho, double rough) 
#endif
{

  int iter = 0;
  double u_t, u_t_n, res, delB, y_plus, hs_plus;
  u_t_n = u; /* initial guess */

  do {

    if(iter > 300) {
      if(res > fabs(u) / 100) {
      	/* printf("log_law did not converge.  final residual %lf\n", res); */
      	return 0;
      } else {
      	return u_t;
      }
      
      break;
    }

    hs_plus = rough * u_t_n * rho / mu;
    delB = (1/kE.vonKarman) * log(1 + 0.3 * hs_plus); 

    y_plus = rho * u_t_n * d / mu;
    if(y_plus < 11.06) {
#ifdef DEBUG
      /*printf("minimum y plus reached - grid resolution too fine near wall\n");*/
#endif
      /*y_plus = 11.06; one option is to set y_plus to a hard limit*/
      
      /* y_plus = u / u_t instead we will assume we are inside the viscous sublayer
         u_t = u / y_plus
         u_t = u * mu / (rho * u_t * d)
         can eliminate u_t_n from y_plus 
         so u_t = sqrt(u * mu / (rho * d))*/
      u_t = sqrt(u * mu / (rho * d));
      break;
    }
    if(y_plus < hs_plus / 2) {
#ifdef DEBUG
      printf("yplus < hs_plus/2\n");
#endif
      y_plus = hs_plus/2;
    }

		/* if(y_plus > 10000) {
			return 0;
		} */

    u_t = fabs(u * kE.vonKarman / ( log(y_plus) + 5.2 - delB )); /* testing fabs */
    res = fabs(u_t_n - u_t);

		if(iter < 100)
    	u_t_n = u_t;
		else if (iter > 100 && iter < 150) 
			u_t = u_t_n + 0.8 * (u_t - u_t_n);
		else if (iter > 200 && iter < 250)
			u_t = u_t_n + 1.2 * (u_t - u_t_n);
			

    iter++;

  } while(res > 0.0001);
  
  return u_t;
}

/*************************************
NEEDS TO BE FIXED - Only iterates on one wall

Currently not used in the solver.  Only obstacles apply shear
*************************************/
int kE_boundary_wall_shear(struct solver_data *solver) {

  /* extension of velocity to add/subtract wall shears from velocity variable
   * wall shear from south wall on U = 2 * nu * u / (dy^2 * AE) * (1 - AS)  if laminar 
   * for turbulent, we can calculate u* based on k (k = u* ^2 / C)
   * and then tau from u*.  Wall stress = 1/rho * tau/dy 
   * unfortunately, this is cell centered and so will need to be averaged for cells i and i+1 */
  long int i,j,k;
  double ws_u, ws_d, tau_u, tau_d;


  i = 0;
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {


        if( ( ((FV(i,j,k) > solver->emf)   && (FV(i,j,k) < 1-solver->emf)) ||
              ((FV(i+1,j,k) > solver->emf) && (FV(i+1,j,k) < 1-solver->emf)) ) &&
              AE(i,j,k) > solver->emf) { 
          

          tau_u = tau_x(i,j,k);
          tau_d = tau_x(i+1,j,k);
          
          ws_u  = tau_u  / (solver->rho * DELY) * fabs((1-AN(i,j,k)) - (1-AN(i,j-1,k)));
          ws_u += tau_u / (solver->rho * DELZ) * fabs((1-AT(i,j,k)) - (1-AT(i,j,k-1)));
          
          ws_d  = tau_d / (solver->rho * DELY) * fabs((1-AN(i+1,j,k)) - (1-AN(i+1,j-1,k)));
          ws_d += tau_d / (solver->rho * DELZ) * fabs((1-AT(i+1,j,k)) - (1-AT(i+1,j,k-1)));
        
          U(i,j,k) += 0.5 * (ws_u + ws_d) * solver->delt;
          
          if(isnan(U(i,j,k))) {
          #ifdef DEBUG
            printf("nan at cell %ld %ld %ld in wall_shear\n",i,j,k);
          #endif
          }
        }

        if( ( ((FV(i,j,k) > solver->emf && FV(i,j,k) < 1-solver->emf)) ||
              ((FV(i,j+1,k) > solver->emf && FV(i,j+1,k) < 1-solver->emf)) ) &&
                AN(i,j,k) > solver->emf) { 
          

          tau_u = tau_y(i,j,k);
          tau_d = tau_y(i,j+1,k);
          
          ws_u  = tau_u / (solver->rho * DELX) * fabs((1-AE(i,j,k)) - (1-AE(i-1,j,k)));
          ws_u += tau_u / (solver->rho * DELZ) * fabs((1-AT(i,j,k)) - (1-AT(i,j,k-1)));
          
          ws_d  = tau_d / (solver->rho * DELX) * fabs((1-AE(i,j+1,k)) - (1-AE(i-1,j+1,k)));
          ws_d += tau_d / (solver->rho * DELZ) * fabs((1-AT(i,j+1,k)) - (1-AT(i,j+1,k-1)));
        
          V(i,j,k) += 0.5 * (ws_u + ws_d) * solver->delt;
        }

        if( ( ((FV(i,j,k) > solver->emf && FV(i,j,k) < 1-solver->emf)) ||
              ((FV(i,j,k+1) > solver->emf && FV(i,j,k+1) < 1-solver->emf)) ) &&
                AT(i,j,k) > solver->emf) { 
          

          tau_u = tau_z(i,j,k);
          tau_d = tau_z(i,j,k+1);
          
          ws_u  = tau_u / (solver->rho * DELX) * fabs((1-AE(i,j,k)) - (1-AE(i-1,j,k)));
          ws_u += tau_u / (solver->rho * DELY) * fabs((1-AN(i,j,k)) - (1-AN(i,j-1,k)));
          
          ws_d  = tau_d / (solver->rho * DELX) * fabs((1-AE(i,j,k+1)) - (1-AE(i-1,j,k+1)));
          ws_d += tau_d / (solver->rho * DELY) * fabs((1-AN(i,j,k+1)) - (1-AN(i,j-1,k+1)));
        
          W(i,j,k) += 0.5 * (ws_u + ws_d) * solver->delt;
        }


      }
    }
  

  return 0;
}

int kE_tau(struct solver_data *solver) {
  long int i,j,k, im1, jm1, km1, l, m, n;
  double d, u_t, wall_n[3], u_c[3], mag, tau;
  double u_perp_n[3], u_perp_c, u_parr_c, u_parr[3];
  double E_limit;
  const double del[3] = { DELX, DELY, DELZ };
  
#pragma omp parallel for shared(solver, del) private(i,j,k,im1, jm1, km1, l, m, n, \
              d, u_t, wall_n, u_c, mag, tau, u_perp_n, u_perp_c, u_parr_c, u_parr, E_limit)  
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
      
        if(FV(i,j,k) > (1-solver->emf) || FV(i,j,k) < solver->emf || VOF(i,j,k) < solver->emf)
          continue;
          
      	if(N_VOF(i,j,k) > 6) continue;

        im1 = i-1;
        jm1 = j-1;
        km1 = k-1;
  
        /* first adjust delx,dely,delz based on space occupied in cell */
        wall_n[0] = AE(i,j,k) - AE(im1,j,k);
        wall_n[1] = AN(i,j,k) - AN(i,jm1,k);
        wall_n[2] = AT(i,j,k) - AT(i,j,km1);

        mag = vector_magnitude(wall_n);
        wall_n[0] /= mag;
        wall_n[1] /= mag;
        wall_n[2] /= mag;

        d = fabs(wall_n[0]) * del[0] + fabs(wall_n[1]) * del[1] + fabs(wall_n[2]) * del[2];
        d /= 2;
        /* d *= FV(i,j,k); */

        u_c[0] = (U(im1,j,k) + U(i,j,k))/2;
        u_c[1] = (V(i,jm1,k) + V(i,j,k))/2;
        u_c[2] = (W(i,j,km1) + W(i,j,k))/2;

        /* find velocity component perpendicular to wall normal vector */
        /* first find parallel component */
        u_parr_c = inner_product(u_c, wall_n);
        
        /* now find parallel vector by multiplying scalar by wall normal */
        vector_multiply(u_parr, wall_n, u_parr_c); /* u_parr[]   = wall_n[] * u_parr_c */
        /* now find perpendicular vector by subtracting u_parr from u_c */
        vector_subtract(u_perp_n, u_c, u_parr);    /* u_perp_n[] = u_c[]    - u_parr[] */
        
        /* now normalize u_perp_n */
        u_perp_c = vector_magnitude(u_perp_n);
        if(u_perp_c > 0.0001) {
          u_perp_n[0] /= u_perp_c;
          u_perp_n[1] /= u_perp_c;
          u_perp_n[2] /= u_perp_c;
        
          u_t = log_law(u_perp_c, d, solver->nu * solver->rho, solver->rho, kE.rough);
        } else u_t = 0;        

        tau = pow(u_t,2) * solver->rho;
        
        tau_x(i,j,k) = tau * u_perp_n[0] * -1.0;
        tau_y(i,j,k) = tau * u_perp_n[1] * -1.0;
        tau_z(i,j,k) = tau * u_perp_n[2] * -1.0;

				if(isnan(u_t)) {
					printf("u_t nan at cell %ld %ld %ld   skipping...\n", i,j,k);
					continue;
				}

				/* testing fabs( ... ) */
        k(i,j,k) = fabs(pow(u_t,2) / sqrt(kE.C_mu));
                
        /* testing fabs( ... ) */        
        E_limit = fabs(kE.C_mu * pow(k(i,j,k), 1.5) / kE.length);
        E(i,j,k) = max(E_limit, fabs(pow(u_t,3) / (d * kE.vonKarman)));
        
        nu_t(i,j,k) = max(kE.C_mu * pow(k(i,j,k),2) / E(i,j,k),0);

      }
    }
  }
  
  return 0;
}

int kE_wall_shear(struct solver_data *solver) {

  /* extension of velocity to add/subtract wall shears from velocity variable
   * wall shear from south wall on U = 2 * nu * u / (dy^2 * AE) * (1 - AS)  if laminar 
   * for turbulent, we can calculate u* based on k (k = u* ^2 / C)
   * and then tau from u*.  Wall stress = 1/rho * tau/dy 
   * unfortunately, this is cell centered and so will need to be averaged for cells i and i+1 */
  long int i,j,k;
  double ws_u, ws_d, tau_u, tau_d, delv;
  
#pragma omp parallel for shared(solver) private(i,j,k,ws_u,ws_d,tau_u,tau_d,delv)  
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {


        if( ( ((FV(i,j,k) > solver->emf)   && (FV(i,j,k) < 1-solver->emf)) ||
              ((FV(i+1,j,k) > solver->emf) && (FV(i+1,j,k) < 1-solver->emf)) ) &&
              AE(i,j,k) > solver->emf) { 
          

          tau_u = tau_x(i,j,k);
          tau_d = tau_x(i+1,j,k);
          
          ws_u  = tau_u  / (solver->rho * DELY) * fabs((1-AN(i,j,k)) - (1-AN(i,j-1,k)));
          ws_u += tau_u / (solver->rho * DELZ) * fabs((1-AT(i,j,k)) - (1-AT(i,j,k-1)));
          
          ws_d  = tau_d / (solver->rho * DELY) * fabs((1-AN(i+1,j,k)) - (1-AN(i+1,j-1,k)));
          ws_d += tau_d / (solver->rho * DELZ) * fabs((1-AT(i+1,j,k)) - (1-AT(i+1,j,k-1)));
        
          delv = 0.5 * (ws_u + ws_d) * solver->delt;
          U(i,j,k) += delv;
          
          /* if(delv * U(i,j,k) > solver->emf) {
          	printf("warning: negative U wall shear. i j k delv U(i,j,k) %ld %ld %ld %lf %lf\n", i,j,k, delv, U(i,j,k));
          } else {
            if(fabs(delv) > fabs(U(i,j,k)) && !isnan(delv)) U(i,j,k) = 0;
            else if(fabs(delv) > solver->emf && !isnan(delv)) U(i,j,k) += delv;
          } */

        }

        if( ( ((FV(i,j,k) > solver->emf && FV(i,j,k) < 1-solver->emf)) ||
              ((FV(i,j+1,k) > solver->emf && FV(i,j+1,k) < 1-solver->emf)) ) &&
                AN(i,j,k) > solver->emf) { 
          

          tau_u = tau_y(i,j,k);
          tau_d = tau_y(i,j+1,k);
          
          ws_u  = tau_u / (solver->rho * DELX) * fabs((1-AE(i,j,k)) - (1-AE(i-1,j,k)));
          ws_u += tau_u / (solver->rho * DELZ) * fabs((1-AT(i,j,k)) - (1-AT(i,j,k-1)));
          
          ws_d  = tau_d / (solver->rho * DELX) * fabs((1-AE(i,j+1,k)) - (1-AE(i-1,j+1,k)));
          ws_d += tau_d / (solver->rho * DELZ) * fabs((1-AT(i,j+1,k)) - (1-AT(i,j+1,k-1)));
        
          delv = 0.5 * (ws_u + ws_d) * solver->delt;
          V(i,j,k) += delv;
          
          /*if(delv * V(i,j,k) > solver->emf) {
          	printf("warning: negative U wall shear. i j k delv V(i,j,k) %ld %ld %ld %lf %lf\n", i,j,k, delv, V(i,j,k));
          } else {
            if(fabs(delv) > fabs(V(i,j,k)) && !isnan(delv)) V(i,j,k) = 0;
          	else if(fabs(delv) > solver->emf && !isnan(delv)) V(i,j,k) += delv;
          }*/
        }

        if( ( ((FV(i,j,k) > solver->emf && FV(i,j,k) < 1-solver->emf)) ||
              ((FV(i,j,k+1) > solver->emf && FV(i,j,k+1) < 1-solver->emf)) ) &&
                AT(i,j,k) > solver->emf) { 
          

          tau_u = tau_z(i,j,k);
          tau_d = tau_z(i,j,k+1);
          
          ws_u  = tau_u / (solver->rho * DELX) * fabs((1-AE(i,j,k)) - (1-AE(i-1,j,k)));
          ws_u += tau_u / (solver->rho * DELY) * fabs((1-AN(i,j,k)) - (1-AN(i,j-1,k)));
          
          ws_d  = tau_d / (solver->rho * DELX) * fabs((1-AE(i,j,k+1)) - (1-AE(i-1,j,k+1)));
          ws_d += tau_d / (solver->rho * DELY) * fabs((1-AN(i,j,k+1)) - (1-AN(i,j-1,k+1)));
        
          delv = 0.5 * (ws_u + ws_d) * solver->delt;
          W(i,j,k) += delv;
          
          /*
          if(delv * W(i,j,k) > solver->emf) {
          	printf("warning: negative U wall shear. i j k delv W(i,j,k) %ld %ld %ld %lf %lf\n", i,j,k, delv, W(i,j,k));
          } else {
            if(fabs(delv) > fabs(W(i,j,k)) && !isnan(delv)) W(i,j,k) = 0;
          	else if(fabs(delv) > solver->emf && !isnan(delv)) W(i,j,k) += delv;
          }*/          
                    
        }


      }
    }
  }

  return 0;
}

int kE_loop_explicit(struct solver_data *solver) {
  /* calculate k, E and nut at each timestep */

  double delk, delE, Production, Diffusion_k, Diffusion_E, nu_t, nu_eff, E_limit;
  double AEdkdx, ANdkdy, ATdkdz, AEdEdx, ANdEdy, ATdEdz;
  double dvdx, dudy, dudz, dwdx, dvdz, dwdy;
  long int i,j,k;
  
#pragma omp parallel for shared(solver) private(i,j,k, \
              delk, delE, Production, Diffusion_k, Diffusion_E, nu_t, nu_eff, E_limit, \
              AEdkdx, ANdkdy, ATdkdz, AEdEdx, ANdEdy, ATdEdz, \
              dvdx, dudy, dudz, dwdx, dvdz, dwdy)  
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        /* exit conditions */
        if(FV(i,j,k) <= (1-solver->emf) || VOF(i,j,k) < solver->emf)
          continue;
      	if(N_VOF(i,j,k) != 0) continue;

/*    K central difference - Area fractions may need to be fixed *
        AEdkdx = 0.5 * (AE(i,j,k) + AE(i-1,j,k)) * RDX * 
                 ( (k_N(i+1,j,k) + k_N(i,j,k)) - (k_N(i,j,k) + k_N(i-1,j,k)) ) / 2;
        ANdkdy = 0.5 * (AN(i,j,k) + AN(i,j-1,k)) * RDY *
                 ( (k_N(i,j+1,k) + k_N(i,j,k)) - (k_N(i,j,k) + k_N(i,j-1,k)) ) / 2;
        ATdkdz = 0.5 * (AT(i,j,k) + AT(i,j,k-1)) * RDZ *
                 ( (k_N(i,j,k+1) + k_N(i,j,k)) - (k_N(i,j,k) + k_N(i,j,k-1)) ) / 2; */

/*    K upwind */
        if(U(i-1,j,k) >= 0) 
          AEdkdx = AE(i-1,j,k) * RDX * 
                 (k_N(i,j,k) - k_N(i-1,j,k));
        else
          AEdkdx = AE(i,j,k) * RDX * 
                 (k_N(i,j,k) - k_N(i+1,j,k));
                 
        if(V(i,j-1,k) >= 0)                       
          ANdkdy = AN(i,j-1,k) * RDY *
                 (k_N(i,j,k) - k_N(i,j-1,k));
        else
          ANdkdy = AN(i,j,k) * RDY * 
                 (k_N(i,j,k) - k_N(i,j+1,k));
        
        if(W(i,j,k-1) >= 0)                     
          ATdkdz = AT(i,j,k-1) * RDZ *
                 (k_N(i,j,k) - k_N(i,j,k-1)); 
        else
          ATdkdz = AT(i,j,k) * RDZ *
                 (k_N(i,j,k) - k_N(i,j,k+1));                 
                 
        nu_t = nu_t_N(i,j,k);
        nu_eff = nu_t + solver->nu / kE.sigma_k;
        nu_eff = max(nu_eff, solver->nu);

        /* Diffusion_k = (1.0 / FV(i,j,k)) * 
                    ( RDX * (nu_eff * AEdkdx) +
                      RDY * (nu_eff * ANdkdy) +
                      RDZ * (nu_eff * ATdkdz)); */
        Diffusion_k = pow(RDX,2) * ( AE(i,  j,k) * (k_N(i+1,j,k) - k_N(i,j,k)) - 
                                     AE(i-1,j,k) * (k_N(i  ,j,k) - k_N(i-1,j,k)) );
        Diffusion_k +=pow(RDY,2) * ( AN(i,  j,k) * (k_N(i,j+1,k) - k_N(i,j,k)) - 
                                     AN(i,j-1,k) * (k_N(i  ,j,k) - k_N(i,j-1,k)) );
        Diffusion_k +=pow(RDZ,2) * ( AT(i,  j,k) * (k_N(i,j,k+1) - k_N(i,j,k)) - 
                                     AT(i,j,k-1) * (k_N(i  ,j,k) - k_N(i,j,k-1)) );
        Diffusion_k *= nu_eff * (1.0 / FV(i,j,k));

        dvdx = (V(i+1,j,k) + V(i+1,j-1,k) + V(i,j,k) + V(i,j-1,k))/4 - 
               (V(i-1,j,k) + V(i-1,j-1,k) + V(i,j,k) + V(i,j-1,k))/4;
        dvdx *= RDX;

        dudy = (U(i,j+1,k) + U(i-1,j+1,k) + U(i,j,k) + U(i-1,j,k))/4 -
               (U(i,j-1,k) + U(i-1,j-1,k) + U(i,j,k) + U(i-1,j,k))/4;
        dudy *= RDY;

        dudz = (U(i,j,k+1) + U(i-1,j,k+1) + U(i,j,k) + U(i-1,j,k))/4 -
               (U(i,j,k-1) + U(i-1,j,k-1) + U(i,j,k) + U(i-1,j,k))/4;
        dudz *= RDZ;

        dwdx = (W(i+1,j,k) + W(i+1,j,k-1) + W(i,j,k) + W(i,j,k-1))/4 -
               (W(i-1,j,k) + W(i-1,j,k-1) + W(i,j,k) + W(i,j,k-1))/4;
        dwdx *= RDX;

        dvdz = (V(i,j,k+1) + V(i,j-1,k+1) + V(i,j,k) + V(i,j-1,k))/4 -
               (V(i,j,k-1) + V(i,j-1,k-1) + V(i,j,k) + V(i,j-1,k))/4;
        dvdz *= RDZ;

        dwdy = (W(i,j+1,k) + W(i,j+1,k-1) + W(i,j,k) + W(i,j,k-1))/4 -
               (W(i,j-1,k) + W(i,j-1,k-1) + W(i,j,k) + W(i,j,k-1))/4;
        dwdy *= RDY;
        
        Production = nu_t * (1.0 / FV(i,j,k)) * 
                     ( 2 * 0.5 * (AE(i,j,k) + AE(i-1,j,k)) * pow((U(i,j,k) - U(i-1,j,k)) * RDX, 2) + /* 10/26 fixed AE factor */
                       2 * 0.5 * (AN(i,j,k) + AN(i,j-1,k)) * pow((V(i,j,k) - V(i,j-1,k)) * RDY, 2) +
                       2 * 0.5 * (AT(i,j,k) + AT(i,j,k-1)) * pow((W(i,j,k) - W(i,j,k-1)) * RDZ, 2) +
                       (dvdx + dudy) * (AE(i,j,k) * dvdx + AN(i,j,k) * dudy) +
                       (dudz + dwdx) * (AT(i,j,k) * dudz + AE(i,j,k) * dwdx) +
                       (dvdz + dwdy) * (AT(i,j,k) * dvdz + AN(i,j,k) * dwdy) );

        delk = ( -1.0 / FV(i,j,k)) * 
                ( fabs((U(i,j,k) + U(i-1,j,k)) / 2) * AEdkdx + /* TESTING FABS 03/07/16 */
                  fabs((V(i,j,k) + V(i,j-1,k)) / 2) * ANdkdy + 
                  fabs((W(i,j,k) + W(i,j,k-1)) / 2) * ATdkdz ) +
                Production + Diffusion_k - E_N(i,j,k);

/* E central difference *
        AEdEdx = 0.5 * (AE(i,j,k) + AE(i-1,j,k)) * RDX * 
                 ( (E_N(i+1,j,k) + E_N(i,j,k)) - (E_N(i,j,k) + E_N(i-1,j,k)) ) / 2;
        ANdEdy = 0.5 * (AN(i,j,k) + AN(i,j-1,k)) * RDY * 
                 ( (E_N(i,j+1,k) + E_N(i,j,k)) - (E_N(i,j,k) + E_N(i,j-1,k)) ) / 2;
        ATdEdz = 0.5 * (AT(i,j,k) + AT(i,j,k-1)) * RDZ *
                 ( (E_N(i,j,k+1) + E_N(i,j,k)) - (E_N(i,j,k) + E_N(i,j,k-1)) ) / 2; */
                 
/*    E upwind */
        if(U(i-1,j,k) >= 0) 
          AEdEdx = AE(i-1,j,k) * RDX * 
                 (E_N(i,j,k) - E_N(i-1,j,k));
        else
          AEdEdx = AE(i,j,k) * RDX * 
                 (E_N(i,j,k) - E_N(i+1,j,k));
                 
        if(V(i,j-1,k) >= 0)                       
          ANdEdy = AN(i,j-1,k) * RDY *
                 (E_N(i,j,k) - E_N(i,j-1,k));
        else
          ANdEdy = AN(i,j,k) * RDY * 
                 (E_N(i,j,k) - E_N(i,j+1,k));
        
        if(W(i,j,k-1) >= 0)                     
          ATdEdz = AT(i,j,k-1) * RDZ *
                 (E_N(i,j,k) - E_N(i,j,k-1)); 
        else
          ATdEdz = AT(i,j,k) * RDZ *
                 (E_N(i,j,k) - E_N(i,j,k+1));  
                 
        nu_eff = nu_t + solver->nu / kE.sigma_E;
        nu_eff = max(nu_eff, solver->nu);


        Diffusion_E = pow(RDX,2) * ( AE(i,  j,k) * (E_N(i+1,j,k) - E_N(i,j,k)) - 
                                     AE(i-1,j,k) * (E_N(i  ,j,k) - E_N(i-1,j,k)) );
        Diffusion_E +=pow(RDY,2) * ( AN(i,  j,k) * (E_N(i,j+1,k) - E_N(i,j,k)) - 
                                     AN(i,j-1,k) * (E_N(i  ,j,k) - E_N(i,j-1,k)) );
        Diffusion_E +=pow(RDZ,2) * ( AT(i,  j,k) * (E_N(i,j,k+1) - E_N(i,j,k)) - 
                                     AT(i,j,k-1) * (E_N(i  ,j,k) - E_N(i,j,k-1)) );
        Diffusion_E *= nu_eff * (1.0 / FV(i,j,k));
                      
        delE = ( -1.0 / FV(i,j,k)) *
                ( fabs((U(i,j,k) + U(i-1,j,k)) / 2) * AEdEdx +  /* TESTING FABS 03//07/16 */
                  fabs((V(i,j,k) + V(i,j-1,k)) / 2) * ANdEdy + 
                  fabs((W(i,j,k) + W(i,j,k-1)) / 2) * ATdEdz ) + 
               ( E_N(i,j,k) / k_N(i,j,k) ) * 
                  ( kE.C1E * Production - kE.C2E * E_N(i,j,k)) +
               Diffusion_E;
        
        /* TESTING: rescuing k and E from negative / nan values */
               
        if(isnan(delk) || (k_N(i,j,k) + delk * solver->delt) < 0.0000001) delk = 0;
        k(i,j,k) = k_N(i,j,k) + delk * solver->delt;
        E_limit = kE.C_mu * pow(k(i,j,k), 1.5) / kE.length;
                
        if(isnan(delE) || (E_N(i,j,k) + delE * solver->delt) < 0.0000001) delE = 0;        
        E(i,j,k) = max(E_limit, E_N(i,j,k) + delE * solver->delt);
        nu_t(i,j,k) = max(kE.C_mu * pow(k(i,j,k),2) / E(i,j,k),0);

#ifdef DEBUG
        if(E(i,j,k) < 0) {
          /* printf("warning: negative E at cell %ld %ld %ld\n",i,j,k);
          printf("k %e   E %e\n",k(i,j,k),E(i,j,k)); */
        }
        if(nu_t(i,j,k) > 0.1) {
          /* printf("warning: very high viscocity at cell %ld %ld %ld\n",i,j,k);
          printf("k %e   E %e\n",k(i,j,k),E(i,j,k)); */
        } 
#endif

      }
    }
  }

  kE_boundaries(solver);
  kE_copy(solver);

  return 0;
}

int kE_loop(struct solver_data *solver) {
  /* calculate k, E and nut at each timestep */
  /* implicit method */

  double delk, delE, Production, Diffusion_k, Diffusion_E, nu_t, nu_eff;
  double AEdkdx, ANdkdy, ATdkdz, AEdEdx, ANdEdy, ATdEdz;
  double dvdx, dudy, dudz, dwdx, dvdz, dwdy;
  long int i,j,k;
  double res;
  int iter = 0;
  double k_i, E_i;
  double relax; 
  const int niter = 500;
  
  do {
   
  if(iter > niter) {
    printf("turbulence did not converge. final residual %lf\n", res);
    break;
  }
  res = 0.0; 
  
  for(i=1; i<IMAX-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        /* exit conditions */
        if(FV(i,j,k) <= (1-solver->emf) || VOF(i,j,k) < solver->emf)
          continue;
          
        if(N_VOF(i,j,k) != 0) continue;

/*    K central difference - Area fractions may need to be fixed */
        AEdkdx = 0.5 * (AE(i,j,k) + AE(i-1,j,k)) * RDX * 
                 ( (k(i+1,j,k) + k(i,j,k)) - (k(i,j,k) + k(i-1,j,k)) ) / 2;
        ANdkdy = 0.5 * (AN(i,j,k) + AN(i,j-1,k)) * RDY *
                 ( (k(i,j+1,k) + k(i,j,k)) - (k(i,j,k) + k(i,j-1,k)) ) / 2;
        ATdkdz = 0.5 * (AT(i,j,k) + AT(i,j,k-1)) * RDZ *
                 ( (k(i,j,k+1) + k(i,j,k)) - (k(i,j,k) + k(i,j,k-1)) ) / 2; 

/*    K upstream *
        AEdkdx = AE(i,j,k) * RDX * 
                 (k(i,j,k) - k(i-1,j,k));
        ANdkdy = AN(i,j,k) * RDY *
                 ( (k(i,j+1,k) + k(i,j,k)) - (k(i,j,k) + k(i,j-1,k)) ) / 2;
        ATdkdz = AT(i,j,k) * RDZ *
                 ( (k(i,j,k+1) + k(i,j,k)) - (k(i,j,k) + k(i,j,k-1)) ) / 2; */
                 
        nu_t = kE.C_mu * pow(k(i,j,k),2) / E(i,j,k);
        nu_eff = nu_t + solver->nu / kE.sigma_k;
        nu_eff = max(nu_eff, solver->nu);

        /* Diffusion_k = (1.0 / FV(i,j,k)) * 
                    ( RDX * (nu_eff * AEdkdx) +
                      RDY * (nu_eff * ANdkdy) +
                      RDZ * (nu_eff * ATdkdz)); */
        Diffusion_k = pow(RDX,2) * ( AE(i,  j,k) * (k(i+1,j,k) - k(i,j,k)) - 
                                     AE(i-1,j,k) * (k(i  ,j,k) - k(i-1,j,k)) );
        Diffusion_k +=pow(RDY,2) * ( AN(i,  j,k) * (k(i,j+1,k) - k(i,j,k)) - 
                                     AN(i,j-1,k) * (k(i  ,j,k) - k(i,j-1,k)) );
        Diffusion_k +=pow(RDZ,2) * ( AT(i,  j,k) * (k(i,j,k+1) - k(i,j,k)) - 
                                     AT(i,j,k-1) * (k(i  ,j,k) - k(i,j,k-1)) );
        Diffusion_k *= nu_eff * (1.0 / FV(i,j,k));

        dvdx = (V(i+1,j,k) + V(i+1,j-1,k) + V(i,j,k) + V(i,j-1,k))/4 - 
               (V(i-1,j,k) + V(i-1,j-1,k) + V(i,j,k) + V(i,j-1,k))/4;
        dvdx *= RDX;

        dudy = (U(i,j+1,k) + U(i-1,j+1,k) + U(i,j,k) + U(i-1,j,k))/4 -
               (U(i,j-1,k) + U(i-1,j-1,k) + U(i,j,k) + U(i-1,j,k))/4;
        dudy *= RDY;

        dudz = (U(i,j,k+1) + U(i-1,j,k+1) + U(i,j,k) + U(i-1,j,k))/4 -
               (U(i,j,k-1) + U(i-1,j,k-1) + U(i,j,k) + U(i-1,j,k))/4;
        dudz *= RDZ;

        dwdx = (W(i+1,j,k) + W(i+1,j,k-1) + W(i,j,k) + W(i,j,k-1))/4 -
               (W(i-1,j,k) + W(i-1,j,k-1) + W(i,j,k) + W(i,j,k-1))/4;
        dwdx *= RDX;

        dvdz = (V(i,j,k+1) + V(i,j-1,k+1) + V(i,j,k) + V(i,j-1,k))/4 -
               (V(i,j,k-1) + V(i,j-1,k-1) + V(i,j,k) + V(i,j-1,k))/4;
        dvdz *= RDZ;

        dwdy = (W(i,j+1,k) + W(i,j+1,k-1) + W(i,j,k) + W(i,j,k-1))/4 -
               (W(i,j-1,k) + W(i,j-1,k-1) + W(i,j,k) + W(i,j,k-1))/4;
        dwdy *= RDY;
        
        Production = nu_t * (1.0 / FV(i,j,k)) * 
                     ( 2 * AE(i,j,k) * pow((U(i,j,k) - U(i-1,j,k)) * RDX, 2) +
                       2 * AN(i,j,k) * pow((V(i,j,k) - V(i,j-1,k)) * RDY, 2) +
                       2 * AT(i,j,k) * pow((W(i,j,k) - W(i,j,k-1)) * RDZ, 2) +
                       (dvdx + dudy) * (AE(i,j,k) * dvdx + AN(i,j,k) * dudy) +
                       (dudz + dwdx) * (AT(i,j,k) * dudz + AE(i,j,k) * dwdx) +
                       (dvdz + dwdy) * (AT(i,j,k) * dvdz + AN(i,j,k) * dwdy) );

        delk = ( -1.0 / FV(i,j,k)) * 
                ( ((U(i,j,k) + U(i-1,j,k)) / 2) * AEdkdx + 
                  ((V(i,j,k) + V(i,j-1,k)) / 2) * ANdkdy + 
                  ((W(i,j,k) + W(i,j,k-1)) / 2) * ATdkdz ) +
                Production + Diffusion_k - E(i,j,k);

/* E central difference */
        AEdEdx = 0.5 * (AE(i,j,k) + AE(i-1,j,k)) * RDX * 
                 ( (E(i+1,j,k) + E(i,j,k)) - (E(i,j,k) + E(i-1,j,k)) ) / 2;
        ANdEdy = 0.5 * (AN(i,j,k) + AN(i,j-1,k)) * RDY * 
                 ( (E(i,j+1,k) + E(i,j,k)) - (E(i,j,k) + E(i,j-1,k)) ) / 2;
        ATdEdz = 0.5 * (AT(i,j,k) + AT(i,j,k-1)) * RDZ *
                 ( (E(i,j,k+1) + E(i,j,k)) - (E(i,j,k) + E(i,j,k-1)) ) / 2; 
                 
/* E upstream 
        AEdEdx = AE(i,j,k) * RDX * 
                 ( E(i,j,k) - E(i-1,j,k) ) ;
        ANdEdy = AN(i,j,k) * RDY * 
                 ( (E(i,j+1,k) + E(i,j,k)) - (E(i,j,k) + E(i,j-1,k)) ) / 2;
        ATdEdz = AT(i,j,k) * RDZ *
                 ( (E(i,j,k+1) + E(i,j,k)) - (E(i,j,k) + E(i,j,k-1)) ) / 2; */
                 
        nu_eff = nu_t + solver->nu / kE.sigma_E;
        nu_eff = max(nu_eff, solver->nu);

        /************************************************ 
          Diffusion term is wrong.  Should match viscocity term of momentum equation. FIX 
           Should be the second derivative, and work out to something like
           us - 2C + ds 
           Also fix for k *
        Diffusion_E = (1.0 / FV(i,j,k)) * 
                    ( RDX * (nu_eff * AEdEdx) +
                      RDY * (nu_eff * ANdEdy) +
                      RDZ * (nu_eff * ATdEdz)); */

        Diffusion_E = pow(RDX,2) * ( AE(i,  j,k) * (E(i+1,j,k) - E(i,j,k)) - 
                                     AE(i-1,j,k) * (E(i  ,j,k) - E(i-1,j,k)) );
        Diffusion_E +=pow(RDY,2) * ( AN(i,  j,k) * (E(i,j+1,k) - E(i,j,k)) - 
                                     AN(i,j-1,k) * (E(i  ,j,k) - E(i,j-1,k)) );
        Diffusion_E +=pow(RDZ,2) * ( AT(i,  j,k) * (E(i,j,k+1) - E(i,j,k)) - 
                                     AT(i,j,k-1) * (E(i  ,j,k) - E(i,j,k-1)) );
        Diffusion_E *= nu_eff * (1.0 / FV(i,j,k));
                      
        delE = ( -1.0 / FV(i,j,k)) *
                ( ((U(i,j,k) + U(i-1,j,k)) / 2) * AEdEdx + 
                  ((V(i,j,k) + V(i,j-1,k)) / 2) * ANdEdy + 
                  ((W(i,j,k) + W(i,j,k-1)) / 2) * ATdEdz ) +
               ( E(i,j,k) / k(i,j,k) ) * 
                  ( kE.C1E * Production - kE.C2E * E(i,j,k)) +
               Diffusion_E;

        k_i = k(i,j,k);
        E_i = E(i,j,k);

        k(i,j,k) = k_N(i,j,k) + delk * solver->delt;
        E(i,j,k) = E_N(i,j,k) + delE * solver->delt;

        res = max(fabs(k_i-k(i,j,k)),res);
        res = max(fabs(E_i-E(i,j,k)),res);
        
        /* under relaxation */
        /* relax = 0.5 * (pow(DELX,2) * pow(DELY,2) * pow(DELZ,2)) / 
                       (nu_eff * (pow(DELX,2) + pow(DELY,2) + pow(DELZ,2)) ); */
        relax = 0.5;
        k(i,j,k) = k_i + (k(i,j,k) - k_i) * relax;
        E(i,j,k) = E_i + (E(i,j,k) - E_i) * relax;
        nu_t(i,j,k) = max(kE.C_mu * pow(k(i,j,k),2) / E(i,j,k),0);

        /* if(E(i,j,k) < 0) {
          printf("warning: negative E at cell %ld %ld %ld\n",i,j,k);
        } */

      }
    }
  }


  iter++;
  } while(res > 1e-9);

  if(iter < niter) {
    printf("turbulence converged in %d iterations. final residual: %e\n", iter, res);
  }

  kE_boundaries(solver);
  kE_copy(solver);

  return 0;
}
int kE_kill(struct solver_data *solver) {

  free(kE.k);
  free(kE.E);
  free(kE.nu_t);
  free(kE_n.k);
  free(kE_n.E);
  free(kE_n.nu_t);
  
  solver->mesh->turbulence_model = NULL;

  return 0;
}
