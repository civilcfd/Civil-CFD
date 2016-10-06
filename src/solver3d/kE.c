/*
 * kE.c 
 *
 * Implements a k-epsilon turbulence model */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <petscksp.h>
#include <libxml/xmlwriter.h>
#include <libxml/xpath.h>

#include "vtk.h"
#include "solver.h"
#include "readfile.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "vof_boundary.h"
#include "solver_mpi.h"
 
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

int kE_length() {
  kE.length = kE.length_scale * kE.raw_length;
  kE_n.length = kE.length;

  return 0;
}

int kE_set_value(char *param, int dims, 
                   double *vector) {


  if(param == NULL || vector == NULL) {
    printf("error: null values passed to kE_set_value\n");
    return(1);
  }

  if(strcmp(param, "length_scale")==0) {
    if(dims != 1) {
      printf("error in source file: length_scale requires 1 arguments\n");
      return(1);
    }

    kE.length_scale = vector[0];
    kE_n.length_scale = kE.length_scale;
    kE_length();
  }
  else if (strcmp(param, "length")==0) {
    if(dims != 1) {
      printf("error in source file: length requires 1 arguments\n");
      return(1);
    }

    kE.raw_length = vector[0];
    kE_n.raw_length = kE.raw_length;
    kE_length();

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

int kE_read_xml(char *filename) {
  const char *kE_properties_double[] = {"length_scale", "length", "rough", "end" };
  xmlXPathContext *xpathCtx;
  xmlDoc *doc;
  double vector[3];
  char buf[256];
  int i;

  xmlInitParser();
  LIBXML_TEST_VERSION
 
  doc = xmlParseFile(filename);
  if(doc == NULL) {
    printf("Could not open %s\n",filename);
    return 1;
  }

  xpathCtx = xmlXPathNewContext(doc);
  if(doc == NULL) {
    printf("Could not create xpath context\n");
    return 1;
  }

  i = 0;
  while(strcmp("end", kE_properties_double[i]) != 0) {
    sprintf(buf, "/Case/Solver/Turbulence/%s", kE_properties_double[i]);

    if(!read_xmlpath_double(&vector[0], buf, xpathCtx)) {
      printf("Could not evaluate %s\n", kE_properties_double[i]);
    }
    else {
      printf("Read %s: %lf\n", buf, vector[0]);
      kE_set_value(kE_properties_double[i], 1, vector);
    }
 
    i++;
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

int kE_write_xml(void *writer_ptr) {
  xmlTextWriterPtr writer = (xmlTextWriterPtr) writer_ptr;


  xmlTextWriterStartElement(writer, BAD_CAST "Turbulence");

  xmlTextWriterWriteFormatElement(writer, BAD_CAST "length_scale", "%e", kE.length_scale);  
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "length", "%e", kE.raw_length);  
  xmlTextWriterWriteFormatElement(writer, BAD_CAST "rough", "%e", kE.rough);  

  xmlTextWriterEndElement(writer);

  return 0;
}

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
  fprintf(fp,"length %le\n",kE.raw_length);
  fprintf(fp,"rough %le\n",kE.rough);

  fclose(fp);
  return 0;
}

int kE_load_values(struct solver_data *solver) {
  
  kE_set_internal(solver, 0.0, 0.0);
  csv_read_k(solver->mesh,solver->t);
  csv_read_E(solver->mesh,solver->t);

  return 0;
}

int kE_broadcast(struct solver_data *solver) {
	
	MPI_Bcast(&kE.raw_length, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&kE.length_scale, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&kE.rough, 1, MPI_LONG_INT, 0, MPI_COMM_WORLD);
  kE_length();

	return 0;
}

int kE_setup(struct solver_data *solver) {

  solver->turbulence_init = kE_init;
  solver->turbulence_loop = kE_loop_explicit;
  solver->turbulence_kill = kE_kill;
  solver->turbulence_read = kE_read_xml;
  solver->turbulence_write = kE_write;
  solver->wall_shear = kE_wall_shear;
  solver->turbulence_nu   = kE_nu;
  solver->turbulence_load_values = kE_load_values;
    
  solver->mesh->turbulence_model = &kE;
  
  kE.C1E = kE_n.C1E = 1.44;
  kE.C2E = kE_n.C2E = 1.92;
  kE.C3E = kE_n.C3E = -0.33;
  kE.C_mu = kE_n.C_mu = 0.09;
  kE.sigma_k = kE_n.sigma_k = 1.0;
  kE.sigma_E = kE_n.sigma_E = 1.3;
  kE.vonKarman = kE_n.vonKarman = 0.41;
  kE.length_scale = 0.038;
  kE.raw_length = min(IMAX * DELX, JMAX * DELY);
  kE.raw_length = min(kE.length, KMAX * DELZ);
  kE_length();

  kE.rough = kE_n.rough = 0.00161; /* concrete */
  
  return 0;
}

int kE_init(struct solver_data *solver) {
  long int size;
  
  size = mesh_mpi_space(solver->mesh);

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
  
  if(!solver->rank) {
    kE_set_internal(solver, 0.0001, 0);
  }
  kE_length();
 
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
  long int l,m,n,irange;
  double E_limit;
  
  E_limit = kE.C_mu * pow(k, 1.5) / kE.length;

  if(!solver->rank) irange = IMAX;
  else irange = IRANGE;

  for(l=0; l<irange; l++) {
    for(m=0; m<JMAX; m++) {
      for(n=0; n<KMAX; n++) {
        k(l,m,n) = k;
        E(l,m,n) = max(E,E_limit);
        
        /* set shear stresses to zero upon initialization */
        tau_x(l,m,n) = 0;
        tau_y(l,m,n) = 0;
        tau_z(l,m,n) = 0;
        
        //nu_t(l,m,n) = 0;
        
        nu_t(l,m,n) = max(kE.C_mu * pow(k(l,m,n),2) / E(l,m,n),0);
      }
    }
  }
      
  //kE_boundaries(solver);
  kE_copy(solver);
  
  return 0;
}

int kE_copy(struct solver_data *solver) {
  long int l,m,n;
 
  for(l=0; l<IRANGE; l++) {
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

  /* sets neumann boundaries on k and epsilon */
  long int i,j,k, l, m, n;

  for(i=0; i<IRANGE; i++) {
    for(j=0; j<JMAX; j++) {
      k(i,j,0) = k(i,j,1);
      k(i,j,KMAX-1) = k(i,j,KMAX-2);

      E(i,j,0) = E(i,j,1);
      E(i,j,KMAX-1) = E(i,j,KMAX-2);

      nu_t(i,j,0) = nu_t(i,j,1);
      nu_t(i,j,KMAX-1) = nu_t(i,j,KMAX-2);

    }
  }
  for(i=0; i<IRANGE; i++) {
    for(k=0; k<KMAX; k++) {
      k(i,0,k) = k(i,1,k);
      k(i,JMAX-1,k) = k(i,JMAX-2,k);

      E(i,0,k) = E(i,1,k);
      E(i,JMAX-1,k) = E(i,JMAX-2,k);

      nu_t(i,0,k) = nu_t(i,1,k);
      nu_t(i,JMAX-1,k) = nu_t(i,JMAX-2,k);

    }
  }
  for(j=0; j<JMAX; j++) {
    for(k=0; k<KMAX; k++) {
      if(ISTART == 0) {
        k(0,j,k) = k(1,j,k);
        E(0,j,k) = E(1,j,k);
        nu_t(0,j,k) = nu_t(1,j,k);
      } else if(ISTART + IRANGE == IMAX) {
        k(IRANGE-1,j,k) = k(IRANGE-2,j,k);
        E(IRANGE-1,j,k) = E(IRANGE-2,j,k);
        nu_t(IRANGE-1,j,k) = nu_t(IRANGE-2,j,k);
      }

    }
  }

  for(i=0; i<IRANGE; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
      
      	tau_x(i,j,k) = 0;
      	tau_y(i,j,k) = 0;
      	tau_z(i,j,k) = 0;
      }
    }
  }

  kE_special_boundaries(solver);
  
  kE_tau(solver);

	/* free surface boundaries - the free surface acts as a von Neumann boundary */
	/* also takes care of cells outside of the solution space */
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
      
      	if(FV(i,j,k) < solver->emf || VOF(i,j,k) < solver->emf || N_VOF(i,j,k) > 6) { 
      		k(i,j,k) = 0;
      		E(i,j,k) = 0;
      		nu_t(i,j,k) = 0;
      		
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
  /* Calculate u_t based on method in Wilcox text */

  int iter = 0;
  double u_t, u_t_n, res, B, y_plus, kr_plus;

  /* initial guess */
  y_plus = 11;
  B = 9;
  res = 1;
  iter = 0;
  u_t_n = fabs(u * kE.vonKarman / (log(y_plus) + B));

  while(res > fabs(u) / 100000 && iter < 50) {
    y_plus  = max(rho * u_t_n * d / mu, 11);
    kr_plus = u_t_n * rough * rho / mu;
    B       = min(8.5 + (1/kE.vonKarman) * log(1/kr_plus),9);
    B       = max(B, log(y_plus) / 2 * -1.0);
    u_t     = fabs(u * kE.vonKarman / (log(y_plus) + B));

    res = fabs(u_t - u_t_n);
    u_t_n = u_t;  

    iter++;
  }

  if(res > fabs(u) / 1000) return 0;
  else return u_t;
}

int kE_tau(struct solver_data *solver) {
  long int i,j,k, im1, jm1, km1;
  double d, u_t, wall_n[3], u_c[3], mag, tau;
  double u_perp_n[3], u_perp_c, u_parr_c, u_parr[3];
  double E_limit;
  const double del[3] = { DELX, DELY, DELZ };
  
  for(i=1; i<IRANGE; i++) {
    for(j=1; j<JMAX; j++) {
      for(k=1; k<KMAX; k++) {
      
        if(FV(i,j,k) < solver->emf || VOF(i,j,k) < solver->emf)
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
        
        if(mag < solver->emf) continue;
        
        wall_n[0] /= mag;
        wall_n[1] /= mag;
        wall_n[2] /= mag;

        d = fabs(wall_n[0]) * del[0] + fabs(wall_n[1]) * del[1] + fabs(wall_n[2]) * del[2];
        d /= 2;

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

        k(i,j,k) = fabs(pow(u_t,2) / sqrt(kE.C_mu));
                  
        E_limit = fabs(kE.C_mu * pow(k(i,j,k), 1.5) / kE.length);
        E(i,j,k) = max(E_limit, fabs(pow(u_t,3) / (d * kE.vonKarman)));
        
        nu_t(i,j,k) = max(tau * d / (u_perp_c*solver->rho),0);

      }
    }
  }
  
  return 0;
}

int kE_wall_shear(struct solver_data *solver) {

  /* extension of velocity to add/subtract wall shears from velocity variable
   * for turbulent, we can calculate u* based on k (k = u* ^2 / C)
   * and then tau from u*.  Wall stress = 1/rho * tau/dy 
   * unfortunately, this is cell centered and so will need to be averaged for cells i and i+1 */
  long int i,j,k;
  double ws_u, ws_d, tau_u, tau_d, delv;
   
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        tau_u = tau_x(i,j,k);
        tau_d = tau_x(i+1,j,k);
        
        if(fabs(tau_u) + fabs(tau_d) > solver->emf) {
          
          ws_u  = tau_u  / (solver->rho * DELY) * fabs((1-AN(i,j,k)) - (1-AN(i,j-1,k)));
          ws_u += tau_u / (solver->rho * DELZ) * fabs((1-AT(i,j,k)) - (1-AT(i,j,k-1)));
          
          ws_d  = tau_d / (solver->rho * DELY) * fabs((1-AN(i+1,j,k)) - (1-AN(i+1,j-1,k)));
          ws_d += tau_d / (solver->rho * DELZ) * fabs((1-AT(i+1,j,k)) - (1-AT(i+1,j,k-1)));
        
          delv = 0.5 * (ws_u + ws_d) * solver->delt;
          if(!isnan(delv)) U(i,j,k) += delv;
        }

        tau_u = tau_y(i,j,k);
        tau_d = tau_y(i,j+1,k);
        
        if(fabs(tau_u) + fabs(tau_d) > solver->emf) {
          
          ws_u  = tau_u / (solver->rho * DELX) * fabs((1-AE(i,j,k)) - (1-AE(i-1,j,k)));
          ws_u += tau_u / (solver->rho * DELZ) * fabs((1-AT(i,j,k)) - (1-AT(i,j,k-1)));
          
          ws_d  = tau_d / (solver->rho * DELX) * fabs((1-AE(i,j+1,k)) - (1-AE(i-1,j+1,k)));
          ws_d += tau_d / (solver->rho * DELZ) * fabs((1-AT(i,j+1,k)) - (1-AT(i,j+1,k-1)));
        
          delv = 0.5 * (ws_u + ws_d) * solver->delt;
          if(!isnan(delv)) V(i,j,k) += delv;
        }

        tau_u = tau_z(i,j,k);
        tau_d = tau_z(i,j,k+1);
        
        if(fabs(tau_u) + fabs(tau_d) > solver->emf) {
          
          ws_u  = tau_u / (solver->rho * DELX) * fabs((1-AE(i,j,k)) - (1-AE(i-1,j,k)));
          ws_u += tau_u / (solver->rho * DELY) * fabs((1-AN(i,j,k)) - (1-AN(i,j-1,k)));
          
          ws_d  = tau_d / (solver->rho * DELX) * fabs((1-AE(i,j,k+1)) - (1-AE(i-1,j,k+1)));
          ws_d += tau_d / (solver->rho * DELY) * fabs((1-AN(i,j,k+1)) - (1-AN(i,j-1,k+1)));
        
          delv = 0.5 * (ws_u + ws_d) * solver->delt;
          if(!isnan(delv)) W(i,j,k) += delv;
        } 


      }
    }
  }

  return 0;
}

int kE_loop_explicit(struct solver_data *solver) {
  /* calculate k, E and nut at each timestep */

  double delk, delE, Production, Diffusion_k, Diffusion_E, nu_t, nu_eff, E_limit;
  double dkdx, dkdy, dkdz, dEdx, dEdy, dEdz;
  double dvdx, dudy, dudz, dwdx, dvdz, dwdy;
  long int i,j,k;
  
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        /* exit conditions */
        if(FV(i,j,k) <= (1-solver->emf) || VOF(i,j,k) < solver->emf)
          continue;
      	if(N_VOF(i,j,k) != 0) continue;

/*    K upwind */
        if(U(i-1,j,k) >= 0) 
          dkdx = RDX * 
                 (k_N(i,j,k) - k_N(i-1,j,k));
        else
          dkdx = RDX * 
                 (k_N(i,j,k) - k_N(i+1,j,k));
                 
        if(V(i,j-1,k) >= 0)                       
          dkdy = RDY *
                 (k_N(i,j,k) - k_N(i,j-1,k));
        else
          dkdy = RDY * 
                 (k_N(i,j,k) - k_N(i,j+1,k));
        
        if(W(i,j,k-1) >= 0)                     
          dkdz = RDZ *
                 (k_N(i,j,k) - k_N(i,j,k-1)); 
        else
          dkdz = RDZ *
                 (k_N(i,j,k) - k_N(i,j,k+1));                 
                 
        nu_t = nu_t_N(i,j,k);
        nu_eff = nu_t + solver->nu / kE.sigma_k;
        nu_eff = max(nu_eff, solver->nu);

        Diffusion_k = pow(RDX,2) * ( (k_N(i+1,j,k) - k_N(i,j,k)) - 
                                     (k_N(i  ,j,k) - k_N(i-1,j,k)) );
        Diffusion_k +=pow(RDY,2) * ( (k_N(i,j+1,k) - k_N(i,j,k)) - 
                                     (k_N(i  ,j,k) - k_N(i,j-1,k)) );
        Diffusion_k +=pow(RDZ,2) * ( (k_N(i,j,k+1) - k_N(i,j,k)) - 
                                     (k_N(i  ,j,k) - k_N(i,j,k-1)) );
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
                     ( pow((U(i,j,k) - U(i-1,j,k)) * RDX, 2) + 
                       pow((V(i,j,k) - V(i,j-1,k)) * RDY, 2) +
                       pow((W(i,j,k) - W(i,j,k-1)) * RDZ, 2) +
                       (dvdx + dudy) * (dvdx + dudy) +
                       (dudz + dwdx) * (dudz + dwdx) +
                       (dvdz + dwdy) * (dvdz + dwdy) );

        delk = ( -1.0 / FV(i,j,k)) * 
                ( fabs((U(i,j,k) + U(i-1,j,k)) / 2) * dkdx + /* TESTING FABS 03/07/16 */
                  fabs((V(i,j,k) + V(i,j-1,k)) / 2) * dkdy + 
                  fabs((W(i,j,k) + W(i,j,k-1)) / 2) * dkdz ) +
                Production + Diffusion_k - E_N(i,j,k);

/*    E upwind */
        if(U(i-1,j,k) >= 0) 
          dEdx = RDX * 
                 (E_N(i,j,k) - E_N(i-1,j,k));
        else
          dEdx = RDX * 
                 (E_N(i,j,k) - E_N(i+1,j,k));
                 
        if(V(i,j-1,k) >= 0)                       
          dEdy = RDY *
                 (E_N(i,j,k) - E_N(i,j-1,k));
        else
          dEdy = RDY * 
                 (E_N(i,j,k) - E_N(i,j+1,k));
        
        if(W(i,j,k-1) >= 0)                     
          dEdz = RDZ *
                 (E_N(i,j,k) - E_N(i,j,k-1)); 
        else
          dEdz = RDZ *
                 (E_N(i,j,k) - E_N(i,j,k+1));  
                 
        nu_eff = nu_t + solver->nu / kE.sigma_E;
        nu_eff = max(nu_eff, solver->nu);


        Diffusion_E = pow(RDX,2) * ( (E_N(i+1,j,k) - E_N(i,j,k)) - 
                                     (E_N(i  ,j,k) - E_N(i-1,j,k)) );
        Diffusion_E +=pow(RDY,2) * ( (E_N(i,j+1,k) - E_N(i,j,k)) - 
                                     (E_N(i  ,j,k) - E_N(i,j-1,k)) );
        Diffusion_E +=pow(RDZ,2) * ( (E_N(i,j,k+1) - E_N(i,j,k)) - 
                                     (E_N(i  ,j,k) - E_N(i,j,k-1)) );
        Diffusion_E *= nu_eff * (1.0 / FV(i,j,k));
                      
        delE = ( -1.0 / FV(i,j,k)) *
                ( fabs((U(i,j,k) + U(i-1,j,k)) / 2) * dEdx +  
                  fabs((V(i,j,k) + V(i,j-1,k)) / 2) * dEdy + 
                  fabs((W(i,j,k) + W(i,j,k-1)) / 2) * dEdz ) + 
               ( E_N(i,j,k) / k_N(i,j,k) ) * 
                  ( kE.C1E * Production - kE.C2E * E_N(i,j,k)) +
               Diffusion_E;
        
               
        if(isnan(delk)) delk = 0;
        k(i,j,k) = max(k_N(i,j,k) + delk * solver->delt, 0);
        E_limit = kE.C_mu * pow(k(i,j,k), 1.5) / kE.length;
                
        if(isnan(delE) || (E_N(i,j,k) + delE * solver->delt) < 0.0000001) delE = 0;        
        E(i,j,k) = max(E_limit, E_N(i,j,k) + delE * solver->delt);
        nu_t(i,j,k) = max(kE.C_mu * pow(k(i,j,k),2) / E(i,j,k),0);

      }
    }
  }

  kE_boundaries(solver);

  solver_sendrecv_edge(solver, kE.k);
  solver_sendrecv_edge(solver, kE.E);
  solver_sendrecv_edge(solver, kE.nu_t);

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
