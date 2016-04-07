/* csv.c
 *
 * routines to read/write csv files */

#include <stdio.h>
#include <math.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <io.h>
#endif

#include "csv.h"
#include "../solver3d/kE.h"

int csv_read_U(struct mesh_data *mesh, double timestep)
{
  char filename[256];


  sprintf(filename, "%4.3lf/U.csv", timestep);

  return csv_read_vector_grid(filename, 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->u, mesh->v, mesh->w);

}
int csv_write_U(struct mesh_data *mesh, double timestep)
{
  char filename[256];
  #ifndef _WIN32
  mode_t process_mask = umask(0);
  #endif
  
  sprintf(filename, "%4.3lf", timestep);

  #ifdef _WIN32  
  mkdir(filename);
  #else
  mkdir(filename, S_IRWXU | S_IRWXG | S_IRWXO);
  umask(process_mask);
  #endif
  
  sprintf(filename, "%4.3lf/U.csv", timestep);

  if(mesh->compress) return csv_compressed_write_vector_grid(filename, "u, v, w", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->u, mesh->v, mesh->w);
  else return csv_write_vector_grid(filename, "u, v, w", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->u, mesh->v, mesh->w);

}
int csv_write_vorticity(struct mesh_data *mesh, double timestep)
{
  char filename[256];
  #ifndef _WIN32
  mode_t process_mask = umask(0);
  #endif
  
  sprintf(filename, "%4.3lf", timestep);

  #ifdef _WIN32  
  mkdir(filename);
  #else
  mkdir(filename, S_IRWXU | S_IRWXG | S_IRWXO);
  umask(process_mask);
  #endif
  
  sprintf(filename, "%4.3lf/vorticity.csv", timestep);

  if(mesh->compress)  csv_compressed_write_vector_grid(filename, "u-vorticity, v-vorticity, w-vorticity", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->u_omega, mesh->v_omega, mesh->w_omega);
  else return csv_write_vector_grid(filename, "u-vorticity, v-vorticity, w-vorticity", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->u_omega, mesh->v_omega, mesh->w_omega);

}
int csv_read_P(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/P.csv", timestep);

  return csv_read_scalar_grid(filename, 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->P);

}
int csv_write_P(struct mesh_data *mesh, double timestep)
{
  char filename[256];
	#ifndef _WIN32
  mode_t process_mask = umask(0);
	#endif
	
  sprintf(filename, "%4.3lf", timestep);
  
	#ifdef _WIN32
	mkdir(filename);
	#else
  mkdir(filename, S_IRWXU | S_IRWXG | S_IRWXO);
  umask(process_mask);
  #endif
	
  sprintf(filename, "%4.3lf/P.csv", timestep);

  if(mesh->compress) return csv_compressed_write_scalar_grid(filename, "P", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->P);
  else return csv_write_scalar_grid(filename, "P", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->P);

}
int csv_read_vof(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/vof.csv", timestep);

  return csv_read_scalar_grid(filename, 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->vof);

}
int csv_read_n_vof(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/n_vof.csv", timestep);

  return csv_read_integer_grid(filename, 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->n_vof);

}
int csv_write_n_vof(struct mesh_data *mesh, double timestep)
{
  char filename[256];
  #ifndef _WIN32
	mode_t process_mask = umask(0);
	#endif
	
  sprintf(filename, "%4.3lf", timestep);
  
	#ifdef _WIN32
	mkdir(filename);
	#else
  mkdir(filename, S_IRWXU | S_IRWXG | S_IRWXO);
  umask(process_mask);
  #endif
	
  sprintf(filename, "%4.3lf/n_vof.csv", timestep);

  if(mesh->compress) return csv_compressed_write_integer_grid(filename, "n_vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->n_vof);
  else return csv_write_integer_grid(filename, "n_vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->n_vof);

}
int csv_write_vof(struct mesh_data *mesh, double timestep)
{
  char filename[256];
  #ifndef _WIN32
	mode_t process_mask = umask(0);
	#endif
	
  sprintf(filename, "%4.3lf", timestep);
  
	#ifdef _WIN32
	mkdir(filename);
	#else
  mkdir(filename, S_IRWXU | S_IRWXG | S_IRWXO);
  umask(process_mask);
  #endif
	
  sprintf(filename, "%4.3lf/vof.csv", timestep);

  if(mesh->compress) return csv_compressed_write_scalar_grid(filename, "vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->vof);
  else return csv_write_scalar_grid(filename, "vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->vof);

}
int csv_read_af(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/af.csv", timestep);

  return csv_read_vector_grid(filename, 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->ae, mesh->an, mesh->at);

}
int csv_write_af(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/af.csv", timestep);

  if(mesh->compress)  return csv_compressed_write_vector_grid(filename, "ae, an, at", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->ae, mesh->an, mesh->at);
  else return csv_write_vector_grid(filename, "ae, an, at", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->ae, mesh->an, mesh->at);

}
int csv_read_fv(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/fv.csv", timestep);

  return csv_read_scalar_grid(filename, 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->fv);

}
int csv_write_fv(struct mesh_data *mesh, double timestep)
{
  char filename[256];

  sprintf(filename, "%4.3lf/fv.csv", timestep);

  if(mesh->compress)  return csv_compressed_write_scalar_grid(filename, "fv", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->fv); 
  else return csv_write_scalar_grid(filename, "fv", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->fv);

}
int csv_write_k(struct mesh_data *mesh, double timestep)
{
  char filename[256];
  struct kE_data *turb;
  
  turb = mesh->turbulence_model;

  sprintf(filename, "%4.3lf/k.csv", timestep);

  if(mesh->compress)  return csv_compressed_write_scalar_grid(filename, "k", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        turb->k);
  else return csv_write_scalar_grid(filename, "k", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        turb->k);

}
int csv_write_E(struct mesh_data *mesh, double timestep)
{
  char filename[256];
  struct kE_data *turb;
  
  turb = mesh->turbulence_model;

  sprintf(filename, "%4.3lf/E.csv", timestep);

  if(mesh->compress) return csv_compressed_write_scalar_grid(filename, "E", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        turb->E); 
  else return csv_write_scalar_grid(filename, "E", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        turb->E);

}
long int csv_read_scalar_grid(char *filename,  
                          long int ni, long int nj, long int nk,
                          double *scalars) {
  FILE *fp;
  long int i, j, k, count, n;
  char text[1024];
  double f;

  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_read_scalar_grid\n");
    return -1;
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    n = csv_compressed_read_scalar_grid(filename, ni, nj, nk, scalars);
    if(n == -1) {
      printf("error: csv_read_scalar_grid cannot open %s to read\n", filename);
      return -1;
    }
    return n;
  }

  if(!fgets(text, sizeof(text), fp))
  { 
    if(!feof(fp)) {
      printf("error: fgets in csv_read_scalar_grid\n");
      return(-1);
    }
    return 0;
  }

  count=0;
  while(!feof(fp))
  {
    fscanf(fp, "%ld%*c %ld%*c %ld%*c %lf", &i, &j, &k, &f);
    
    if(i>ni || j>nj || k>nk) {
      printf("error: data out of bounds in csv_read_scalar_grid\n");
      break;
    }

    scalars[i + ni * (j + k * nj)] = f;
    count++;
  }

  fclose(fp);

  return count;
}
long int csv_read_integer_grid(char *filename,  
                          long int ni, long int nj, long int nk,
                          int *scalars) {
  FILE *fp;
  long int i, j, k, count, n;
  char text[1024];
  int f;

  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_read_scalar_grid\n");
    return -1;
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    n = csv_compressed_read_integer_grid(filename, ni, nj, nk, scalars);
    if(n == -1) {
      printf("error: csv_read_scalar_grid cannot open %s to read\n", filename);
      return -1;
    }
    return n;
  }

  if(!fgets(text, sizeof(text), fp))
  { 
    if(!feof(fp)) {
      printf("error: fgets in csv_read_scalar_grid\n");
      return(-1);
    }
    return 0;
  }

  count=0;
  while(!feof(fp))
  {
    fscanf(fp, "%ld%*c %ld%*c %ld%*c %d", &i, &j, &k, &f);
    
    if(i>ni || j>nj || k>nk) {
      printf("error: data out of bounds in csv_read_scalar_grid\n");
      break;
    }

    scalars[i + ni * (j + k * nj)] = f;
    count++;
  }

  fclose(fp);

  return count;
}

long int csv_read_vector_grid(char *filename,  
                          long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2) {
  FILE *fp;
  long int i, j, k, count, n;
  char text[1024];
  double f, g, h;

  if(filename == NULL || v0 == NULL || v1 == NULL || v2 == NULL) {
    printf("error: passed null arguments to csv_read_vector_grid\n");
    return -1;
  }

  fp = fopen(filename, "r");

  if(fp == NULL) {
    n = csv_compressed_read_vector_grid(filename, ni, nj, nk, v0, v1, v2);
    if(n == -1) {
      printf("error: csv_read_vector_grid cannot open %s to read\n", filename);
      return -1;
    }
    return n;
  }

  if(!fgets(text, sizeof(text), fp))
  { 
    if(!feof(fp)) {
      printf("error: fgets in csv_read_vector_grid\n");
      return(-1);
    }
    return 0;
  }

  count=0;
  while(!feof(fp))
  {
    fscanf(fp, "%ld%*c %ld%*c %ld%*c %lf%*c %lf%*c %lf%*c", &i, &j, &k, &f, &g, &h);
    
    if(i>ni || j>nj || k>nk) {
      printf("error: data out of bounds in csv_read_vector_grid\n");
      break;
    }

    v0[i + ni * (j + k * nj)] = f;
    v1[i + ni * (j + k * nj)] = g;
    v2[i + ni * (j + k * nj)] = h;
    count++;
  }

  fclose(fp);

  return count;
}

int csv_write_vector_grid(char *filename, char *dataset_name, long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2) {
  FILE *fp;
  long int i, j, k;

  const double emf = 0.000001;

  if(filename == NULL || v0 == NULL || v1 == NULL || v2 == NULL) {
    printf("error: passed null arguments to csv_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename);
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_write_vector_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp,"x, y, z, %s\n",dataset_name);

  for(i=0; i<ni; i++) {
    for(j=0; j<nj; j++) {
      for(k=0; k<nk; k++) {
        
        if(fabs(v0[i+ni * (j+ k*nj)]) > emf || fabs(v1[i+ni * (j+ k*nj)]) > emf ||  
           fabs(v2[i+ni * (j+ k*nj)]) > emf )
          fprintf(fp, "%ld, %ld, %ld, %lf, %lf, %lf\n", i, j, k, 
                v0[i + ni * (j + k * nj)], v1[i + ni * (j + k * nj)], 
                v2[i + ni * (j + k * nj)]); 
     
      }
    }
  }

  fclose(fp);

  return 0;
}

int csv_write_scalar_grid(char *filename, char *dataset_name, long int ni, long int nj, long int nk,
                          double *scalars) {
  FILE *fp;
  long int i, j, k;

  const double emf = 0.00000001;

  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename);
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp,"x, y, z, %s\n",dataset_name);

  for(i=0; i<ni; i++) {
    for(j=0; j<nj; j++) {
      for(k=0; k<nk; k++) {
        
        if(fabs(scalars[i+ni * (j+ k*nj)]) > emf) 
          fprintf(fp, "%ld, %ld, %ld, %10.8lf\n", i, j, k, 
                scalars[i + ni * (j + k * nj)]); 
      
      }
    }
  }

  fclose(fp);

  return 0;
}

int csv_write_integer_grid(char *filename, char *dataset_name, long int ni, long int nj, long int nk,
                          int *scalars) {
  FILE *fp;
  long int i, j, k;


  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename);
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp,"x, y, z, %s\n",dataset_name);

  for(i=0; i<ni; i++) {
    for(j=0; j<nj; j++) {
      for(k=0; k<nk; k++) {
        
          fprintf(fp, "%ld, %ld, %ld, %d\n", i, j, k, 
                scalars[i + ni * (j + k * nj)]); 
      
      }
    }
  }

  fclose(fp);

  return 0;
}

int csv_write_scalar_grid_paraview(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars) {
  FILE *fp;
  long int i, j, k;

  double ci, cj, ck;
  const double emf = 0.000001;

  ci = oi + di/2;
  cj = oj + dj/2;
  ck = ok + dk/2;

  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename);
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp,"x,y,z,%s\n",dataset_name);

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
        
        if(scalars[i+ni * (j+ k*nj)] > emf) 
          fprintf(fp, "%lf,%lf,%lf,%lf\n", di * i + ci,
                dj * j + cj, dk * k + ck,
                scalars[i + ni * (j + k * nj)]); 
      
      }
    }
  }

  fclose(fp);

  return 0;
}

void csv_remove(char *filename) {
  char filename_gz[1024];
  
  strncpy(filename_gz, filename, strlen(filename) + 1);
  strncat(filename_gz, ".gz", 3);
  remove(filename);
  remove(filename_gz);
}

