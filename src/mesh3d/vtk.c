/* vtk.c
 *
 * routines to write to vtk files */

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "vtk.h"
#include "kE.h"

int vtk_write_fv(struct mesh_data *mesh, int timestep)
{
  char filename[256];

  sprintf(filename, "vtk/fv_%d.vtk", timestep);

  // not compressing fv because it's opened so often in the gui
  /* if(mesh->compress) return vtk_compressed_write_scalar_grid(filename, "fv", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->fv);
  else */ 
  return vtk_write_scalar_grid(filename, "fv", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->fv);

}

int vtk_write_n_vof(struct mesh_data *mesh, int timestep)
{
  char filename[256];

  sprintf(filename, "vtk/n_vof_%d.vtk", timestep);

  if(mesh->compress) return vtk_compressed_write_integer_grid(filename, "n_vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->n_vof);
  else return vtk_write_integer_grid(filename, "n_vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->n_vof);

}

int vtk_write_U(struct mesh_data *mesh, int timestep)
{
  char filename[256];
  
  sprintf(filename, "vtk/U_%d.vtk", timestep);

  if(mesh->compress) return vtk_compressed_write_vector_grid(filename, "U", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, 
                        mesh->u, mesh->v, mesh->w);
  else return vtk_write_vector_grid(filename, "U", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, 
                        mesh->u, mesh->v, mesh->w);

}

int vtk_write_P(struct mesh_data *mesh, int timestep)
{
  char filename[256];

  sprintf(filename, "vtk/P_%d.vtk", timestep);

  if(mesh->compress) return vtk_compressed_write_scalar_grid(filename, "P", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->P);
  else return vtk_write_scalar_grid(filename, "P", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->P);

}

int vtk_write_vorticity(struct mesh_data *mesh, int timestep)
{
  char filename[256];
  
  /* only writes z-vorticity because why would you want to visualize any other vorticity */
  /* the csv has all vorticities written */
 
 /*
  sprintf(filename, "vtk/z_vorticity_%d.vtk", timestep);

  return vtk_write_scalar_magnitude_grid(filename, "z_vorticity", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->w_omega); */

 /* changed my mind */
 sprintf(filename, "vtk/vorticity_%d.vtk", timestep);
 
 if(mesh->compress)  return vtk_compressed_write_vector_grid(filename, "vorticity", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, 
                        mesh->u_omega, mesh->v_omega, mesh->w_omega);
 else return vtk_write_vector_grid(filename, "vorticity", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, 
                        mesh->u_omega, mesh->v_omega, mesh->w_omega);

}

int vtk_write_vof(struct mesh_data *mesh, int timestep)
{
  char filename[256];

  sprintf(filename, "vtk/vof_%d.vtk", timestep);

  if(mesh->compress) return vtk_compressed_write_scalar_grid(filename, "vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->vof);
  else return vtk_write_scalar_grid(filename, "vof", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, mesh->vof);

}

int vtk_write_k(struct mesh_data *mesh, int timestep)
{
  char filename[256];
  struct kE_data *turb;
  
  turb = mesh->turbulence_model;

  sprintf(filename, "vtk/k_%d.vtk", timestep);

  if(mesh->compress) return vtk_compressed_write_scalar_grid(filename, "k", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, turb->k); 
  else return vtk_write_scalar_grid(filename, "k", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, turb->k);

}

int vtk_write_E(struct mesh_data *mesh, int timestep)
{
  char filename[256];
  struct kE_data *turb;
  
  turb = mesh->turbulence_model;

  sprintf(filename, "vtk/E_%d.vtk", timestep);

  if(mesh->compress)  return vtk_compressed_write_scalar_grid(filename, "E", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, turb->E);
  else return vtk_write_scalar_grid(filename, "E", 
                        mesh->imax, mesh->jmax, mesh->kmax,
                        mesh->origin[0], mesh->origin[1], mesh->origin[2],
                        mesh->delx, mesh->dely, mesh->delz, turb->E);

}

int vtk_write_scalar_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars) {
  FILE *fp;
  long int i, j, k;
  double d;

  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to vtk_write_scalar_grid\n");
    return 1;
  }

  vtk_remove(filename);
  
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  fprintf(fp, "BINARY\n");
#else
  fprintf(fp, "ASCII\n");
#endif

  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni, nj, nk);
  fprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  fprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  fprintf(fp, "POINT_DATA %ld\n", ni*nj*nk);
  fprintf(fp, "SCALARS %s_data double 1\n", dataset_name);
  fprintf(fp, "LOOKUP_TABLE default\n");

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
#ifdef VTK_BINARY
        d = double_swap(scalars[i + ni * (j + k * nj)]);
        fwrite(&d, sizeof(double), 1, fp)
#else
        fprintf(fp, "%4.6lf\n", scalars[i + ni * (j + k * nj)]); 
#endif 
      }
    }
  }

  fclose(fp);

  return 0;
}

int vtk_write_scalar_magnitude_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars) {
  FILE *fp;
  long int i, j, k;
  double d;

  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to vtk_write_scalar_grid\n");
    return 1;
  }

  vtk_remove(filename);
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  fprintf(fp, "BINARY\n");
#else
  fprintf(fp, "ASCII\n");
#endif

  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni, nj, nk);
  fprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  fprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  fprintf(fp, "POINT_DATA %ld\n", ni*nj*nk);
  fprintf(fp, "SCALARS %s_data double 1\n", dataset_name);
  fprintf(fp, "LOOKUP_TABLE default\n");

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
#ifdef VTK_BINARY
        d = double_swap(scalars[i + ni * (j + k * nj)]);
        fwrite(&d, sizeof(double), 1, fp)
#else
        fprintf(fp, "%4.6lf\n", fabs(scalars[i + ni * (j + k * nj)])); 
#endif 
      }
    }
  }

  fclose(fp);

  return 0;
}


int vtk_write_integer_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          int *scalars) {
  FILE *fp;
  long int i, j, k;
  uint32_t n;
  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to vtk_write_integer_grid\n");
    return 1;
  }

  vtk_remove(filename);  
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  fprintf(fp, "BINARY\n");
#else
  fprintf(fp, "ASCII\n");
#endif

  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni, nj, nk);
  fprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  fprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  fprintf(fp, "POINT_DATA %ld\n", ni*nj*nk);
  fprintf(fp, "SCALARS %s_data integer 1\n", dataset_name);
  fprintf(fp, "LOOKUP_TABLE default\n");

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
#ifdef VTK_BINARY
        n = int_swap(scalars[i + ni * (j + k * nj)]);
        fwrite(&n, sizeof(uint32_t), 1, fp)
#else
        fprintf(fp, "%d\n", scalars[i + ni * (j + k * nj)]); 
#endif 
      }
    }
  }

  fclose(fp);

  return 0;
}


int vtk_write_vector_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *v1, double *v2, double *v3) {
  FILE *fp;
  long int i, j, k;
  const double emf = 0.000001, d;
  
  if(filename == NULL || v1 == NULL || v2 == NULL || v3 == NULL) {
    printf("error: passed null arguments to vtk_write_vector_grid\n");
    return 1;
  }

  vtk_remove(filename);
  
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_write_vector_grid cannot open %s to write\n", filename);
    return 1;
  }

  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  fprintf(fp, "BINARY\n");
#else
  fprintf(fp, "ASCII\n");
#endif

  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  fprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni-1, nj-1, nk-1);
  fprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  fprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  fprintf(fp, "POINT_DATA %ld\n", (ni-1)*(nj-1)*(nk-1));
  fprintf(fp, "VECTORS %s_data double\n", dataset_name);

  for(k=0; k<nk-1; k++) {
    for(j=0; j<nj-1; j++) {
      for(i=0; i<ni-1; i++) {
#ifdef VTK_BINARY
        d=double_swap(v1[i + ni * (j + k * nj)]);
        fwrite(&d, sizeof(double), 1, fp)
        
        d=double_swap(v2[i + ni * (j + k * nj)]);
        fwrite(&d, sizeof(double), 1, fp)
        
        d=double_swap(v3[i + ni * (j + k * nj)]);
        fwrite(&d, sizeof(double), 1, fp)
#else
        if( (fabs(v1[i+ni * (j+ k*nj)]) > emf && fabs(v1[i+1 +ni * (j+ k*nj)]) > emf) || 
            (fabs(v2[i+ni * (j+ k*nj)]) > emf && fabs(v2[i+ni * ((j+1)+ k*nj)]) > emf) ||  
            (fabs(v3[i+ni * (j+ k*nj)]) > emf && fabs(v3[i+ni * (j+ (k+1)*nj)]) > emf) ) 
          fprintf(fp, "%4.6lf %4.6lf %4.6lf\n", 
            (v1[i + ni * (j + k * nj)] + v1[i+1 + ni * (j + k * nj)])/2,
            (v2[i + ni * (j + k * nj)] + v2[i + ni * ((j+1) + k * nj)])/2,
            (v3[i + ni * (j + k * nj)] + v3[i + ni * (j + (k+1) * nj)])/2); 
        else  
          fprintf(fp, "0.0 0.0 0.0\n");
          
#endif 
      }
    }
  }

  fclose(fp);

  return 0;
}
uint32_t int_swap(uint32_t n) {
   union
   {
      uint32_t n;
      uint8_t b[4];
   } dat1, dat2;

   dat1.n = n;
   dat2.b[0] = dat1.b[3];
   dat2.b[1] = dat1.b[2];
   dat2.b[2] = dat1.b[1];
   dat2.b[3] = dat1.b[0];
   return dat2.n;
}
double double_swap(double d) {
   union
   {
      double d;
      uint8_t b[8];
   } dat1, dat2;

   dat1.d = d;
   dat2.b[0] = dat1.b[7];
   dat2.b[1] = dat1.b[6];
   dat2.b[2] = dat1.b[5];
   dat2.b[3] = dat1.b[4];
   dat2.b[4] = dat1.b[3];
   dat2.b[5] = dat1.b[2];
   dat2.b[6] = dat1.b[1];
   dat2.b[7] = dat1.b[0];
   return dat2.d;
}

void vtk_remove(char *filename) {
  char filename_gz[1024];
  
  strncpy(filename_gz, filename, strlen(filename) + 1);
  strncat(filename_gz, ".gz", 3);
  remove(filename);
  remove(filename_gz);
}
