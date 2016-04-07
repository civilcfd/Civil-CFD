/* vtk_compressed.c
 *
 * routines to write to compressed vtk files */

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <zlib.h>

#include "vtk.h"
#include "kE.h"

int vtk_compressed_write_scalar_grid(char *filename_vtk, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars) {
  gzFile *fp;
  long int i, j, k;
  double d;
  char filename[1024];

  if(filename_vtk == NULL || scalars == NULL) {
    printf("error: passed null arguments to vtk_write_scalar_grid\n");
    return 1;
  }

  strncpy(filename, filename_vtk, strlen(filename_vtk) + 1);
  strncat(filename, ".gz", 3);
  
  vtk_remove(filename_vtk);
  
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_compressed_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  gzprintf(fp, "# vtk DataFile Version 2.0\n");
  gzprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  gzprintf(fp, "BINARY\n");
#else
  gzprintf(fp, "ASCII\n");
#endif

  gzprintf(fp, "DATASET STRUCTURED_POINTS\n");
  gzprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni, nj, nk);
  gzprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  gzprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  gzprintf(fp, "POINT_DATA %ld\n", ni*nj*nk);
  gzprintf(fp, "SCALARS %s_data double 1\n", dataset_name);
  gzprintf(fp, "LOOKUP_TABLE default\n");

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
#ifdef VTK_BINARY
        d = double_swap(scalars[i + ni * (j + k * nj)]);
        gzwrite(&d, sizeof(double), 1, fp)
#else
        gzprintf(fp, "%4.6lf\n", scalars[i + ni * (j + k * nj)]); 
#endif 
      }
    }
  }

  gzclose(fp);

  return 0;
}

int vtk_compressed_write_integer_grid(char *filename_vtk, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          int *scalars) {
  gzFile *fp;
  long int i, j, k;
  uint32_t n;
  char filename[1024];
  
  if(filename_vtk == NULL || scalars == NULL) {
    printf("error: passed null arguments to vtk_write_integer_grid\n");
    return 1;
  }

  strncpy(filename, filename_vtk, strlen(filename_vtk) + 1);
  strncat(filename, ".gz", 3);
  
  vtk_remove(filename_vtk);
  
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_compressed_write_integer_grid cannot open %s to write\n", filename);
    return 1;
  }

  gzprintf(fp, "# vtk DataFile Version 2.0\n");
  gzprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  gzprintf(fp, "BINARY\n");
#else
  gzprintf(fp, "ASCII\n");
#endif

  gzprintf(fp, "DATASET STRUCTURED_POINTS\n");
  gzprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni, nj, nk);
  gzprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  gzprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  gzprintf(fp, "POINT_DATA %ld\n", ni*nj*nk);
  gzprintf(fp, "SCALARS %s_data integer 1\n", dataset_name);
  gzprintf(fp, "LOOKUP_TABLE default\n");

  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
#ifdef VTK_BINARY
        n = int_swap(scalars[i + ni * (j + k * nj)]);
        gzwrite(&n, sizeof(uint32_t), 1, fp)
#else
        gzprintf(fp, "%d\n", scalars[i + ni * (j + k * nj)]); 
#endif 
      }
    }
  }

  gzclose(fp);

  return 0;
}

int vtk_compressed_write_vector_grid(char *filename_vtk, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *v1, double *v2, double *v3) {
  gzFile *fp;
  long int i, j, k;
  const double emf = 0.000001, d;
  char filename[1024];
  
  if(filename_vtk == NULL || v1 == NULL || v2 == NULL || v3 == NULL) {
    printf("error: passed null arguments to vtk_compressed_write_vector_grid\n");
    return 1;
  }

  strncpy(filename, filename_vtk, strlen(filename_vtk) + 1);
  strncat(filename, ".gz", 3);
  
  vtk_remove(filename_vtk);
  
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_compressed_write_vector_grid cannot open %s to write\n", filename);
    return 1;
  }

  gzprintf(fp, "# vtk DataFile Version 2.0\n");
  gzprintf(fp, "%s\n", dataset_name);
#ifdef VTK_BINARY
  gzprintf(fp, "BINARY\n");
#else
  gzprintf(fp, "ASCII\n");
#endif

  gzprintf(fp, "DATASET STRUCTURED_POINTS\n");
  gzprintf(fp, "DIMENSIONS %ld %ld %ld\n", ni-1, nj-1, nk-1);
  gzprintf(fp, "ORIGIN %lf %lf %lf\n", oi+di/2, oj+dj/2, ok+dk/2);
  gzprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  gzprintf(fp, "POINT_DATA %ld\n", (ni-1)*(nj-1)*(nk-1));
  gzprintf(fp, "VECTORS %s_data double\n", dataset_name);

  for(k=0; k<nk-1; k++) {
    for(j=0; j<nj-1; j++) {
      for(i=0; i<ni-1; i++) {
#ifdef VTK_BINARY
        d=double_swap(v1[i + ni * (j + k * nj)]);
        gzwrite(&d, sizeof(double), 1, fp)
        
        d=double_swap(v2[i + ni * (j + k * nj)]);
        gzwrite(&d, sizeof(double), 1, fp)
        
        d=double_swap(v3[i + ni * (j + k * nj)]);
        gzwrite(&d, sizeof(double), 1, fp)
#else
        if( (fabs(v1[i+ni * (j+ k*nj)]) > emf && fabs(v1[i+1 +ni * (j+ k*nj)]) > emf) || 
            (fabs(v2[i+ni * (j+ k*nj)]) > emf && fabs(v2[i+ni * ((j+1)+ k*nj)]) > emf) ||  
            (fabs(v3[i+ni * (j+ k*nj)]) > emf && fabs(v3[i+ni * (j+ (k+1)*nj)]) > emf) ) 
          gzprintf(fp, "%4.6lf %4.6lf %4.6lf\n", 
            (v1[i + ni * (j + k * nj)] + v1[i+1 + ni * (j + k * nj)])/2,
            (v2[i + ni * (j + k * nj)] + v2[i + ni * ((j+1) + k * nj)])/2,
            (v3[i + ni * (j + k * nj)] + v3[i + ni * (j + (k+1) * nj)])/2); 
        else  
          gzprintf(fp, "0.0 0.0 0.0\n");
          
#endif 
      }
    }
  }

  gzclose(fp);

  return 0;
}
