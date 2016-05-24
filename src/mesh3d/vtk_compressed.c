/* vtk_compressed.c
 *
 * routines to write to compressed vtk files */

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <zlib.h>

#include "vtk.h"
#include "kE.h"

#define CELL_INDEX(i,j,k) (k + nk * (j + i * nj))

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
        d = double_swap(scalars[CELL_INDEX(i,j,k)]);
        gzwrite(&d, sizeof(double), 1, fp)
#else
        gzprintf(fp, "%4.6lf\n", scalars[CELL_INDEX(i,j,k)]); 
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
        n = int_swap(scalars[CELL_INDEX(i,j,k)]);
        gzwrite(&n, sizeof(uint32_t), 1, fp)
#else
        gzprintf(fp, "%d\n", scalars[CELL_INDEX(i,j,k)]); 
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
  gzprintf(fp, "ORIGIN %lf %lf %lf\n", oi, oj, ok);
  gzprintf(fp, "SPACING %lf %lf %lf\n", di, dj, dk);

  gzprintf(fp, "POINT_DATA %ld\n", (ni-1)*(nj-1)*(nk-1));
  gzprintf(fp, "VECTORS %s_data double\n", dataset_name);

  for(k=0; k<nk-1; k++) {
    for(j=0; j<nj-1; j++) {
      for(i=0; i<ni-1; i++) {
#ifdef VTK_BINARY
        d=double_swap(v1[CELL_INDEX(i,j,k)]);
        gzwrite(&d, sizeof(double), 1, fp)
        
        d=double_swap(v2[CELL_INDEX(i,j,k)]);
        gzwrite(&d, sizeof(double), 1, fp)
        
        d=double_swap(v3[CELL_INDEX(i,j,k)]);
        gzwrite(&d, sizeof(double), 1, fp)
#else
        if( (fabs(v1[CELL_INDEX(i,j,k)]) > emf && fabs(v1[CELL_INDEX(i+1,j,k)]) > emf) || 
            (fabs(v2[CELL_INDEX(i,j,k)]) > emf && fabs(v2[CELL_INDEX(i,j+1,k)]) > emf) ||  
            (fabs(v3[CELL_INDEX(i,j,k)]) > emf && fabs(v3[CELL_INDEX(i,j,k+1)]) > emf) ) 
          gzprintf(fp, "%4.6lf %4.6lf %4.6lf\n", 
            (v1[CELL_INDEX(i,j,k)] + v1[CELL_INDEX(i+1,j,k)])/2,
            (v2[CELL_INDEX(i,j,k)] + v2[CELL_INDEX(i,j+1,k)])/2,
            (v3[CELL_INDEX(i,j,k)] + v3[CELL_INDEX(i,j,k+1)])/2); 
        else  
          gzprintf(fp, "0.0 0.0 0.0\n");
          
#endif 
      }
    }
  }

  gzclose(fp);

  return 0;
}

int vtk_decompress(const char *cstr) {
  char buf[1024*1024*16];
  char filename[1024];
  char filename_gz[1024];
  int len, wlen;
  gzFile *fi;
  FILE *fp;
  
  strncpy(filename, cstr, strlen(cstr) + 1);
    
  strncpy(filename_gz, filename, strlen(filename) + 1);
  strncat(filename_gz, ".gz", 3);
  
  fi = gzopen(filename_gz,"r");
  if(fi == NULL) return 0;
  
  fp = fopen(filename,"w");
  if(fp == NULL) return 0;
  
  gzrewind(fi);
  while(!gzeof(fi))
  {
      len  = gzread(fi,buf,sizeof(buf));
      wlen = fwrite(buf, len, 1, fp);
      if(len != wlen) {
        printf("vtk_decompress: error writing to file\n");
        return 0;
      }
  }
  gzclose(fi);  
  fclose(fp);
  
  return 1;
}
