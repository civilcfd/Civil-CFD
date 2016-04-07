/* csv_compressed.c
 *
 * routines to read/write compressed csv files */

#include <stdio.h>
#include <math.h>
#include <sys/stat.h>
#include <string.h>
#include <zlib.h>

#ifdef _WIN32
#include <io.h>
#endif

#include "csv.h"
#include "../solver3d/kE.h"

int csv_compressed_write_scalar_grid(char *filename_csv, char *dataset_name, long int ni, long int nj, long int nk,
                          double *scalars) {
  gzFile *fp;
  long int i, j, k;
  char filename[1024];

  const double emf = 0.00000001;
  

  if(filename_csv == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_compressed_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename_csv);
  strncpy(filename, filename_csv, strlen(filename_csv) + 1);
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_compressed_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  gzprintf(fp,"x, y, z, %s\n",dataset_name);

  for(i=0; i<ni; i++) {
    for(j=0; j<nj; j++) {
      for(k=0; k<nk; k++) {
        
        if(fabs(scalars[i+ni * (j+ k*nj)]) > emf) 
          gzprintf(fp, "%ld, %ld, %ld, %10.8lf\n", i, j, k, 
                scalars[i + ni * (j + k * nj)]); 
      
      }
    }
  }

  gzclose(fp);

  return 0;
}

int csv_compressed_write_vector_grid(char *filename_csv, char *dataset_name, long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2) {
  gzFile *fp;
  long int i, j, k;
  char filename[1024];

  const double emf = 0.000001;

  if(filename_csv == NULL || v0 == NULL || v1 == NULL || v2 == NULL) {
    printf("error: passed null arguments to csv_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename_csv);
  strncpy(filename, filename_csv, strlen(filename_csv) + 1);
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_compressed_write_vector_grid cannot open %s to write\n", filename);
    return 1;
  }

  gzprintf(fp,"x, y, z, %s\n",dataset_name);

  for(i=0; i<ni; i++) {
    for(j=0; j<nj; j++) {
      for(k=0; k<nk; k++) {
        
        if(fabs(v0[i+ni * (j+ k*nj)]) > emf || fabs(v1[i+ni * (j+ k*nj)]) > emf ||  
           fabs(v2[i+ni * (j+ k*nj)]) > emf )
          gzprintf(fp, "%ld, %ld, %ld, %lf, %lf, %lf\n", i, j, k, 
                v0[i + ni * (j + k * nj)], v1[i + ni * (j + k * nj)], 
                v2[i + ni * (j + k * nj)]); 
     
      }
    }
  }

  gzclose(fp);

  return 0;
}

int csv_compressed_write_integer_grid(char *filename_csv, char *dataset_name, long int ni, long int nj, long int nk,
                          int *scalars) {
  gzFile *fp;
  long int i, j, k;
  char filename[1024];


  if(filename_csv == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_write_scalar_grid\n");
    return 1;
  }

  csv_remove(filename_csv);
  strncpy(filename, filename_csv, strlen(filename_csv) + 1);
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_compressed_write_integer_grid cannot open %s to write\n", filename);
    return 1;
  }

  gzprintf(fp,"x, y, z, %s\n",dataset_name);

  for(i=0; i<ni; i++) {
    for(j=0; j<nj; j++) {
      for(k=0; k<nk; k++) {
        
          gzprintf(fp, "%ld, %ld, %ld, %d\n", i, j, k, 
                scalars[i + ni * (j + k * nj)]); 
      
      }
    }
  }

  gzclose(fp);

  return 0;
}

long int csv_compressed_read_scalar_grid(char *filename_csv,  
                          long int ni, long int nj, long int nk,
                          double *scalars) {
  gzFile *fp;
  long int i, j, k, count;
  char text[1024];
  char filename[1024];
  double f;

  if(filename_csv == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_read_scalar_grid\n");
    return -1;
  }

  strncpy(filename, filename_csv, strlen(filename_csv) + 1);
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "r");

  if(fp == NULL) {
    printf("error: csv_compressed_read_scalar_grid cannot open %s to read\n", filename);
    return -1;
  }

  if(!gzgets(fp, text, sizeof(text)))
  { 
    if(!gzeof(fp)) {
      printf("error: gzgets in csv_compressed_read_scalar_grid\n");
      return(-1);
    }
    return 0;
  }

  count=0;
  while(!gzeof(fp))
  {
    if(!gzgets(fp, text, sizeof(text)))
    { 
      if(!gzeof(fp)) {
        printf("error: gzgets in csv_compressed_read_scalar_grid\n");
        return(-1);
      }
      break;
    }
    sscanf(text, "%ld%*c %ld%*c %ld%*c %lf", &i, &j, &k, &f);
    
    if(i>ni || j>nj || k>nk) {
      printf("error: data out of bounds in csv_compressed_read_scalar_grid\n");
      break;
    }

    scalars[i + ni * (j + k * nj)] = f;
    count++;
  }

  gzclose(fp);

  return count;
}

long int csv_compressed_read_integer_grid(char *filename_csv,  
                          long int ni, long int nj, long int nk,
                          int *scalars) {
  FILE *fp;
  long int i, j, k, count;
  char text[1024];
  char filename[1024];
  int f;

  if(filename_csv == NULL || scalars == NULL) {
    printf("error: passed null arguments to csv_read_scalar_grid\n");
    return -1;
  }

  strncpy(filename, filename_csv, strlen(filename_csv) + 1);
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "r");

  if(fp == NULL) {
    printf("error: csv_compressed_read_scalar_grid cannot open %s to read\n", filename);
    return -1;
  }

  if(!gzgets(fp, text, sizeof(text)))
  { 
    if(!gzeof(fp)) {
      printf("error: fgets in csv_compressed_read_scalar_grid\n");
      return(-1);
    }
    return 0;
  }

  count=0;
  while(!feof(fp))
  {
    if(!gzgets(fp, text, sizeof(text)))
    { 
      if(!gzeof(fp)) {
        printf("error: gzgets in csv_compressed_read_integer_grid\n");
        return(-1);
      }
      break;
    }
    sscanf(text, "%ld%*c %ld%*c %ld%*c %d", &i, &j, &k, &f);
    
    if(i>ni || j>nj || k>nk) {
      printf("error: data out of bounds in csv_compressed_read_integer_grid\n");
      break;
    }

    scalars[i + ni * (j + k * nj)] = f;
    count++;
  }

  gzclose(fp);

  return count;
}

long int csv_compressed_read_vector_grid(char *filename_csv,  
                          long int ni, long int nj, long int nk,
                          double *v0, double *v1, double *v2) {
  gzFile *fp;
  long int i, j, k, count;
  char text[1024];
  char filename[1024];
  double f, g, h;

  if(filename_csv == NULL || v0 == NULL || v1 == NULL || v2 == NULL) {
    printf("error: passed null arguments to csv_read_scalar_grid\n");
    return -1;
  }

  strncpy(filename, filename_csv, strlen(filename_csv) + 1);
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "r");

  if(fp == NULL) {
    printf("error: csv_compressed_read_vector_grid cannot open %s to read\n", filename);
    return -1;
  }

  if(!gzgets(fp, text, sizeof(text)))
  { 
    if(!gzeof(fp)) {
      printf("error: gzgets in csv_compressed_read_vector_grid\n");
      return(-1);
    }
    return 0;
  }

  count=0;
  while(!gzeof(fp))
  {
    if(!gzgets(fp, text, sizeof(text)))
    { 
      if(!gzeof(fp)) {
        printf("error: gzgets in csv_compressed_read_vector_grid\n");
        return(-1);
      }
      break;
    }
    sscanf(text, "%ld%*c %ld%*c %ld%*c %lf%*c %lf%*c %lf%*c", &i, &j, &k, &f, &g, &h);
    
    if(i>ni || j>nj || k>nk) {
      printf("error: data out of bounds in csv_compressed_read_vector_grid\n");
      break;
    }

    v0[i + ni * (j + k * nj)] = f;
    v1[i + ni * (j + k * nj)] = g;
    v2[i + ni * (j + k * nj)] = h;
    count++;
  }

  gzclose(fp);

  return count;
}
