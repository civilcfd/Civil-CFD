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

  strncpy(filename, filename_csv, strlen(filename_csv));
  strncat(filename, ".gz", 3);
  fp = gzopen(filename, "w");

  if(fp == NULL) {
    printf("error: csv_compressd_write_scalar_grid cannot open %s to write\n", filename);
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

  strncpy(filename, filename_csv, strlen(filename_csv));
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
      return 0;
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