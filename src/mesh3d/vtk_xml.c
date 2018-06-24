/* vtk_xml_xml.c
 *
 * routines to write to vtk files using the xml format */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <string.h>
#include <zlib.h>

#include "vtk_xml.h"
#include "kE.h"

#define CELL_INDEX(i,j,k) ((k) + nk * ((j) + (i) * nj))

unsigned char * base64_encode(const unsigned char *src, size_t len,
			      size_t *out_len);

int vtk_xml_write_scalar_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *scalars) {
  FILE *fp;
  long int i, j, k, n;
  size_t len, out_len;
  char *enc;
  uint8_t arr[24];
  long int count;
  double dd[3];


  int ret, flush;
  unsigned have;
  z_stream strm;
  double *scalars_reorder;
  unsigned char *in;
  unsigned char *out;

  uint32_t header[4];
  unsigned char header_char[16];


  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  ret = deflateInit(&strm, Z_DEFAULT_COMPRESSION);
  if (ret != Z_OK) {
    printf("error: could not initialize zlib\n");
    return 1;
  }


  if(filename == NULL || scalars == NULL) {
    printf("error: passed null arguments to vtk_xml_write_scalar_grid\n");
    return 1;
  }

  vtk_xml_remove(filename);
  
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_xml_write_scalar_grid cannot open %s to write\n", filename);
    return 1;
  }

  scalars_reorder = malloc(sizeof(double) * ni * nj * nk);
  in = malloc(sizeof(unsigned char) * ni * nj * nk * 8);
  out = malloc(sizeof(unsigned char) * ni * nj * nk * 8);

  if(scalars_reorder == NULL || in == NULL || out == NULL) {
    printf("error: could not malloc in vtk_xml_write_scalar_grid\n");
    return 1;
  }

  fprintf(fp, "<?xml version=\"1.0\"?>\n");

  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\" > \n");
  
  fprintf(fp, "<ImageData WholeExtent=\"%ld %ld %ld %ld %ld %ld\" Origin=\"%lf %lf %lf\" Spacing=\"%lf %lf %lf\">\n",
              0, ni-1, 0, nj-1, 0, nk-1, oi, oj, ok, di, dj, dk);
  fprintf(fp,"<Piece Extent=\"%ld %ld %ld %ld %ld %ld\">\n",
              0, ni-1, 0, nj-1, 0, nk-1);
  
  fprintf(fp, "<PointData Scalars=\"%s\">\n",dataset_name);
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"0\"  />\n",dataset_name);
  
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");


  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</ImageData>\n");

  fprintf(fp, "<AppendedData encoding=\"base64\">\n_");

  n = 0;
  for(k=0; k<nk; k++) {
    for(j=0; j<nj; j++) {
      for(i=0; i<ni; i++) {
        scalars_reorder[n] = scalars[CELL_INDEX(i,j,k)];
        n++;
      }
    }
  }

  memcpy(in, scalars_reorder, ni * nj * nk * 8);
  strm.avail_in = ni * nj * nk * 8;
  strm.next_in  = in;

  have = 0;

  strm.avail_out = ni * nj * nk * 8;
  strm.next_out = out;
  ret = deflate(&strm, Z_FINISH);    
  if (ret == Z_STREAM_ERROR) {
    printf("error: could not deflate scalars\n");
    free(scalars_reorder);
    free(in);
    free(out);
    return 1;
  }
  have = ni * nj * nk * 8 - strm.avail_out;

  if(ret != Z_STREAM_END)    
  {
    printf("Could not deflate all data.  Deflated %d from %ld\n",have,ni*nj*nk*8);
  }
  else {
  #ifdef DEBUG
    printf("Deflated %d from %ld\n",have,ni*nj*nk*8);
  #endif
  }
  deflateEnd(&strm);

  /* first write header */
  header[0] = 1; // number of blocks
  header[1] = ni * nj * nk * 8;
  header[2] = ni * nj * nk * 8;
  header[3] = have;
  memcpy(header_char, header, 16);
  enc = base64_encode(header_char, 16, &out_len);
  fwrite(enc, 1, out_len, fp);
  free(enc);
  

  /* now write data */
  enc = base64_encode(out, have, &out_len);
  fwrite(enc, 1, out_len, fp);
  free(enc);

  fprintf(fp, "</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);

  free(scalars_reorder);
  free(in);
  free(out);

  return 0;
}


int vtk_xml_write_vector_grid(char *filename, char *dataset_name, 
                          long int ni, long int nj, long int nk,
                          double oi, double oj, double ok,
                          double di, double dj, double dk,
                          double *v1, double *v2, double *v3) {
  FILE *fp;
  long int i, j, k, n;
  const double emf = 0.000001, d;
  size_t len, out_len;
  char *enc;
  uint8_t arr[24];
  long int count;
  
  int ret, flush;
  unsigned have;
  z_stream strm;
  double *v_reorder;
  unsigned char *in;
  unsigned char *out;

  uint32_t header[4];
  unsigned char header_char[16];


  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  ret = deflateInit(&strm, Z_DEFAULT_COMPRESSION);
  if (ret != Z_OK) {
    printf("error: could not initialize zlib\n");
    return 1;
  }

  if(filename == NULL || v1 == NULL || v2 == NULL || v3 == NULL) {
    printf("error: passed null arguments to vtk_xml_write_vector_grid\n");
    return 1;
  }

  vtk_xml_remove(filename);
  
  fp = fopen(filename, "w");

  if(fp == NULL) {
    printf("error: vtk_xml_write_vector_grid cannot open %s to write\n", filename);
    return 1;
  }

  v_reorder = malloc(sizeof(double) * (ni-1) * (nj-1) * (nk-1) * 3);
  in = malloc(sizeof(unsigned char) * (ni-1) * (nj-1) * (nk-1) * 3 * 8);
  out = malloc(sizeof(unsigned char) * (ni-1) * (nj-1) * (nk-1) * 3 * 8);

  if(v_reorder == NULL || in == NULL || out == NULL) {
    printf("error: could not malloc in vtk_xml_write_scalar_grid\n");
    return 1;
  }

  fprintf(fp, "<?xml version=\"1.0\"?>\n");

  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\" compressor=\"vtkZLibDataCompressor\" > \n");
  
  fprintf(fp, "<ImageData WholeExtent=\"%ld %ld %ld %ld %ld %ld\" Origin=\"%lf %lf %lf\" Spacing=\"%lf %lf %lf\">\n",
              0, ni-2, 0, nj-2, 0, nk-2, oi, oj, ok, di, dj, dk);
  fprintf(fp,"<Piece Extent=\"%ld %ld %ld %ld %ld %ld\">\n",
              0, ni-2, 0, nj-2, 0, nk-2);
  
  fprintf(fp, "<PointData Vectors=\"%s\">\n",dataset_name);
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"  />\n",dataset_name);
  
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");


  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</ImageData>\n");

  fprintf(fp, "<AppendedData encoding=\"base64\">\n_");

  n = 0;
  for(k=0; k<nk-1; k++) {
    for(j=0; j<nj-1; j++) {
      for(i=0; i<ni-1; i++) {

        if( (fabs(v1[CELL_INDEX(i,j,k)]) > emf && fabs(v1[CELL_INDEX(i+1,j,k)]) > emf) || 
            (fabs(v2[CELL_INDEX(i,j,k)]) > emf && fabs(v2[CELL_INDEX(i,j+1,k)]) > emf) ||  
            (fabs(v3[CELL_INDEX(i,j,k)]) > emf && fabs(v3[CELL_INDEX(i,j,k+1)]) > emf) ) {

          v_reorder[n] = (v1[CELL_INDEX(i,j,k)] + v1[CELL_INDEX(i+1,j,k)])/2; n++;
          v_reorder[n] = (v2[CELL_INDEX(i,j,k)] + v2[CELL_INDEX(i,j+1,k)])/2; n++;
          v_reorder[n] = (v3[CELL_INDEX(i,j,k)] + v3[CELL_INDEX(i,j,k+1)])/2; n++;
        }
        else  {
          v_reorder[n] = 0; n++;
          v_reorder[n] = 0; n++;
          v_reorder[n] = 0; n++;
        }
          

      }
    }
  }

  memcpy(in, v_reorder, (ni-1) * (nj-1) * (nk-1) * 8 * 3);
  strm.avail_in = (ni-1) * (nj-1) * (nk-1) * 8 * 3;
  strm.next_in  = in;

  have = 0;

  strm.avail_out = (ni-1) * (nj-1) * (nk-1) * 8 * 3;
  strm.next_out = out;
  ret = deflate(&strm, Z_FINISH);    
  if (ret == Z_STREAM_ERROR) {
    printf("error: could not deflate vectors\n");
    free(v_reorder);
    free(in);
    free(out);
    return 1;
  }
  have = (ni-1) * (nj-1) * (nk-1) * 8 * 3 - strm.avail_out;

  if(ret != Z_STREAM_END)    
  {
    printf("Could not deflate all data.  Deflated %d from %ld\n",have,(ni-1) * (nj-1) * (nk-1)*8*3);
  }
  else {
  #ifdef DEBUG
    printf("Deflated %d from %ld\n",have,(ni-1) * (nj-1) * (nk-1)*8*3);
  #endif
  }
  deflateEnd(&strm);

  /* first write header */
  header[0] = 1; // number of blocks
  header[1] = (ni-1) * (nj-1) * (nk-1) * 8 * 3;
  header[2] = (ni-1) * (nj-1) * (nk-1) * 8 * 3;
  header[3] = have;
  memcpy(header_char, header, 16);
  enc = base64_encode(header_char, 16, &out_len);
  fwrite(enc, 1, out_len, fp);
  free(enc);
  

  /* now write data */
  enc = base64_encode(out, have, &out_len);
  fwrite(enc, 1, out_len, fp);
  free(enc);

  fprintf(fp, "</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);

  free(v_reorder);
  free(in);
  free(out);

  return 0;
}

void vtk_xml_remove(char *filename) {
  char filename_gz[1024];
  
  strncpy(filename_gz, filename, strlen(filename) + 1);
  strncat(filename_gz, ".gz", 3);
  remove(filename);
  remove(filename_gz);
}

/*
 * Base64 encoding/decoding (RFC1341)
 * Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
 *
 * This software may be distributed under the terms of the BSD license.
 * See README for more details.
 */

static const unsigned char base64_table[65] =
	"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/**
 * base64_encode - Base64 encode
 * @src: Data to be encoded
 * @len: Length of the data to be encoded
 * @out_len: Pointer to output length variable, or %NULL if not used
 * Returns: Allocated buffer of out_len bytes of encoded data,
 * or %NULL on failure
 *
 * Caller is responsible for freeing the returned buffer. Returned buffer is
 * nul terminated to make it easier to use as a C string. The nul terminator is
 * not included in out_len.
 */
unsigned char * base64_encode(const unsigned char *src, size_t len,
			      size_t *out_len)
{
	unsigned char *out, *pos;
	const unsigned char *end, *in;
	size_t olen;
	int line_len;

	olen = len * 4 / 3 + 4; /* 3-byte blocks to 4-byte */
	olen += olen / 72; /* line feeds */
	olen++; /* nul termination */
	if (olen < len)
		return NULL; /* integer overflow */
	out = malloc(olen);
	if (out == NULL)
		return NULL;

	end = src + len;
	in = src;
	pos = out;
	line_len = 0;
	while (end - in >= 3) {
		*pos++ = base64_table[in[0] >> 2];
		*pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
		*pos++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
		*pos++ = base64_table[in[2] & 0x3f];
		in += 3;
		line_len += 4;
		if (line_len >= 72) {
			//*pos++ = '\n';
			line_len = 0;
		}
	}

	if (end - in) {
		*pos++ = base64_table[in[0] >> 2];
		if (end - in == 1) {
			*pos++ = base64_table[(in[0] & 0x03) << 4];
			*pos++ = '=';
		} else {
			*pos++ = base64_table[((in[0] & 0x03) << 4) |
					      (in[1] >> 4)];
			*pos++ = base64_table[(in[1] & 0x0f) << 2];
		}
		*pos++ = '=';
		line_len += 4;
	}

	//if (line_len)
		//*pos++ = '\n';

	*pos = '\0';
	if (out_len)
		*out_len = pos - out;
	return out;
}