#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "intersections.h"
#include "mesh.h"
#include "readfile.h"
#include "stl.h"
#include "volfract.h"
#include "vtk.h"
#include "csv.h"
#include "markcells.h"

int main(int argc, char *argv[])
{
  struct mesh_data *mesh;
  struct stl_data *stl; 
  double limits[3];
  char filename[1024];
  int i;

  printf("mesh3d: fractional area volume mesh generator\n");

  if (argc<3) {
    printf("usage: mesh3d <source file> <stl file>\n");
    return(1);
  }

  #ifdef DEBUG
    printf("%s %s %s\n", argv[0], argv[1], argv[2]);
  #endif

  mesh = mesh_init_empty();


  if(read_mesh(mesh, argv[1])==1) return 1;

  if(mesh_init_complete(mesh) == 1) {
    printf("error: failed to initialize mesh in main()\n");
    return(1);
  }


  printf("Mesh successfully initialized\n");

  stl = stl_init_empty(); 

	limits[0] = mesh->origin[0] + mesh->delx * (mesh->imax + 1) + 0.0001;
	limits[1] = mesh->origin[1] + mesh->dely * (mesh->jmax + 1) + 0.0001;
	limits[2] = mesh->origin[2] + mesh->delz * (mesh->kmax + 1) + 0.0001;
	
	filename[0] = 0;
	for(i=2;i<argc;i++) {
		strncat(filename, argv[i], strlen(argv[i]));
		if(i+1<argc) strncat(filename, " ", 1);
	}
	printf("Reading geometry from stl file: %s\n",filename); 
	
  if(read_stl(stl, filename, limits)==1) return 1;  

  if(!stl_check(stl)) return 1;
  
  printf("Marking cells with no intersections\n");
  
  if(markcells_initialize(mesh,stl)==1) return 1;
  /* exit(0); */

  printf("Calculating area fractions\n");

  if(intersect_area_fractions(mesh, stl)==1) return 1; 

  printf("Calculating volume fractions\n");

  if(volume_fractions(mesh, stl)==1) return 1;

  printf("Filling mesh cells around obstacles\n");

  mesh_fill(mesh, stl);  

  printf("Eliminating very small mesh cells\n");

  mesh_normalize(mesh);
	/* */
	printf("Checking area / velocity ratios\n");
	
	mesh_avratio(mesh, 2.8); /* */
	mesh_avratio(mesh, 2.8); /* */
	mesh_avratio(mesh, 2.8); /* */
  mesh_normalize(mesh);

  printf("Writing mesh to file\n");

  if(vtk_write_fv(mesh, 0) == 1) return 1;
  if(csv_write_fv(mesh, 0) == 1) return 1;
  if(csv_write_af(mesh, 0) == 1) return 1;

  mesh_free(mesh);
  stl_free(stl);

  return(0);
}
