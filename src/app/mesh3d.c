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

  printf("mesh3d: fractional area volume mesh generator\n");fflush(stdout);

  if (argc<3) {
    printf("usage: mesh3d <source file> <stl file>\n");
    return(1);
  }

  #ifdef DEBUG
    printf("%s %s %s\n", argv[0], argv[1], argv[2]);
  #endif

  mesh = mesh_init_empty();


  if(read_mesh_xml(mesh, argv[1])==1) return 1;

  if(mesh_init_complete(mesh) == 1) {
    printf("error: failed to initialize mesh in main()\n");
    return(1);
  }


  printf("Mesh successfully initialized\n");fflush(stdout);

  stl = stl_init_empty(); 

	limits[0] = mesh->origin[0] + mesh->delx * (mesh->imax + 1) + 0.0001;
	limits[1] = mesh->origin[1] + mesh->dely * (mesh->jmax + 1) + 0.0001;
	limits[2] = mesh->origin[2] + mesh->delz * (mesh->kmax + 1) + 0.0001;
	
	filename[0] = 0;
	for(i=2;i<argc;i++) {
		strncat(filename, argv[i], strlen(argv[i]));
		if(i+1<argc) strncat(filename, " ", 1);
	}
	printf("Reading geometry from stl file: %s\n",filename);fflush(stdout); 
	
  if(read_stl(stl, filename, limits)==1) return 1;  

  if(!stl_check(stl)) return 1;
  
  printf("\nMarking cells with no intersections\n");fflush(stdout);
  
  if(markcells_initialize(mesh,stl)==1) return 1;
  /* exit(0); */

  printf("\nCalculating area fractions\n");fflush(stdout);

  if(intersect_area_fractions(mesh, stl)==1) return 1; 

  printf("\nCalculating volume fractions\n");fflush(stdout);

  if(volume_fractions(mesh, stl)==1) return 1;

  printf("\nFilling mesh cells around obstacles\n");fflush(stdout);

  mesh_fill(mesh, stl);  

  printf("\nEliminating very small mesh cells\n");fflush(stdout);

  mesh_normalize(mesh);
	/* */
	printf("\nChecking area / velocity ratios\n");fflush(stdout);
	
	mesh_avratio(mesh, 4.0); /* */
	mesh_avratio(mesh, 4.0); /* */
	mesh_avratio(mesh, 4.0); /* */
  mesh_normalize(mesh);

  printf("\nWriting mesh to file\n");fflush(stdout);

  if(vtk_write_fv(mesh, 0) == 1) return 1;
  if(csv_write_fv(mesh, 0) == 1) return 1;
  if(csv_write_af(mesh, 0) == 1) return 1;

  mesh_free(mesh);
  stl_free(stl);

  return(0);
}
