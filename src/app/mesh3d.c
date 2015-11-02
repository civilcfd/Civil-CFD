#include <stdio.h>
#include <stdlib.h>

#include "intersections.h"
#include "mesh.h"
#include "readfile.h"
#include "stl.h"
#include "volfract.h"
#include "vtk.h"
#include "csv.h"

int main(int argc, char *argv[])
{
  struct mesh_data *mesh;
  struct stl_data *stl; 

  printf("mesh3d: fractional area volume mesh generator\n");

  if (argc!=3) {
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

  if(read_stl(stl, argv[2])==1) return 1;  

  if(!stl_check(stl)) return 1;

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
	
	mesh_avratio(mesh, 4.0); /* */
	mesh_avratio(mesh, 4.0); /* */
	mesh_avratio(mesh, 4.0); /* */
  mesh_normalize(mesh);

  printf("Writing mesh to file\n");

  if(vtk_write_fv(mesh, 0) == 1) return 1;
  if(csv_write_fv(mesh, 0) == 1) return 1;
  if(csv_write_af(mesh, 0) == 1) return 1;

  mesh_free(mesh);
  stl_free(stl);

  return(0);
}
