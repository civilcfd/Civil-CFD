#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "solver_data.h"
#include "solver.h"
#include "vof.h"
#include "intersections.h"
#include "mesh.h"
#include "readfile.h"
#include "volfract.h"
#include "vtk.h"
#include "csv.h"
#include "kE.h"
#include "track.h"
#include "vof_macros.h"

#ifdef _WIN32
#include <io.h>
char * const nulFileName = "NUL";
#define CROSS_DUP(fd) _dup(fd)
#define CROSS_DUP2(fd, newfd) _dup2(fd, newfd)
#define STDOUT_FILENO 1
#else
#include <unistd.h>
char * const nulFileName = "/dev/null";
#define CROSS_DUP(fd) dup(fd)
#define CROSS_DUP2(fd, newfd) dup2(fd, newfd)
#endif

int main(int argc, char *argv[])
{
  struct solver_data *solver;
  double timestep;
  long int i,j,k;

  int stdoutBackupFd;
  FILE *nullOut;

  printf("inspect_cell: diagnostic tool that outputs relevant data for a cell\n");

  if (argc<4) {
    printf("usage: inspect_cell <timestep> <i> <j> <k>\n");
    return(1);
  } 

  #ifdef DEBUG
    printf("%s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4]);
  #endif

  /* duplicate stdout and supress stdout */
  stdoutBackupFd = CROSS_DUP(STDOUT_FILENO);
  fflush(stdout);
  nullOut = fopen(nulFileName, "w");
#ifndef _WIN32
  CROSS_DUP2(fileno(nullOut), STDOUT_FILENO);
#else
  CROSS_DUP2(_fileno(nullOut), STDOUT_FILENO);
#endif

  timestep = atof(argv[1]);
  i = atol(argv[2]);
  j = atol(argv[3]);
  k = atol(argv[4]);

  solver = solver_init_empty();
  if(solver == NULL) return 1;

  vof_setup_solver(solver);
  
  if(solver_load(solver, "solverfile", "meshfile", "initials")==1) {
    fflush(stdout);
    fclose(nullOut);  
    CROSS_DUP2(stdoutBackupFd, STDOUT_FILENO);
#ifndef _WIN32
    close(stdoutBackupFd);
#else
	_close(stdoutBackupFd);
#endif
    printf("Error reading solver and mesh data.\n");
  
    return 1;
  }

  if(solver_init_complete(solver)==1) {
    fflush(stdout);
    fclose(nullOut);  
	CROSS_DUP2(stdoutBackupFd, STDOUT_FILENO);
#ifndef _WIN32
    close(stdoutBackupFd);
#else
	_close(stdoutBackupFd);
#endif
    printf("Error reading solver and mesh data.\n");
  
    return 1;
  }

  if(solver->turbulence_read != NULL) solver->turbulence_read("turbulencefile");
  
  solver->init(solver);
  solver->turbulence_init(solver);
  
  if(mesh_load_csv(solver->mesh, 0) == 1) {
    fflush(stdout);
    fclose(nullOut);  
    CROSS_DUP2(stdoutBackupFd, STDOUT_FILENO);
#ifndef _WIN32
	close(stdoutBackupFd);
#else
	_close(stdoutBackupFd);
#endif
    printf("Error reading solver and mesh data.\n");
  
    return 1;
  }
  solver_initial_values(solver);
  solver->petacal(solver);

  csv_read_U(solver->mesh,timestep);
  csv_read_P(solver->mesh,timestep);
  csv_read_vof(solver->mesh,timestep);
  csv_read_n_vof(solver->mesh,timestep);
  solver->t = timestep;

  /* Restore stdout */
  fflush(stdout);
  fclose(nullOut);  
  CROSS_DUP2(stdoutBackupFd, STDOUT_FILENO);
#ifndef _WIN32
  close(stdoutBackupFd);
#else
  _close(stdoutBackupFd);
#endif

  if(i > solver->mesh->imax-1 || j > solver->mesh->jmax-1 || k > solver->mesh->kmax-1 ||
     i < 0 || j < 0 || k < 0) {
    printf("Cannot read cell - out of range\n");
    return 0;
  }
  
  printf("\ntimestep %lf     cell: %ld %ld %ld\n",timestep, i, j, k);
  printf("position: %lf %lf %lf\n\ncell data:\n\n", i * DELX + solver->mesh->origin[0], 
                                                    j * DELY + solver->mesh->origin[1],
                                                    k * DELZ + solver->mesh->origin[2]);
  
  printf("FV: %lf\n",FV(i,j,k));
  printf("P: %lf      VOF: %lf     N_VOF: %d\n\nedge data:\n\n",P(i,j,k),VOF(i,j,k),N_VOF(i,j,k));
  printf("A(e/n/t) %ld %ld %ld: %lf %lf %lf\n",i,j,k,AE(i,j,k),AN(i,j,k),AT(i,j,k));
  printf("A(w/s/b) %ld %ld %ld: %lf %lf %lf\n\n",i,j,k,AE(i-1,j,k),AN(i,j-1,k),AT(i,j,k-1));  
  printf("U(e/n/t) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i,j,k),V(i,j,k),W(i,j,k));
  printf("U(w/s/b) %ld %ld %ld: %lf %lf %lf\n",i,j,k,U(i-1,j,k),V(i,j-1,k),W(i,j,k-1));  

  return(0);
}
