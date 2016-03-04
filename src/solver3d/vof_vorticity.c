/* vof_vorticity.c
 *
 * calculate vorticity
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <omp.h>
 
#include "vtk.h"
#include "vof.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "track.h"

#include "vof_macros.h"


int vof_vorticity(struct solver_data *solver) {
  long int i,j,k;
  
#pragma omp parallel for shared (solver) private(i,j,k) \
            collapse(3) schedule(dynamic, 100)  
  for(i=0; i<IMAX; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
        
        if(i == 0 || i==IMAX-1 || j==0 || j==JMAX-1 || k==0 || k==KMAX-1) {
          U_OMEGA(i,j,k)=0;
          V_OMEGA(i,j,k)=0;
          W_OMEGA(i,j,k)=0;
          continue;
        }
        
        U_OMEGA(i,j,k) = (W(i,j,k+1) - W(i,j,k))/DELY - (V(i,j+1,k) - V(i,j,k))/DELZ;
        V_OMEGA(i,j,k) = (U(i+1,j,k) - U(i,j,k))/DELZ - (W(i,j,k+1) - W(i,j,k))/DELX;
        W_OMEGA(i,j,k) = (V(i,j+1,k) - V(i,j,k))/DELX - (U(i+1,j,k) - U(i,j,k))/DELY;
        
      }
    }
  }
  
  return 0;
}