/* vof_vorticity.c
 *
 * calculate vorticity
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <omp.h>
 
#include "vtk.h"
#include "vof_mpi.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "track.h"

#include "vof_macros.h"


int vof_vorticity(struct solver_data *solver) {
  long int i,j,k;
  
  if(solver->rank > 0) return 0;

  for(i=0; i<IMAX; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
        
        if(i == 0 || i==IMAX-1 || j==0 || j==JMAX-1 || k==0 || k==KMAX-1) {
          U_OMEGA(i,j,k)=0;
          V_OMEGA(i,j,k)=0;
          W_OMEGA(i,j,k)=0;
          continue;
        }
        /* wrong formula fixed 3/28/18 
        U_OMEGA(i,j,k) = (W(i,j,k+1) - W(i,j,k))/DELY - (V(i,j+1,k) - V(i,j,k))/DELZ;
        V_OMEGA(i,j,k) = (U(i+1,j,k) - U(i,j,k))/DELZ - (W(i,j,k+1) - W(i,j,k))/DELX;
        W_OMEGA(i,j,k) = (V(i,j+1,k) - V(i,j,k))/DELX - (U(i+1,j,k) - U(i,j,k))/DELY; */
        if(AT(i,j,k) > solver->emf && AT(i,j-1,k) > solver->emf && 
           AN(i,j,k) > solver->emf && AN(i,j,k-1) > solver->emf)
          U_OMEGA(i,j,k) = (W(i,j,k) - W(i,j-1,k))/DELY - (V(i,j,k) - V(i,j,k-1))/DELZ;
        else
          U_OMEGA(i,j,k) = 0;

        if(AE(i,j,k-1) > solver->emf && AE(i,j,k) > solver->emf && 
           AT(i-1,j,k) > solver->emf && AT(i,j,k) > solver->emf)
          V_OMEGA(i,j,k) = (U(i,j,k) - U(i,j,k-1))/DELZ - (W(i,j,k) - W(i-1,j,k))/DELX;
        else
          V_OMEGA(i,j,k) = 0;

        if(AN(i-1,j,k) > solver->emf && AN(i,j,k) > solver->emf && 
           AE(i,j-1,k) > solver->emf && AE(i,j,k) > solver->emf)
          W_OMEGA(i,j,k) = (V(i,j,k) - V(i-1,j,k))/DELX - (U(i,j,k) - U(i,j-1,k))/DELY;
        else
          W_OMEGA(i,j,k) = 0;
        
      }
    }
  }
  
  return 0;
}