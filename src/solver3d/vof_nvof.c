/*
 * vof_nvof.c: Calculate NVOF and interpolation
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <omp.h>
#include <mpi.h>
#include <petscsys.h>
 
#include "vtk.h"
#include "vof_mpi.h"
#include "vof_boundary.h"
#include "solver.h"
#include "mesh.h"
#include "csv.h"
#include "kE.h"
#include "track.h"
#include "vof_baffles.h"
#include "mesh_mpi.h"
#include "solver_mpi.h"

#include "vof_macros.h"

#define emf solver->emf
#define emf_c solver->emf_c

double surface_interpolate(struct solver_data *solver, long int i, long int j, long int k) {
  double interpolate=0;
  long int l,m,n;

  l = i; m = j; n = k;

  switch (N_VOF(i,j,k))
  {
  case west:
    l=i-1;
    break;

  case east:
    l=i+1;
    break;

  case south:
    m=j-1;
    break;

  case north:
    m=j+1;  
    break;

  case bottom:
    n=k-1;
    break;

  case top:
    n=k+1;
    break;

  case none:
  case empty:
    return 1;

  }

  if(FV(l,m,n) > emf) {
    interpolate = 1 / ((VOF(i,j,k) + VOF(l,m,n) * 0.5));    
  }
  else {
    interpolate = 1 / (0.5 + VOF(i,j,k));
  }

  if ((FV(l,m,n) >= emf) && (N_VOF(l,m,n)!=0.0)) interpolate=1.0;

  interpolate = min(interpolate , 2);
  interpolate = max(interpolate , 0);
            
  return interpolate;
}

int vof_mpi_nvof(struct solver_data *solver) {
  long int i,j,k,l,m,n,ii,jj,kk;
  int mobs, inf, infcr, iobs;
  double vf, fxm, fxp, fym, fyp, fzm, fzp; 
  double vfxm, vfxp, vfym, vfyp, vfzm, vfzp;


  const double nemf = -1.0 * emf; 
 
  for(i=1; i<IRANGE-1; i++) {
    for(j=0; j<JMAX; j++) {
      for(k=0; k<KMAX; k++) {
        N_VOF(i,j,k) = none;
        if(j==0 || k==0 || j==JMAX-1 || k==KMAX-1 || FV(i,j,k) == 0)
          N_VOF(i,j,k) = 0;
        else if(VOF(i+1,j,k) >= emf && VOF(i,j+1,k) >= emf 
          && VOF(i-1,j,k) >= emf && VOF(i,j-1,k) >= emf
          && VOF(i,j,k+1) >= emf && VOF(i,j,k-1) >= emf) {
          N_VOF(i,j,k) = 0;          
        }
      }
    }
  }

  if(ISTART ==0) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        N_VOF(0,j,k) = 0;
      }
    }
  }
  if(IRANGE+ISTART == IMAX) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {
        N_VOF(IRANGE-1,j,k) = 0;
      }
    }
  }
 
  for(i=1; i<IRANGE-1; i++) {
    for(j=1; j<JMAX-1; j++) {
      for(k=1; k<KMAX-1; k++) {

        if(FV(i,j,k) == 0) {
          N_VOF(i,j,k) = 0;
          continue;
        }

        if(VOF(i,j,k) < emf) {
          N_VOF(i,j,k) = 8;
          continue;
        }
  
        if(VOF(i,j,k) > emf_c) {
          N_VOF(i,j,k) = 0;
          continue;
        }

        if(N_VOF(i,j,k) == 0) continue;

        /*# now calculate the partial derivatives of F
        # this code calculates how much F changes in each axis
        # the goal is to determine where the free surface lies
        # and the slope of the free surface */
               
               fxm=VOF(i,j,k);
               fxp=VOF(i,j,k);
               fym=VOF(i,j,k);
               fyp=VOF(i,j,k);
               fzm=VOF(i,j,k);
               fzp=VOF(i,j,k);

               if (AE(i-1,j,k)>emf) fxm=VOF(i-1,j,k);
               if (AE(i,j,k)>emf)   fxp=VOF(i+1,j,k);

               if (AN(i,j-1,k)>emf) fym=VOF(i,j-1,k);
               if (AN(i,j,k)>emf)   fyp=VOF(i,j+1,k);

               if (AT(i,j,k-1)>emf) fzm=VOF(i,j,k-1);
               if (AT(i,j,k)>emf)   fzp=VOF(i,j,k+1);

               mobs=1;
               inf=1;
               iobs=1;
               vf=0.0;

               vfxm=0.0;
               vfxp=0.0;

               for (kk=1;kk<4;kk++)
                {
                  n=k-2+kk;
                  for (jj=1;jj<4;jj++)
                    {
                       m=j-2+jj;
                       vfxm=vfxm+VOF(i-1,m,n);
                       vfxp=vfxp+VOF(i+1,m,n);
                    }
                }

               if (FV(i-1,j,k)==0.0 && solver->gx > nemf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fxm>=emf))
                 {
                    inf=inf+1;
                    if (vfxm>vf) N_VOF(i,j,k)=west;
                    if (N_VOF(i,j,k)==west) vf=vfxm;
                 }

               iobs=1;
               if (FV(i+1,j,k)==0.0 && solver->gx < emf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fxp>=emf))
                 {
                    inf=inf+1;
                    if (vfxp>vf) N_VOF(i,j,k)=east;
                    if (N_VOF(i,j,k)==east) vf=vfxp;
                 }
               iobs=1;
                                /* z-axis **/
               vfzm=0.0;
               vfzp=0.0;
               for (ii=1;ii<4;ii++)
                 {
                   l=i-2+ii;
                   for (jj=1;jj<4;jj++)
                     {
                       m=j-2+jj;
                       vfzm=vfzm+VOF(l,m,k-1);
                       vfzp=vfzp+VOF(l,m,k+1);
                     }
                 }
/*
       if(i==241 && j==48 && k==42) {
          printf("break\n");
        }  */
        
               if (FV(i,j,k-1)==0.0 && solver->gz > nemf) iobs=2; 
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fzm>=emf))
                 {
                   inf=inf+1;
                   if (vfzm>vf) N_VOF(i,j,k)=bottom;
                   if (N_VOF(i,j,k)==bottom) vf=vfzm;
                 }
               iobs=1;

               if (FV(i,j,k+1)==0.0 && solver->gz < emf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fzp>=emf))
                 {
                   inf=inf+1;
                   if (vfzp>vf) N_VOF(i,j,k)=top;
                   if (N_VOF(i,j,k)==top) vf=vfzp;
                 }
               iobs=1;
                                /*y-axis**/

               vfym=0.0;
               vfyp=0.0;

               for (kk=1;kk<4;kk++)
                 {
                   n=k-2+kk;
                   for (ii=1;ii<4;ii++)
                     {
                       l=i-2+ii;
                       vfym=vfym+VOF(l,j-1,n); 
                       vfyp=vfyp+VOF(l,j+1,n);
                     }
                 }

               if (FV(i,j-1,k)==0.0  && solver->gy > nemf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fym>=emf))
                 {
                   inf=inf+1;
                   if (vfym>vf) N_VOF(i,j,k)=south;
                   if (N_VOF(i,j,k)==south) vf=vfym;
                 }
               iobs=1;

               if (FV(i,j+1,k)==0.0  && solver->gy < emf) iobs=2;
               mobs=mobs+(iobs-1);
               if ((iobs!=2)&&(fyp>=emf))
                 {
                   inf=inf+1;
                   if (vfyp>vf) N_VOF(i,j,k)=north;
                   if (N_VOF(i,j,k)==north) vf=vfyp;
                 }
               iobs=1;

              
               /* check if it is not a free surface, but is bounded by an obstacle */
               /* essentially we have fluid in each direction that isn't an obstacle */
               infcr=8-mobs;
               if ((inf==infcr)&&(infcr>1)) N_VOF(i,j,k)=0; 
               

      }
    }
  }


  return 0;

#undef emf
#undef emf_c
}
