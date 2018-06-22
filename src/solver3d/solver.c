/*
 * solver.c
 *
 * creates and defines a CFD solver, including the mesh
 * does not perform the solution
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mesh.h"
#include "readfile.h"
#include "solver.h"
#include "readsolver.h"
#include "laminar.h"
#include "vof_mpi.h"
#include "kE.h"

struct solver_data *solver_init_empty() {

  struct solver_data *solver;

  solver = malloc(sizeof(struct solver_data));
  if (solver == NULL) {
    printf("error: could not allocate solver_data in solver_init_empty\n");
    return NULL;
  }

  /* initialize with default values */
  /* sensible defaults are used so that the user doesn't need
   * to specify standard values
   * based on water at 20 deg. C */

  solver->emf   = 0.00001;
  solver->emf_c = 0.99999;
  

  solver->min_vof = 0.0001;
  solver->max_vof = 0.9999;

  solver->vchgt = 0;
  solver->delt  = 0.001; /* 0.001 based on other solvers, this is generally
                          * a stable delta t for civil engineering
                          * problems */
  solver->delt_n = solver->delt;
  solver->delt_min = 0.000001;
  
  solver->endt = 200;    /* a good steady state value again for
                          * civil engineering problems */

  solver->writet = 1;

  solver->nu   = 1.004e-6;
  solver->rho  = 1000;

  solver->abstol = 50;
  solver->reltol = 1.0e-3;

  solver->con = 0.45;

  solver->gx   = 0;
  solver->gy   = 0;
  solver->gz   = -9.81;

  solver->umax = 0;
  solver->vmax = 0;
  solver->wmax = 0;

  solver->t    = 0.0;

  solver->iter = 0;
  solver->niter= 150;
  solver->vof_flag = 0;

  solver->loop = NULL;
  solver->boundaries = NULL;
  solver->special_boundaries = NULL;
  solver->pressure = NULL;
  solver->velocity = NULL;
  solver->convect   = NULL;

	time(&solver->start_time);

  strcpy(solver->ic[0].param,"end");

  solver->mesh = mesh_init_empty();

  if(solver->mesh == NULL) {
    printf("error: could not allocate mesh in solver_init_empty\n");
    free(solver);
    return NULL;
  }

  laminar_setup(solver);

  solver->ready = 0;

  return solver;
}

int solver_load(struct solver_data *solver, char *solverfile) {

 /* loads the solver and mesh into memory from file */

  if(read_mesh_xml(solver->mesh, solverfile) == 1) return 1;
  if(read_solver_xml(solver, solverfile) == 1) return 1;
  if(read_initial_xml(solver, solverfile) == 1) return 1;

  return 0;
}

int solver_run(struct solver_data *solver) {

  if(solver->ready == 0) {
    printf("error: uninitialized solver passed to solver_run\n");
    return(1);
  }

  if(solver_check(solver) == 1)
    return(1);

  return solver->loop(solver);
}

int solver_init_complete(struct solver_data *solver) {

  if(solver_check(solver) == 1) {
    return(1);
  }

  if(mesh_init_complete(solver->mesh) == 1) {
    return(1);
  }

  solver->ready = 1;

  return 0;

}

int solver_check(struct solver_data *solver) {

  if(solver == NULL) {
    printf("error: null solver passed to solver_init_complete\n");
    return(1);
  }

  if(solver->loop == NULL || solver->boundaries == NULL || 
     solver->pressure == NULL || solver->velocity == NULL) {
    printf("error: solver function pointers are null\n");
    return(1);
  }

  return 0; 
}

int solver_store_initial(struct solver_data *solver, char *param, int dims, 
                   double *vector) {
  static int index=0;

  if(solver == NULL || param == NULL) {
    printf("error: null values passed to solver_store_initial\n");
    return(1);
  }
  if(index == 16) {
    printf("error: too many initial conditions in solver_store_initial\n");
    return(1);
  }

  if(strcmp(param, "velocity")==0) {
    if(dims != 3) {
      printf("error in source file: velocity requires 3 arguments\n");
      return(1);
    }

    strcpy(solver->ic[index].param, "velocity");
    solver->ic[index].value[0] = vector[0];
    solver->ic[index].value[1] = vector[1];
    solver->ic[index].value[2] = vector[2];
    solver->ic[index].dims = dims;
    index++;        
    
  } else if(strcmp(param, "inside")==0) {
    if(dims != 3) {
      printf("error in source file: inside requires 3 arguments\n");
      return(1);
    }

    strcpy(solver->ic[index].param, "inside");
    solver->ic[index].value[0] = vector[0];
    solver->ic[index].value[1] = vector[1];
    solver->ic[index].value[2] = vector[2];
    solver->ic[index].dims = dims;
    index++;        
    
  } else if(strcmp(param, "vof_height")==0) {
    if(dims < 1) {
      printf("error in source file: vof_height requires 1 argument\n");
      return(1);
    }

    solver->ic[index].value[0] = vector[0];
    solver->ic[index].dims = 2;
    strcpy(solver->ic[index].param, "vof_height");
        
    if(dims >= 2) /* second optional dimension is delay */ {
      solver->ic[index].value[1] = vector[1];      
    }
    else {
      solver->ic[index].value[1] = 0.0;      
    }
    
    index++;   
    
  } else if(strcmp(param, "hydrostatic")==0) {
    if(dims < 1) {
      printf("error in source file: hydrostatic requires 1 argument\n");
      return(1);
    }
    
    solver->ic[index].value[0] = vector[0];
    solver->ic[index].dims = dims;
    strcpy(solver->ic[index].param, "hydrostatic");
    index++;   
    
  } else if(strcmp(param, "kE_k")==0) {
    if(dims < 1) {
      printf("error in source file: kE_k requires 1 argument\n");
      return(1);
    }

    solver->ic[index].value[0] = vector[0];
    solver->ic[index].dims = dims;
    strcpy(solver->ic[index].param, "kE_k");
    index++;   
  } else if(strcmp(param, "end")==0) {

    solver->ic[index].value[0] = 0;
    solver->ic[index].dims = dims;
    strcpy(solver->ic[index].param, "end");
    index = 0;   
  }

  return 0;
}                   

double solver_get_initial_scalar(struct solver_data *solver, char *param) {
  int i;
  
  for(i=0; i<16; i++) {
  
    if(strcmp(solver->ic[i].param,"end")==0) break;
    
    if(strcmp(solver->ic[i].param,param)==0) return solver->ic[i].value[0];
    
  }  
  
  return -1;
}

int solver_get_initial_vector(struct solver_data *solver, char *param, double *vector) {
  static int i = 0;

  for(; i<16; i++) {
  
    if(strcmp(solver->ic[i].param,"end")==0) break;
    
    if(strcmp(solver->ic[i].param,param)==0) {
      vector[0] = solver->ic[i].value[0];
      vector[1] = solver->ic[i].value[1];
      vector[2] = solver->ic[i].value[2];
      i++;
      return 0;
    }
  }  
  
  i=0;
  return 1;
}

int solver_initial_values(struct solver_data *solver) {
  int i;
  double vector[3];
  
  for(i=0; i<16; i++) {
  
    if(strcmp(solver->ic[i].param,"end")==0) break;
    
    vector[0] = solver->ic[i].value[0];
    vector[1] = solver->ic[i].value[1];
    vector[2] = solver->ic[i].value[2];
    solver_set_initial(solver, solver->ic[i].param, solver->ic[i].dims, vector);
  
  }
  
  return 0;
}
                   
int solver_set_initial(struct solver_data *solver, char *param, int dims, 
                   double *vector) {
  long int z_height;

  if(solver == NULL || param == NULL || vector == NULL) {
    printf("error: null values passed to solver_set_initial\n");
    return(1);
  }

  if(strcmp(param, "inside")==0) {
    if(dims < 3) {
      printf("error in source file: inside requires 3 arguments\n");
      return(1);
    }  
  
    mesh_fill_vof(solver->mesh, vector);
  
  } else if(strcmp(param, "velocity")==0) {
    if(dims < 3) {
      printf("error in source file: velocity requires 3 arguments\n");
      return(1);
    }
 
    mesh_set_array(solver->mesh, "u", vector[0], -1, 0, 0, 0, 0, 0);   
    mesh_set_array(solver->mesh, "v", vector[1], -1, 0, 0, 0, 0, 0);
    mesh_set_array(solver->mesh, "w", vector[2], -1, 0, 0, 0, 0, 0);          
  
  } else if(strcmp(param, "u")==0) {
    if(dims != 1) {
      printf("error in source file: u requires 1 argument\n");
      return(1);
    }

    mesh_set_array(solver->mesh, "u", vector[0], -1, 0, 0, 0, 0, 0);
  } else if(strcmp(param, "v")==0) {
    if(dims != 1) {
      printf("error in source file: v requires 1 argument\n");
      return(1);
    }

    mesh_set_array(solver->mesh, "v", vector[0], -1, 0, 0, 0, 0, 0);  
  } else if(strcmp(param, "w")==0) {
    if(dims != 1) {
      printf("error in source file: w requires 1 argument\n");
      return(1);
    }

    mesh_set_array(solver->mesh, "w", vector[0], -1, 0, 0, 0, 0, 0);  
  } else if(strcmp(param, "vof_height")==0) {
    if(dims < 1) {
      printf("error in source file: vof_height requires 1 argument\n");
      return(1);
    }
    
    solver->vof_height = vector[0];
    solver->vof_delay = vector[1];
    
    if(fmod(solver->vof_height,solver->mesh->delz) < solver->emf) solver->vof_height -= solver->min_vof * 100;
    
    if(solver->vof_delay < solver->emf) {
    
      z_height = floor(solver->vof_height / solver->mesh->delz);

      mesh_set_array(solver->mesh, "vof", 1.0, 0, solver->mesh->imax, 
                                               0, solver->mesh->jmax, 
                                               0, z_height);  
                                               
      if(solver->vof_height / solver->mesh->delz - z_height > 0) {
        mesh_set_array(solver->mesh, "vof", solver->vof_height / solver->mesh->delz - z_height, 
                                               0, solver->mesh->imax, 
                                               0, solver->mesh->jmax, 
                                               z_height, z_height+1);        
      }
    
    }
                                             
  } else if(strcmp(param, "hydrostatic")==0) {
    if(dims < 1) {
      printf("error in source file: hydrostatic requires 1 argument\n");
      return(1);
    }
    
    if(vector[0]>0) mesh_set_hydrostatic(solver->mesh, fabs(solver->gz), solver->rho);  
    
  } else if(strcmp(param, "kE_k")==0) {
    if(dims < 1) {
      printf("error in source file: kE_k requires 1 argument\n");
      return(1);
    }
    if(!kE_check(solver)) {
      printf("solver_set_initial: kE_k ignored as turbulence model is not k-Epsilon\n");
      return(1);
    }
    
    kE_set_internal(solver, vector[0], 0);
  }

  return 0;
}                   
                   
int solver_set_value(struct solver_data *solver, char *param, int dims, 
                   double *vector) {


  if(solver == NULL || param == NULL || vector == NULL) {
    printf("error: null values passed to solver_set_value\n");
    return(1);
  }
  else if(strcmp(param, "abstol")==0) {
    if(dims != 1) {
      printf("error in source file: abstol requires 1 argument\n");
      return(1);
    }

    solver->abstol = vector[0];
  }
  else if(strcmp(param, "reltol")==0) {
    if(dims != 1) {
      printf("error in source file: reltol requires 1 argument\n");
      return(1);
    }

    solver->reltol = vector[0];
  }
  else if(strcmp(param, "gravity")==0) {
    if(dims != 3) {
      printf("error in source file: gravity requires 3 arguments\n");
      return(1);
    }

    solver->gx = vector[0];
    solver->gy = vector[1];
    solver->gz = vector[2];
  }
  else if (strcmp(param, "nu")==0) {
    if(dims != 1) {
      printf("error in source file: nu requires 1 arguments\n");
      return(1);
    }

    solver->nu = vector[0];

  }
  else if (strcmp(param, "rho")==0) {
    if(dims != 1) {
      printf("error in source file: rho requires 1 arguments\n");
      return(1);
    }

    solver->rho = vector[0];

  }
  else if (strcmp(param, "turbulence")==0) {
    if(dims != 1) {
      printf("error in source file: turbulence requires 1 arguments\n");
      return(1);
    }

    if(vector[0] == 1)  kE_setup(solver);
    else laminar_setup(solver);   

  }
  else if (strcmp(param, "t")==0) {
    if(dims != 1) {
      printf("error in source file: t requires 1 arguments\n");
      return(1);
    }

    solver->t = vector[0];
  }
  else if (strcmp(param, "endt")==0) {
    if(dims != 1) {
      printf("error in source file: endt requires 1 arguments\n");
      return(1);
    }

    solver->endt = vector[0];
  }
  else if (strcmp(param, "writet")==0) {
    if(dims != 1) {
      printf("error in source file: writet requires 1 arguments\n");
      return(1);
    }

    solver->writet = vector[0];
  }
  else if (strcmp(param, "delt")==0) {
    if(dims != 1) {
      printf("error in source file: delt requires 1 arguments\n");
      return(1);
    }

    solver->delt = vector[0];
  }
  else if (strcmp(param, "autot")==0) {
    if(dims != 1) {
      printf("error in source file: autot requires 1 arguments\n");
      return(1);
    }

    if(vector[0] == 0) solver->deltcal=NULL;
  }
  else if(strncmp(param, "end", 3)==0) {
    return(0);
  }
  else {
    printf("error in source file: unrecognized command: %s\n", param);
    return(1);
  }

  return(0);
}
