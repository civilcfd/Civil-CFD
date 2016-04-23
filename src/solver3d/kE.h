/*
 * kE.h 
 *
 * Describes a k-epsilon turbulence model */
 
#ifndef _KE_H
#define _KE_H

#include "solver_data.h"

struct kE_data {

  double *k;
  double *E;
  double *nu_t;
  
  /* wall shear stress */
  double *tau_x;
  double *tau_y;
  double *tau_z;

  /* store constants */
  double C1E;
  double C2E;
  double C3E;
  double C_mu;
  double sigma_k;
  double sigma_E;
  double vonKarman;
  double length_scale;

  double rough;
  double length;

};

double log_law(double u, double d, double mu, double rho, double rough);
int kE_tau(struct solver_data *solver);
int kE_set_value(char *param, int dims, 
                   double *vector);
int kE_init(struct solver_data *solver);
int kE_loop(struct solver_data *solver);
int kE_kill(struct solver_data *solver);
int kE_setup(struct solver_data *solver);
int kE_read(char *filename);
int kE_write(char *filename);
int kE_boundaries(struct solver_data *solver);
int kE_special_boundaries(struct solver_data *solver);
int kE_copy(struct solver_data *solver);
int kE_set_internal(struct solver_data *solver, double k, double E);
int kE_wall_shear(struct solver_data *solver);
double kE_nu(struct solver_data *solver, long int i, long int j, long int k);
int kE_loop_explicit(struct solver_data *solver);
int kE_boundary_wall_shear(struct solver_data *solver);
int kE_check(struct solver_data *solver);
int kE_boundary_fixed_k(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence);

#endif
