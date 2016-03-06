/* vof_boundary.h
*/

#ifndef VOF_BOUNDARY_H
#define VOF_BOUNDARY_H

int sboundary_setup(struct solver_data *solver, int x, long int *imin, long int *jmin, long int *kmin, 
                           long int *imax, long int *jmax, long int *kmax,
                           double min_1, double min_2, double max_1, double max_2);
                           
int boundary_fixed_velocity(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence);

int boundary_mass_outflow(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence);                                                       
/*
int boundary_mass_inflow(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence);   */

int boundary_hgl(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence); 
                            
int vof_vof_height_boundary(struct solver_data *solver);

#endif
