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
int boundary_weir(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence);   
double calc_flow(struct solver_data *solver, int x, long int imin, long int imax, 
                 long int jmin, long int jmax, long int kmin, long int kmax, double *area_ref);
int boundary_hgl(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, double turbulence); 
                            
enum special_boundaries vof_boundaries_check_inside_sb(struct solver_data *solver, long int a, long int b,
                                 int x);
#endif
