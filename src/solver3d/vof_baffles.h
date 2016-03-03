/* vof_baffles.h
*/

#ifndef VOF_BAFFLES_H
#define VOF_BAFFLES_H

int vof_baffles(struct solver_data *solver);

int baffle_k(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, long int pos);
                            
int baffle_slip(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double value, long int pos);

int baffle_flow(struct solver_data *solver, 
                            int x, double min_1, double min_2, double max_1, double max_2, 
                            double *value, long int pos);

int baffle_setup(struct solver_data *solver, int x, long int pos, 
                           long int *imin, long int *jmin, long int *kmin, 
                           long int *imax, long int *jmax, long int *kmax,
                           double min_1, double min_2, double max_1, double max_2);

                           

#endif
