/*
* solver_mpi.h
*
* function definitions for mpi parallelism 
*/

int solver_mpi_range(struct solver_data *solver);
int solver_mpi(struct solver_data *solver, double timestep, double delt);
int solver_mpi_high_rank(struct solver_data *solver, double timestep);
int solver_recv_all(struct solver_data *solver);
int solver_sendrecv_edge(struct solver_data *solver, double *data);
int solver_sendrecv_edge_int(struct solver_data *solver, int *data);
int solver_send_all(struct solver_data *solver);
int solver_mpi_send(struct solver_data *solver, double *data, int to, long int i_start, long int i_range);
int solver_mpi_recv(struct solver_data *solver, double *data, int from, long int i_start, long int i_range);
int solver_mpi_sendrecv(struct solver_data *solver, int to, double *send, long int send_i_start, long int send_i_range, int from, double *recv, long int recv_i_start, long int recv_i_range);
int solver_mpi_sendrecv_int(struct solver_data *solver, int to, int *send, long int send_i_start, long int send_i_range, int from, int *recv, long int recv_i_start, long int recv_i_range);
int solver_broadcast_all(struct solver_data *solver);
double solver_mpi_max(struct solver_data *solver, double x);
double solver_mpi_min(struct solver_data *solver, double x); 
int solver_sum_edge(struct solver_data *solver, double *data);
int solver_mpi_init_comm(struct solver_data *solver);
int solver_mpi_sum(struct solver_data *solver, double *data, MPI_Comm comm, long int i_start, long int i_range);
int solver_sendrecv_delu(struct solver_data *solver);
int solver_mpi_sendrecv_replace(struct solver_data *solver, double *data, long int start, long int range, int to, int from);
int solver_mpi_init_complete(struct solver_data *solver);