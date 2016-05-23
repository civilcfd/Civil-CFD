/*
 * solver_mpi.h
 *
 * function definitions for mpi parallelism 
 */
 
 int solver_mpi_range(struct solver_data *solver);
 int solver_mpi(struct solver_data *solver, double timestep, double delt);
 int solver_mpi_high_rank(struct solver_data *solver);
 int solver_recv_all(struct solver_data *solver);
 int solver_send_recv_edge(struct solver_data *solver);
 int solver_send_all(struct solver_data *solver);
 int solver_mpi_send(struct solver_data *solver, double *data, int to, long int i_start, long int i_range);
 int solver_mpi_recv(struct solver_data *solver, double *data, int from, long int i_start, long int i_range);
 int solver_mpi_sendrecv(struct solver_data *solver, int to, double *send, long int send_i_start, long int send_i_range, int from, double *recv, long int recv_i_start, long int recv_i_range);
 int solver_broadcast_all(struct solver_data *solver);
 