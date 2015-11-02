/* vof_macros.h
 *
 * Macros used in the VOF solver */
 
#ifndef _VOF_MACROS_H
#define _VOF_MACROS_H
 
#ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
    #define min(a,b) ((a) < (b) ? (a) : (b))
#endif

/* macros to reduce code needed for access to array elements */
#define U(i, j, k) solver->mesh->u[mesh_index(solver->mesh, i, j, k)]
#define V(i, j, k) solver->mesh->v[mesh_index(solver->mesh, i, j, k)]
#define W(i, j, k) solver->mesh->w[mesh_index(solver->mesh, i, j, k)]
#define UN(i, j, k) mesh_n->u[mesh_index(mesh_n, i, j, k)]
#define VN(i, j, k) mesh_n->v[mesh_index(mesh_n, i, j, k)]
#define WN(i, j, k) mesh_n->w[mesh_index(mesh_n, i, j, k)]
#define P(i, j, k) solver->mesh->P[mesh_index(solver->mesh, i, j, k)]
#define PN(i, j, k) mesh_n->P[mesh_index(mesh_n, i, j, k)]
#define D(i, j, k) solver->mesh->D[mesh_index(solver->mesh, i, j, k)]
#define DN(i, j, k) mesh_n->D[mesh_index(mesh_n, i, j, k)]

#define VOF(i, j, k) solver->mesh->vof[mesh_index(solver->mesh, i, j, k)]
#define VOF_N(i, j, k) mesh_n->vof[mesh_index(mesh_n, i, j, k)]
#define N_VOF(i, j, k) solver->mesh->n_vof[mesh_index(solver->mesh, i, j, k)]
#define N_VOF_N(i, j, k) mesh_n->n_vof[mesh_index(mesh_n, i, j, k)]


#ifdef FV
#undef FV
#endif
#define FV(i, j, k) solver->mesh->fv[mesh_index(solver->mesh, i, j, k)]

#ifdef AE
#undef AE
#endif
#define AE(i, j, k) solver->mesh->ae[mesh_index(solver->mesh, i, j, k)]

#ifdef AN
#undef AN
#endif
#define AN(i, j, k) solver->mesh->an[mesh_index(solver->mesh, i, j, k)]

#ifdef AT
#undef AT
#endif
#define AT(i, j, k) solver->mesh->at[mesh_index(solver->mesh, i, j, k)]

#define PETA(i, j, k) solver->mesh->peta[mesh_index(solver->mesh, i, j, k)]
#define BETA(i, j, k) solver->mesh->beta[mesh_index(solver->mesh, i, j, k)]
#define TANTH(i, j, k) solver->mesh->tanth[mesh_index(solver->mesh, i, j, k)]

#define NUT(i, j, k) solver->mesh->nut[mesh_index(solver->mesh, i, j, k)]
#define NUT_N(i, j, k)  mesh_n->nut[mesh_index(solver->mesh, i, j, k)]

#define DELX solver->mesh->delx
#define DELY solver->mesh->dely
#define DELZ solver->mesh->delz
#define RDX solver->mesh->rdx
#define RDY solver->mesh->rdy
#define RDZ solver->mesh->rdz
#define IMAX solver->mesh->imax
#define JMAX solver->mesh->jmax
#define KMAX solver->mesh->kmax


/* access functions */

#ifdef DEBUG
float u(struct solver_data *solver, long int i,long int j,long int k);
float v(struct solver_data *solver, long int i,long int j,long int k);
float w(struct solver_data *solver, long int i,long int j,long int k);
float un(struct solver_data *solver, long int i,long int j,long int k);
float vn(struct solver_data *solver, long int i,long int j,long int k);
float wn(struct solver_data *solver, long int i,long int j,long int k);
float p(struct solver_data *solver, long int i,long int j,long int k);
float vof(struct solver_data *solver, long int i,long int j,long int k);
float n_vof(struct solver_data *solver, long int i,long int j,long int k);
float ae(struct solver_data *solver, long int i,long int j,long int k);
float an(struct solver_data *solver, long int i,long int j,long int k);
float at(struct solver_data *solver, long int i,long int j,long int k);
float fv(struct solver_data *solver, long int i,long int j,long int k);
float peta(struct solver_data *solver, long int i,long int j,long int k);
void track_cell(struct solver_data *solver, long int i,long int j,long int k);
#endif

#endif
