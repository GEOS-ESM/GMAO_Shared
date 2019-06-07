/* filename: slustub.c
   The wrapper to SuperLU for factorizing and solving diagonal blocks
   of local problem on distributed environment.
   */
#include <stdlib.h>
#include <stdio.h>

#include "dsp_defs.h"
#include "util.h"
#include "Cnames.h"

#if (F77_CALL_C == ADD_)
#define SLU_fact_ slu_fact_
#define SLU_solve_ slu_solve_
#define SLU_solve_t_ slu_solve_t_
#define SLU_tidy_ slu_tidy_
#elif (F77_CALL_C == UPCASE)
#define SLU_fact_ SLU_FACT
#define SLU_solve_ SLU_SOLVE
#define SLU_solve_t_ SLU_SOLVE_T
#define SLU_tidy_ SLU_TIDY
#elif (F77_CALL_C == NOCHANGE)
#define SLU_fact_ slu_fact
#define SLU_solve_ slu_solve
#define SLU_solve_t_ slu_solve_t
#define SLU_tidy_ slu_tidy
#else
#define SLU_fact_ slu_fact_
#define SLU_solve_ slu_solve_
#define SLU_solve_t_ slu_solve_t_
#define SLU_tidy_ slu_tidy_
#endif

/* prototype of functions used but not defined in the headers */
extern void get_perm_c(int ispec, SuperMatrix *A, int *perm_c);

/* prototypes of the functions defined in this file */
void SLU_tidy_(void);
void SLU_fact_(int *nrow, double *aval, int *ja, int *ia, int *info);
void SLU_solve_(int *nrow, int *ncol, int *ldx, double *x, double *flops);
void SLU_solve_t_(int *nrow, int *ncol, int *ldx, double *x, double *flops);

/* The following four variables are saved so that we can factorize a
   matrix once and used it many times.  However, this structure also
   limits that only one factorization can be saved at any given time. */
/* the L and U factors of a sparse matrix */
static SuperMatrix *L=NULL, *U=NULL;
/* left and right permutation used during factorization */
static int *perm_r=NULL, *perm_c=NULL;

void SLU_tidy_(void)
{
  if (perm_r != NULL) {SUPERLU_FREE(perm_r); perm_r=NULL;}
  if (perm_c != NULL) {SUPERLU_FREE(perm_c); perm_c=NULL;}
  if (L != NULL) {Destroy_SuperNode_Matrix(L); L=NULL;}
  if (U != NULL) {Destroy_CompCol_Matrix(U); U=NULL;}
  return;
} /* end of SLU_tidy */

void SLU_fact_(int *nrow, double *aval, int *ja, int *ia, int *info)
{
  SuperMatrix A, AC;
  mem_usage_t mem_usage;
  int i, nnz, offset;
  int panel_size, relax;
  int *etree;
  char refact;

  /* the default size of 8 is used for panel_size and relax */
  panel_size = 8; relax = 4;
  StatInit(panel_size, relax);

  /* transform the base of array index to zero (0) */
  offset = ia[0];
  nnz = ia[*nrow] - ia[0];
  if (offset != 0) {
    for (i=0; i<=*nrow; ++i) ia[i] -= offset;
  }
  for (i=0; i<nnz; ++i) --ja[i];

  /* Construct SuperMatrix out of aval, ja, ia
     -- note that compressed column format is assumed,
     -- the matrix is also assumed to be square */
  dCreate_CompCol_Matrix(&A, *nrow, *nrow, nnz, aval, ja, ia, NC, _D, GE);
  SLU_tidy_();
  L = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  U = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  i = sizeof(int)* (*nrow);
  perm_c = (int *) malloc(i);
  perm_r = (int *) malloc(i);
  etree = (int *) malloc(i);
  if (L == NULL || U == NULL || perm_r == NULL || perm_c == NULL ||
      etree == NULL) {
    fprintf(stderr, "SLU_fact cannot allocate space for SuperLU.\n");
    exit(-1);
  }

  /* reorder using minimum degree on A^T*A and factor the matrix */
  refact = 'N';
  get_perm_c(1, &A, perm_c);
  sp_preorder(&refact, &A, perm_c, etree, &AC);
  dgstrf(&refact, &AC, 1.0, 0.0, relax, panel_size, etree, NULL, 0,
	 perm_r, perm_c, L, U, info);
  if (*info == 0) {
	dQuerySpace(L, U, panel_size, &mem_usage);
	*info = mem_usage.for_lu/sizeof(double);
  }
  else if (*info < 0) {
    printf("SLU_fact: the %d-th argument to DGSTRF is an illegal value.\n",
	   -(*info));
  }
  else if (*info > 0 && *info <= *nrow) {
    SLU_tidy_();
    printf("SLU_fact: the %d-th diagonal element of U is zero.\n", *info);
    *info = -15;
  }
  else {
    SLU_tidy_();
    printf("SLU_fact: unknown error code %d.\n", *info);
    *info = -16;
  }

  /* reclaim workspace */
  Destroy_CompCol_Permuted(&AC);
  SUPERLU_FREE(A.Store);
  SUPERLU_FREE(etree);

  /* restoring the shifted offset */
  if (offset != 0) {
    for (i=0; i<=*nrow; ++i) ia[i] += offset;
  }
  for (i=0; i<nnz; ++i) ++ja[i];
  return;
}  /* end of SLU_fact_ */

/* using the factorization L/U to solve a linear system 
   X is both right-hand size and the solution   */
void SLU_solve_(int *nrow, int *ncol, int *ldx, double *x, double *flops)
{
  extern SuperLUStat_t SuperLUStat;
  SuperMatrix B;
  char trans;
  int ierr;

  if ((L==NULL) || (U==NULL)) {
    printf("SLU_fact must be called before calling SLU_solve.\n");
    return;
  }

  trans = 'N';
  dCreate_Dense_Matrix(&B, *nrow, *ncol, x, *ldx, DN, _D, GE);
  dgstrs(&trans, L, U, perm_r, perm_c, &B, &ierr);
  SUPERLU_FREE(B.Store);
  if (ierr < 0) {
    printf("SLU_solve: the %d-th argument to dgstrs is an illegal value.\n");
  }
  else {
    *flops += (double)SuperLUStat.ops[SOLVE];
  }
  return;
} /* end of SLU_solve_ */
void SLU_solve_t_(int *nrow, int *ncol, int *ldx, double *x, double *flops)
{
  extern SuperLUStat_t SuperLUStat;
  SuperMatrix B;
  char trans;
  int ierr;

  if ((L==NULL) || (U==NULL)) {
    printf("SLU_fact must be called before calling SLU_solve.\n");
    return;
  }

  trans = 'T';
  dCreate_Dense_Matrix(&B, *nrow, *ncol, x, *ldx, DN, _D, GE);
  dgstrs(&trans, L, U, perm_r, perm_c, &B, &ierr);
  SUPERLU_FREE(B.Store);
  if (ierr < 0) {
    printf("SLU_solve: the %d-th argument to dgstrs is an illegal value.\n");
  }
  else {
    *flops += (double)SuperLUStat.ops[SOLVE];
  }
  return;
} /* end of SLU_solve_t_ */
