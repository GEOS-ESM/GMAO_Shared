/* filename: azlan.c
   This file implements a driver for simple eigenvalue problem using
   pLANSO and Aztec.
   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "az_aztec.h"

/* a macro function to convert between Fortran and C function name */
#if defined(append_)
#define F77_PLANDR plandr_
#define F77_OP op_
#define F77_OPM opm_
#define F77_STORE store_
#define F77_PPURGE ppurge2_
#define F77_DGEMV dgemv_
#elif defined(caps)
#define F77_PLANDR PLANDR
#define F77_OP OP
#define F77_OPM OPM
#define F77_STORE STORE
#define F77_PPURGE PPURGE2
#define F77_DGEMV SGEMV
#endif

/* global variables needed by matrix-vector multiplications routine OP */
static double *az_val, *lanso_vec=0;
static int *az_indx, *az_bindx, *az_rpntr, *az_cpntr, *az_bpntr, *az_data_org,
  az_exchange_flag;
static double mmvtime=0, kmvtime=0, orthtime=0;
static int mmvcnt=0, kmvcnt=0, north=0;

/* function prototypes */
extern void F77_PLANDR(int *, int *, int *, double *, double *, double *,
		       int *, double *, int *, int *, double *, double *,
		       double *, int *, int *, int *, int *);

extern void AZ_matvec_mult(double *val, int *indx, int *bindx, int *rpntr,
			   int *cpntr, int *bpntr, double *b,
			   register double *c, int exchange_flag,
			   int *data_org);

/* matrix vector multiplicaiton routine -- need to copy input vector before
   calling Aztec function.
*/
void F77_OP(int *ndim, double *pp, double *qq, double *rr, int *mpicom)
{
  int i;
  double tmptime;
  register double *r2;
  static int length=0;
  static double *tmp=NULL;

  tmptime = MPI_Wtime();
  if (length < (*ndim)+az_data_org[AZ_N_external]) {
    if (tmp != NULL) free(tmp);
    length = (*ndim) + az_data_org[AZ_N_external];
    tmp = (double *) malloc(length*sizeof(double));
    if (tmp == NULL) {
      fprintf(stderr, "F77_OP: unable to allocate workspace of size %d.\n",
	      length);
      length = 0;
      i=MPI_Abort(MPI_COMM_WORLD, length);
    }
  }

  r2 = rr;
  for (i=0; i<*ndim; ++i) tmp[i] = qq[i];
  AZ_matvec_mult(az_val, az_indx, az_bindx, az_rpntr, az_cpntr, az_bpntr,
		 tmp, r2, az_exchange_flag, az_data_org);

  kmvcnt++;
  kmvtime += MPI_Wtime() - tmptime;
} /* end of OP */

void F77_OPM(int *ndim, double *pp, double *qq, int *mpicom)
{
  int i;
  double tmptime;

  tmptime = MPI_Wtime();
  for (i=0; i< *ndim; ++i) qq[i] = pp[i];

  mmvcnt++;
  mmvtime += MPI_Wtime() - tmptime;
} /* end of OPM */

/* the function store is required by pLANSO to access the Lanczos vectors */
void F77_STORE(int *ndim, int *isw, int *ind, double *vec)
{
  int i;
  double *s, *t;
  switch (*isw) {
  case 1:
    s = vec; t = lanso_vec+(*ndim)*(*ind + 1);
    break;
  case 2:
    s = lanso_vec+(*ndim)*(*ind + 1); t = vec;
    break;
  case 3:
    s = vec; t = lanso_vec+(*ndim)*(*ind - 1);
    break;
  case 4:
    s = lanso_vec+(*ndim)*(*ind - 1); t = vec;
    break;
  default:
    printf("STORE: unknown switching flag %d.\n", *isw);
    return;
  }
  for (i=0; i<*ndim; ++i, ++s, ++t) *t = *s;
  return;
} /* end of store */

/* PPURGE function --
   called to decide whether to perform re-orthogonalization
   */
void F77_PPURGE(int *ndim, int *jnd, double rr[], double qq[],
		double ra[], double qa[], double wrk[], double eta[],
		double oldeta[], double *reps1, double *rnm,
		int *msglvl, int *mpicom)
{
  char trans, notrans;
  int i, k, ierr, myid, loop;
  double one, zero, mone, tmp, tmptime;

  tmptime = MPI_Wtime();
  if (*msglvl > 11) {
    ierr = MPI_Comm_rank((MPI_Comm)(*mpicom), &myid);
    if (ierr != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, ierr);
  }
  else
    myid = -1;

  /* perform reorthogonalization -- first orthogonalize q */
  trans='T';
  notrans='N';
  one = 1.0;
  zero=0.0;
  mone=-1.0;
  i = 1;
  k = (*jnd) - 1;
#if defined(caps)
  F77_DGEMV(_cptofcd(&trans,i), ndim, &k, &one, lanso_vec+(*ndim)+(*ndim),
	    ndim, qa, &i, &zero, eta, &i);
#else
  F77_DGEMV(&trans, ndim, &k, &one, lanso_vec+(*ndim)+(*ndim), ndim,
	    qa, &i, &zero, eta, &i, i);
#endif
  ierr = MPI_Allreduce(eta, oldeta, k, MPI_DOUBLE, MPI_SUM,
		       (MPI_Comm)(*mpicom));
  if (ierr != MPI_SUCCESS) {
    printf("PPURGE failed to perform MPI_Allreduce(error code = %d).\n",
	   ierr);
    MPI_Abort(MPI_COMM_WORLD, ierr);
  }
#if defined(caps)
  F77_DGEMV(_cptofcd(&notrans,i), ndim, &k, &mone, lanso_vec+(*ndim)+(*ndim),
	    ndim, oldeta, &i, &one, qq, &i);
#else
  F77_DGEMV(&notrans, ndim, &k, &mone, lanso_vec+(*ndim)+(*ndim), ndim,
	    oldeta, &i, &one, qq, &i, i);
#endif
  *wrk = 0.0;
  for (i=0; i<*ndim; ++i) *wrk += qq[i]*qq[i];
  ierr = MPI_Allreduce(wrk, &tmp, 1, MPI_DOUBLE, MPI_SUM, (MPI_Comm)(*mpicom));
  if (ierr != MPI_SUCCESS) {
    printf("PPURGE failed to perform MPI_Allreduce(error code = %d).\n",
	   ierr);
    MPI_Abort(MPI_COMM_WORLD, ierr);
  }
  tmp = 1.0/sqrt(tmp);
  for (i=0; i<*ndim; ++i) {
    qq[i] *= tmp;
    qa[i] = qq[i];
  }
  if (myid == 0) {
    tmp = 0.0;
    for (i=0; i<k; ++i) tmp += oldeta[i]*oldeta[i];
    printf("Purge q_%d (error=%e).\n", *jnd, sqrt(tmp));
  }

  /* orthogonalize r */
  loop = 0;
  do {
    i = 1;
#if defined(caps)
    F77_DGEMV(_cptofcd(&trans,i), ndim, &k, &one, lanso_vec+(*ndim)+(*ndim),
	      ndim, ra, &i, &zero, eta, &i);
#else
    F77_DGEMV(&trans, ndim, &k, &one, lanso_vec+(*ndim)+(*ndim), ndim,
	      ra, &i, &zero, eta, &i, i);
#endif
    eta[k]=0.0;
    for (i=0; i<*ndim; ++i) eta[k] += qq[i]*ra[i];
    i = 1;
    ierr = MPI_Allreduce(eta, oldeta, *jnd, MPI_DOUBLE, MPI_SUM,
			 (MPI_Comm)(*mpicom));
    if (ierr != MPI_SUCCESS) {
      printf("PPURGE failed to perform MPI_Allreduce(error code = %d).\n",
	     ierr);
      MPI_Abort(MPI_COMM_WORLD, ierr);
    }
#if defined(caps)
    F77_DGEMV(_cptofcd(&notrans,i), ndim, &k, &mone, lanso_vec+(*ndim)+(*ndim),
	      ndim, oldeta, &i, &one, rr, &i);
#else
    F77_DGEMV(&notrans, ndim, &k, &mone, lanso_vec+(*ndim)+(*ndim), ndim,
	      oldeta, &i, &one, rr, &i, i);
#endif
    *rnm = 0.0;
    for (i=0; i<*ndim; ++i) {
      rr[i] -= oldeta[k]*qq[i];
      *rnm += rr[i]*rr[i];
      ra[i] = rr[i];
    }
    ierr = MPI_Allreduce(rnm, wrk, 1, MPI_DOUBLE, MPI_SUM,
			 (MPI_Comm)(*mpicom));
    if (ierr != MPI_SUCCESS) {
      printf("PPURGE failed to perform MPI_Allreduce(error code = %d).\n",
	     ierr);
      MPI_Abort(MPI_COMM_WORLD, ierr);
    }
    tmp = 0.0;
    for (i=0; i<*jnd; ++i) tmp += oldeta[i]*oldeta[i];
    ++loop;
    if (myid == 0)
      printf("Purge r_%d (Loop #%d) %e/%e.\n", *jnd, loop, sqrt(tmp),
	     sqrt(*wrk));
  } while (tmp > (*reps1)*(*wrk) && loop<5);
  *rnm = sqrt(*wrk);

  north++;
  orthtime = MPI_Wtime() - tmptime;
}  /* end of PPURGE */

/* Plain pLANSO caller:
   It allocate MAX_MEM byte of total space (includes matrix, workspace, ...)
   and call planso to compute as many eigenpairs as possible.

   Without computing eigenvectors, the workspace required by pLANSO is
   5*ndim + 4*lanmax + max(ndim, lanmax+1) + 1

   The space required to store Lanczos vectors are (lanmax+2)*ndim.
   */
void plain_lanso(int *options, int *params, int *indx, int *bindx, int *rpntr,
		 int *cpntr, int *bpntr, double *val, int *data_org,
		 int *status, int *proc_config)
{
  int i, ndim, lanmax, ierr, maxprs, my_id, nw, ev, msglvl, neig, jnd;
  MPI_Comm mpicom;
  double lanso_time, condm, kappa, endl, endr;
  double *lanso_work, *ritz, *bnd;

#if defined(MSGLVL)
  msglvl = MSGLVL;
#elif defined(DEBUG)
  msglvl = 100;
#else
  msglvl = 1;
#endif
  /* Aztec can only use MPI_COMM_WORLD */
  mpicom = MPI_COMM_WORLD;

  /* copy matrix information to global variables */
  az_indx = indx;
  az_bindx = bindx;
  az_rpntr = rpntr;
  az_cpntr = cpntr;
  az_bpntr = bpntr;
  az_val = val;
  az_data_org = data_org;
  lanso_vec = NULL;
  lanso_work = NULL;
  ritz = NULL;
  bnd = NULL;

  /* allocate workspace for pLANSO */
  ndim = data_org[AZ_N_internal] + data_org[AZ_N_border];
  lanmax = (MAX_MEM_SIZE/sizeof(double) - (8*ndim + 1)) / (ndim + 7);
  ierr = MPI_Allreduce(&ndim, &i, 1, MPI_INT, MPI_SUM, mpicom);
  if (lanmax > i) lanmax = i;
  nw = 6*ndim+5*lanmax;
  lanso_vec = (double *)malloc(sizeof(double)*ndim*(lanmax+2));
  lanso_work = (double *)malloc(sizeof(double)*nw);
  ritz = (double *)malloc(sizeof(double)*lanmax);
  bnd = (double *)malloc(sizeof(double)*lanmax);
  if (lanso_vec == NULL || lanso_work == NULL || ritz == NULL ||
      bnd == NULL) {
    fprintf(stderr, "plain_lanso: unable allocate desired workspace.\n");
    fprintf(stderr, "       Reduce MAX_MEM_SIZE defined in Makefile.\n");
    exit(-1);
  }

  /* initialize the workspace to set starting vector to [1, 1, ..., 1]^T */
  for (i=0; i<ndim; ++i) lanso_work[i] = 1.0;
  ierr = MPI_Comm_rank(mpicom, &my_id);

/*
  if (msglvl>1 && my_id==0) printf("azl: testing maxvec ...\n");
  lanso_time = MPI_Wtime();
  for (i=0; i<100; ++i)
    F77_OP(&ndim, lanso_work, lanso_work, lanso_work+ndim, (int *)&mpicom);
  lanso_time = 0.01*(MPI_Wtime() - lanso_time);
  ierr = MPI_Allreduce(&lanso_time, &endr, 1, MPI_DOUBLE, MPI_SUM, mpicom);
  ierr = MPI_Allreduce(&lanso_time, &endl, 1, MPI_DOUBLE, MPI_MAX, mpicom);
  if (my_id == 0)
    printf("azl matvec time: (%i PE) %f, (average) %f, (max) %f (sec).\n",
	   proc_config[AZ_N_procs], endr, endr/proc_config[AZ_N_procs], endl);
*/
  if (msglvl>1 && my_id==0) printf("azl: calling PLANDR ...\n");

  jnd = 0;
  ierr = 0; 
  neig = 0;
  ev = 0;
  condm = 1.0;
  kappa = 1.0e-8;
  endl = -1.0e-30;
  endr = 1.0e-30;
  maxprs = lanmax;
  /* call the driver for LANSO */
  lanso_time = MPI_Wtime();
  F77_PLANDR(&ndim, &lanmax, &maxprs, &condm, &endl, &endr, &ev, &kappa,
	     &jnd, &neig, ritz, bnd, lanso_work, &nw, &ierr, &msglvl,
	     (int *)&mpicom);
  lanso_time = MPI_Wtime() - lanso_time;
  if (ierr != 0) {
    fprintf(stderr, "PLANDR returned with error code %i on PE #%i.\n",
	    ierr, my_id);
  }
  lanso_work[0] = lanso_time;
  lanso_work[1] = kmvtime;
  lanso_work[2] = mmvtime;
  lanso_work[3] = orthtime;
  ierr = MPI_Allreduce(lanso_work, lanso_vec, 4, MPI_DOUBLE, MPI_SUM, mpicom);
  ierr = MPI_Allreduce(&lanso_time, &endl, 1, MPI_DOUBLE, MPI_MAX, mpicom);
  lanso_time = endl;

  if (jnd > 0)
    condm = (fabs(ritz[0])>fabs(ritz[jnd-1]))?fabs(ritz[0]):fabs(ritz[jnd-1]);
  else
    condm = -1.0;
  endl = kappa*condm;
  /* print out some information about the solution process */
  if (my_id == 0) {
    printf("PLANSO steps = %i (lanmax=%i),  ierr=%i\n", jnd, lanmax, ierr);
    printf("PLANSO time: (%d PE) %f, (average) %f, (max) %f (sec).\n",
	   proc_config[AZ_N_procs], lanso_vec[0],
	   lanso_vec[0]/proc_config[AZ_N_procs], lanso_time);
    printf("PLANSO MATVEC (K %d %f) (M %d %f).\n", kmvcnt,
	   lanso_vec[1]/proc_config[AZ_N_procs], mmvcnt,
	   lanso_vec[2]/proc_config[AZ_N_procs]);
    printf("PLANSO %d re-orthogonalization used %f sec.\n", north,
	   lanso_vec[3]/proc_config[AZ_N_procs]);
    printf("PLANSO computed %i eigenvalues, est. spectral radius=%10.2e .\n",
	   neig, condm);
    if (jnd > 0 && msglvl > 3) {
      printf("Ritz[%i] = %25.17g (%10.2e)\n", 1, ritz[0], bnd[0]);
      for (i=1; i<jnd-1; ++i) {
	if (bnd[i] < endl)
	  printf("Ritz[%i] = %25.17g (%10.2e)\n", i+1, ritz[i], bnd[i]);
      }
      i = jnd-1;
      printf("Ritz[%i] = %25.17g (%10.2e)\n", i+1, ritz[i], bnd[i]);
    }
  }

  free(ritz);
  free(bnd);
  free(lanso_work);
  free(lanso_vec);
  return;
} /* end of plain LANSO */
