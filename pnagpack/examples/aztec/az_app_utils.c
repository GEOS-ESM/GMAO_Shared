/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id$";
#endif


/******************************************************************************
 *Copyright 1995, Sandia Corporation.  The United States Government retains a *
 *nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 *Export of this program may require a license from the United States         *
 *Government.                                                                 *
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
extern int AZ_using_fortran;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_options(int options[], double params[])

{

  /*
   * Choose among AZTEC options (see User's Guide).
   */

  AZ_defaults(options, params);

  options[AZ_solver]   = AZ_cgs;
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_ls;
  options[AZ_conv]     = AZ_r0;
  options[AZ_output]   = 1;
  options[AZ_pre_calc] = AZ_calc;
  options[AZ_max_iter] = 1550;
  options[AZ_poly_ord] = 5;
  options[AZ_overlap]  = AZ_none;
  options[AZ_kspace]   = 130;
  options[AZ_orthog]   = AZ_modified;
  options[AZ_aux_vec]  = AZ_resid;

  params[AZ_tol]       = 4.00e-9;
  params[AZ_drop]      = 0.0;

} /* init_options */

void my_read_update(int *N_update, int *update[], int proc_config[],
                    int N, int chunk)

/******************************************************************************

  This is modified from AZ_read_update to partition the given grid into
  sub-cube if the problem size and NPE are both cubic, else linear partition
  is used.

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        On Output, number of unknowns updated by this processor.

  update:          On Output, list of unknowns updated by this processor in
                   ascending order.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

  N:               Total number of chunks to be distributed.

  chunk:           Size of each chunk to be treated as a single unit.
                   The unknowns contained in the kth chunk are given
                   by {k*chunk, k*chunk + 1, ..... , (k+1)*chunk - 1}
                   and 'N*chunk' is the total number of unknowns to be
                   distributed.
******************************************************************************/

{

  /* local variables */

  int   t1, t2, i;
  int   ii, j;
  int   allocated, length;
  int   cflag;
  int   partner;
  int   proc_x, proc_y, proc_z,m;
  int   pts_x, pts_y, pts_z;
  int   px, py, pz, k;
  int   start_x, start_y, start_z;
  int   end_x, end_y, end_z;
  int   pt_number;
  int   count, check;
  int   proc, nprocs;
  int   type, type2;
  FILE *fp = NULL;
  MPI_Request request;  /* Message handle */


  /**************************** execution begins *****************************/

  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  /*
   * Figure out which chunks should be assigned to this processor using a box
   * decomposition.  That is, it is assumed that all the chunks are ordered
   * naturally corresponding to an m x m x m box where m = N^(1/3).  Boxes of
   * chunks are assigned to processors.
   */

  m = (int) (pow( (double) N, 1.0/3.0)+0.5);
  if (m*m*m != N) {
    printf("Number of grid points must be an integer cubed.\n");
    exit(-1);
  }
  proc_x = (int) (pow( (double) nprocs, 1.0/3.0)+0.5);

  if (proc_x*proc_x*proc_x == nprocs) {
    pts_x = m; pts_y = m; pts_z = m;
    proc_y = proc_x; proc_z = proc_x;

    if ( proc_x*proc_y*proc_z != nprocs) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: %d x %d x %d != %d ",
 			 proc_x, proc_y, proc_z, nprocs);
          (void) fprintf(stdout," (total number of processors)\n");
        }
	exit(1);
    }

    if ( pts_x * pts_y * pts_z != N ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: %d x %d x %d != %d ",
 			 pts_x, pts_y, pts_z, N);
          (void) fprintf(stdout," (total number of grid points)\n");
        }
	exit(1);
    }
    if ( pts_x%proc_x != 0 ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: grid points along x axis are not an ");
          (void) fprintf(stdout,"even multiple of processors\n");
	  (void) fprintf(stdout,"       along x axis.");
        }
	exit(1);
    }
    if ( pts_y%proc_y != 0 ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: grid points along y axis is not an ");
          (void) fprintf(stdout,"even multiple of processors\n");
	  (void) fprintf(stdout,"       along y axis.");
        }
	exit(1);
    }
    if ( pts_z%proc_z != 0 ) {
        if (proc == 0) {
          (void) fprintf(stdout,"Error: grid points along z axis is not an ");
          (void) fprintf(stdout,"even multiple of processors\n");
	  (void) fprintf(stdout,"       along z axis.");
        }
	exit(1);
    }
    pts_x = pts_x/proc_x;
    pts_y = pts_y/proc_y;
    pts_z = pts_z/proc_z;

    *N_update = pts_x * pts_y * pts_z * chunk;
    if (!AZ_using_fortran) *update  = (int *) calloc(*N_update, sizeof(int));

    /* compute the lower left corner and the upper right corner */

    px = proc % proc_x;
    pz = (proc-px) / proc_x;
    py = pz % proc_y;
    pz = (pz-py) / proc_y;

    start_x = px * pts_x;
    end_x   = start_x + pts_x;
    start_y = py * pts_y;
    end_y   = start_y + pts_y;
    start_z = pz * pts_z;
    end_z   = start_z + pts_z;

    /* set update[] */

    count = 0;
    for (k = start_z; k < end_z; k++ ) {
      for (j = start_y; j < end_y; j++ ) {
        for (i = start_x; i < end_x; i++ ) {
          for (ii = 0; ii < chunk; ii++ ) {
            pt_number = (i + j * m + k * m * m) * chunk + ii;
            (*update)[count++] = pt_number;
          }
        }
      }
    }
  }

  else {

    /*
     * Figure out which chunks should be assigned to this processor for linear
     * partitioning.  This means that processor 0 is assigned the chunks
     * approximately corresponding to 0, ... , N/nprocs and processor 1 is
     * approximately assigned the chunks 1+N/nprocs to 2*N/nprocs.
     */

    t1 = N/nprocs;
    t2 = N - t1 * nprocs;

    if ( proc >= t2) t2 += (proc * t1);
    else {
      t1++;
      t2    = proc*t1;
    }
    *N_update = t1*chunk;
    t2   *= chunk;

    if (!AZ_using_fortran) *update = (int *) calloc(*N_update,sizeof(int));
    if ( (*update == NULL) && (*N_update != 0)) {
      (void) fprintf (stderr, "Not enough space to allocate 'update'\n");
      exit(-1);
    }

    for (i = 0; i < *N_update; i++) (*update)[i] = i + t2;
  }

} /* my_read_update */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_guess_and_rhs(int update_index[], int update[], double *x[],double
                        *ax[],int data_org[], double val[], int indx[], int
                        bindx[], int rpntr[], int cpntr[], int bpntr[], int
                        proc_config[])

/*
 * Set the initial guess and the right hand side where the right hand side
 * is obtained by doing a matrix-vector multiplication.
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date :  3/15/95
 *
 * Parameters
 *
 *    update_index   ==      On input, ordering of update and external
 *                           locally on this processor. For example
 *                           'update_index[i]' gives the index location
 *                           of the block which has the global index
 *                           'update[i]'.
 *    update         ==      On input, list of pts to be updated on this node
 *    data_org       ==      On input, indicates how data is set on this node.
 *                           For example, data_org[] contains information on
 *                           how many unknowns are internal, external and
 *                           border unknowns as well as which points need
 *                           to be communicated. See User's Guide for more
 *                           details.
 *    val, indx,     ==      On input, holds matrix nonzeros. See User's Guide
 *    bindx, rpntr,          for more details.
 *    cpntr, bpntr
 *    x              ==      On output, 'x' is allocated and set to all zeros.
 *    ax             ==      On output, 'ax' is allocated and is set to the
 *                           result of a matrix-vector product.
 */

{

  int i,j;
  int temp,num;
  extern int num_PDE_eqns;   /* number of PDEs being solved */
  extern int application;    /* problem being solved        */
  double sum = 0.0;

  temp = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  num  = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* allocate vectors */

  i       = num + data_org[AZ_N_external];
  *x      = (double *) calloc(i, sizeof(double));
  *ax     = (double *) calloc(i, sizeof(double));
  if ((*ax == NULL) && (i != 0)) {
    (void) fprintf(stderr, "Not enough space in init_guess_and_rhs() for ax\n");
    exit(1);
  }

  /* initialize 'x' to a function which will be used in matrix-vector product*/

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
    for (i = 0; i < temp; i++) {
      for (j = rpntr[i]; j < rpntr[i+1]; j++) {
        (*x)[j] = (double) (update[i]) + (double)(j-rpntr[i]) /
          (double)(num_PDE_eqns);
      }
    }
  }
  else {
    for (i = 0; i < temp; i++) {
      (*x)[i] = (double) (update[i]) / (double) (num_PDE_eqns);
    }
  }

  /* Reorder 'x' so that it conforms to the transformed matrix */
 
  AZ_reorder_vec(*x,data_org,update_index,rpntr);

  if (application == 2) {

    /* take out the constant vector. Used for the */
    /* finite element problem because it is singular */

    sum = AZ_gsum_double(sum, proc_config);
    i   = AZ_gsum_int(num, proc_config);
    if (i != 0) sum = sum / ((double) i);
    for (i = 0; i < num; i++) (*x)[i] -= sum;
  }

  AZ_matvec_mult(val, indx, bindx, rpntr, cpntr, bpntr, *x, *ax, 1, data_org);

  for (i = 0; i < num; i++) (*x)[i] = 0.0;

} /* init_guess_and_rhs */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_matrix_vector_structures(int proc_config[], int *update_index[], int
                                   *update[], int  *data_org[], int *external[],
                                   int *extern_index[], int input_option, double
                                   *val[], int *bindx[], int *indx[], int
                                   *bpntr[], int *rpntr[], int *cpntr[])


/*
 * Read in the points to be updated on this processor, create the global
 * distributed form of the application matrix, and then convert it to a
 * local distributed form for AZTEC kernels. Along the way, initialize the
 * following quantities:
 *     update_index[], update[], data_org[], a[], bindx[], bpntr[], cpntr[],
 *     rpntr[], indx[], external[], extern_index[].
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date:   3/15/95
 *
 * Parameters
 *
 *    proc_config    ==      On input, processor information:
 *                              proc_config[AZ_node] = name of this processor
 *                              proc_config[AZ_N_procs] = # of processors used
 *    update         ==      On output, list of pts to be updated on this node
 *    val,bindx      ==      On output, local distributed form of arrays
 *                           holding matrix values
 *    external       ==      On output, list of external vector elements
 *    update_index   ==      On output, ordering of update and external
 *    extern_index   ==      locally on this processor. For example
 *                           'update_index[i]' gives the index location
 *                           of the block which has the global index
 *                           'update[i]'.
 *    data_org       ==      On output, indicates how the data is set out on
 *                           this node. For example, data_org[] contains
 *                           information on how many unknowns are internal,
 *                           external, and border unknowns as well as which
 *                           points need to be communicated. See User's Guide
 *                           for more details.
 *    input_option   ==      Indicates how update[] will be initialized.
 *                           = 0, linear decomposition
 *                           = 1, points read from file 'update'.
 *                           = 2, box decomposition
 *                           See AZ_read_update() comments for more details.
 *
 *      The default finite difference MSR problem corresponds to a setting up
 *      a series of uncoupled 3D Poisson equations on a cube.
 *      To solve other problems, the call 'add_row_3D(...)' in
 *      'create_msr_matrix()' can be changed to 'add_row_5pt()' or
 *      'add_row_9pt()'.
 */

{

  int    N_update;            /* Number of pts updated on this processor     */
  int    MSRorVBR;
  int    chunks;
int blk_size, num_blk_cols,num_blk_rows,size,kk, convert_to_vbr = 0;
double *val2;
int    *bindx2;

  extern int N_grid_pts,      /* N_grid_pts is the total number of grid
                                 points used in the simulation.              */
    num_PDE_eqns;             /* num_PDE_eqns refers to the number of PDEs
                                 that are being solved.                      */

  /* ----------------- external function declarations ------------------- */

  extern void create_msr_matrix(int update[], double **, int **bindx,
                                int N_update);
  extern void AZ_transform(int proc_config[], int *external[], int bindx[],
                           double a[], int update[], int *update_index[],
                           int *extern_index[], int *data_org[], int, int *,
                           int *, int *, int **, int mat_type);
  extern void create_fe_matrix(int update[], int proc, int **bindx, double **,
                               int N_update);
  extern void create_vbr_matrix(int update[], double **, int **indx,
                                int N_update, int **rpntr, int **bpntr,
                                int **bindx);
  extern int application;

  MSRorVBR = AZ_MSR_MATRIX;
  if (application == 1) MSRorVBR = AZ_VBR_MATRIX;

  chunks = num_PDE_eqns;
  if (MSRorVBR == AZ_VBR_MATRIX) chunks = 1;

  /* initialize the list of global indices. NOTE: the list of global */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly. */

  my_read_update(&N_update, update, proc_config, N_grid_pts, chunks);

  /* create the matrix: each processor creates only the      */
  /* rows appearing in update[] ... however this row is      */
  /* created as if it were on a serial machine (i.e. using   */
  /* the global column numbers)                              */

  if (application == 1)
    create_vbr_matrix(*update, val, indx, N_update, rpntr, bpntr, bindx);
  else {
    *indx = NULL; *bpntr = NULL; *rpntr = NULL; *cpntr = NULL;

    if (application == 0) create_msr_matrix(*update, val, bindx, N_update);
    if (application == 2) create_fe_matrix(*update, proc_config[AZ_node],
                                           bindx, val, N_update);
    if (application == 3) { 
        AZ_read_msr_matrix(*update, val, bindx, N_update, proc_config);
    }
  }

  /* convert matrix to a distributed parallel matrix */

  AZ_transform(proc_config, external, *bindx, *val, *update, update_index,
               extern_index, data_org, N_update, *indx, *bpntr, *rpntr, cpntr,
               MSRorVBR);

  if ( (convert_to_vbr == 1) && (application == 3) ) {
     if (proc_config[AZ_node] == 0 ) {
	 printf("enter the block size\n");
	 scanf("%d",&blk_size);
     }
     AZ_broadcast((char *) &blk_size,  sizeof(int), proc_config, AZ_PACK);
     AZ_broadcast((char *) NULL         , 0          , proc_config, AZ_SEND);

     if ( N_update%blk_size != 0 ) {
        (void) fprintf(stderr," The block size must be a multiple of the number of rows per processor.\n");
        exit(-1);
     }

     num_blk_rows = N_update/blk_size;
     num_blk_cols = ( (*data_org)[AZ_N_external] + N_update)/blk_size;
     *cpntr = (int *) calloc(num_blk_cols+2, sizeof(int));
     *rpntr = (int *) calloc(num_blk_cols+2, sizeof(int));
     *bpntr = (int *) calloc(num_blk_cols+2, sizeof(int));
     size   = 20*(num_blk_cols+2);
     *indx  =  (int *) calloc(size, sizeof(int));
     bindx2 = *bindx;
     val2   = *val;
     *bindx = (int *) calloc(size, sizeof(int));
     *val   =  (double *) calloc(size*blk_size*blk_size, sizeof(double));

     for (kk = 0 ; kk < num_blk_cols ; kk++ ) (*cpntr)[kk] = blk_size;
     AZ_msr2vbr(*val,*indx,*rpntr,*cpntr,*bpntr,*bindx,bindx2,val2,
		num_blk_rows,num_blk_cols,size,size*blk_size*blk_size,blk_size);
     MSRorVBR = AZ_VBR_MATRIX;
     N_update /= blk_size;
     num_PDE_eqns = blk_size; 
     for (kk = 0 ; kk < N_update ; kk++ )
           (*update)[kk] = (*update)[blk_size*kk]/blk_size;
     for (kk = 0 ; kk < (*data_org)[AZ_N_external] ; kk++ ) 
           (*external)[kk] = (*external)[blk_size*kk]/blk_size;

     (*data_org)[AZ_matrix_type] = AZ_VBR_MATRIX;
     (*data_org)[AZ_N_int_blk ] /= blk_size;
     (*data_org)[AZ_N_bord_blk] /= blk_size;
     (*data_org)[AZ_N_ext_blk ] /= blk_size;
     free(bindx2);  free(val2);
  }


} /* init_matrix_vector_structures */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_global_element(int element,int update[],int data_org[],
           int update_index[], int rpntr[], double vector[],int proc_config[])
{
/*
 * Print out the vector element corresponding to the global numbering
 * 'element'. Note: if the VBR format is used, this routine will print
 * out all the vector elements corresponding to this block.
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date:   6/15/96
 *
 * Parameters
 *
 *    element        ==      On input, global number of vector element that
 *                           will be printed.
 *    update         ==      On input, list of pts updated on this node
 *    data_org       ==      On input, indicates how the data is set out on
 *                           this node. For example, data_org[] contains
 *                           information on how many unknowns are internal,
 *                           external, and border unknowns as well as which
 *                           points need to be communicated. See User's Guide
 *                           for more details.
 *    update_index   ==      On input, ordering of update locally on this
 *                           processor. For example, 'update_index[i]' gives 
 *                           the index location of the block which has the 
 *                           global index 'update[i]'.
 *    rpntr          ==      On input, rpntr[i+1]-rpntr[i] gives the block
 *                           size of the ith local block.
 *    vector         ==      On input, vector to be printed (just one element).
 *    proc_config    ==      On input, processor information:
 *                              proc_config[AZ_node] = name of this processor
 *                              proc_config[AZ_N_procs] = # of processors used
 */
   int i,k,count;

   /* synchronize things */

   i = AZ_gsum_int(1,proc_config);
 


   i = AZ_find_index(element,update,
                     data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk]);
   if (i !=-1) {
      i = update_index[i];
      if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
         fprintf(stdout,"(%d) = %e\n",element,vector[i]);
      else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
        for (k = rpntr[i]; k < rpntr[i+1]; k++ ) 
           fprintf(stdout,"(%d,%d) = %e\n",element,k-rpntr[i],vector[k]);
      }
      fflush(stdout);
   }

   /* synchronize things */
   i = AZ_gsum_int(i,proc_config);

}

