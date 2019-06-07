/* converted from Aztec example az_main.c to run azlan.c. - K. Wu, Aug, 1997 */

/******************************************************************************
 * Sample driver for pLANSO using AZTEC package for performing matrix-vector
 * multplicaitons.  The test problem is a matrix generated from finite
 * difference descritization of Poisson equation on 3-D cube.  The matrix is
 * is generated using AZTEC code as well.
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "az_aztec.h"

int N_grid_pts;                 /* number of grid points in application      */

int num_PDE_eqns;               /* number of equations (usually PDEs)        */
                                /* associated with each grid point.          */

int application;                /*  0: Finite Difference Poisson (MSR format)*/

void main(int argc, char *argv[])
{

  int    i,input_option,c;
  extern char *optarg;
  extern int optind;
  double total_time, average_time;

  /* See Aztec User's Guide for more information   */
  /* on the variables that follow.                 */

  int    proc_config[AZ_PROC_SIZE];/* Processor information:                 */
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.   */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.        */
  int    *data_org;                /* Array to specify data layout */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve()
                                      indicating success or failure.         */

  int    *update,                  /* vector elements (global index) updated
                                      on this processor.                     */
    *external;                     /* vector elements needed by this node.   */

  int    *update_index;            /* ordering of update[] and external[]    */
  int    *extern_index;            /* locally on this processor. For example
                                      update_index[i] gives the index
                                      location of the vector element which
                                      has the global index 'update[i]'.      */

                                   /* Sparse matrix to be solved is stored
                                      in these arrays.                       */
  int    *rpntr,*cpntr,*indx, *bpntr, *bindx;
  double *val;

  /* -------------  external function declarations ------------------------- */

  extern void init_matrix_vector_structures(int proc_config[],
                                            int *update_index[], int *update[],
                                            int  *data_org[], int *external[],
                                            int *extern_index[], int
                                            input_option, double *a[],
                                            int *bindx[], int *indx[],
                                            int *bpntr[], int *rpntr[],
                                            int *cpntr[]);

  extern void init_options(int options[], double params[]);

  extern void fill_fe_matrix(double val[],int bindx[], int update[],
                             int update_index[], int external[],
                             int extern_index[], int data_org[]);

  /* ----------------------- execution begins ------------------------------*/

  /* get number of processors and the name of this processor */
  MPI_Init(&argc,&argv);
  AZ_processor_info(proc_config);
  total_time = MPI_Wtime();

  /*
   * Read and broadcast: problem choice, problem size, equations per grid point
   * and how we wish to initialize 'update'.
   */
  application = 0;
  N_grid_pts = 64;
  num_PDE_eqns = 1;
  input_option = 0;
  if (proc_config[AZ_node] == 0) {
    while ((c = getopt(argc, argv, "n:d:?")) != EOF) {
      switch(c) {
      case 'n':
	i = sscanf(optarg, "%d", &N_grid_pts);
	if (i == 1) {
	  i = (int) (pow((double)N_grid_pts, 1.0/3.0)+0.5);
	  N_grid_pts = i*i*i;
	}
	else
	  N_grid_pts = 64;
	break;
      case 'd':
	i = sscanf(optarg, "%d", &num_PDE_eqns);
	if (i == 1) num_PDE_eqns = (num_PDE_eqns>0)?num_PDE_eqns:1;
	else num_PDE_eqns = 1;
	break;
      case '?':
	fprintf(stderr,
		"Usage:\n %s [-n num_grid_points] [-d variables_per_node]\n",
		argv[0]);
	exit(0);
	break;
      }
    }
  }
  AZ_broadcast((char *) &N_grid_pts  , sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &num_PDE_eqns, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &input_option, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &application , sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) NULL         , 0          , proc_config, AZ_SEND);

  /* create an application matrix for AZTEC */

  init_matrix_vector_structures(proc_config, &update_index, &update, &data_org,
                                &external, &extern_index, input_option, &val,
                                &bindx, &indx, &bpntr, &rpntr, &cpntr);

  /* initialize AZTEC options */

  init_options(options,params);
#ifdef DEBUG
  if ( (i = AZ_check_input(data_org, options, params, proc_config) ) < 0) {
    AZ_print_error(i);
    exit(-1);
  }
#endif
  /* Matrix fill for finite element example (see Aztec User's Guide). */

  if (application == 2)
    fill_fe_matrix(val, bindx, update, update_index, external, extern_index,
                   data_org);

  /* update[], update_index[], external[], extern_index[] are used to map
   * between Aztec's ordering of the equations and the user's ordering
   * (see the User's guide for more details). If these mapping arrays 
   * are not needed by the user, they can be deallocated as they are not 
   * used by AZ_solve().
   */

  free((void *) update);   free((void *) update_index);
  free((void *) external); free((void *) extern_index);

  plain_lanso(options, params, indx, bindx, rpntr, cpntr, bpntr, val,
	      data_org, status, proc_config);

  average_time = MPI_Wtime() - total_time;
  MPI_Allreduce(&average_time, &total_time, 1, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);

  if (proc_config[AZ_node] == 0) {
    average_time = total_time / proc_config[AZ_N_procs];
    printf("azl time: (%d PE) %f sec, (average) %f sec.\n",
	   proc_config[AZ_N_procs], total_time, average_time);
  }

  /* Free allocated memory */
  free((void *) indx);
  free((void *) bindx);    free((void *) rpntr);       free((void *) cpntr);
  free((void *) bpntr);    free((void *) val);         free((void *) data_org);

  MPI_Finalize();

}
