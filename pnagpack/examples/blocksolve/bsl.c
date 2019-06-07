/* bsl.c --
   Converted from master0.c and worker0.c of BlockSolve95.
   It construct a matrix as in grid0 example of BlockSolve95,
   then call pLANSO to find as many eigenvalues as possible.
   The limitation is the memory available, see MAX_MEM_SIZE in
   Makefile.
   */
#include "pargrid.h"

extern int plain_lanso(BSpar_mat *, BScomm *, BSprocinfo *);

/*+ worker - Solve a sparse eigenvalue problem associated with a grid

     Input Parameters:
     grid - the given grid
     procinfo - the processor information (in BlockSolve format)

 +*/

static void worker(par_grid *grid, BSprocinfo *procinfo)
{
  int	i, global_nnz, write_option = 0;
  BSspmat *A;
  BSpar_mat *pA;
  BScomm *Acomm;

  /* number grid to use in matrix assembly */
  num_grid3d(grid,procinfo);

  /* now call the routines to set up the matrix */
  A = get_mat3d(grid,procinfo);

  /* Set symmetry and storage scheme to be used */
  BSset_mat_symmetric(A,grid->symmetric);
  BSset_mat_icc_storage(A,grid->icc_storage);

	/* write out matrix */
  if(write_option) {
    write_mat_matlab("MAT.m",A,procinfo);
  }

  /* permute the matrix */
  pA = BSmain_perm(procinfo,A); CHKERR(0);

  /* count nnzs for display */
  global_nnz = 2*pA->local_nnz - pA->num_rows;
  GISUM(&global_nnz,1,&i,procinfo->procset);
  if(procinfo->my_id==0) {
    printf("o  ");
    printf("Number of nonzeros = %d\n",global_nnz);
  }

  Acomm = BSsetup_forward(pA,procinfo); CHKERR(0);
  /* solve for eigenvalues */
  i=plain_lanso(pA,Acomm,procinfo); CHKERR(0);

  /* free the grid */
  free_grid(grid);

  /* free the spmat */
  BSfree_easymat(A);

  /* free the par mat, etc. */
  BSfree_par_mat(pA);
  BSfree_comm(Acomm);
}

/*+ main - the main routine for this grid program

     Input Parameters:
     argc, argv - the usual
     argv[1] - the number of processors in the x direction
     argv[2] - the number of processors in the y direction
     argv[3] - the number of processors in the z direction
     argv[4] - the number of points in the x direction on each processor
               do not use less than 3 points
     argv[5] - the number of points in the y direction on each processor
               do not use less than 3 points
     argv[6] - the number of points in the z direction on each processor
               do not use less than 3 points

     To Run, see tools.

     Notes: The grid program solves a linear system associated with
            a 3-D grid distributed across the processors.  The 3-D
            grid is partitioned in all three dimensions amongst the
            processors.  A 7pt stencil is used.

 +*/

int main(int argc, char **argv)
{
  double total_time, average_time;
  par_grid grid;
  BSprocinfo *procinfo;

  /* Call BSinit() to initialize BlocklSolve and MPI */
  BSinit(&argc,&argv);
  total_time = MPI_Wtime();

  /* set up the context for BlockSolve */
  procinfo = BScreate_ctx(); CHKERRN(0);
  /* tell it that this matrix has no i-nodes or cliques */
  BSctx_set_si(procinfo,TRUE); CHKERRN(0);
  /* tell it to print out some information on the reordering */
  BSctx_set_pr(procinfo,TRUE); CHKERRN(0);
  BSctx_set_scaling(procinfo,FALSE); CHKERRN(0);
  BSctx_set_num_rhs(procinfo,1); CHKERRN(0);
  BSctx_set_pre(procinfo,PRE_ICC); CHKERRN(0);
  grid.icc_storage = TRUE;
  grid.symmetric = TRUE;
  grid.ncomp = 1;
  grid.positive = FALSE;

  if(procinfo->my_id==0) {
    printf("\n");
    printf("************** pLANSO example using Blocksolve **************\n");
  }

  /* read in grid parameters */
  if (argc < 7) {
    SETERRC(ARG_ERR,"Argument list too small\n");
    return 0;
  }
  sscanf(argv[1],"%d",&grid.worker_x);
  sscanf(argv[2],"%d",&grid.worker_y);
  sscanf(argv[3],"%d",&grid.worker_z);
  if (procinfo->my_id == 0) {
    printf("o  Number of workers (x,y,z): %d %d %d\n",
	   grid.worker_x,grid.worker_y,grid.worker_z);
  }
  if (procinfo->nprocs != grid.worker_x*grid.worker_y*grid.worker_z) {
    SETERRC(ARG_ERR,"Number of processors is not correct\n");
    return 0;
  }
  sscanf(argv[4],"%d",&grid.l_num_x);
  sscanf(argv[5],"%d",&grid.l_num_y);
  sscanf(argv[6],"%d",&grid.l_num_z);

  if (procinfo->my_id == 0) {
    printf("o  Local discretizations (x,y,z): %d %d %d\n",
	   grid.l_num_x,grid.l_num_y,grid.l_num_z);
  }

  /* call the worker */
  worker(&grid,procinfo); 

  average_time = MPI_Wtime() - total_time;
  MPI_Allreduce(&average_time, &total_time, 1, MPI_DOUBLE, MPI_SUM,
		procinfo->procset);

  if (procinfo->my_id==0) {
    average_time = total_time / procinfo->nprocs;
    printf("bsl time: (%d PE) %f sec, (average) %f sec.\n",
	   procinfo->nprocs, total_time, average_time);
    printf("************ End Blocksolve Example *****************\n");
    printf("\n");
  }

  /* print logging if enabled */
  BSprint_log(procinfo); CHKERRN(0);

  /* free the context */
  BSfree_ctx(procinfo); CHKERRN(0);

  /* finalize BlockSolve and MPI */
  BSfinalize();

  return(0);
}
