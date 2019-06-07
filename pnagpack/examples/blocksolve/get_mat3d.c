#include "pargrid.h"

/*+ get_mat - Generate a sparse matrix from the given grid

     Input Parameters:
     grid - the given grid
     ncomp - the number of components per grid point
     procinfo - the usual procinfo stuff

     Returns:
     the sparse matrix

 +*/

BSspmat *get_mat3d(par_grid *grid, BSprocinfo *procinfo)
{
  BSspmat *A;
  int	i, j, k, l, m;
  int	count, ncomp, nonsym, positive;
  point ***points = grid->points;
  int	*rp, *cval;
  FLOAT *aval, fncomp, off_diag, sign;

  ncomp = grid->ncomp;
  nonsym = (!grid->symmetric);
  positive = grid->positive;

  /* allocate space for the sparse matrix */
  rp = (int *) MALLOC((ncomp*grid->local_total+1)*sizeof(int));
  cval = (int *) MALLOC((7*ncomp*ncomp*grid->local_total)*sizeof(int));
  aval = (FLOAT *) MALLOC((7*ncomp*ncomp*grid->local_total)*sizeof(FLOAT));
  grid->rp = rp;
  grid->cval = cval;
  grid->aval = aval;

  /* ****************************************************** */
  /* now put values in the matrix                           */
  /* ****************************************************** */
  /* fncomp = (FLOAT) ncomp;
  off_diag = -(1.0/(6.0*fncomp)); */
  /* the following modification is made so that a trivial matrix will be generated - KW, Sept. 97 */
  fncomp = 6.0;
  off_diag = -1.0;
  count = 1;
  rp[0] = 0;
  sign = 1.0;
  srand48((long)(11311));
  for (i=1;i<grid->l_num_x+1;i++) {
    for (j=1;j<grid->l_num_y+1;j++) {
      for (k=1;k<grid->l_num_z+1;k++) {
	for (l=0;l<ncomp;l++) {
	  rp[count] = rp[count-1];

	  /* diagonal, set the column number, the nonzero values */
	  /* and increment the length of the row */
	  for (m=0;m<ncomp;m++) {
	    cval[rp[count]] = ncomp*points[i][j][k].num+m;
	    if(m==l) {
	      /* various diag possibilities for testing */
	      /* aval[rp[count]] = sign*(m+l+1)*fncomp;
	      aval[rp[count]] = sign*fncomp;
	      aval[rp[count]] = fncomp + drand48(); */
	      aval[rp[count]] = fncomp;
	      sign *= -1.0;
	    } else {
	      if(nonsym) {
		if(l<m)
		  aval[rp[count]] = -1.5;
		else
		  aval[rp[count]] = -0.5;
	      } else
		aval[rp[count]] = -1.0;
	    }
	    rp[count]++;
	  }

	  /* east */
	  if (points[i+1][j][k].num >= 0) {
	    for (m=0;m<ncomp;m++) {
	      cval[rp[count]] = ncomp*points[i+1][j][k].num+m;
	      if(nonsym)
		aval[rp[count]] = (1.0/(12.0*fncomp));
	      else
		aval[rp[count]] = off_diag;
	      rp[count]++;
	    }
	  }

	  /* west */
	  if (points[i-1][j][k].num >= 0) {
	    for (m=0;m<ncomp;m++) {
	      cval[rp[count]] = ncomp*points[i-1][j][k].num+m;
	      if(nonsym)
		aval[rp[count]] = -(1.0/(3.0*fncomp));
	      else
		aval[rp[count]] = off_diag;
	      rp[count]++;
	    }
	  }

	  /* north */
	  if (points[i][j+1][k].num >= 0) {
	    for (m=0;m<ncomp;m++) {
	      cval[rp[count]] = ncomp*points[i][j+1][k].num+m;
	      aval[rp[count]] = off_diag;
	      rp[count]++;
	    }
	  }

	  /* south */
	  if (points[i][j-1][k].num >= 0) {
	    for (m=0;m<ncomp;m++) {
	      cval[rp[count]] = ncomp*points[i][j-1][k].num+m;
	      aval[rp[count]] = off_diag;
	      rp[count]++;
	    }
	  }

	  /* up */
	  if (points[i][j][k+1].num >= 0) {
	    for (m=0;m<ncomp;m++) {
	      cval[rp[count]] = ncomp*points[i][j][k+1].num+m;
	      aval[rp[count]] = off_diag;
	      rp[count]++;
	    }
	  }

	  /* down */
	  if (points[i][j][k-1].num >= 0) {
	    for (m=0;m<ncomp;m++) {
	      cval[rp[count]] = ncomp*points[i][j][k-1].num+m;
	      aval[rp[count]] = off_diag;
	      rp[count]++;
	    }
	  }
	  /* take abs value if positive matrix requested */
	  if(positive) {
	    for (m=rp[count-1];m<rp[count];m++) {
	      aval[m]=fabs(aval[m]);
	    }
	  }
	  /* sort the values in this row */
	  BSheap_sort1d(rp[count]-rp[count-1],&(cval[rp[count-1]]),
			&(aval[rp[count-1]])); CHKERRN(0);
			count++;
	}
      }
    }
  }
  A = BSeasy_A(ncomp*grid->offset,ncomp*grid->local_total,
	       rp,cval,aval,procinfo); 
  CHKERRN(0);
  return(A);
}
