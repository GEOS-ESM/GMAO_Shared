#include "pargrid.h"

/*+ free_grid - Free the grid data structure

     Input Parameters:
     grid - the given grid

 +*/

void free_grid(par_grid *grid)
{
	int i, j;

	/* free the grid */
	for (i=0;i<grid->l_num_x+2;i++) {
		for (j=0;j<grid->l_num_y+2;j++) {
			FREE(grid->points[i][j]);
		}
		FREE(grid->points[i]);
	}
	FREE(grid->points);
	if(grid->rp != NULL) FREE(grid->rp);
	if(grid->cval != NULL) FREE(grid->cval);
	if(grid->aval != NULL) FREE(grid->aval);

}
