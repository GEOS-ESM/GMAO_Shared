#include "pargrid.h"

/*+ num_grid - Number grid points on a 3d mesh

     Input Parameters:
     grid - the given grid
     procinfo - the processor information (in BlockSolve format)

 +*/

void num_grid3d(par_grid *grid, BSprocinfo *procinfo)
{
	int	i, j, k, ncomp;
	int	*offset;
	point *msg;
	int	count, max_msg_size;
	BSmsg_list	*msg_list;

	msg_list = NULL;

	/* the number of components per grid point */
	ncomp = grid->ncomp;

	/* find the number of local grid points */
	grid->local_total = grid->l_num_x*grid->l_num_y*grid->l_num_z;

	/* determine the beginning number of my grid points in the global */
	/* grid point numbering */
	BSoffset(1,&(grid->local_total),&offset,procinfo); CHKERR(0);
	grid->offset = (*offset);
	FREE(offset);

	/* determine the number of global grid points */
	grid->global_total = grid->local_total;
	GISUM(&(grid->global_total),1,&i,procinfo->procset);

	/* determine the maximum number of grid points on a single processor */
	j = grid->local_total;
	GIMAX(&j,1,&i,procinfo->procset);

	/* print out the number of unknows and the max unknowns */
	/* on a single processor */
	if (procinfo->my_id == 0) {
		printf("o  ");
		printf("Global total unknowns = %d, Local max unknowns = %d\n",
			ncomp*grid->global_total,ncomp*j);
	}

	/******************************************************************/
	/* set up grid with ghost points */
	/******************************************************************/
	grid->points = (point ***) MALLOC(sizeof(point **)*(grid->l_num_x+2));
	for (i=0;i<grid->l_num_x+2;i++) {
		grid->points[i] = (point **) MALLOC(sizeof(point *)*(grid->l_num_y+2));
		for (j=0;j<grid->l_num_y+2;j++) {
			grid->points[i][j] = (point *) 
				MALLOC(sizeof(point)*(grid->l_num_z+2));
			for (k=0;k<grid->l_num_z+2;k++) {
				grid->points[i][j][k].num = -1;
				grid->points[i][j][k].type = -1;
			}
		}
	}

	/* number local part of grid */
	count = 0;
	for (i=1;i<grid->l_num_x+1;i++) {
		for (j=1;j<grid->l_num_y+1;j++) {
			for (k=1;k<grid->l_num_z+1;k++) {
				grid->points[i][j][k].num = count + grid->offset;
				grid->points[i][j][k].type = grid->type;
				count++;
			}
		}
	}

	/******************************************************************/
	/* exchange edge information with other processors                 */
	/******************************************************************/

	/* allocate a message */
	max_msg_size = grid->l_num_x;
	if (max_msg_size < grid->l_num_y) {
		max_msg_size =  grid->l_num_y;
	}
	if (max_msg_size < grid->l_num_z) {
		max_msg_size =  grid->l_num_z;
	}
	max_msg_size += 2;
	max_msg_size *= (max_msg_size*sizeof(point));
	msg = (point *) MALLOC(max_msg_size);

	/* now send my east grid edge to the east */
	if (Mxpos(grid,procinfo) != grid->worker_x-1) {
		Msend_border_msg(msg_list,grid->points,msg,EAST_MSG,
			Meast(grid,procinfo),grid->l_num_x,grid->l_num_x,0,grid->l_num_y+1,
			0,grid->l_num_z+1,procinfo);
	}

	/* receive from the west */
	if (Mxpos(grid,procinfo) != 0) {
		Mrecv_border_msg(grid->points,EAST_MSG,0,0,0,grid->l_num_y+1,
			0,grid->l_num_z+1,procinfo);
	}

	/* now send my west grid edge to the west */
	if (Mxpos(grid,procinfo) != 0) {
		Msend_border_msg(msg_list,grid->points,msg,WEST_MSG,
			Mwest(grid,procinfo),1,1,0,grid->l_num_y+1,
			0,grid->l_num_z+1,procinfo);
	}

	/* receive from the east */
	if (Mxpos(grid,procinfo) != grid->worker_x-1) {
		Mrecv_border_msg(grid->points,WEST_MSG,grid->l_num_x+1,grid->l_num_x+1,
			0,grid->l_num_y+1,0,grid->l_num_z+1,procinfo);
	}

	/* now send my north grid edge to the north */
	if (Mypos(grid,procinfo) != grid->worker_y-1) {
		Msend_border_msg(msg_list,grid->points,msg,NORTH_MSG,
			Mnorth(grid,procinfo),0,grid->l_num_x+1,grid->l_num_y,grid->l_num_y,
			0,grid->l_num_z+1,procinfo);
	}

	/* receive from the south */
	if (Mypos(grid,procinfo) != 0) {
		Mrecv_border_msg(grid->points,NORTH_MSG,0,grid->l_num_x+1,0,0,
			0,grid->l_num_z+1,procinfo);
	}

	/* now send my south grid edge to the south */
	if (Mypos(grid,procinfo) != 0) {
		Msend_border_msg(msg_list,grid->points,msg,SOUTH_MSG,
			Msouth(grid,procinfo),0,grid->l_num_x+1,1,1,
			0,grid->l_num_z+1,procinfo);
	}

	/* receive from the north */
	if (Mypos(grid,procinfo) != grid->worker_y-1) {
		Mrecv_border_msg(grid->points,SOUTH_MSG,0,grid->l_num_x+1,
		grid->l_num_y+1,grid->l_num_y+1,0,grid->l_num_z+1,procinfo);
	}

	/* now send my upper grid edge up */
	if (Mzpos(grid,procinfo) != grid->worker_z-1) {
		Msend_border_msg(msg_list,grid->points,msg,UP_MSG,Mup(grid,procinfo),
			0,grid->l_num_x+1,0,grid->l_num_y+1,grid->l_num_z,grid->l_num_z,
			procinfo);
	}

	/* receive from below */
	if (Mzpos(grid,procinfo) != 0) {
		Mrecv_border_msg(grid->points,UP_MSG,0,grid->l_num_x+1,
			0,grid->l_num_y+1,0,0,procinfo);
	}

	/* now send my lower grid edge down */
	if (Mzpos(grid,procinfo) != 0) {
		Msend_border_msg(msg_list,grid->points,msg,DOWN_MSG,
			Mdown(grid,procinfo),0,grid->l_num_x+1,0,grid->l_num_y+1,
			1,1,procinfo);
	}

	/* receive from above */
	if (Mzpos(grid,procinfo) != grid->worker_z-1) {
		Mrecv_border_msg(grid->points,DOWN_MSG,0,grid->l_num_x+1,
			0,grid->l_num_y+1,grid->l_num_z+1,grid->l_num_z+1,procinfo);
	}

	FINISH_SEND_LIST(msg_list);

	FREE(msg);

	/* check to make sure that the local grid is numbered */
	for (i=1;i<grid->l_num_x+1;i++) {
		for (j=1;j<grid->l_num_y+1;j++) {
			for (k=1;k<grid->l_num_z+1;k++) {
				if (grid->points[i][j][k].num == -1) {
					printf("Proc %d: bad point %d %d %d\n",
						procinfo->my_id,i,j,k);
				}
			}
		}
	}
}
