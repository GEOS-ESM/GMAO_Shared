/* get the structure and function definitions for BlockSolve */
#include "BSprivate.h"

/* set up max stencil sizes, etc */
#define MAX_LEN 27

/* point on the grid */
typedef struct __point {
	int	num;
	int	type;
} point;

/* grid structure for each worker */
typedef struct __par_grid {
	int	worker_x;
	int	worker_y;
	int	worker_z;
	int	num_x;
	int	num_y;
	int	num_z;
	int	l_num_x;
	int	l_num_y;
	int	l_num_z;
	int	split_line;
	int	type;
	int	local_total;
	int	global_total;
	int	offset;
	int	ncomp;
	int	symmetric;
	int	icc_storage;
	int	positive;
	int	bmatrix;
	FLOAT shift;
	point	***points;
	int *rp;
	int *cval;
	FLOAT *aval;
} par_grid;

/* subroutine types */
/*
extern int sscanf(const char *,...);
extern int fprintf(FILE *,...);
*/
extern int fclose(FILE *);
extern double drand48();
extern void srand48(long);
extern void write_mat_matlab(char *, BSspmat *, BSprocinfo *);
extern void write_vec_matlab(char *,FLOAT *,BSspmat *,BSprocinfo *);
extern void write_vec_matlab2(char *,FLOAT *,BSpar_mat *,BSprocinfo *);
extern void	worker(par_grid *, BSprocinfo *);
extern void num_grid3d(par_grid *,BSprocinfo *);
extern void free_grid(par_grid *);
extern BSspmat *get_mat(par_grid *, BSprocinfo *);
extern BSspmat *get_mat3d(par_grid *,BSprocinfo *);
extern void ex3_freel2g(int *);
extern void ex3_loc2glob(int ,int *,int *, BSprocinfo *, BSmapping *);
extern void ex3_glob2proc(int , int *, int *, BSprocinfo *, BSmapping *);
extern void ex3_freeg2l(int *);
extern void ex3_glob2loc(int , int *, int *, BSprocinfo *, BSmapping *);

#define TRUE 1
#define FALSE 0

/* error defns */
#define ARG_ERR -100

/* message types */
#define EAST_MSG 50
#define WEST_MSG 51
#define NORTH_MSG 52
#define SOUTH_MSG 53
#define UP_MSG 54
#define DOWN_MSG 55

/* macros */
#define Mxpos(grid,procinfo) (procinfo->my_id % grid->worker_x)
#define Mypos(grid,procinfo) ((procinfo->my_id / grid->worker_x) % grid->worker_y)
#define Mzpos(grid,procinfo) (procinfo->my_id / (grid->worker_x*grid->worker_y))
#define Meast(grid,procinfo) (Mzpos(grid,procinfo)*(grid->worker_x*grid->worker_y) \
+ Mypos(grid,procinfo)*grid->worker_x \
+ (procinfo->my_id+1) % grid->worker_x)
#define Mwest(grid,procinfo) (Mzpos(grid,procinfo)*(grid->worker_x*grid->worker_y) \
+ Mypos(grid,procinfo)*grid->worker_x \
+ (procinfo->my_id+grid->worker_x-1) % grid->worker_x)
#define Mnorth(grid,procinfo) (Mzpos(grid,procinfo)*(grid->worker_x*grid->worker_y) \
+ ((Mypos(grid,procinfo)+1) % grid->worker_y)*grid->worker_x \
+ Mxpos(grid,procinfo))
#define Msouth(grid,procinfo) (Mzpos(grid,procinfo)*(grid->worker_x*grid->worker_y) \
+ ((Mypos(grid,procinfo)+grid->worker_y-1) % grid->worker_y)*grid->worker_x \
+ Mxpos(grid,procinfo))
#define Mup(grid,procinfo) (((Mzpos(grid,procinfo)+1) % grid->worker_z)*(grid->worker_x*grid->worker_y) \
+ Mypos(grid,procinfo)*grid->worker_x \
+ Mxpos(grid,procinfo))
#define Mdown(grid,procinfo) (((Mzpos(grid,procinfo)+grid->worker_z-1) % grid->worker_z)*(grid->worker_x*grid->worker_y) \
+ Mypos(grid,procinfo)*grid->worker_x \
+ Mxpos(grid,procinfo))

#define Msend_border_msg(msg_list,points,msg,msg_type,to,x1,x2,y1,y2,z1,z2,pi) \
{ \
	int	count99, i99, j99, k99; \
	count99 = 0; \
	for (i99=x1;i99<=x2;i99++) { \
		for (j99=y1;j99<=y2;j99++) { \
			for (k99=z1;k99<=z2;k99++) { \
				msg[count99].num = points[i99][j99][k99].num; \
				msg[count99].type = points[i99][j99][k99].type; \
				count99++; \
			} \
		} \
	} \
	MY_SEND_SYNC(msg_list,msg_type,msg,(count99*sizeof(point)/sizeof(int)), \
		(to),MPI_INT,pi); \
}

#define Mrecv_border_msg(points,intype,x1,x2,y1,y2,z1,z2,procinfo) \
{ \
	int	count99, i99, j99, k99, in_msg_size99; \
	point	*in_msg99; \
	MPI_Status	stat99; \
	in_msg_size99 = sizeof(point)*((x2)-(x1)+1)* \
		((y2)-(y1)+1)*((z2)-(z1)+1); \
	in_msg99 = (point *) MALLOC(in_msg_size99); \
	RECVSYNCNOMEM(intype,in_msg99,(in_msg_size99/sizeof(int)),MPI_INT, \
		procinfo,stat99); \
	count99 = 0; \
	for (i99=x1;i99<=x2;i99++) { \
		for (j99=y1;j99<=y2;j99++) { \
			for (k99=z1;k99<=z2;k99++) { \
				points[i99][j99][k99].num = in_msg99[count99].num; \
				points[i99][j99][k99].type = in_msg99[count99].type; \
				count99++; \
			} \
		} \
	} \
	FREE(in_msg99); \
}

void ex3_freel2g();
void ex3_loc2glob();
void ex3_glob2proc();
void ex3_freeg2l();
void ex3_glob2loc();
