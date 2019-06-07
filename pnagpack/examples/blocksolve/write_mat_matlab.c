#include	"pargrid.h"
#include	"stdio.h"

void	write_mat_matlab(char *str, BSspmat *A, BSprocinfo  *procinfo)
{
	FILE	*fp;
	int	i, j, row, count = 1, nnz;

	/* try to open the file */
	if ((fp = fopen(str,"w")) == NULL) {
		fprintf(stderr,"Cannot open %s\n",str);
		exit(-1);
	}

	nnz = 0;
	for (i=0;i<A->num_rows;i++) {
		for (j=0;j<A->rows[i]->length;j++) {
			nnz++;
		}
	}

	/* print out the header */
	fprintf(fp,"n = %d;number_nz=%d;t=zeros(number_nz,3);\n",A->num_rows,nnz);

	/* print out the triplets */
	for (i=0;i<A->num_rows;i++) {
		A->map->flocal2global(1,&i,&row,procinfo,A->map);
		for (j=0;j<A->rows[i]->length;j++) {
			fprintf(fp,"t(%d,1:3) = [%d %d %4.16e];\n",count,row+1,A->rows[i]->col[j]+1,A->rows[i]->nz[j]);
			count++;
		}
	}
	fprintf(fp,"a = sparse(t(:,1),t(:,2),t(:,3));\n");

	fclose(fp);
}

void	write_vec_matlab(char *str,FLOAT *x,BSspmat *A,BSprocinfo *procinfo)
{
	FILE	*fp;
	int	i, j, row, cnt, num_rhs;

	/* try to open the file */
	if ((fp = fopen(str,"w")) == NULL) {
		fprintf(stderr,"Cannot open %s\n",str);
		exit(-1);
	}

	num_rhs = procinfo->num_rhs;

	/* print out the header */
	fprintf(fp,"v=zeros(%d,%d);\n",A->num_rows,num_rhs);

	/* print out the vectors */
	cnt = 0;
	for (j=0;j<num_rhs;j++) {
		for (i=0;i<A->num_rows;i++) {
			A->map->flocal2global(1,&i,&row,procinfo,A->map);
			fprintf(fp,"v(%d,%d) = %4.16e;\n",row+1,j+1,x[cnt]);
			cnt++;
		}
	}

	fclose(fp);
}

void	write_vec_matlab2(char *str,FLOAT *x,BSpar_mat *A,BSprocinfo *procinfo)
{
	FILE	*fp;
	int	i, j, row, cnt, num_rhs;

	/* try to open the file */
	if ((fp = fopen(str,"w")) == NULL) {
		fprintf(stderr,"Cannot open %s\n",str);
		exit(-1);
	}

	num_rhs = procinfo->num_rhs;

	/* print out the header */
	fprintf(fp,"v=zeros(%d,%d);\n",A->num_rows,num_rhs);
	printf("here 1\n");

	/* print out the vectors */
	cnt = 0;
	for (j=0;j<num_rhs;j++) {
		for (i=0;i<A->num_rows;i++) {
	printf("i=%d\n",i);
			A->map->flocal2global(1,&i,&row,procinfo,A->map);
			fprintf(fp,"v(%d,%d) = %4.16e;\n",row+1,j+1,x[cnt]);
			cnt++;
		}
	}

	fclose(fp);
}

