/* filename: dmem.c
 
   This is an equivalent of dmem.F, please see the comments in dmem.F
   for more detail infromation.

   Please notice the naming convertion used here. Most unix supports
   this type mix use of Fortran and C code.
*/

#include <stdio.h>
#include <malloc.h>
/* define names of the timer functions */
#ifdef _IBM
#define dmalloc dmalloc
#define r8alloc r8alloc
#define dmfree dmfree
#else
#ifdef _SGI
#define dmalloc dmalloc_
#define r8alloc r8alloc_
#define dmfree dmfree_
#else
#ifdef _CRAY
#define dmalloc DMALLOC 
#define r8alloc R8ALLOC
#define dmfree DMFREE
#else
#define dmalloc dmalloc
#define r8alloc r8alloc
#define dmfree dmfree
#endif
#endif
#endif

#ifdef _CRAY
void dmalloc(ptr, nelm, type)
void **ptr; unsigned int *nelm; int *type;
{
    int nsize;
     
     nsize=sizeof(ptr)* *nelm;
    *ptr = (void *) malloc(nsize);
    if (*ptr == NULL) {
        fprintf(stderr,
	"dmalloc: malloc failed to allocate %d bytes from memory.\n",
           nsize);
	*ptr = 0;
    }
}

void r8alloc(ptr, nelm)
void **ptr; unsigned int *nelm;
{
    *ptr = (void *) malloc(8 * *nelm);
    if (*ptr == NULL) {
        fprintf(stderr,
	"r8alloc: malloc failed to allocate %d 8-byte words from heap.\n",
		*nelm);
	*ptr = 0;
    }
}

void dmfree(ptr)
void **ptr;
{
    free((char *) *ptr);
}
#else
void dmalloc(ptr, nelm, type)
void **ptr; unsigned int *nelm; int *type;
{
    /* this size information only correct for SUN */
    static int size[8] = {4,4,4,8,8,1,16,8};

    *ptr = (void *) malloc(size[*type - 1] * *nelm);
    if (*ptr == NULL) {
        fprintf(stderr,
        "dmalloc: malloc failed to allocate %d bytes from memory.\n",
                size[*type - 1] * *nelm);
        *ptr = 0;
    }
}

void r8alloc(ptr, nelm)
void **ptr; unsigned int *nelm;
{
    *ptr = (void *) malloc(8 * *nelm);
    if (*ptr == NULL) {
        fprintf(stderr,
        "r8alloc: malloc failed to allocate %d 8-byte words from heap.\n",
                *nelm);
	*ptr = 0;
    }
}
void dmfree(ptr)
void **ptr;
{
    free((char *) *ptr);
}
#endif

