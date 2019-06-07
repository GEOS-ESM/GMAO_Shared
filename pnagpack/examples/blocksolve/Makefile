#
# To use:
# make PETSC_ARCH=ARCH BOPT=O
#
PETSC_ARCH=solaris
BOPT=O

all:  bsl

BS_DIR = ../../../BlockSolve95
LIBBASE =

#Compiler flags:
#MLOG turns on logging and flop counting
#DEBUG_ALL turns on error code checking and is set in BSdepend.h

CFLAGS   = $(BS_INCLUDE) $(MPI_INCLUDE) 
SYSLIB    = $(BS_LIB) $(LAPACK_LIB) $(BLAS_LIB) \
$(MPI_LIB) $(FC_LIB) $(SYS_LIB) -lm \
-Wl"-D equiv(DCOPY)=SCOPY;equiv(DDOT)=SDOT;equiv(DSCAL)=SSCAL;equiv(DAXPY)=SAXPY;equiv(IDAMAX)=ISAMAX;equiv(DSWAP)=SSWAP"
LINCLUDE  = $(SOURCEH)

bslan.o: bslan.c
	-$(CC) -c $(CFLAGS) $(BASEOPT) -DMSGLVL=20 -DMAX_MEM_SIZE=16731136 $<

bsl: bsl.o bslan.o get_mat3d.o write_mat_matlab.o num_grid3d.o \
	free_grid.o ppurge.o chkopts
	-$(CLINKER) -o $@ bsl.o bslan.o get_mat3d.o write_mat_matlab.o \
	num_grid3d.o free_grid.o ppurge.o -L../.. -lplan -llanso $(SYSLIB) 

include $(BS_DIR)/bmake/$(PETSC_ARCH)/$(PETSC_ARCH)
 
#clean:
#	/bin/rm -f *.o *~ core *.trace

