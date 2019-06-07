SRC_C    = azl.o azlan.o az_examples.o az_app_utils.o
SRC_F    = ppurge.o

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# Need to modify the following flags and commands according
# to the rules of Aztec
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

AZTECDIR        = ../../../Aztec
CC              = mpicc
FC              = mpif77
CFORT           = -Dappend_
#CFORT           = -Dcaps
CFLAGS2         = -O
FFLAGS2         = -O
#LDFLAGS2        = -Wl"-D equiv(DCOPY)=SCOPY;equiv(DDOT)=SDOT;equiv(DSQRT)=SQRT;equiv(DNRM2)=SNRM2;equiv(DSCAL)=SSCAL;equiv(DAXPY)=SAXPY;equiv(DASUM)=SASUM;equiv(IDAMAX)=ISAMAX;equiv(DSWAP)=SSWAP;equiv(DGEMV)=SGEMV;equiv(DGEMM)=SGEMM;equiv(DTRSM)=STRSM;equiv(DGER)=SGER;equiv(DTRMM)=STRMM;equiv(DTRMV)=STRMV;equiv(DSYRK)=SSYRK"
LIB_MACH_DEP    = -lF77 -lV77 -lM77 -lsunmath
MAIN_COMPILER   = $(CC)


MAX_MEM = 16731136
GEN_CFLAGS   = -I${AZTECDIR}/lib -DMAX_MEM_SIZE=$(MAX_MEM) -DMSGLVL=20
GEN_LIBS     = -L../.. -lplan -llanso -L${AZTECDIR}/lib -laztec  -lm

LDFLAGS      = $(LDFLAGS2) $(GEN_LIBS) $(LIB_MACH_DEP)
CFLAGS       = $(CFLAGS2)  $(CFORT)    $(GEN_CFLAGS)
FFLAGS       = $(FFLAGS2)  $(CFORT)    $(GEN_CFLAGS)

OBJS   = $(SRC_C) $(SRC_F)

#####################################
# Define implicit compilation rules #
#####################################
all: azl

%.o:%.F
	$(FC) $(FFLAGS) -c $<
%.o:%.c
	$(CC) $(CFLAGS) -c $<

# application code

azl: $(OBJS)
	$(MAIN_COMPILER) $(OBJS) $(LDFLAGS) -o azl
clean:
	@echo "cleaning ..."
	/bin/rm -f *.o core *~ azl
