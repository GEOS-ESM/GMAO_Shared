# include top level Make.inc file for machine specific options
include ../../Make.inc
#
# definitions for SuperLU
SLUINC=../../../SuperLU/SRC
SLULIB=../../../SuperLU/superlu_solaris.a ../../../blas_solaris.a
#
LIBALL=-L../.. -lplan -llanso $(LIBS)
COMMON=hbread1.o psp_setup.o dmem.o getarg.o
SIMPHB=simphb.o shbstore.o shbpurge.o ${COMMON}
SHBHI=shbhi.o shbstore.o shbpurge.o hbread3.o ${COMMON}
SHBCG=shbcg.o hbread3.o slustub.o ${COMMON}
GRID=grid.o shbstore.o shbpurge.o hbread3.o ${COMMON}
SHBOOC=shbooc.o dskstore.o ${COMMON}
SHBSTD=shbstd.o shbstore.o shbpurge.o hbread2.o ${COMMON}
STDHI=stdhi.o shbstore.o shbpurge.o hbread2.o ${COMMON}
DYNSRC=simphb.o shbstore.o shbpurge.o hbread1.o psp_setup.o grid.o shbcg.o \
	shbooc.o dskstore.o shbstd.o hbread2.o shbhi.o stdhi.o hbread3.o
ALL=simphb shbhi shbcg grid shbooc shbstd stdhi
#
.DEFAULT:
	@$(ECHO) 'Unknown target $@.'
	@$(ECHO) 'The following targets are defined in this Makefile.'
	@$(ECHO) '${ALL} ooc std clean'
#
all: ${ALL}
#
simphb: $(SIMPHB)
	$(MPIFC) -o simphb $(SIMPHB) $(LIBALL) $(LDFLAGS)
#
shbhi: $(SHBHI)
	$(MPIFC) -o shbhi $(SHBHI) $(LIBALL) $(LDFLAGS)
#
shbcg: $(SHBCG)
	$(MPIFC) -o shbcg $(SHBCG) $(SLULIB) $(LIBALL) $(LDFLAGS)
#
grid: $(GRID)
	$(MPIFC) -o grid $(GRID) $(LIBALL) $(LDFLAGS)
#
ooc shbooc: $(SHBOOC)
	$(MPIFC) -o shbooc $(SHBOOC) $(LIBALL) $(LDFLAGS)
#
std shbstd: $(SHBSTD)
	$(MPIFC) -o shbstd $(SHBSTD) $(LIBALL) $(LDFLAGS)
#
stdhi: $(STDHI)
	$(MPIFC) -o stdhi $(STDHI) $(LIBALL) $(LDFLAGS)
#
clean:
	$(RM) -f *.o *trace core *~ ${ALL} _LANVEC_* *.T *.lst fort.*
#
# compile all files that uses dynamic memory allocation with
# lower optimization level
#
$(DYNSRC):
	$(MPIFC) $(NOOPT) -c $<
.f.o:
	$(MPIFC) $(FFLAGS) -c $<
#
simphb.o shbstore.o shbpurge.o shbooc.o dskstore.o shbstd.o shbhi.o stdhi.o grid.o: simphb.h
shbcg.o: shbcg.h
#
# special target dmem.o, getarg.o and getarg.o
# on a CRAY machine
#dmem.o: dmem.F
#	$(FC) $(NOOPT) -DCRAY -c dmem.F
#getarg.o: getarg.F
#	$(FC) $(FFLAGS) -DCRAY -c getarg.F
#slustub.o: slustub.c
#	$(CC) -DUpCase $(CFLAGS) -I ${SLUINC} -c $<
#
# on Most workstations
#
dmem.o: dmem.F
	$(FC) $(NOOPT) -DSUN -c dmem.F
slustub.o: slustub.c
	$(CC) -DAdd_ $(CFLAGS) -I ${SLUINC} -c $<
#
# In some special cases where there is no fortran interface to dynamic
# memory allocation routine malloc, the following line is for SGI r10000
#
#dmem.o: dmem.c
#	cc $(NOOPT) -D_SGI -c dmem.c
#
#
