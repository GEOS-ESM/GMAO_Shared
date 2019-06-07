# include ../Make.inc for machine specific commands and flags
include ../Make.inc
.SUFFIXES: .f .o
#
OBJS=plandr2.o ddotmpi.o plandr.o pstartv.o planso.o ppurge.o pstpone.o
#
.DEFALUT:
	@$(ECHO) "Unknown target: $@"
	@$(ECHO) "This makefile has defined the following targets:"
	@$(ECHO) "lib clean"
#
lib: $(OBJS)
	$(AR) ruv ../libplan.a $(OBJS)
	$(RANLIB) ../libplan.a
#
clean:
	$(RM) -f *.o *~ core *.T *.lst
#
#  Command to build .o files from .f files.
#
.f.o:
	$(MPIFC) -c $(FFLAGS) $<
#
