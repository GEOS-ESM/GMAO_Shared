# include machine specific definitions from top level Make.inc
include ../../Make.inc
#
# For Cray T3E
#
OBJS=diag.o convck.o
LIBALL=../../libplan.a ../../liblanso.a $(LIBS)
#
diag: ${OBJS}
	${MPIFC} -o $@ ${OBJS} ${LIBALL} ${LDFLAGS}
#
clean:
	rm -f *.o *.trace core *~ fort.* diag *.T *.lst
#
diag.o: diag.f
	${MPIFC} ${FFLAGS} -c $<
convck.o: convck.f
	${FC} ${FFLAGS} -c $<
#
#  Default command.
#
.DEFAULT:
	@$(ECHO) "Unknown target $@"
	@$(ECHO) "Available targets are: diag, clean"
#
