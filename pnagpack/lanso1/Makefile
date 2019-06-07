#
# Modified from the original form to produce a library file
#
include ../Make.inc
#
# System V machines might be using ${F77} and ${F77FLAGS} instead of
# ${FC} and ${FFLAGS} - at least the Ardent Titan is.
# Uncomment the next line if your "make" understands only ${F77}.
#FC=${F77}
#
# For use with Sun FORTRAN 1.3.1 or 1.2.
# Note that -dalign is needed whenever -O2 or -O3 is in effect
# to work around a code generation bug (cf. Sun Bug Id 1046392).
#
F77FLAGS=${FFLAGS}
LDFLAGS=
LIBS=

NOSUNOS4=./nosunos4
NODEBUG=./nodebug
FILES=	INSTALL README CHANGELOG loadev.f \
	doit std201.5 sni207.1 sni207.5 hs.5 diag100.5

#
# objects with embedded debugging code and SunOS4.x-dependent code
#
D4OBJ=	landr.o lanso.o purge.o ritvec.o startv.o
D4SRC=	landr.f lanso.f purge.f ritvec.f startv.f

#
# BLAS level 1
#
BLASOBJ=blas.o
BLASSRC=blas.f

#
# debug utility routines
#
DBGOBJ=	#debug.o
DBGSRC=	#debug.f

#
# eXtra supporting objects:
#
#	MACHAR is written by W. C. Cody,
#	TQL2 and PYTHAG are part of EISPACK.
#
# Note that TQLB needs PYTHAG also.
#
X_OBJ=	machar.o tql2.o pythag.o
X_SRC=	machar.f tql2.f pythag.f

#
# common objects
#
C_OBJ=	${X_OBJ} dsort2.o tqlb.o ortbnd.o random.o store.o stpone.o
C_SRC=	${X_SRC} dsort2.f tqlb.f ortbnd.f random.f store.f stpone.f

OBJS=	${BLASOBJ} ${DBGOBJ} ${C_OBJ} ${D4OBJ}

DIAG=	diag.o
STD=	std.o
SNI=	main.o op.o
HS=	decomp.o opsimon.o

C_MAIN=	op.f opsimon.f
D4MAIN=	diag.f std.f main.f decomp.f

C_SRCS=	${BLASSRC} ${DBGSRC} ${C_SRC} ${C_MAIN}
D4SRCS=	${D4SRC} ${D4MAIN}

ALL=	diag std sni hs
#
#  Default command.
#
.DEFAULT:
	@$(ECHO) "Unknown target $@"
	@$(ECHO) "Available targets are: lib, clean, all, diag, std, sni, hs"
#
# library liblanso.a
#
lib: ${C_OBJ} ${D4OBJ} datx.o
	$(AR) ruv ../liblanso.a ${C_OBJ} ${D4OBJ} datx.o
	${RANLIB} ../liblanso.a
#
all: ${ALL}

diag: ${DIAG} ${OBJS}
	${FC} -o diag ${LDFLAGS} ${FFLAGS} ${DIAG} ${OBJS} ${LIBS}

std: ${STD} ${OBJS}
	${FC} -o std ${LDFLAGS} ${FFLAGS} ${STD} ${OBJS} ${LIBS}

sni: ${SNI} ${OBJS}
	${FC} -o sni ${LDFLAGS} ${FFLAGS} ${SNI} ${OBJS} ${LIBS}

hs: ${HS} ${OBJS}
	${FC} -o hs ${LDFLAGS} ${FFLAGS} ${HS} ${OBJS} ${LIBS}

clean: FRC
	rm -f ${OBJS} ${STD} ${SNI} ${HS} ${DIAG} datx.o

realclean: clean
	rm -f ${ALL} errs

dist_s: FRC
	-@rm -rf ../dist_s
	-@mkdir ../dist_s
	-@for i in ${D4SRCS}; do\
		${NOSUNOS4} $$i | ${NODEBUG} > ../dist_s/$$i; done
	-@sed -e 's/\(debug.[fo]\)$$/#\1/' Makefile > ../dist_s/Makefile
	-@cp -p ${FILES} ${C_SRCS} ../dist_s
	-@rm -f ../dist_s/#debug.f

dist_d: FRC
	-@rm -rf ../dist_d
	-@mkdir ../dist_d
	-@for i in ${D4SRCS}; do\
		${NOSUNOS4} $$i > ../dist_d/$$i; done
	-@cp -p Makefile nodebug ${FILES} ${C_SRCS} ../dist_d

dist: dist_s dist_d

ilpdist: FRC
	-@rm -rf ../ilpdist/lanso.1.0
	-@mkdir ../ilpdist/lanso.1.0
	-@for i in ${D4SRCS}; do\
		${NOSUNOS4} $$i | ${NODEBUG} > ../ilpdist/lanso.1.0/$$i; done
	-@sed -e 's/\(debug.[fo]\)$$/#\1/' Makefile > ../ilpdist/lanso.1.0/Makefile
	-@cp -p ${FILES} ${C_SRCS} ../ilpdist/lanso.1.0
	-@rm -f ../ilpdist/lanso.1.0/#debug.f

eispac: eispac.o
	${FC} -o eispac ${LDFLAGS} ${FFLAGS} eispac.o -leispackd-g

FRC:
#
machar.o: machar.f
	${FC} ${NOOPT} -c machar.f
.f.o:
	${FC} ${FFLAGS} -c $<
#
