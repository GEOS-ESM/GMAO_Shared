# Include build configuration:
include GNUmakefile.config

FSON_DEPS+=nc_diag_fson/GNUmakefile \
    nc_diag_fson/src/nc_diag_fson.f90 \
    nc_diag_fson/src/ncdf_path_m.f90 \
    nc_diag_fson/src/ncdf_string_m.f90 \
    nc_diag_fson/src/ncdf_value_m.f90

FSON_OBJS+= \
    nc_diag_fson/build/nc_diag_fson.o \
    nc_diag_fson/build/ncdf_path_m.o \
    nc_diag_fson/build/ncdf_string_m.o \
    nc_diag_fson/build/ncdf_value_m.o

NC_DIAG_RES_DEPS+=GNUmakefile \
    ncdres_climsg.o

NC_DIAG_RES_OBJS+= \
    $(NC_DIAG_RES_DEPS) \
    nc_diag_res.o

.PHONY: all clean tests check

all: libnc_diag_res.a nc_diag_res.o tests

tests:
	make -C tests

check:
	make -C tests check

libnc_diag_res.a: $(NC_DIAG_RES_OBJS)
	$(NCDRES_MAKE_FANCY_BIG_BEGIN)
	$(NCDRES_MAKE_QUIET)ar rcs libnc_diag_res.a $(NC_DIAG_RES_OBJS) $(FSON_OBJS)
	$(NCDRES_MAKE_FANCY_BIG_END)

nc_diag_res.o: nc_diag_res.f90 nc_diag_fson/lib/libnc_diag_fson.a $(NC_DIAG_RES_DEPS)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) $(NCDRES_CONFIG_FLAGS) \
		-c nc_diag_res.f90 \
		-Inc_diag_fson/lib -Lnc_diag_fson/lib -lnc_diag_fson
	$(NCDRES_MAKE_FANCY_END)

ncdres_climsg.o: ncdres_climsg.F90 GNUmakefile.config GNUmakefile.libconfig
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) $(NCDRES_CONFIG_FLAGS) \
		-c ncdres_climsg.F90
	$(NCDRES_MAKE_FANCY_END)

nc_diag_fson/lib/libnc_diag_fson.a: $(FSON_DEPS)
	make -C nc_diag_fson lib

clean:
	rm -f *.mod *.x *.o *.a
	make -C nc_diag_fson clean clobber
	make -C tests clean
