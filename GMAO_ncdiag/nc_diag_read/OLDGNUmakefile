# Include build configuration:
include GNUmakefile.config

NC_DIAG_READ_OBJS+= \
    kinds.o \
    ncdr_climsg.o \
    ncdr_state.o ncdr_types.o \
    netcdf_unlimdims.o \
    ncdr_realloc_mod.o \
    ncdr_dims.o ncdr_check.o \
    ncdr_vars.o \
    ncdr_attrs.o \
    ncdr_global_attrs.o \
    ncdr_alloc_assert.o \
    ncdr_vars_fetch.o \
    ncdr_attrs_fetch.o \
    ncdr_global_attrs_fetch.o

NC_DIAG_READ_DEPS+=GNUmakefile \
    $(NC_DIAG_READ_OBJS)

.PHONY: all clean tests check

all: libnc_diag_read.a nc_diag_read.o tests

tests:
	make -C tests

check:
	make -C tests check

libnc_diag_read.a: nc_diag_read.o $(NC_DIAG_READ_DEPS)
	$(NCDR_MAKE_FANCY_BIG_BEGIN)
	$(NCDR_MAKE_QUIET)ar rcs libnc_diag_read.a nc_diag_read.o $(NC_DIAG_READ_OBJS)
	$(NCDR_MAKE_FANCY_BIG_END)

kinds.o: GNUmakefile GNUmakefile.config GNUmakefile.libconfig kinds.F90
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) -c kinds.F90 $(NCDR_CONFIG_FLAGS)
	$(NCDR_MAKE_FANCY_END)

nc_diag_read.o: $(NC_DIAG_READ_DEPS) nc_diag_read.F90
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c nc_diag_read.F90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_climsg.o: ncdr_climsg.F90
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_climsg.F90
	$(NCDR_MAKE_FANCY_END)

ncdr_types.o: ncdr_types.f90 kinds.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_types.f90
	$(NCDR_MAKE_FANCY_END)

ncdr_state.o: ncdr_state.f90 kinds.o ncdr_types.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_state.f90
	$(NCDR_MAKE_FANCY_END)

ncdr_vars.o: ncdr_vars.f90 kinds.o ncdr_state.o ncdr_climsg.o ncdr_check.o ncdr_alloc_assert.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_vars.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_vars_fetch.o: ncdr_vars_fetch.f90 kinds.o ncdr_state.o ncdr_check.o ncdr_alloc_assert.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_vars_fetch.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_attrs.o: ncdr_attrs.f90 kinds.o ncdr_state.o ncdr_climsg.o ncdr_check.o ncdr_alloc_assert.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_attrs.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_attrs_fetch.o: ncdr_attrs_fetch.f90 kinds.o ncdr_state.o ncdr_check.o ncdr_alloc_assert.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_attrs_fetch.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_global_attrs.o: ncdr_global_attrs.f90 kinds.o ncdr_state.o ncdr_climsg.o ncdr_check.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_global_attrs.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_global_attrs_fetch.o: ncdr_global_attrs_fetch.f90 kinds.o ncdr_state.o ncdr_check.o ncdr_alloc_assert.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_global_attrs_fetch.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_alloc_assert.o: ncdr_alloc_assert.f90 kinds.o ncdr_state.o ncdr_climsg.o ncdr_check.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_alloc_assert.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_dims.o: ncdr_dims.f90 kinds.o ncdr_state.o ncdr_climsg.o ncdr_check.o netcdf_unlimdims.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_dims.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_check.o: ncdr_check.f90 kinds.o ncdr_climsg.o ncdr_state.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_check.f90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

ncdr_realloc_mod.o: ncdr_realloc_mod.F90 kinds.o ncdr_types.o
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c ncdr_realloc_mod.F90
	$(NCDR_MAKE_FANCY_END)

netcdf_unlimdims.o: netcdf_unlimdims.F90 GNUmakefile
	$(NCDR_MAKE_FANCY_BEGIN)
	$(NCDR_MAKE_QUIET)$(NCDR_FORT90C) $(NCDR_F90FLAGS) $(NCDR_CONFIG_FLAGS) \
		-c netcdf_unlimdims.F90 \
		$(NCDR_NETCDF_FLAGS)
	$(NCDR_MAKE_FANCY_END)

clean:
	rm -f *.mod *.x *.o *.a
	make -C tests clean
