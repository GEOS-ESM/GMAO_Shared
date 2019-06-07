# Include build configuration
include GNUmakefile.config

NC_DIAG_CAT_OBJS+=\
    ncdc_realloc.o ncdc_state.o ncdc_climsg.o ncdc_types.o \
    ncdc_cli_process.o \
    ncdc_dims.o ncdc_vars.o \
    ncdc_metadata.o \
    netcdf_unlimdims.o

NC_DIAG_CAT_DEPS+=Makefile \
    $(NC_DIAG_CAT_OBJS)

all: nc_diag_cat.x test_nc_unlimdims.x

kinds.o: kinds.F90 Makefile Makefile.config
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) $(NCDC_CONFIG_FLAGS) -c kinds.F90
	$(NCDC_MAKE_FANCY_END)

ncdc_types.o: ncdc_types.f90 kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_types.f90
	$(NCDC_MAKE_FANCY_END)

ncdc_climsg.o: ncdc_climsg.F90 ncdc_state.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_climsg.F90 $(NCDC_NETCDF_FLAGS)
	$(NCDC_MAKE_FANCY_END)

ncdc_cli_process.o: ncdc_cli_process.F90 ncdc_state.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_cli_process.F90
	$(NCDC_MAKE_FANCY_END)

ncdc_realloc.o: ncdc_realloc.F90 ncdc_types.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_realloc.F90
	$(NCDC_MAKE_FANCY_END)

ncdc_state.o: ncdc_state.F90 ncdc_types.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_state.F90
	$(NCDC_MAKE_FANCY_END)

ncdc_dims.o: ncdc_dims.F90 ncdc_state.o ncdc_realloc.o ncdc_climsg.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_dims.F90
	$(NCDC_MAKE_FANCY_END)

ncdc_vars.o: ncdc_vars.F90 ncdc_state.o ncdc_dims.o ncdc_realloc.o ncdc_climsg.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_vars.F90 $(NCDC_NETCDF_FLAGS)
	$(NCDC_MAKE_FANCY_END)

ncdc_metadata.o: ncdc_metadata.F90 ncdc_state.o ncdc_dims.o ncdc_vars.o ncdc_types.o ncdc_climsg.o ncdc_cli_process.o netcdf_unlimdims.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_metadata.F90 $(NCDC_NETCDF_FLAGS)
	$(NCDC_MAKE_FANCY_END)

ncdc_data.o: ncdc_data.F90 ncdc_state.o ncdc_climsg.o ncdc_cli_process.o ncdc_dims.o ncdc_vars.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_data.F90 $(NCDC_NETCDF_FLAGS)
	$(NCDC_MAKE_FANCY_END)

ncdc_data_MPI.o: ncdc_data_MPI.F90 ncdc_state.o ncdc_dims.o ncdc_vars.o ncdc_climsg.o ncdc_cli_process.o kinds.o
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c ncdc_data_MPI.F90 $(NCDC_NETCDF_FLAGS)
	$(NCDC_MAKE_FANCY_END)

nc_diag_cat.x: nc_diag_cat.F90 $(NC_DIAG_CAT_DEPS)
	$(NCDC_MAKE_FANCY_BIG_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) nc_diag_cat.F90 $(NC_DIAG_CAT_OBJS) $(NCDC_NETCDF_FLAGS) -o nc_diag_cat.x
	$(NCDC_MAKE_FANCY_BIG_END)

test_nc_unlimdims.x: test_nc_unlimdims.F90 netcdf_unlimdims.o ncdc_state.o ncdc_climsg.o
	$(NCDC_MAKE_FANCY_BIG_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) test_nc_unlimdims.F90 ncdc_climsg.o ncdc_state.o netcdf_unlimdims.o $(NCDC_NETCDF_FLAGS) -o test_nc_unlimdims.x
	$(NCDC_MAKE_FANCY_BIG_END)

netcdf_unlimdims.o: netcdf_unlimdims.F90 Makefile
	$(NCDC_MAKE_FANCY_BEGIN)
	$(NCDC_MAKE_QUIET)$(NCDC_FORT90C) $(NCDC_F90FLAGS) $(NCDC_CONFIG_FLAGS) -c netcdf_unlimdims.F90 $(NCDC_NETCDF_FLAGS)
	$(NCDC_MAKE_FANCY_END)

clean:
	rm -f *.mod *.x *.o
