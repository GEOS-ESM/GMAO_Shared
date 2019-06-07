#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 09Jun2003  da Silva  First crack. 
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # Generic stuff

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help help:
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE)"

THIS := $(shell basename `pwd`)
LIB  = lib$(THIS).a

BINS = test_nc_unlimdims.x nc_diag_cat.x

esma_install install: $(LIB) $(BINS)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a        $(ESMALIB)
	$(CP) -p *.mod      $(ESMAINC)/$(THIS)
	$(CP) -p $(BINS)    $(ESMABIN)

esma_clean clean:
	-$(RM) *~ *.[aox] *.mod *.x

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.mod *.x

esma_doc doc:
	@$(PROTEX) $(PROTEX_FLAGS) *GridComp*.[fF]* > $(ESMADOC)/$(THIS).tex


#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS := gen_diff_size_netcdf.f90 kinds.F90 ncdc_climsg.F90 ncdc_cli_process.F90 \
        ncdc_data_MPI.F90 \
        ncdc_dims.F90 ncdc_metadata.F90 ncdc_realloc.F90 ncdc_state.F90 \
        ncdc_types.f90 ncdc_vars.F90 netcdf_unlimdims.F90

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) \
        $(addsuffix .d, $(basename $(BINS))) 

LIBS  = $(LIB) $(LIB_SDF) $(LIB_MPI) $(LIB_SYS)

INC_DIRS = . $(INC_GMAO_SHARED) $(ESMAINC)/nc_diag_write
MOD_DIRS = . $(INC_DIRS) $(INC_ESMF) \
             $(foreach dir,$(SUBDIRS),$(ESMAINC)/$(dir))

USER_FDEFS  = $(D)_REAL8_ $(D)ERROR_TRACEBACK $(D)QUIET $(D)USE_MPI
USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir))

FREAL = $(FREAL4) # for now, require 32 bit reals (R4)

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

test_nc_unlimdims.x : test_nc_unlimdims.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB) $(LIB_SDF) $(LIB_MPI) $(LIB_SYS)

nc_diag_cat.x : nc_diag_cat.o $(LIB)
	$(LD) -o $@ $+ $(LDFLAGS) $(LIB) $(LIB_SDF) $(LIB_MPI) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.


  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

