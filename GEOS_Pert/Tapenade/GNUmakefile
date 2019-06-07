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
       ESMADIR := $(PWD)/../../../..
endif

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

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

THIS = $(shell basename `pwd`)
LIB  = lib$(THIS).a

esma_install install: $(LIB)
	$(MKDIR) $(ESMALIB) $(ESMAINC)/$(THIS)
	$(CP) -p *.a    $(ESMALIB)
	$(CP) -p *.mod  $(ESMAINC)/$(THIS)

esma_clean clean:
	-$(RM) *~ *.[aox] *.mod 

esma_distclean distclean:
	-$(RM) *~ *.[aoxd] *.mod

esma_doc doc:
	@echo no documentation here

#                  --------------------
#                  User Defined Targets
#                  --------------------

SRCS = $(wildcard \
         adBuffer.f  \
         adStack.c \
         tapenade_iter.F90 )

DEPS := $(addsuffix .d, $(basename $(SRCS)))
OBJS := $(addsuffix .o, $(basename $(SRCS)))

FREAL = $(FREAL4)

INC_DIRS = . $(INC_MPI)
MOD_DIRS = . $(INC_DIRS) 

USER_FINCS = $(foreach dir,$(INC_DIRS),$(I)$(dir))
USER_FMODS = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FDEFS   = $(D)SPMD
 #$(D)BPOINT

vpath % $(INC_DIRS) $(MOD_DIRS)

THIS_CFIO = MAPL_cfio_r4

$(LIB) lib : $(DEPS) $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

# Make sure dep files are not remade during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

#.
