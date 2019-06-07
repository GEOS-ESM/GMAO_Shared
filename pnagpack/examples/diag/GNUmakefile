#
# Makefile for ESMA components.
#
# REVISION HISTORY:
#
# 3mar2004  Zaslavsky  Initial imlementation.
# 20Oct2004  da Silva  Standardization
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
  include $(ESMADIR)/Config/GMAO_base.mk  # System dependencies

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------


THIS = 
LIB  = 
BIN  = diag.x

esma_install install: $(LIB) $(BIN)
	$(MKDIR) $(ESMALIB) $(ESMAETC) $(ESMAINC)/$(THIS) $(ESMABIN)

esma_clean esma_distclean clean:
	$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *.x

esma_doc doc:
	@echo "Target $@ not implemented yet in `pwd`"


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
	@echo "        FREAL = $(FREAL)"

#                  --------------------
#                  User Defined Targets
#                  --------------------

LIB_LANSO  = $(ESMALIB)/liblanso.a 
LIB_PLANSO = $(ESMALIB)/libplanso.a 

SRCS = convck.f diag.f 
 
OBJS := $(addsuffix .o, $(basename $(SRCS)))
DEPS := $(addsuffix .d, $(basename $(SRCS)))

FOPT = $(FOPT3)

USER_FFLAGS = $(BIG_ENDIAN) -I$(INC_MPI)
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
FDP =

vpath % $(MOD_DIRS)

$(LIB) lib : $(DEPS) $(OBJS)
	$(RM) $(LIB)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

diag.x : $(LIB) diag.o convck.o 
	$(LD) $(LDFLAGS) -o $@ diag.o convck.o $(LIB_PLANSO) $(LIB_LANSO) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)
#%.x : $(LIB) %.o 
#	$(LD) $(LDFLAGS) -o $@ $*.o convck.o $(LIB_PLANSO) $(LIB_LANSO) $(LIB_MPI) $(LIB_SCI) $(LIB_SYS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.
