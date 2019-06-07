#
# Makefile for ESMA components.
#
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

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
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
	@echo "         SITE = $(SITE) "


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------
UTIL    =

SUBDIRS = lanso1 plan

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          help install clean doc

export ESMADIR BASEDIR ARCH SITE

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done

local_esma_install local_install: $(LIB)
	@echo No local install in here

#                  --------------------
#                  User Defined Targets
#                  --------------------

dist:
	$(MAKE) distclean
	( cd ..  ;\
          $(TAR)  cvf evac-`date '+%d%b%Y'`.tar src ;\
          $(GZIP) evac-`date '+%d%b%Y'`.tar )

distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd]
	$(MAKE) esma_distclean

binclean:
	-$(RM) -r $(ESMADIR)/Config 
	-$(RM) -r $(ESMADIR)/$(ARCH)

realclean:
	$(MAKE) distclean
	$(MAKE) binclean

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros
#.
