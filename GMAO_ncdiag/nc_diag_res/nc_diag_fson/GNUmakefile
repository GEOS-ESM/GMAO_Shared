# This file has been set up to build with nc_diag_res.
# Do NOT change this file!
include ../GNUmakefile.config

SRC=src
DIST=lib
BUILD=build

F95=.f90
OBJ=.o
SO=.so
MOD=.mod
EXE=

BUILD_SHARED=0
BUILD_STATIC=1

LIBTARGET=$(DIST)/libnc_diag_fson$(SO)
LIBTARGET_STATIC=$(DIST)/libnc_diag_fson.a

ifeq ($(BUILD_SHARED), 1)
    BUILD_LIBS += $(LIBTARGET)
endif

ifeq ($(BUILD_STATIC), 1)
    BUILD_LIBS += $(LIBTARGET_STATIC)
endif

FC = ifort
FCFLAGS = $(NCDRES_F90FLAGS)
LDFLAGS=

# "make" builds all
all: lib

JSON = $(shell find src -name '*.json')
json: $(patsubst $(SRC)%, $(DIST)%, $(JSON))

FSON = \
    ncdf_string_m \
    ncdf_value_m \
    ncdf_path_m \
    nc_diag_fson
OBJECTS = $(patsubst %, $(BUILD)/%$(OBJ), $(FSON))

$(LIBTARGET) : $(OBJECTS)
	$(NCDRES_MAKE_FANCY_BIG_BEGIN)
	$(NCDRES_MAKE_QUIET)mkdir -p `dirname $@`
	$(NCDRES_MAKE_QUIET)$(FC) -shared -o $(LIBTARGET) $^
	$(NCDRES_MAKE_FANCY_BIG_END)

$(LIBTARGET_STATIC) : $(OBJECTS)
	$(NCDRES_MAKE_FANCY_BIG_BEGIN)
	$(NCDRES_MAKE_QUIET)mkdir -p `dirname $@`
	$(NCDRES_MAKE_QUIET)ar rcs $(LIBTARGET_STATIC) $^
	$(NCDRES_MAKE_FANCY_BIG_END)

lib: $(BUILD_LIBS)

$(DIST)/%.json : $(SRC)/%.json
	mkdir -p $(@D)
	cp -f $< $@

# build fson objects
$(BUILD)/%$(OBJ): $(SRC)/%$(F95)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)mkdir -p `dirname $@`
	$(NCDRES_MAKE_QUIET)$(FC) $(FCFLAGS) -fpic -module $(DIST) -c $< -o $@
	$(NCDRES_MAKE_FANCY_END)

clean:
	rm -rf $(BUILD)
	rm -rf *$(MOD)

clobber: clean
	rm -rf $(DIST)/*
