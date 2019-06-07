include ../GNUmakefile.config
include ../GNUmakefile.libconfig

NC_DIAG_RES_TESTS=test_ncdres_simple.x test_ncdres_invalid.x
NC_DIAG_RES_TESTS_FAIL=test_ncdres_fopen_crash.x \
    test_ncdres_fclose_crash_1.x \
    test_ncdres_fclose_crash_2.x

# Coverage build flags
ifndef COVERAGE
    COVERAGE=0
endif

ifeq ($(COVERAGE), 1)
    NCDRES_F90FLAGS+=-prof-gen=srcpos
    STATIC_COMPILE=1
endif

# Static compilation flags
ifndef STATIC_COMPILE
    STATIC_COMPILE=0
endif

ifeq ($(STATIC_COMPILE), 1)
    NC_DIAG_RES_FLAGS=
    NCDRES_F90FLAGS+=-I.. \
        ../fson/src/fson_string_m.f90 ../fson/src/fson_value_m.f90 \
        ../fson/src/fson_path_m.f90 ../fson/src/fson.f90 ../nc_diag_res.f90 ../ncdres_*.F90 
endif

.PHONY: all clean check precoverage coverage

all: $(NC_DIAG_RES_TESTS) $(NC_DIAG_RES_TESTS_FAIL)

check: $(NC_DIAG_RES_TESTS) $(NC_DIAG_RES_TESTS_FAIL)
	@echo "Running regular tests..."
	@$(foreach ncdr_test,$(NC_DIAG_RES_TESTS),echo "Running: $(ncdr_test)";./$(ncdr_test) || exit 1;echo "Pass: $(ncdr_test)";)
	@echo "Running tests that should fail..."
	@$(foreach ncdr_test,$(NC_DIAG_RES_TESTS_FAIL),echo "Running FAIL test: $(ncdr_test)";./$(ncdr_test) && exit 1 || true;echo "Pass: $(ncdr_test)";)

test_ncdres_simple.x: test_ncdres_simple.f90 $(NCDRES_LIB_PATH)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) \
		test_ncdres_simple.f90 -o test_ncdres_simple.x \
		$(NC_DIAG_RES_FLAGS)
	$(NCDRES_MAKE_FANCY_END)

test_ncdres_invalid.x: test_ncdres_invalid.f90 $(NCDRES_LIB_PATH)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) \
		test_ncdres_invalid.f90 -o test_ncdres_invalid.x \
		$(NC_DIAG_RES_FLAGS)
	$(NCDRES_MAKE_FANCY_END)

test_ncdres_fopen_crash.x: test_ncdres_fopen_crash.f90 $(NCDRES_LIB_PATH)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) \
		test_ncdres_fopen_crash.f90 -o test_ncdres_fopen_crash.x \
		$(NC_DIAG_RES_FLAGS)
	$(NCDRES_MAKE_FANCY_END)

test_ncdres_fclose_crash_1.x: test_ncdres_fclose_crash_1.f90 $(NCDRES_LIB_PATH)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) \
		test_ncdres_fclose_crash_1.f90 -o test_ncdres_fclose_crash_1.x \
		$(NC_DIAG_RES_FLAGS)
	$(NCDRES_MAKE_FANCY_END)

test_ncdres_fclose_crash_2.x: test_ncdres_fclose_crash_2.f90 $(NCDRES_LIB_PATH)
	$(NCDRES_MAKE_FANCY_BEGIN)
	$(NCDRES_MAKE_QUIET)$(NCDRES_FORT90C) $(NCDRES_F90FLAGS) \
		test_ncdres_fclose_crash_2.f90 -o test_ncdres_fclose_crash_2.x \
		$(NC_DIAG_RES_FLAGS)
	$(NCDRES_MAKE_FANCY_END)

$(NCDRES_LIB_PATH): $(NCDRES_LIB_DEPS)
	@echo " *** nc_diag_read library not built yet (or outdated)."
	@echo " *** (Re-)building it now..."
	$(NCDRES_BUILD_CMD)

precoverage:
	@[ ! "$$COVERAGE" = "1" ] && \
		echo "ERROR: COVERAGE must be set to 1 to set build flags and use" && \
		echo "       use coverage tools." && exit 1 || true
	@[ ! -f ~/launch_browser ] && \
		echo "ERROR: You must have a script called launch_browser in your" && \
		echo "       home directory to display code coverage results. It" && \
		echo "       should accept a file as an argument to display the" && \
		echo "       code coverage result page." && exit 1 || true
	@echo " ** Code coverage enabled for this build."
	@echo " ** Tests will be run after build and results will be displayed."

coverage: precoverage clean check
	profmerge
	codecov -srcroot ..
	~/launch_browser CODE_COVERAGE.HTML &
	rm -f *.mod *.o

clean:
	rm -f *.mod *.x *.o
	rm -rf pgopti.* CodeCoverage CODE_COVERAGE.HTML *.dyn
