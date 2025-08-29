# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

### Removed

### Deprecated


## [2.1.4] - 2025-08-29

### Changed

- Updated `esma_mpirun` for `openmpi` to run with `-map-by node -bind-to core` which is a good default.


## [2.1.3] - 2025-07-07

### Fixed

- Fix issue with installing `gmao2ioda` scripts as the `gmao2ioda/` directory might not exist in all fixtures
- Fixed `cmpdir.pl` script and add extra functionality

## [2.1.2] - 2025-06-30

### Added

- Added `gmao2ioda` scripts

## [2.1.1] - 2025-05-13

### Fixed

- Fixed issue with unprotected f2py call in `GMAO_gfio`. Caused failures if build system did not have f2py available.

### Added

- Added `.editorconfig`

## [2.1.0] - 2025-04-23

### Added

- Added Github Actions for CI

### Changed

- Bring in JEDI supporting code
- Updates for `esma_mpirun` and MVAPICH 4
  - Always bind to core
- Refactored helfsurface subroutines to improve the performance

### Fixed

- fix vED
- SLES15 stuff
- fix to reading of table in time ave program

## [2.0.0] - 2025-03-14

### Fixed

- Corrected perl scripts for left/right brace error introduced with SLES-15:
  - `Perl_Config.pm`
  - `vED`

### Removed

- Remove `GMAO_perllib` into its own repository

## [1.9.9] - 2024-10-03

### Changed

- Update ESMF CMake target to `ESMF::ESMF`
- Update `n4zip` with long option for non-alphabetical sorting

### Fixed

- Declare `main` in `makdep.c` as `int` to satisfy gcc 14

## [1.9.8] - 2024-07-31

### Fixed

- Fixed invalid Fortran in `GMAO_ods/ods_diagnc4.f90`
- Added `stdlib.h` in `GMAO_pyobs/VegType_io.c` to fix compilation error with `icx`

## [1.9.7] - 2024-01-30

### Added
- Add new `CHOOSEZ0` option that triggers the use of Charnock coefficient

### Changed
- Move to circleci-tools v2 orb
- Update CI to v2 orb

## [1.9.6] - 2023-11-01

### Fixed

- Fixes for Perl 5.26 on SLES15

## [1.9.5] - 2023-10-18

### Changed

- Updated `esma_mpirun` for MVAPICH to use `mpiexec`
- add info on LEO-GEO: obs table upadate (zero-diff)

## [1.9.4] - 2023-09-15

### Changed

- add Planet IQ to odsmatch.rc

## [1.9.3] - 2023-09-15

### Added

- `dyn_hydro`
  - This program adds the hydrometeors to bkg.eta files that don't have them (such as those coming from MERRA-2). (Mistakenly not
    brought over during CVS-to-Git transition)

### Changed

- Syncs obsys with that used in FP as of 14 Sep 2023
- dyn2dyn updated to allow for user-specified IM x JM when regrinding

## [1.9.2] - 2023-08-24

### Changed

- Update CI to use Baselibs and BCs from CircleCI Orb

## [1.9.1] - 2023-06-06

### Added

- Added protection against specific error conditions to reset the ice to a default profile instead of crashing. This is valid only for MITgcm

### Changed

- Change CICE4 build to shared library, in preparation for CICE6 integration

## [1.9.0] - 2023-05-09

### Changed

- Updated 181 and 91L ak/bk to remove kinks in 1st derivative of DZ. Modified ak at top four levels of L137.
- Converted `pyrob` and `pyrob_CF` to Python 3. Requires MAPL 2.36.0 and ESMA_cmake v3.28.0 for f2py support

## [1.8.0] - 2023-03-10

### Changed

- Merging in changes associated with refactored physics and `gmap4`.

## [1.7.2] - 2023-03-17

### Changed

- Improved trace gas ODS support; incorporates most recent changes to ODS in rt_g530_1.6.0_merge (AMSR, more OMPS)

## [1.7.1] - 2023-02-10

### Changed

- Updated files to be consistent with das tags 5.30
- Updated `CODEOWNERS` as GEOS_Util has been moved

## [1.7.0] - 2023-01-20

### Removed

- Removed GEOS_Util as it is its own repo in GEOS-ESM now

## [1.6.5] - 2023-01-20

### Added

- New coupled diagnostics package based on python3

## [1.6.4] - 2023-01-17

### Changed

- Updates for WSTAR Calculation using Model Omega. Also, several basic QUICKPLOT updates.

### Removed

- Remove `hasw` _constraint_ choice from remapping questions since they are being [decommissioned by NCCS](https://www.nccs.nasa.gov/nccs-users/SCU10%2613-Decommismion)

## [1.6.3] - 2022-12-08

### Changed

- Moved to GitHub Actions for label enforcement
- Updated CircleCI to Baselibs 7.7.0
- Set default data ocean to be `CS` at C90+ resolution in `remap_restarts.py`
- Updated `remap_restarts.py` to use argparse

### Fixed

- Fix an issue with `regrid.pl` and `remap_restarts.py` regridding from MERRA-2 from 2021-06 to 2021-09
- Fix an issue with `remap_restarts.py` for post 2011 regridding from MERRA2

## [1.6.2] - 2022-11-17

### Fixed

- Fixed `gcmpost.script` logic to generate `xdf.tabl` files for collections with `monthly: 1` attribute
- Update EASE grid tile file name to contain pfafxxx
- Add a dependency on GMAO_perllib in GMAO_etc

## [1.6.1] - 2022-09-22

### Changed

Apply latest updates from L. Takacs to the plot package

## [1.6.0] - 2022-09-12

### Changed

- Updates the plot package with the renamed carbon species.

## [1.5.8] - 2022-09-08

### Added

- Add YAML validator

### Changed
- Converted obsys_rc.py, check_obsysrc.py, and associated unittests to Python 3
- Added features to cmpdir.pl (option to toggle between diff and xxdiff; added -B to -bwi flags)

- Generated command line to run the program
- changed the location of the temporary folder for remap_restarts MERRA-2 case
- Added new remap tests
  - amip_c180Toc90
  - c180Toc360
  - c360Toc24
  - s2sv3Toc12MOM6

### Fixed

- Add `target` attribute to variables in `nc_diag_cat`. Needed to build with Intel 2022.1

## [1.5.7] - 2022-08-01

### Changed

Updates consistent with CVS Jason-4_0_p1_sky. Releases non-needed memory from zcmp/progz. Otherwise, minor plot updates.

## [1.5.6] - 2022-07-22

### Added

- Introduce series of remap_restarts Python scripts to GEOS_Util/post
- Add a series of tests for remap_restarts

### Changed

- Updated CI to modern v1 orb
- Updated CI to use Baselibs 7.5.0

### Fixed

- Implemented a more robust check to test if a history's monthly attribute is 1 (on) and not commented out.
- Fixed setting of `ks` in `m_set_eta.F90` to be based on `bk`
- Added check for infinity in time_ave.F and replace with undef

## [1.5.5] - 2022-06-08

### Changed

- Changed restart file geosachem to achem in regrid.pl
- More updates to CMake to more canonical CMake style (NetCDF, ESMF, etc.). These were missed in previous go-arounds as they are only built with ADAS. (Requires ESMA_cmake v3.15.1)

## [1.5.4] - 2022-05-03

### Added

- Added statsNx.rc for screen level variable fstats

### Changed

- Added options for land-only and screen-level variables in `fstats.x` and `g5fcst_stats.pl`
- Added a few new tags to `regrid.pl`
- Changes allow `gcmpost.script` and `gcmclim.script` to handle `collection.monthly: 1` output and stop automatic archiving.

### Fixed

- Fixed a minor CMake issue to keep all mod files in `build/include`

## [1.5.3] - 2022-03-18

### Changed

- Modified chckhist.new to handle OPS HISTORY.rc,  fixed minor bugs in 3CH.F90 and 3CH.j

## [1.5.2] - 2022-03-18

### Added

- Added preprocessing team as CODEOWNER for the GEOS_Util/pre directory
- added a way to process Reynolds only for producing SST and Ice Concentration data, using a land-sea mask.
  The _new_ file is: `proc_SST_FRACI_reynolds_quart.F90` and modified: `read_Reynolds.F90`, `CMakeLists.txt`
  **Note**: the contents of this directory will be defunct once `ExtData` mechanics are implemented, WIP with @bena-nasa

### Changed

- bugfix to token_resolve() routine in GMAO_etc/Manipulate_time.pm

### Fixed

- Updates to CMake to support Spack

### Added

- Added capability to produce netcdf ocean pre-processing datasets, with doc and notebooks to demo.

## [1.5.1] - 2022-02-04

### Changed

- Compress CircleCI artifacts
- Updated CircleCI to use Orb

## [1.5.0] - 2021-12-16

### Added

- Add support for GOCART2G restarts in `regrid.pl`

### Removed

- Moved lightning files to Chem_Shared

## [1.4.13] - 2021-12-15

### Fixed

- Quickplot and quickstat bugs

### Added

- Quickplot now supports plotting GOCART-2G collections
- Support for Three Corner Hat (3CH) Analysis

## [1.4.12] - 2021-12-09

### Changed

- Update `regrid.pl`
  - Add options for MOM5 and MOM6 tile files
  - Add ability to use git tags for "tagin" and "tagout"

### Fixed

- Update to `idcheck.pl` to allow fvsetup to check whether the expid already exists in the SemperPy databases

- Moved Lightning_mod.F90 and lightning_toolbox_mod.F90 from GEOS_Shared to a different repository (GEOSchem_GridComp)

## [1.4.11] - 2021-11-03

### Changed

- add Cascade knob to g5fcst_stats.pl and regrid.pl
- revised dyn_blob: more general on the blobs
- make sure echorc.x exits w/ success code when applicable
- Updated CI to use Baselibs 6.2.8
- Updated `pyrob` to work with GEOS-IT files
- Changed the Intel MPI and MVAPICH2 flags in `regrid.pl` to be modern

## [1.4.10] - 2021-10-08

### Added

- Added a new `GMAO_eu` target to create a `libGMAO_eu.a` like in GNU Make days
- Added `parallel_untar.py` script
- Added `dyn_blob.x` and `dyn_fsens_conv.x` to Hermes

### Fixed

- CMake fix for non-Intel compilers in GMAO_ods

### Changed

- Updates to documentation tables for KX values
- Updates to `obsys-nccs.rc`
- Updates to `g5fcst_stats.pl`
- Updated the CI for GMAO_Shared to do Intel and GNU

## [1.4.9] - 2021-10-06

### Fixed

- Fixed issue with `regrid.pl` and regridding `catch_internal_rst`

## [1.4.8] - 2021-10-05

### Fixed

- Fixed issue with CICE4 by compiling with old non-vectorized `Release` flags when compiling with Intel. Requires ESMA_cmake v3.6.1

## [1.4.7] - 2021-10-04

### Added

- Added `pyrob_CF` script

### Changed

- Updates to support Catchment-CN.4.5 in addition to Catchment-CN.4.0

## [1.4.6] - 2021-09-15

### Changed

- Updates to `regrid.pl` to allow processing of restarts for two different versions of Catchment-CN land model
- change is zero-diff for Catchment and Catchment-CN4.0

## [1.4.6] - 2021-07-21

### Changed

- Updates to plots package from L. Takacs

## [1.4.5] - 2021-06-25

### Fixed

- Fix for IASI (#202)

## [1.4.4] - 2021-06-25

### Added

- Add new `echorc.pl` script (alternative to `echorc.x`)
- Added 181 levs to `GMAO_hermes/dyn2dyn.f90` and a frequency change in `GMAO_etc/obsys-nccs.rc`

### Changed

- Changed `esma_mpirun` for MVAPICH2

## [1.4.3] - 2021-06-11

### Added

- Add changes consistent with what is in GEOSadas 5.28

## [1.4.2] - 2021-05-25

### Added

Add ability to write out energy components to file.

### Fixed

Bugfix to prevent a seg-fault when calculating the lightning flash rate implemented in HEMCO/GEOS-Chem.

## [1.4.1] - 2021-05-14

- Renamed `LANL_Shared/LANL_cice` to `LANL_Shared/CICE4`

### Fixed

- In `regrid.pl`: Fixed the -wemin and -wemout options so that they will accept integer values; Also added Jason-NL BCS tag choice
- Fixes for PSAS code and Intel MPI

### Changed

- Multiple updates brought over from GEOSadas work (see #166)
- Update F2PY module calls to support both Python2 and Python3 loaded at same time

## [1.4.0] - 2021-04-15

### Fixed

- Sync atmOcnIntlayer with that in GEOS-FP GEOS-5.27.1 (02/2021) GEOSadas-5_27_1_p3
- Fixed build for directories that are built as part of the GEOSdas

## [1.3.10] - 2021-04-02

### Changed

- Multiple updates brought over from GEOSadas work (see #166)

## [1.3.9] - 2021-03-17

### Fixed

- Fix for out-of-bounds error in lightning module (#99)

### Changed

- Stats plot updates

### Added

- Extend binary tile to supprot future river route component development

### Removed

- Remove `CMIP_1977_1982` directory in `GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend`

## [1.3.8] - 2021-02-12

### Changed

- Plot and stats updates

### Removed

- Eliminated references to MAPL_COMM - these are no longer used and will
be deleted in future releases of MAPL

## [1.3.7] - 2021-01-29

### Changed

- Plot and stats updates

### Fixed

- Fix flags for `zonal.f` compilation

## [1.3.6] - 2021-01-12

### Changed

- Added `-partition` option to `regrid.pl`

### Fixed

- Fixed compilation flags for `zonal.f` to match that of CVS
- Fixed bug in `res/zonal.gs` for `zonal.x` location
- Added flag to regrid_forcing_esmf.x to force a 0 to 1 range when regridding files that should use fractions

### Added

- Allows PRs with "0-diff trivial" labels to skip updating `CHANGELOG.md`

## [1.3.5] - 2020-12-10

### Fixed

- Use `CONFIGURE_DEPENDS` with `file(GLOB)` calls
- Fix OpenMP in GMAO_stoch

### Added

- Add support for `Aggressive` build type

## [1.3.4] - 2020-11-25

### Fixed

- Updates for DSO work

### Changed

- Update plots consistent with latest cvs tags

## [1.3.3] - 2020-10-28

### Changed

- Added Docker authentication for CI
- Update CI images for Baselibs 6.0.22

### Fixed

- Updates to `plots/configure` so that it is run at install time. This
  should allow others to then run plots on another person's build

### Removed
### Added

- Add warning-level loggers to `GEOS_Utilities` (D)QSAT code for NaN detection

## [1.3.2] - 2020-10-14

### Fixed

- Fixed regrid_forcing.x and regrid_forcing_esmf.x so they will work with MAPL2.2 and beyond

### Changed

- Updated the NCPUs detection in various post scripts to be SLURM-aware

## [1.3.1] - 2020-10-13

### Changed

- Updated `changelog-enforcer.yml` to v1.4.0

### Added

- Modify plots to display new constraint diagnostics
- Update `atmOceanIntLayer.F90` add AVOIL_v0 which wraps the update of surface skin variables (temperature, water mass, salinity, etc) into a single call.

## [1.3.0] - 2020-09-28

### Fixed

- Undo the change to `GEOS_Utilities.F90` in v1.2.0. This has a bug at the end of the table (#123)

## [1.2.0] - 2020-09-25

### Changed

- Update `GEOS_Utilities.F90` to match `GEOSadas-5_27_0` (#115)

## [1.1.10] - 2020-09-25

### Changed

- Updated the CircleCI Image

## [1.1.9] - 2020-09-15

### Changed

- All the subroutines that are called from within OPENWATERCORE of GEOS_OpenWaterGridComp.F90 have been moved to GMAO_Shared so they can be shared across applications and components.

## [1.1.8] - 2020-08-12

### Fixed

- Fixes to allow Intel 19.1.2 to use CICE

## [1.1.7] - 2020-07-23

### Changed

- Updates to plots
- Modify mpirun flags for Open MPI 4.0.4

### Fixed

- Allows the JCAP functionality to work again in mkiau gridcomp as that needs the r4 version of ncep_sp in gmao_transf

## [1.1.6] - 2020-07-06

### Fixed

- CMake update for building with GCC 10 in Release mode
- Fix for incorrect OpenMP usage in `lightning_toolbox_mod.F90` as detected by Intel Fortran 19

## [1.1.5] - 2020-06-26

### Added

- Support for new vertical resolutions

### Changed

- Updates to coupled plotting package (coupled_diagnostics)

## [1.1.4] - 2020-06-05

### Changed

- Rolls back the constraint on gcmpost.script to only operate on pressure-level collections.

### Fixed

- Enables correct post proccessing of MAPL monthly collections.
- Added ignore_nan option for `time_ave.F` (off by default).
- Update allowing both ifort and gfortran compilations to read CICE binary grid files.

## [1.1.3] - 2020-04-09

### Added

- Adds lightning module to compute flash rate (requires: GEOSchem_GridComp v1.3.3+ and GEOSgcm_GridComp v1.8.3+)

### Fixed

- Allows regrid.pl to run on SLES-12
- Enabled compilation of convert_aerosols.x as R8

## [1.1.2] - 2020-03-26

### Fixed

- Modified gcmpost.script to only operate on pressure-level collections.
- Fixes issues where systems do not have either ImageMagick or F2Py. If F2Py is not found, then F2Py targets are not built.

### Changed

- Adding an option to pass in the model (default: hasw) to quickstat.
- Changes to GMAO_hermes:
  - Split off the independent sections of module m_topo_remap used by FV core into shared_topo_remap
  - write_eta.F90 added to prepare for a config file which eventually will replace m_set_eta module
  - Add option to build HERMES_LIGHT

## [1.1.1] - 2020-03-04

- Enable additional upper levels for forecast stats plots.
- Add QITOT & QLTOT to horizontal plots.
- Add aerosols to time series plots.
- Use sbatch at NCCS, qsub at NAS.

## [1.1.0] - 2020-02-13

### Changed

- Many [updates to support MAPL 2.0](https://github.com/GEOS-ESM/GMAO_Shared/releases/tag/v1.1.0)

## [1.0.17] - 2020-02-04

### Changed

- Updates to produce stratospheric forecast statistics up to 1mb.

## [1.0.16] - 2020-01-28

### Fixed

- Better representation of statistical significance in stats plots accounting for line thickness
- Updates to perl scripts in post to accommodate SLES12

## [1.0.15] - 2020-01-16

### Changed

- Updates to catchment restart regridding.
- Add WEMIN input and output and an attempt to make in/out restarts consistent.
- Clearer names for WEmin in/out.
- Provide info about WEmin before prompting so the user is informed about their choices.
- Update post/plot to CVS tag Jason-3_3_aoil.
- Updated file in pre directory to be consistent with ops.

## [1.0.13] - 2019-11-07

### Changed

- Build stats.x as big-endian

### Fixed

- Fixes undefined references to GEOSSRC, GEOSBIN, and GEOSAPP in quickstat

## [1.0.12] - 2019-09-27

### Changed

- Updates for s2s

### Fixed

- Fixed bug in utils

## [1.0.11] - 2019-09-13

### Fixed

- Moved getco2 for ldas

### Changed

- Updates for IAU_Error_plots

## [1.0.10] - 2019-09-06

### Fixed

- CMake fix for GFIO existence test

## [1.0.9] - 2019-09-04

### Fixed

- Make GMAO_gfio build contingent on the subdir being present.

## [1.0.8] - 2019-08-28

### Changed

- Updates from CVS GEOSadas-5_25_2 tag

## [1.0.7] - 2019-08-27

### Changed

- Modified color shading based on 99.99% confidence

## [1.0.6] - 2019-08-27

### Removed

- Extracted FMS to separate repo

## [1.0.5] - 2019-08-27

### Fixed

- Yet another file missing from installation for plots to function

## [1.0.3] - 2019-08-05

### Changed

- Montage Plots Updates from `Jason-3_2`

## [1.0.2] - 2019-07-26

### Changed

- Updates from Jason-UNSTABLE as of 2019-Jul-26

## [1.0.1] - 2019-07-26

### Changed

- Updates from Jason-3_1

## [1.0.0] - 2019-07-25

### Changed

- Initial Release with Semantic Versioning
  - Repository split from CVS ESMA Repository.
  - Equivalent to cvs/GEOSadas-5_25_0 release.


