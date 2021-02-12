# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Fixed
### Changed
### Added
### Removed

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


