# GMAO_etc: architecture

**Audience:** coding agents and developers working on the GMAO_etc scripts and utilities in GMAO_Shared.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMAO_Shared` @ `595999b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GMAO_etc/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

> **Build status on this line:** not built, not ported. On `release/MAPL-v3`, `CMakeLists.txt` (repository root) calls `esma_add_subdirectories` only for `GMAO_mpeu`, `GMAO_hermes` and `GEOS_Shared` (comment: "not yet ported to MAPL3 - only build GMAO_mpeu, GMAO_hermes (light), and GEOS_Shared for FVdycoreCubed_GridComp"); `GMAO_etc` appears only in the commented-out `if ( NOT HERMES_LIGHT )` block. The code in this directory is the MAPL2 code, kept for a later port.

## 1. Purpose
GMAO_etc is not a Fortran library. It is a bag of run-time and data-assimilation
helper scripts (Perl, Python, csh/sh), a few small command-line programs
(Fortran, C), GrADS scripts and a documentation extractor (`protex`). All of it
is installed into `bin/`, `etc/` or `lib/grads/`.

## 2. Build (`GMAO_etc/CMakeLists.txt`)
- Target `GMAO_etc` is an empty custom target (`add_custom_target ... ALL`)
  that depends on `GMAO_perllib` (the Perl modules `Manipulate_time.pm` etc.
  live in that separate repository).
- Compiled executables: `argopt` (C), `fresolv.x`, `gethdftinc.x` (also links
  NetCDF Fortran), `zeit_ci.x`, `zeit_co.x`, `zeit_pr.x` (with `zeit_usage.f`,
  `zeit_read.f`, `zeit_jday.f`), `mkdrstdate.x`, `rst_date`. All Fortran ones use
  `${FREAL8}`; all but the last two link `GMAO_mpeu`.
- Installed to `bin`: all `*.pl`, `*.pm`, `*.py`; a list of extension-less Perl
  scripts (`esma_mpirun`, `tick`, `token_resolve`, `pasta`, `pesto`,
  `require_obsys`, `now`, `num_days_in_month`, `n4zip`, `r_dist`, ...);
  `Functions.csh`, `pyrob`, `pyrob_CF`, `RmShmKeys_sshmpi.csh`,
  `rmshmkeyhere.sh`, `protex`.
- Installed to `bin` through `configure_file` (`@DASPERL@` -> `/usr/bin/perl`):
  `bacon`, `fvgetf`, `grads2grib`, `rget`, `rput`, `rst4ncepf`, `rst2rcyc`.
- Installed to `etc`: `obsys*.rc`, `*.tmpl`, `pesto.arc`, `PUBLICTAG`.
  Installed to `lib/grads`: `lats4d.gs`, `lats4d`, `info.gs`.
- Checked by comparing the directory listing with the CMake rules. Covered by
  the `*.pl`/`*.pm`/`*.py` globs: all Perl and Python files. Covered by the
  `obsys*.rc` glob: `obsys*.rc`. No compile or install rule: `pasta.c` (the
  installed `pasta` is the Perl script), `m_zeits.F90`, `fvtimes`, `datank`,
  `fln`, `lnlist`, `mkfake`, `acquire.rc`, `bacon.rc`, `lnlist.rc`,
  `Functions.ppt`, `utFunctions.csh`, `utFunctions.csm`, `unused_dummy.H`,
  `zeit.README`, `unittests/`.

## 3. Contents
| Group | Files | Role |
|---|---|---|
| Job / run-time | `esma_mpirun` (mpirun wrapper that filters flags per MPI), `prund.pl`, `jobIDfilter`, `lsfpL`, `subset_nodefile.pl`, `RmShmKeys_sshmpi.csh`, `rmshmkeyhere.sh`, `Run_parallel.pm`, `g5_modules_perl_wrapper`, `parallel-untar.py` | Launch and clean up jobs |
| csh helpers | `Functions.csh` (alias-based function emulation; example in `utFunctions.csh`), `show_children.pl`, `stripname`, `fln`, `lnlist` | Script support |
| Date / time | `tick` and `token_resolve` (thin wrappers over `Manipulate_time.pm` in GMAO_perllib), `now`, `num_days_in_month`, `dateutils.py`, `diffdates.py`, `randates.py`, `mkdrstdate.f`, `rst_date.f` | Date arithmetic, `%y4`-style token expansion, restart date files |
| Observing-system (obsys) | `obsys_check.pl` (is a data class available for a range per `obsys.rc`), `require_obsys[.pl]`, `verify_obsys.pl`, `obsclass_filter.pl`, `OScheck.pl`, `obsys_rc.py` (class `ObsysRc`), `check_obsysrc.py`, `obsys*.rc`, `setup_modis_l2data.pl` | Query and validate observation-resource files |
| Remote copy / file transfer | `rget`, `rput`, `pasta`/`Pasta.pm`, `pesto`/`pesto.arc`, `bbscp.pl`, `mkdir_remote`, `acquire.rc`, `r_dist`, `fvgetf`, `n4zip` | Get/put files across hosts and archives |
| Timing (zeit) | `zeit_ci`, `zeit_co`, `zeit_pr` (check in / out / print), `zeit_summary.pl`, `m_zeits.F90`, `zeit.README` | Command-line multi-timer for scripts |
| Small programs | `fresolv.f90` (file-name resolution with `m_StrTemplate`, `m_FileResolv`), `argopt.c`, `gethdftinc.f` | Command-line tools |
| Perl modules | `Arch_utils`, `Conv_utils`, `GMAO_utils`, `Perl_Config`, `Verbose`, `WriteLog`, `Query`, `Recd_State`, `Round`, `getsponsor` | Shared by the scripts above |
| Data conversion / misc | `grads2grib`, `rst4ncepf`, `rst2rcyc`, `ana4dupd.pl`, `extract_*.pl`, `mptfix.pl`, `pyrob`, `pyrob_CF`, `bacon`, `subst_path.pl` | DAS and post-processing helpers |
| GrADS | `lats4d.gs`, `lats4d`, `info.gs` | Visualization / subsetting |
| Documentation | `protex` (ProTeX 2.00, turns `!BOP`/`!EOP` prologues into LaTeX) | Doc extractor |
| Tests | `unittests/` (`test.csh`, `test_check_obsysrc.py`, `input`, `outexp`) | Test for `check_obsysrc.py` only; run by hand (`unittests/README`) |

## 4. Who uses it
GEOSgcm fixture (`src`), greps excluding GMAO_etc itself, output limited:
- `esma_mpirun` (`$GEOSBIN/esma_mpirun`, e.g. `set RUN_CMD = '$GEOSBIN/esma_mpirun -np '`):
  `Applications/@GEOSgcm_App/gcm_setup`, `scm_setup`, `gmichem_setup`,
  `geoschemchem_setup`, `stratchem_setup`; FVdycoreCubed
  `scripts/fv3.j` and `scripts/fv3_setup` (under
  `Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSsuperdyn_GridComp/@FVdycoreCubed_GridComp`);
  `Shared/@GMAO_Shared/@GEOS_Util/post/3CH.j`;
  `@GEOS_Util/pre/remap_restart/remap_upper.py`, `remap_other_restarts.py`,
  `remap_catchANDcn.py`.
- `tick` (`$GEOSBIN/tick`): `Applications/@GEOSgcm_App/gcm_forecast.tmpl`,
  `gcm_run_CPLFCST360S2S.j`, `gcm_run_CPLFCST360S2Sall.j`;
  `@GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/make_ostia_quart.csh`,
  `make_ostia_eigth.csh`, `SST_ICE_bcs.csh`. Other calls use
  `$GEOSUTIL/post/tick`, a copy in GEOS_Util (`gcm_forecast.setup`,
  `gcm_forecast.tmpl`, `@GEOS_Util/post/gcmpost_CPLFCSTfull.script`): see
  Known traps.
- `token_resolve`: `@GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/BLEND-OSTIA-REYNOLDS.pl`;
  also named in `@GMAO_perllib` (`acquire`, `Manipulate_time.pm`). No hit in
  `Applications`.
- `mkdrstdate`: `@GEOS_Util/pre/remap_restart/remap_analysis.py`.
- `rst_date`: `@GEOS_Util/plots_v2/stats_v2/plots_corcmp.py`, `scorecard.py`,
  `plots_syscmp.py`.
- `lats4d`: `@GEOS_Util/pre/prepare_ocnExtData/sst_ice_helpers.F90`,
  `@GEOS_Util/pre/NSIDC-OSTIA_SST-ICE_blend/make_ostia_quart.csh`;
  `GMAO_hermes/README.geos2fv`.
- `protex`: no caller; `GMAO_hermes/m_vdc2vdc.F90` and `@MAPL/base/TeX`
  mention it only in prologue text / build files.
- `zeit_*`: named in `GMAO_mpeu/doc/Timing.tex`, `GMAO_mpeu/m_zeit.F90`,
  `GMAO_hermes/ana5sfc.F90`.
- No user found in the GEOSgcm fixture for `obsys_*`, `require_obsys`,
  `Functions.csh`, `protex`, and none in `src/Applications` for `token_resolve`.
  *Unverified:* they may be used by GEOSadas or operational job scripts.

## 5. Known traps
- `m_zeits.F90` is not compiled by this directory's CMake, though the zeit
  programs read like its users (no `use m_zeits` found in the fixture).
- Several scripts are Perl wrappers that need `Manipulate_time.pm` or
  `Pasta.pm` on the Perl path; the former is in `@GMAO_perllib`, not here.
- `bacon`, `rget`, `rput` etc. carry an `@DASPERL@` shebang and run only
  after `configure_file`.
- `tick` exists twice: `GMAO_etc/tick` (installed as `$GEOSBIN/tick`) and
  `@GEOS_Util/post/tick` (`$GEOSUTIL/post/tick`); both are called by scripts
  in `Applications/@GEOSgcm_App`.
- `GMAO_etc` build target has no sources; linking it provides nothing.

## 6. Observations (noted, not fixed)
- `*.pl`, `*.pm`, `*.py` are picked by glob: any new file there is installed.
- `unittests/` covers one script; no CTest registration seen in CMake.
- Several scripts have legacy headers (DAO, `/usr/bin/perl` paths, Python 2.3
  era notes in `dateutils.py`).
