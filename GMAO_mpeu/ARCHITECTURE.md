# GMAO_mpeu: architecture

**Audience:** coding agents and developers working on the GMAO_mpeu library (message passing and environment utilities) in GMAO_Shared.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMAO_Shared` @ `595999b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`; this directory is `GMAO_mpeu/`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Library-level document: purpose, contents and users. Read the source for routine detail.

> **Build status on this line:** built. `CMakeLists.txt` (repository root) adds `GMAO_mpeu` unconditionally (`esma_add_subdirectories(GMAO_mpeu)`), in the restricted build "only build GMAO_mpeu, GMAO_hermes (light), and GEOS_Shared". Nothing in this directory differs from the MAPL2 line `main` (`git diff origin/main...origin/release/MAPL-v3` has no entry under `GMAO_mpeu`).

## 1. Purpose

"mpeu" (message-passing environment utilities) is a Fortran utility library of
small, mostly independent `m_*` modules: error exit, standard I/O units,
parallel-safe log output, resource-file parsing (`inpak90`), string and
file-name templates, timing, sorting, and kind definitions. It is the oldest
low-level layer under the other GMAO_Shared libraries and the GOCART/Chem code.

## 2. Build

`GMAO_mpeu/CMakeLists.txt` compiles one source list (`srcs`: `get_zeits.c`, 39
Fortran files) into two libraries:

| Target | Depends on | Notes |
|---|---|---|
| `GMAO_mpeu` (`${this}`) | `MPI::MPI_Fortran` | real MPI; the target most CMake users link |
| `GMAO_eu` | `mpi0` | same sources, with the serial MPI stub; also creates `include/GMAO_eu` and `etc/GMAO_eu` directories |

- `mpi0` is a subdirectory target (`mpi0/CMakeLists.txt`): "no-MPI" stubs
  (`mpi_init`, `mpi_send`, `mpi_allreduce`, `mpi_wtime`, ... one file per call,
  `m_mpi0.F`, `m_SHMEM.F90`, `Const_MPItype.H`, plus `../get_zeits.c`). Its
  `mpi0/README` says the interfaces are deliberately plain implicit F77 style.
  `jic/` (headers only: a fake `mpif.h`, plus `README`) is on its include path
  (`target_include_directories ... ../jic`), so it is included by `mpi0`.
- `GMAO_eu` is linked by `GMAO_transf` executables (`dyn_divo.x`, `ec2fv.x`, ...).
  Other users in the fixture link `GMAO_mpeu`.
- `assert.H` is copied to the include directories of both targets.
- `m_FileResolv.f90` is lower-case suffix but gets the preprocessor flag `${PP}`.
- Both targets get `${MISMATCH}` Fortran options (argument type/rank mismatch
  tolerance; the CMake comment cites NAG).
- Present but not compiled: `m_StringLinkedList.F90`; `unused_dummy.H`;
  `doc/` (LaTeX manual: `mpeu.tex` and chapters); `README` (2005 change log by Jing Guo); `About..` (empty).

## 3. Contents

| Group | Modules | Role |
|---|---|---|
| Kinds | `m_realkinds`, `m_intkinds`, `m_sysclocks` | real/integer KIND definitions; binding to system clocks in C (`get_zeits.c`) |
| Standard I/O, errors | `m_stdio`, `m_ioutil`, `m_mpout`, `m_die`, `m_dropdead`, `m_flow`, `m_fpe`, `m_ior4` | standard unit numbers; convenient I/O helpers; mergeable parallel log output; abort with `mpout` flushed; styled abort; call-tree tracing; floating-point-exception control; read/write a real*4 record |
| Resource files | `m_inpak90` | "Inpak90" resource-file reader (table/label parsing) |
| Strings and names | `m_chars`, `m_String`, `m_StrTemplate`, `m_FileResolv`, `m_Filename` | character-class operations; string pointers; GrADS-style `%y4%m2...` templates; template-based file resolution; file-name manipulation |
| Dates | `m_date`, `m_datetime` | Julian-day arithmetic; (date,time) algebra |
| Lists, sorting | `m_List`, `m_mergedList`, `m_SortingTools`, `m_MergeSorts`, `m_IndexBin_char`, `m_IndexBin_integer`, `m_IndexBin_logical`, `m_Group`, `m_Permuter`, `m_rankMerge` | list manager; lists merged across PEs; sorting collection; merge sorts; indexed bin-sorting templates; group maker; permute/unpermute; merge by ranking |
| MPI helpers | `m_mpif90`, `m_mpif`, `m_parDOT`, `m_showDistrib` | portable MPI interface module; `mpif.h` COMMON interface; parallel dot product; show distributions |
| Memory, misc | `m_mall`, `m_random` | bookkeeper of user-allocated memory; ZUFALL random numbers (F77) |
| GrADS files | `m_GrADS`, `m_rGrADS`, `m_wGrADS` | read/write GrADS data files (`m_GrADS` wraps the reader and writer) |
| Timing | `m_zeit` | check in/out time stamps; uses `get_zeits.c` |

## 4. Who uses it

Counts are `use m_<name>` files in `GEOSgcm/src` (fixture), approximate, carried over from the MAPL2-line document and not recounted on this line (many users are in components not built here).

| Module | Files | Examples |
|---|---|---|
| `m_die` | ~149 | widely used, including `GMAO_hermes`, `GMAO_transf` |
| `m_stdio` | ~69 | |
| `m_mpif90` | ~58 | |
| `m_realkinds` | ~40 | |
| `m_ioutil` | ~34 | |
| `m_inpak90` | ~29 | `GOCART_GridComp/CO2_GridComp/CO2_GridCompMod.F90`; `Chem_Shared/Chem_SettlingCalculator.F90`, `SmokeEmissionMod.F90` |
| `m_mpout`, `m_StrTemplate` | ~27, ~26 | |
| `m_mall`, `m_chars`, `m_List` | 22, 21, 17 | |
| `m_zeit`, `m_SortingTools`, `m_FileResolv`, `m_fpe` | 7, 5, 5, 4 | |
| `m_GrADS` | 1 | |

- GMAO_Shared siblings: `GMAO_hermes` (31 files use `m_die`/`m_stdio`/`m_inpak90`/`m_mpout`/`m_StrTemplate`), `GMAO_transf` (35), `GMAO_etc` (2), `GMAO_stoch` (3); `GMAO_transf/CMakeLists.txt` links `GMAO_eu`.
- Chemistry: `CMakeLists.txt` of `Chem_Base` (`DEPENDENCIES GMAO_mpeu MAPL ...`, `Chem_Base/CMakeLists.txt`, `esma_add_library`), `Chem_Shared`, `@GAAS`, `@ACHEM`, `@TR`, `@StratChem`, `@CARMA`, `@MAM`, GOCART and its `CO`, `CO2`, `CH4`, `CFC`, `Rn`, `O3` gridded components.
- Other: `GEOSphysics_GridComp`, `GEOSmoist_GridComp`, `GEOSsurface_GridComp`, `GEOSmkiau_GridComp` CMakeLists reference `GMAO_mpeu`/`GMAO_eu`.
- MAPL: `MAPL_cfio/ESMF_CFIOMod.F90` and `ESMF_CFIOUtilMod.F90` use mpeu modules; `MAPL/base/tests/utDistIO.F90` and `MAPL/Tests/ExtDataRoot_GridComp.F90` use `m_*` modules.
- WAVEWATCH III (`@ww3`): the ~55 `use` hits are all under `@ww3/regtests/ww3_tp2.14/input/oasis3-mct/lib/mct/` (22 in `mpeu/`, a vendored copy of mpeu with its own `m_die.F90` etc.; 28 in the MCT sources in `mct/`; 4 in `examples/`), not this library; no ww3 CMake reference to `GMAO_mpeu` was found.
- *Unverified:* use by GEOSadas or other fixtures.

## 5. Known traps

- Two targets, same sources: a program must not link both `GMAO_mpeu` and `GMAO_eu` (the latter provides `mpi_*` stub symbols). *Inferred* from the build; not tested.
- Duplicate module names: the vendored OASIS3-MCT `mpeu` under `@ww3` defines `m_die`, `m_stdio`, etc. with the same names.
- Files named `m_IndexBin_*.F90` are templates for three types, not one generic module.
- `m_mpif.F` and `m_random.F` are fixed-form `.F`; `m_FileResolv` is `.f90` but preprocessed.

## 6. Observations (noted, not fixed)

- `m_StringLinkedList.F90` is in the directory but not in `srcs`; no user was checked.
- `About..` is empty; `README` is a 2005 change log, not a description.
- `m_datetime` has no `use` in the fixture (grep); `m_GrADS` has one.
- `GMAO_eu` is built from `GMAO_mpeu` sources, so `m_mpif90` and `m_mpif` compile in both configurations (stub `mpi0` vs real MPI); how they differ at link time was not examined.
