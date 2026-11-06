# GMAO_Shared: architecture

**Audience:** coding agents and developers working on the GMAO_Shared repository.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `GMAO_Shared` @ `b0ebdca` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Repository-level document: how the libraries fit together. Each library has
its own library-level document; read the source for routine detail.

## 1. Placement
- Parent: the GEOSgcm fixture ([`AGENTS.md`](../../../../AGENTS.md),
  [`docs/ARCHITECTURE.md`](../../../../docs/ARCHITECTURE.md)); the fixture lists
  it with `sparse: config/GMAO_Shared.sparse` (fixture `docs/ARCHITECTURE.md`, repository table).
- Nested repositories inside this directory: `@GEOS_Util` and `@GMAO_perllib`
  (separate repositories, not documented here).

## 2. Directory map
| Directory | Kind | Document |
|---|---|---|
| `GEOS_Shared` | shared physics utilities (Fortran library) | [link](../GEOS_Shared/ARCHITECTURE.md) |
| `GMAO_mpeu` | message passing and environment utilities | [link](../GMAO_mpeu/ARCHITECTURE.md) |
| `GMAO_gfio` | gridded file I/O, built as r4 and r8 | [link](../GMAO_gfio/ARCHITECTURE.md) |
| `GMAO_hermes` | dynamics-state library and executables | [link](../GMAO_hermes/ARCHITECTURE.md) |
| `GMAO_transf` | grid transforms | [link](../GMAO_transf/ARCHITECTURE.md) |
| `GMAO_etc` | scripts and small utilities | [link](../GMAO_etc/ARCHITECTURE.md) |
| `GMAO_stoch` | stochastic physics patterns | [link](../GMAO_stoch/ARCHITECTURE.md) |
| `LANL_Shared` | CICE4 sea-ice library | [link](../LANL_Shared/ARCHITECTURE.md) |
| `@GEOS_Util`, `@GMAO_perllib` | nested repositories | separate repositories |
| `.circleci`, `.github` | continuous-integration configuration | none |

Top-level files: `CMakeLists.txt`, `README.md`, `CHANGELOG.md`, `LICENSE*`, `COPYRIGHT`.

## 3. Build structure (`CMakeLists.txt`)
- Option `HERMES_LIGHT` (default `OFF`, "Build light-weighted GMAO_shared library") (`CMakeLists.txt`, `option`).
- Always added: `GEOS_Shared`, `GMAO_hermes`.
- Only if `HERMES_LIGHT` is `OFF`, `esma_add_subdirectories`: `GMAO_mpeu`,
  `GMAO_pilgrim`, `GMAO_etc`, `GEOS_Util`, `LANL_Shared`, `GMAO_perllib`,
  `GMAO_transf`, `GMAO_stoch`, `GEOS_Pert`, `GMAO_mfhdf3`, `GMAO_gfioeos`,
  `GMAO_ods`, `GMAO_psas`, `GMAO_gems`, `GMAO_iret`, `GMAO_pyobs`,
  `GMAO_radmon`, `GMAO_ncdiag`, `arpack`, `pnagpack`.
- `GMAO_gfio` is a special case: added twice, as `GMAO_gfio_r4` and
  `GMAO_gfio_r8`, with `GMAO_gfio_r4` depending on `GMAO_gfio_r8`; guarded by `EXISTS`.
- `gmao2ioda` (scripts) is installed to `bin` if the directory exists.
- `esma_add_subdirectories` adds only directories that exist; a missing one is
  skipped with a debug message "possibly sparse checkout"
  (`@cmake/esma_support/esma_add_subdirectories.cmake` in the fixture).
- With `HERMES_LIGHT` on, `GMAO_hermes` builds a reduced library (`m_set_eta.F90` and `shared_topo_remap.F90`, no
  dependencies; `write_eta.x` is built in both cases); `GMAO_gfio_r8`, which `GMAO_hermes` otherwise links, is not added.
- No other build file in the fixture sets `HERMES_LIGHT` (grep of fixture `CMakeLists.txt`/`*.cmake`).

### Directories listed but absent from this checkout
Absent: `GMAO_pilgrim`, `GEOS_Pert`, `GMAO_mfhdf3`, `GMAO_gfioeos`, `GMAO_ods`,
`GMAO_psas`, `GMAO_gems`, `GMAO_iret`, `GMAO_pyobs`, `GMAO_radmon`,
`GMAO_ncdiag`, `arpack`, `pnagpack`, `gmao2ioda`.
- `config/GMAO_Shared.sparse` (fixture) starts with `/*` and then excludes
  `GMAO_ncdiag`, `GMAO_ods`, `GMAO_psas`, `GMAO_pyobs`, `GEOS_Pert`, `arpack`,
  `pnagpack`, `GMAO_pilgrim`, `gmao2ioda` (patterns `!/<dir>`). These nine are
  absent because of the sparse checkout.
- `GMAO_mfhdf3`, `GMAO_gfioeos`, `GMAO_gems`, `GMAO_iret`, `GMAO_radmon` are
  not in the sparse file and are not in the repository tree at `b0ebdca`
  (`git ls-tree HEAD`); `GMAO_pyobs` is excluded by the sparse file and is also
  not in the tree.

## 4. Dependency graph (CMake `DEPENDENCIES`)
```
GMAO_mpeu  -> MPI::MPI_Fortran          (GMAO_eu: same sources -> mpi0)
GMAO_gfio  -> GMAO_mpeu, esmf, NetCDF   (targets GMAO_gfio_r4, GMAO_gfio_r8)
GMAO_hermes-> GMAO_mpeu, GMAO_gfio_r8   (also include dirs of GEOS_Shared;
                                         executables link GEOS_Shared, MAPL)
GMAO_transf-> GMAO_mpeu, GMAO_hermes, NCEP_sp_r4i4, NCEP_bufr_r4i4,
              NCEP_w3_r8i4, NCEP_sfcio, NCEP_sigio
GMAO_stoch -> GMAO_mpeu, GMAO_gfio_r8, GMAO_transf, MAPL, NCEP_sfcio,
              NCEP_sp_r4i4, ESMF, NetCDF
GEOS_Shared-> MAPL, ESMF, NetCDF
CICE4 (LANL_Shared) -> MAPL
GMAO_etc   -> add_dependencies on GMAO_perllib (no library link)
```
(`<lib>/CMakeLists.txt`, `esma_add_library ... DEPENDENCIES`; `GMAO_transf` has a
commented older line and a live line at the end; `LANL_Shared/CMakeLists.txt`
only calls `add_subdirectory (CICE4)`; the `MAPL` dependency is in `LANL_Shared/CICE4/CMakeLists.txt`.)
`NCEP_*` come from the separate `NCEP_Shared` repository.

## 5. Who uses which library (GEOSgcm fixture)
| Library | Users (from the library documents) |
|---|---|
| `GEOS_Shared` | Moist, Surface, Turbulence, Radiation, Chemistry components (`GEOS_UtilsMod`); `GEOS_GcmGridComp`; `GMAO_hermes` |
| `GMAO_mpeu` | widely (`m_die` in about 150 files, count from the `GMAO_mpeu` document); Chemistry (`m_inpak90`); `GMAO_hermes`, `GMAO_transf` |
| `GMAO_gfio` | `GMAO_hermes`, `GMAO_stoch` (link), `GMAO_transf`; Chemistry (`r4`); `@GEOS_Util/post` |
| `GMAO_hermes` | `@GEOS_Util/post` programs; `GMAO_transf`; `GMAO_stoch`; `FVdycoreCubed_GridComp` (`m_set_eta`, `m_topo_remap`) |
| `GMAO_transf` | `GMAO_stoch` (link only); no other user found in the fixture |
| `GMAO_etc` | `gcm_setup`, `scm_setup` and other setup scripts (`esma_mpirun`, `tick`); FVdycoreCubed scripts; `@GEOS_Util` scripts |
| `GMAO_stoch` | `GEOSphysics_GridComp` (`use stoch_module`) |
| `LANL_Shared` (CICE4) | `GEOSogcm`, `GEOSseaice`, `GEOSdataseaice`, `GEOSCICEDyna`, `GEOSMITDyna`, `GEOSsaltwater` CMake files; `GEOSCICEDyna` Fortran |

Use from GEOSadas or other fixtures: *Unverified:*.

## 6. Cross-library traps
- `mpeu` is vendored a second time under `@ww3` (OASIS3-MCT) with the same module names (`m_die`, `m_stdio`, ...) ([GMAO_mpeu](../GMAO_mpeu/ARCHITECTURE.md)).
- `GMAO_mpeu` and `GMAO_eu` are two targets from the same sources; do not link both ([GMAO_mpeu](../GMAO_mpeu/ARCHITECTURE.md)).
- `mersenne_twister` exists in `GMAO_stoch`, in the MOM6 plug's `@stochastic_physics`, and in `@NCEP_Shared/NCEP_w3` ([GMAO_stoch](../GMAO_stoch/ARCHITECTURE.md)).
- `CICE4` (`LANL_Shared`) and `cice6` (`CICE_GEOSPlug`) both provide `ice_*` modules; name-based greps also hit `@cice6` and `@sis2` ([LANL_Shared](../LANL_Shared/ARCHITECTURE.md)).
- `dynamics_lattice_module` exists in `GMAO_transf` and `@GEOS_Util/post`; `util` exists in `GMAO_transf` and MOM `gotm-4.0` ([GMAO_transf](../GMAO_transf/ARCHITECTURE.md)).
- `tick` exists twice: `GMAO_etc/tick` and `@GEOS_Util/post/tick`; `m_tick` module is in `GMAO_hermes` ([GMAO_etc](../GMAO_etc/ARCHITECTURE.md), [GMAO_hermes](../GMAO_hermes/ARCHITECTURE.md)).
- Module names `m_dyn`, `m_const` and similar occur in several libraries; check which one defines the module in use ([GMAO_hermes](../GMAO_hermes/ARCHITECTURE.md)).
- `GMAO_gfio` precision: consumers choose `_r4` or `_r8`; the unsuffixed include directory is a symlink to `_r8`; `GMAO_gfioeos` is a different name ([GMAO_gfio](../GMAO_gfio/ARCHITECTURE.md)).
- `GMAO_transf` library links `GMAO_mpeu` but its executables link `GMAO_eu`; NCEP precision differs between library and executables ([GMAO_transf](../GMAO_transf/ARCHITECTURE.md)).
- `GEOS_Shared` file names differ from module names (`GEOS_Utilities.F90` is `GEOS_UtilsMod`); `.code` files are included text ([GEOS_Shared](../GEOS_Shared/ARCHITECTURE.md)).

## 7. Children and parent
- Children: [GEOS_Shared](../GEOS_Shared/ARCHITECTURE.md), [GMAO_mpeu](../GMAO_mpeu/ARCHITECTURE.md), [GMAO_gfio](../GMAO_gfio/ARCHITECTURE.md), [GMAO_hermes](../GMAO_hermes/ARCHITECTURE.md), [GMAO_transf](../GMAO_transf/ARCHITECTURE.md), [GMAO_etc](../GMAO_etc/ARCHITECTURE.md), [GMAO_stoch](../GMAO_stoch/ARCHITECTURE.md), [LANL_Shared](../LANL_Shared/ARCHITECTURE.md).
- Parent: [GEOSgcm](../../../../docs/ARCHITECTURE.md).
- Decisions: [`decisions.md`](decisions.md).

## 8. Observations (noted, not fixed)
- Fixture `docs/ARCHITECTURE.md` (repository table) lists only `GEOS_Shared`, `GMAO_mpeu`, `GMAO_etc`, `LANL_Shared/CICE4` for this repository; the checkout also has `GMAO_gfio`, `GMAO_hermes`, `GMAO_transf`, `GMAO_stoch`.
- Six directories named in `CMakeLists.txt` are not in the repository tree (`GMAO_mfhdf3`, `GMAO_gfioeos`, `GMAO_gems`, `GMAO_iret`, `GMAO_radmon`, `GMAO_pyobs`; section 3).
