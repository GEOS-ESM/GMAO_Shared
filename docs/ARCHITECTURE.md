# GMAO_Shared: architecture

**Audience:** coding agents and developers working on the GMAO_Shared repository.
**Describes:** the MAPL3 line, branch `release/MAPL-v3`.
**Verified against:** `GMAO_Shared` @ `595999b` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `f50622d`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

Repository-level document: how the libraries fit together. Each library has
its own library-level document; read the source for routine detail.

> **Build status on this line:** restricted build. On `release/MAPL-v3`, `CMakeLists.txt` sets `HERMES_LIGHT ON`
> and adds only `GMAO_mpeu`, `GMAO_hermes` (light) and `GEOS_Shared` (comment: "not yet ported to MAPL3 - only build
> GMAO_mpeu, GMAO_hermes (light), and GEOS_Shared for FVdycoreCubed_GridComp"); the parent `src/Shared/CMakeLists.txt` says the
> same ("GMAO_Shared is included but restricted to HERMES_LIGHT"). Not built, not ported: `GMAO_gfio`, `GMAO_transf`,
> `GMAO_etc`, `GMAO_stoch`, `LANL_Shared` and every other library in the commented-out list; their code is the MAPL2 code,
> kept for a later port. Each library document states its own status.

## 1. Placement
- Parent: the GEOSgcm fixture (`AGENTS.md` (`../../../../AGENTS.md` in the fixture),
  `docs/ARCHITECTURE.md` (`../../../../docs/ARCHITECTURE.md` in the fixture)); the fixture lists
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
- On this line the file then executes `set (HERMES_LIGHT ON)` and calls `esma_add_subdirectories` for `GMAO_mpeu`,
  `GMAO_hermes`, `GEOS_Shared` (in that order), unconditionally. The `if ( NOT HERMES_LIGHT )` block is commented out
  and never runs. The text below describes that commented block, as on the MAPL2 line where it is live.
- MAPL2 line: always added `GEOS_Shared`, `GMAO_hermes`.
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
- With `HERMES_LIGHT` on (always, on this line), `GMAO_hermes` builds a reduced library (`m_set_eta.F90`, `shared_topo_remap.F90`,
  `m_topo_remap.F90`, `m_dyn.f90`, `m_const.f90`; dependency `GMAO_mpeu`; `write_eta.x` is built in both cases; the MAPL2 line's reduced
  library has only the first two and no dependencies); `GMAO_gfio_r8`, which `GMAO_hermes` otherwise links, is not added.
- Other fixtures: *Unverified:* whether any build file outside this repository sets `HERMES_LIGHT`; on this line it is set in `CMakeLists.txt` here.

### Directories listed but absent from this checkout
Absent: `GMAO_pilgrim`, `GEOS_Pert`, `GMAO_mfhdf3`, `GMAO_gfioeos`, `GMAO_ods`,
`GMAO_psas`, `GMAO_gems`, `GMAO_iret`, `GMAO_pyobs`, `GMAO_radmon`,
`GMAO_ncdiag`, `arpack`, `pnagpack`, `gmao2ioda`.
- `config/GMAO_Shared.sparse` (fixture) starts with `/*` and then excludes
  `GMAO_ncdiag`, `GMAO_ods`, `GMAO_psas`, `GMAO_pyobs`, `GEOS_Pert`, `arpack`,
  `pnagpack`, `GMAO_pilgrim`, `gmao2ioda` (patterns `!/<dir>`). These nine are
  absent because of the sparse checkout.
- `GMAO_mfhdf3`, `GMAO_gfioeos`, `GMAO_gems`, `GMAO_iret`, `GMAO_radmon` are
  not in the sparse file and are not in the repository tree at `595999b`
  (`git ls-tree HEAD`); `GMAO_pyobs` is excluded by the sparse file and is also
  not in the tree.

## 4. Dependency graph (CMake `DEPENDENCIES`)
On this line only `GMAO_mpeu`, `GMAO_hermes` (light: `GMAO_mpeu` only) and `GEOS_Shared` are built; the rest of the graph is the MAPL2 code's.
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
The table is carried over from the MAPL2-line documents. On this line the built consumer in the `GMAO_Shared` `CMakeLists.txt` comment is `FVdycoreCubed_GridComp`, which links `GEOS_Shared` and `GMAO_hermes`; `GMAO_stoch`, `GMAO_transf`, `GMAO_gfio`, `GMAO_etc`, `LANL_Shared` are not built.
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
- Parent: GEOSgcm (`../../../../docs/ARCHITECTURE.md` in the fixture).
- Decisions: [`decisions.md`](decisions.md).

## 8. Differences from the MAPL2 line
- Own changes (`git diff origin/main...origin/release/MAPL-v3`): `CMakeLists.txt` (restricted build), `GMAO_hermes/CMakeLists.txt` (light source list),
  nine `GEOS_Shared` files (MAPL3 API; see the `GEOS_Shared` document), `GMAO_hermes/lcv2prs.F90`, `LANL_Shared/CICE4/drivers/cice4/ice_constants.F90`
  (not built; `MAPL_ConstantsMod` to `MAPL_Constants`), and `.circleci`/`.github` configuration.
- Only on `main`: `GEOS_Shared` (`compute_radar_reflectivity.F90`, `surfacelayer.F90`), `GMAO_etc` (`esma_mpirun`, `pyrob`, `rmshmkeyhere.sh`, `RmShmKeys_sshmpi.csh`), `GMAO_stoch` (`mod_param.F90`, `glats.f`), `GMAO_ods/CMakeLists.txt`.

## 9. Observations (noted, not fixed)
- Fixture `docs/ARCHITECTURE.md` (repository table) lists only `GEOS_Shared`, `GMAO_mpeu`, `GMAO_etc`, `LANL_Shared/CICE4` for this repository; the checkout also has `GMAO_gfio`, `GMAO_hermes`, `GMAO_transf`, `GMAO_stoch`.
- Six directories named in `CMakeLists.txt` are not in the repository tree (`GMAO_mfhdf3`, `GMAO_gfioeos`, `GMAO_gems`, `GMAO_iret`, `GMAO_radmon`, `GMAO_pyobs`; section 3).
