# GMAO_Shared

**Audience:** coding agents and developers working on the GMAO_Shared repository.
**Describes:** the MAPL2 line, branch `main`.
**Verified against:** `GMAO_Shared` @ `b0ebdca` (2026-11-06).
**Placement:** as in the GEOSgcm fixture (`components.yaml` @ `0f2e0fb`): key `GMAO_Shared`, checked out at `src/Shared/@GMAO_Shared`.
Other fixtures (for example GEOSadas, GEOSldas) may place this repository
elsewhere; paths inside this repository do not change.

## What this repository is
A set of Fortran libraries, scripts and tools shared by GEOS components
(message passing utilities, gridded file I/O, grid transforms, shared physics
utilities, the CICE4 sea-ice library, and others). It holds no gridded
components. Eight directories are documented at library level; the rest of
the library list in `CMakeLists.txt` is absent from the GEOSgcm checkout
(sparse checkout).

## Place in the GEOSgcm fixture
Child of the GEOSgcm fixture: see [`../../../AGENTS.md`](../../../AGENTS.md) and
[`../../../docs/ARCHITECTURE.md`](../../../docs/ARCHITECTURE.md). The fixture
checks it out with the sparse pattern `config/GMAO_Shared.sparse`.
Nested repositories (not documented here):
- `@GEOS_Util`: separate repository (documents come later).
- `@GMAO_perllib`: separate repository.

## Documents
| Document | Content |
|---|---|
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | Directory map, build structure, library dependency graph, users, cross-library traps |
| [`docs/decisions.md`](docs/decisions.md) | Visible structural decisions |
| [`GEOS_Shared/ARCHITECTURE.md`](GEOS_Shared/ARCHITECTURE.md) | Shared physics utilities |
| [`GMAO_mpeu/ARCHITECTURE.md`](GMAO_mpeu/ARCHITECTURE.md) | Message passing and environment utilities |
| [`GMAO_gfio/ARCHITECTURE.md`](GMAO_gfio/ARCHITECTURE.md) | Gridded file I/O (built as r4 and r8) |
| [`GMAO_hermes/ARCHITECTURE.md`](GMAO_hermes/ARCHITECTURE.md) | Dynamics-state library and tools |
| [`GMAO_transf/ARCHITECTURE.md`](GMAO_transf/ARCHITECTURE.md) | Grid transforms |
| [`GMAO_etc/ARCHITECTURE.md`](GMAO_etc/ARCHITECTURE.md) | Scripts and utilities (`esma_mpirun`, `tick`, ...) |
| [`GMAO_stoch/ARCHITECTURE.md`](GMAO_stoch/ARCHITECTURE.md) | Stochastic physics patterns |
| [`LANL_Shared/ARCHITECTURE.md`](LANL_Shared/ARCHITECTURE.md) | CICE4 sea-ice library |

## Rules for agents
- Note bugs and oddities; do not fix them unless asked.
- Keep the documents current when you change the code they describe.
- The documents are library-level: purpose, contents, users. Read the source
  for routine detail.
- Check which library defines a module before editing: several module names
  occur in more than one library (see `docs/ARCHITECTURE.md`, cross-library traps).

## Build and test
- `CMakeLists.txt` adds `GEOS_Shared` and `GMAO_hermes` always; with option
  `HERMES_LIGHT` OFF (default) it also adds the other subdirectories.
  `GMAO_gfio` is added twice (`GMAO_gfio_r4`, `GMAO_gfio_r8`).
- Listed subdirectories that do not exist in a checkout are skipped by
  `esma_add_subdirectories`; `GMAO_gfio` and `gmao2ioda` are guarded by `EXISTS`.
- Tests: `GMAO_etc/unittests` exists; `GMAO_gfio` has unbuilt test sources
  (`gfio_test.F`, `testcall.c`). No repository-level test target in `CMakeLists.txt`.
- `.circleci` and `.github` hold continuous-integration configuration.
