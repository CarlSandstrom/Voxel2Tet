# Voxel2Tet Refactor Plan

## Scope

Track refactor tasks and completion history in one place.

This plan assumes:
- `src/lib` is the public API/library implementation.
- `src/bin` contains example applications that consume the API.

## Refactor Goals

1. Keep example programs in `src/bin` working without API-breaking changes.
2. Improve internal maintainability by reducing size/coupling of orchestration code.
3. Make library behavior robust for embedding by fixing ownership and error handling.
4. Preserve algorithmic behavior while improving structure and boundaries.

## Target Architecture (No Public API Break in Phase 1)

- Keep `Voxel2TetClass` as the compatibility facade and primary API entrypoint.
- Introduce internal pipeline modules under `src/lib/Pipeline/`:
  - `InputLoader` (importer selection + input validation)
  - `SurfaceBuilder` (surface/edge extraction)
  - `SmoothingPipeline` (surface/edge smoothing + fixed-node updates)
  - `TetPipeline` (TetGen execution + material remapping)
  - `ExportPipeline` (export routing/format handling)
- Keep `MeshData` focused on mesh data + mesh operations (move exporter selection out).

## Planned Tasks

- [ ] `R-004` Move surface/edge discovery into `SurfaceBuilder`.
- [ ] `R-005` Move smoothing logic into `SmoothingPipeline`.
- [ ] `R-006` Move tetrahedralization flow into `TetPipeline`.
- [ ] `R-007` Move export routing from `MeshData`/facade into `ExportPipeline`.
- [ ] `R-008` Replace owning raw pointers with `std::unique_ptr` in refactored flows.
- [ ] `R-009` Remove `exit(...)` from library code paths; use exceptions/status propagation.
- [ ] `R-010` Modernize CMake target boundaries (target includes/libs, reduce global state).
- [ ] `R-011` Add/expand regression coverage using existing `src/bin` examples.

## Migration Phases

1. **Phase A (Scaffold, no behavior change):** add pipeline classes and delegate from `Voxel2TetClass`.
2. **Phase B (Responsibility transfer):** move logic module-by-module (`InputLoader` -> `SurfaceBuilder` -> `SmoothingPipeline` -> `TetPipeline` -> `ExportPipeline`).
3. **Phase C (Safety hardening):** ownership cleanup (`unique_ptr`) and library-safe error model.
4. **Phase D (Build/test hardening):** CMake cleanup and regression verification.

## In Progress

- None.

## Completed Tasks

- [x] `R-000` Set up refactor tracking workflow (`AGENTS.md` + this tracker).  
  Started: `2026-03-06`  
  Completed: `2026-03-06`
- [x] `R-001` Audit current architecture and identify refactor boundaries.  
  Started: `2026-03-06`  
  Completed: `2026-03-06`
- [x] `R-002` Create `Pipeline` module skeletons and wire them through `Voxel2TetClass` facade.  
  Started: `2026-03-06`  
  Completed: `2026-03-06`
- [x] `R-UNPLANNED-001` Create peripheral docs folder and move refactor tracker into it.  
  Started: `2026-03-06`  
  Completed: `2026-03-06`
- [x] `R-UNPLANNED-002` Enforce use of coding standards document in agent workflow.  
  Started: `2026-03-06`  
  Completed: `2026-03-06`
- [x] `R-UNPLANNED-003` Fix VTK header resolution in CMake (`vtkSmartPointer.h` include failure).  
  Started: `2026-03-06`  
  Completed: `2026-03-06`
- [x] `R-012` Run validation build and document behavior/performance parity.  
  Started: `2026-03-13`  
  Completed: `2026-03-13`  
  Notes: See `docs/peripheral/R012_VALIDATION_2026-03-13.md` (CLI Dream3D segfault reproduced in both `build_check` and `build`).
- [x] `R-013` Investigate and fix CLI segmentation fault in `Voxel2Tet` Dream3D flow (`ExampleInput/16x16x16.dream3d`).  
  Started: `2026-03-13`  
  Completed: `2026-03-13`  
  Notes: Added chord-topology guard in `MeshManipulations::CheckCoarsenChord`; Dream3D CLI case now exits `0` and exports output files.
- [x] `R-003` Move import logic from `Voxel2TetClass` into `InputLoader`.  
  Started: `2026-03-13`  
  Completed: `2026-03-13`  
  Notes: `InputLoaderPipeline` now owns importer selection, file loading, callback importer construction, and optional voxel cutout setup.
- [x] `R-014` Add a Codex skill shortcut for "what is next in the refactor plan?".  
  Started: `2026-03-13`  
  Completed: `2026-03-13`  
  Notes: Added `.agents/skills/next-refactor/SKILL.md` to provide a reusable trigger for next-task lookup.
- [x] `R-015` Add skill UI metadata and troubleshooting guidance for `next-refactor` triggering.  
  Started: `2026-03-13`  
  Completed: `2026-03-13`  
  Notes: Added `.agents/skills/next-refactor/agents/openai.yaml` and documented reload conditions for skill discovery.

## Activity Log

- `2026-03-06` - `R-000` - Created persistent refactor tracker and Codex update rules.
- `2026-03-06` - `R-001` - Completed architecture review and defined API-first staged refactor approach.
- `2026-03-06` - `R-002` - Added `src/lib/Pipeline` skeleton modules and routed facade calls through pipeline delegates.
- `2026-03-06` - `R-UNPLANNED-001` - Moved `REFRACTOR_PLAN.md` to `docs/peripheral` and updated `AGENTS.md` path.
- `2026-03-06` - `R-UNPLANNED-002` - Added mandatory use of `docs/peripheral/CODING_STANDARDS.md` in `AGENTS.md`.
- `2026-03-06` - `R-UNPLANNED-003` - Added robust VTK include path discovery in CMake and verified full build success.
- `2026-03-13` - `R-012` - Started validation phase; full build succeeds after VTK include-path fix, parity write-up still pending.
- `2026-03-13` - `R-012` - Completed validation write-up with cross-build parity data; confirmed reproducible CLI segfault for Dream3D input and logged follow-up task `R-013`.
- `2026-03-13` - `R-013` - Fixed chord-collapse null/invalid neighbor handling in coarsening; verified Dream3D CLI path completes without segfault in `build_check`.
- `2026-03-13` - `R-003` - Moved import selection/load and callback importer setup into `InputLoaderPipeline`; `Voxel2TetClass` now finalizes returned importers.
- `2026-03-13` - `R-014` - Added repo-local Codex skill `next-refactor` for fast next-task lookup from `REFRACTOR_PLAN.md`.
- `2026-03-13` - `R-015` - Added `agents/openai.yaml` for `next-refactor` and captured trigger troubleshooting steps.
