# Codex Workflow For Refactor Tracking

This repository uses a single source of truth for refactor planning:
- `/home/carl/dev/Voxel2Tet/docs/peripheral/REFRACTOR_PLAN.md`

This repository also defines coding conventions in:
- `/home/carl/dev/Voxel2Tet/docs/peripheral/CODING_STANDARDS.md`

When working in this repo, Codex must keep that file current in the same turn as code changes.

Before implementing code changes, Codex should review and follow `CODING_STANDARDS.md`.

## Required update rules

1. If a new task is requested, add it to `Planned Tasks` with a stable task ID (`R-001`, `R-002`, ...).
2. If work starts on a task, move it to `In Progress` and set/update `Started`.
3. If a task is completed, mark it complete in `Completed Tasks` with `Completed` date.
4. For every completion, append one line in `Activity Log` (date, task ID, short result).
5. If work was completed that was not pre-planned, add it under `Completed Tasks` with ID `R-UNPLANNED-###`.
6. Keep entries concise and do not delete history unless the user explicitly asks.

## Date format

Use ISO date format: `YYYY-MM-DD`.
