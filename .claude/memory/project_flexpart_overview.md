---
name: FLEXPART project overview
description: High-level facts about the FLEXPART codebase structure, key modules, and known architectural issues
type: project
---

FLEXPART v11.0 is a Lagrangian Particle Dispersion Model written in Fortran 90, located in src/. 46 source files. Supports ECMWF and GFS/NCEP meteorology, OpenMP parallelisation, NetCDF and binary output.

**Key structural facts:**
- `par_mod` = compile-time parameters only
- `com_mod` = ~500-line global state store, used by almost every module (main refactoring target)
- `particle_mod` = derived type `particle` with coordinates/velocities subtypes + separate mass arrays
- `windfields_mod` = raw eta-level met fields; `getfields_mod` = height-level processed fields
- Main loop in `timemanager_mod::timemanager`; initialisation in `FLEXPART.f90` (not yet a module)
- Compiled with `-DETA` for eta-coordinate mode, otherwise metre-height mode

**Why:** ARCHITECTURE.md written at user request (2026-04-03) as full reference for program flow, variable meanings, dependency graph, and refactoring roadmap.

**How to apply:** When suggesting changes, refer user to ARCHITECTURE.md. Refactoring suggestions should follow the 8-point plan in that document.
