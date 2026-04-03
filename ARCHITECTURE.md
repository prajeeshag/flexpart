# FLEXPART Architecture Reference

FLEXPART is a Lagrangian Particle Dispersion Model (LPDM) that computes the transport and dispersion of atmospheric tracers by tracking large numbers of computational particles through time-varying 3D wind fields. This document describes the program flow, module structure, key variables, dependency graph, and recommendations for refactoring toward a more modular and testable design.

---

## Table of Contents

1. [High-Level Program Flow](#1-high-level-program-flow)
2. [Module Inventory and Variables](#2-module-inventory-and-variables)
3. [Dependency Graph](#3-dependency-graph)
4. [Data Structures](#4-data-structures)
5. [Physics Subsystems](#5-physics-subsystems)
6. [I/O Subsystems](#6-io-subsystems)
7. [Parallelisation Model](#7-parallelisation-model)
8. [Refactoring Recommendations](#8-refactoring-recommendations)

---

## 1. High-Level Program Flow

```
FLEXPART (FLEXPART.f90)
│
├── read_options_and_initialise_flexpart
│   ├── readpaths              – parse pathnames file
│   ├── readcommand            – parse COMMAND namelist (run specs)
│   ├── alloc_random           – allocate per-thread random number arrays
│   ├── gasdev1 (×maxrand)     – fill random number pool
│   ├── readageclasses         – age-spectrum bin edges
│   ├── readavailable          – list of available wind field files
│   ├── detectformat           – ECMWF vs NCEP/GFS GRIB detection
│   ├── gridcheck_ecmwf|gfs    – read met-grid dimensions, hybrid A/B coefficients
│   ├── set_conv_top            – uppermost convection level (50 hPa)
│   ├── gridcheck_nest         – nested-domain grid check
│   ├── readoutgrid [+nest]    – output grid specification
│   ├── read_satellite_info    – satellite retrieval geometry
│   ├── readreceptors          – fixed receptor coordinates
│   ├── readlanduse            – 1200×600 global land-use inventory
│   ├── readreagents           – chemical reagent paths (LCM)
│   ├── initialise_particles   – dispatch to readreleases / readrestart / readpartpositions
│   ├── coordtrafo             – geo → grid coordinates for release points
│   ├── readdepo               – surface resistance tables
│   ├── alloc_drydepo / alloc_convect / alloc_getfields / alloc_interpol
│   ├── assignland             – map land-use to model grid (bilinear, ×10 refinement)
│   └── outgrid_init [+nest]   – output-cell volumes and areas
│
└── timemanager                              ← main computation loop
    │
    ├── [initialize output timing]
    │
    ├── LOOP over simulation time (step = lsynctime)
    │   │
    │   ├── getfields            – load/interpolate met fields into memory slots
    │   │   ├── readwind_ecmwf|gfs  – read GRIB → 3D arrays
    │   │   └── verttransform_ecmwf – eta→meter coordinate transform, cloud detection
    │   │
    │   ├── convmix              – Emanuel convective mixing (if lconvection=1)
    │   │
    │   ├── releaseparticles     – inject new particles from release points
    │   │   └── spawn_particle   – initialise particle_mod derived type
    │   │
    │   ├── [OMP PARALLEL DO over particles]
    │   │   └── advance          – integrate trajectory by one lsynctime
    │   │       ├── init_interpol      – find grid indices and weights
    │   │       ├── interpol_wind_*    – bilinear+vertical wind interpolation
    │   │       ├── turbulence_pbl     – Langevin equations in PBL
    │   │       ├── cbl (optional)     – skewed-PDF CBL turbulence
    │   │       ├── settling           – gravitational settling velocity
    │   │       ├── get_vdep_prob      – dry-deposition probability
    │   │       └── wetdepo            – below-cloud / in-cloud scavenging
    │   │
    │   ├── [concentration accumulation onto output grid]
    │   │   ├── drydepokernel    – accumulate dry deposition
    │   │   ├── calcfluxes       – accumulate gross mass fluxes
    │   │   └── receptor_conc    – accumulate receptor concentrations
    │   │
    │   ├── [at loutstep interval]
    │   │   ├── output_conc      – write gridded concentrations
    │   │   ├── receptoroutput   – write receptor file
    │   │   └── output_particles – particle dump (if ipout>0)
    │   │
    │   └── [at loutrestart interval]
    │       └── output_restart   – write restart.bin
    │
    └── [deallocate all]
```

### Key timing variables in `timemanager`

| Variable | Meaning |
|---|---|
| `itime` | Current model time in seconds relative to `bdate` |
| `lsynctime` | Synchronisation interval — all particles step by this amount |
| `loutstep` | Interval between concentration output files |
| `loutaver` | Averaging window for each output file |
| `loutsample` | Sampling interval within the averaging window |
| `loutrestart` | Interval for restart file writing |
| `loutnext` | Next scheduled output time |

---

## 2. Module Inventory and Variables

### 2.1 `par_mod` — Compile-time Parameters

All entries are `parameter` (compile-time constants). Nothing is allocatable.

| Parameter | Value / Type | Description |
|---|---|---|
| `dp`, `sp` | `integer` | Double / single precision kind selectors |
| `dep_prec` | `= dp` | Precision used for deposition arrays |
| `lusekerneloutput` | `.true.` | Enable kernel smoother for concentration output |
| `lparticlecountoutput` | `.false.` | Output particle counts rather than concentrations |
| `numpath` | `4` | Number of I/O path entries |
| `pi`, `r_earth`, `r_air`, `ga` | `real` | π, Earth radius (m), gas const for dry air (J/kg/K), gravity (m/s²) |
| `cpa`, `kappa` | `real` | Specific heat (J/kg/K), Poisson exponent |
| `vonkarman`, `karman` | `real` | von Kármán constant (0.4) |
| `rgas`, `r_water` | `real` | Universal gas constant; specific gas constant for H₂O |
| `href` | `15.` m | Reference height for dry deposition |
| `hmixmin`, `hmixmax` | `100.`, `4500.` m | PBL height bounds |
| `d_trop`, `d_strat` | `50.`, `0.1` m²/s | Horizontal and stratospheric diffusivities |
| `fturbmeso` | `0.16` | Mesoscale turbulence scaling factor |
| `idiffnorm`, `idiffmax` | `10800`, `21600` s | Normal/maximum met-field interval |
| `minstep` | `1` s | Minimum integration time step |
| `maxnests` | `5` | Maximum number of nested domains |
| `maxcolumn` | `3000` | Max particles per atmospheric column (domain-fill) |
| `maxrand` | `6 000 000` | Size of Gaussian random number pool |
| `nclassunc` | `1` | Number of uncertainty classes |
| `maxtable` | `1000` | Max chemical species in lookup table |
| `numclass` | `13` | Number of land-use classes |
| `numpf` | `1` | Number of precipitation fields |
| `numwfmem` | `2` | Number of wind field time slots kept in memory |
| `maxndia` | `1` | Max number of particle diameter bins |
| `ncluster` | `5` | Number of trajectory clusters for plume output |
| `maxreagent` | `5` | Max chemical reagents (LCM) |
| `maxrecsample` | `2000` | Max receptors per sampling interval |
| `maxxOH`, `maxyOH`, `maxzOH` | `72`, `46`, `7` | OH-field grid dimensions |
| `bclr_a/b/c/e`, `bcls_a/b/c/e` | `real(:)` | Below-cloud scavenging polynomial coefficients (Wang et al. 2014) |
| `max_cloudthck`, `min_cloudthck` | `19000`, `50` m | Cloud thickness bounds |
| `switchnorth`, `switchsouth` | `75.`, `-75.` ° | Latitude thresholds for polar stereo projection |
| `icmv` | `-9999` | Missing value integer code |
| `ispeed` | `1` | Memory-vs-speed trade-off for terminated particles |
| Unit parameters | `unitpath`, `unitcommand`, … | Fortran unit numbers for all I/O files |

---

### 2.2 `com_mod` — Global Run-State Variables

This is the largest module; it acts as the global state store. All allocatable arrays depend on `nspec`, `numpoint`, `numreceptor`, or `numpart`.

#### Run identity and paths

| Variable | Type | Description |
|---|---|---|
| `path(numpath+2*maxnests)` | `character*120` | Directory paths for input/output |
| `length(numpath+2*maxnests)` | `integer` | Length of each path string |
| `pathfile` | `character(256)` | File containing the path names |
| `flexversion`, `flexversion_major` | `character(256)` | Full and major version strings |
| `gitversion` | `character(256)` | Git SHA + timestamp |
| `arg1`, `arg2` | `character(256)` | Command-line arguments |

#### Simulation time control

| Variable | Type | Description |
|---|---|---|
| `ibdate`, `ibtime` | `integer` | Begin date (YYYYMMDD) and time (HHMMSS) |
| `iedate`, `ietime` | `integer` | End date and time |
| `bdate`, `edate` | `real(dp)` | Julian begin/end dates |
| `ldirect` | `integer` | +1 forward, −1 backward |
| `ideltas` | `integer` | Total simulation length (s) |
| `lsynctime` | `integer` | Synchronisation step (s) |
| `loutstep`, `loutaver`, `loutsample` | `integer` | Output timing (s) |
| `lrecoutstep`, `lrecoutaver`, `lrecoutsample` | `integer` | Receptor output timing (s) |
| `loutrestart` | `integer` | Restart file interval (s) |

#### Model physics switches

| Variable | Type | Description |
|---|---|---|
| `lsubgrid` | `integer` | 1 = sub-grid topography on |
| `lconvection` | `integer` | 1 = Emanuel convection on |
| `lturbulence` | `integer` | 1 = turbulence on |
| `lagespectra` | `integer` | 1 = age spectra on |
| `lmesoscale_turb` | `logical` | Mesoscale turbulence (default off) |
| `cblflag` | `integer` | 1 = CBL skewed-PDF scheme |
| `lsettling` | `logical` | Gravitational settling |
| `lcw`, `lcwsum` | `logical` | Cloud water found in GRIB; sum clwc+ciwc |
| `lprecint` | `logical` | New precipitation interpolation scheme (#295) |
| `turbswitch` | `logical` | Markov chain formulation selector |

#### Output control

| Variable | Type | Description |
|---|---|---|
| `iout` | `integer` | 1=conc (ng/m³), 2=mix ratio (pptv), 3=both, 4=plume, 5=plume+conc |
| `ipout` | `integer` | Particle dump: 0=none, 1=each interval, 2=end only |
| `ipin` | `integer` | Restart mode: 0=none, 1=restart.bin, 2=netCDF, 3=IC, 4=both |
| `iflux` | `integer` | 1 = compute gross mass fluxes |
| `lnetcdfout` | `integer` | 1 = NetCDF output |
| `linversionout` | `integer` | 1 = one output file per release |
| `nested_output` | `integer` | 1 = also output nested grid |
| `sfc_only` | `integer` | 1 = surface layer only in grid_time files |
| `ind_source`, `ind_receptor` | `integer` | 1=mass units, 2=mixing ratio units |
| `linit_cond` | `integer` | Sensitivity to initial conditions (backward) |
| `mquasilag` | `integer` | 1 = quasi-Lagrangian condensed particle output |
| `mdomainfill` | `integer` | 1 = domain-filling mode |

#### Species and release data (all allocatable, size `nspec`)

| Variable | Type | Description |
|---|---|---|
| `nspec` | `integer` | Number of chemical species per release |
| `maxspec` | `integer` | Max species per release point |
| `species(nspec)` | `character(10)` | Species names |
| `numpoint` | `integer` | Number of release locations |
| `compoint(1001)` | `character*45` | Release-point comment/label |
| `decay(nspec)` | `real` | Radioactive decay constant (1/s) |
| `weta_gas`, `wetb_gas` | `real(nspec)` | Below-cloud gas scavenging A/B coefficients |
| `crain_aero`, `csnow_aero` | `real(nspec)` | Below-cloud aerosol scavenging (rain/snow) |
| `ccn_aero`, `in_aero` | `real(nspec)` | In-cloud CCN/IN scavenging |
| `reldiff`, `henry`, `f0` | `real(nspec)` | Dry deposition: relative diffusivity, Henry const, O₃ reactivity |
| `ri`, `rcl`, `rgs`, `rlu`, `rm` | `real(:,:)` | Surface resistances (species×landclass) |
| `dryvel(nspec)` | `real` | Constant dry deposition velocity (m/s) |
| `density`, `dquer`, `dsigma` | `real(nspec)` | Particle: density (kg/m³), mean diameter (m), log-std |
| `vset`, `schmi`, `fract` | `real(:,:)` | Settling velocity, Schmidt number, mass fraction per diameter bin |
| `vsetaver`, `cunningham`, `weightmolar` | `real(nspec)` | Averages and Cunningham slip factor |
| `Fn`, `Fs`, `ks1`, `ks2`, `kn2` | `real(nspec)` | Non-spherical shape drag parameters (Bagheri & Bonadonna 2016) |
| `ishape`, `orient` | `integer(nspec)` | Shape and orientation flags |
| `DEP`, `DRYDEP`, `WETDEP` | `logical` | Global deposition switches |
| `DRYDEPSPEC`, `WETDEPSPEC` | `logical(nspec)` | Per-species deposition switches |
| `LEMIS`, `LDECAY`, `CLREA` | `logical` | Emission, decay, chemistry switches |

#### Met-field memory management

| Variable | Type | Description |
|---|---|---|
| `memtime(numwfmem)` | `integer` | Validation time of each slot in memory (s) |
| `memind(3)` | `integer` | Slot pointer (avoids data shuffling) |
| `lwindinterv` | `integer` | Interval between wind fields currently in memory (s) |
| `numbnests` | `integer` | Number of nested input domains |
| `xglobal`, `sglobal`, `nglobal` | `logical` | Global domain flags |
| `southpolemap`, `northpolemap` | `real(9)` | Polar stereographic map parameters |

#### Grid and receptor

| Variable | Type | Description |
|---|---|---|
| `numxgrid`, `numygrid`, `numzgrid` | `integer` | Output grid dimensions |
| `dxout`, `dyout` | `real` | Output grid spacing (°) |
| `outlon0`, `outlat0` | `real` | Output grid lower-left corner |
| `numxgridn`, `numygridn` | `integer` | Nested output grid dimensions |
| `xreceptor`, `yreceptor`, `zreceptor` | `real(:)` | Receptor positions |
| `creceptor(:,:)` | `real` | Receptor concentrations |
| `numreceptor` | `integer` | Number of receptors |
| `landinvent(1200,600,6)` | `integer(1)` | Global land-use inventory (fixed 1°/12 resolution) |
| `z0(numclass)` | `real` | Roughness length per land-use class |

---

### 2.3 `windfields_mod` — Meteorological Input Fields

These arrays hold the 3D wind, thermodynamic, and cloud fields read from GRIB files, stored in `numwfmem` time slots.

| Variable | Dimensions | Description |
|---|---|---|
| `numbwf` | scalar | Total number of available wind field files |
| `wftime(:)` | `integer` | Time of each wind field (s from simulation start) |
| `wfname(:)`, `wfnamen(:)` | `character` | GRIB filenames for mother/nested domains |
| `uueta(:,:,:,:)` | `nx,ny,nuvz,numwfmem` | U-wind on eta levels (m/s) |
| `vveta(:,:,:,:)` | `nx,ny,nuvz,numwfmem` | V-wind on eta levels (m/s) |
| `wweta(:,:,:,:)` | `nx,ny,nwz,numwfmem` | Vertical wind (omega or w) on eta levels |
| `tteta(:,:,:,:)` | `nx,ny,nuvz,numwfmem` | Temperature on eta levels (K) |
| `pveta(:,:,:,:)` | `nx,ny,nuvz,numwfmem` | Potential vorticity on eta levels |
| `rhoeta(:,:,:,:)` | `nx,ny,nuvz,numwfmem` | Air density on eta levels (kg/m³) |
| `prseta(:,:,:,:)` | `nx,ny,nuvz,numwfmem` | Pressure on eta levels (Pa) |
| `drhodzeta(:,:,:,:)` | — | Vertical density gradient |
| `clwc`, `ciwc` | `nx,ny,nuvz,numwfmem` | Cloud liquid/ice water content (kg/kg) |
| `ctwc` | `nx,ny,numwfmem` | Total cloud water content |
| `icloudbot`, `icloudtop` | `nx,ny,numwfmem` | Cloud bottom/top level indices |
| `oro(:,:)` | `nx,ny` | Surface orography (m), time-invariant |
| `lsm(:,:)` | `nx,ny` | Land-sea mask (0–1), time-invariant |
| `excessoro(:,:)` | `nx,ny` | Sub-grid orography standard deviation |
| `etauvheight`, `etawheight` | `nx,ny,nuvz|nwz,numwfmem` | Heights of eta levels above sea level (m) |
| All above with `n` suffix | same with nest dimension | Nested-domain equivalents |

---

### 2.4 `particle_mod` — Particle Derived Types

```fortran
type :: coordinates
  real :: xlon   ! Longitude (grid units)
  real :: ylat   ! Latitude (grid units)
  real :: z      ! Height (m above sea level or AGL depending on kindz)
  real :: zeta   ! Eta coordinate (only used when compiled with -DETA)
end type

type :: velocities
  real :: u, v, w     ! Cartesian velocity components (m/s)
  real :: weta        ! Vertical velocity in eta coordinates
end type

type :: particle
  type(coordinates) :: pos   ! Current position
  type(velocities)  :: vel   ! Current turbulent velocity (Langevin state)
  integer :: npoint          ! Release point index
  integer :: nclass          ! Uncertainty class
  integer :: idt             ! Current integration time step (s)
  integer :: itra1           ! Current time (s from simulation start)
  integer :: itramem         ! Release time
  integer :: itrasplit       ! Next particle-split time
  logical :: alive           ! Whether particle is active
end type

type :: particlecount
  integer :: alive      ! Currently active particles
  integer :: spawned    ! Total particles ever created
  integer :: terminated ! Total particles removed
  integer :: allocated  ! Array size
end type
```

Per-particle allocatable arrays (separate from the derived type for memory layout reasons):

| Array | Dimensions | Description |
|---|---|---|
| `part(:)` | `maxpart` | Main particle array |
| `mass(:,:)` | `maxpart, nspec` | Current mass per species (kg) |
| `mass_init(:,:)` | `maxpart, nspec` | Initial mass (for normalisation) |
| `wetdeposit(:,:)` | `maxpart, nspec` | Cumulative wet deposition (kg) |
| `drydeposit(:,:)` | `maxpart, nspec` | Cumulative dry deposition (kg) |
| `prob(:,:)` | `maxpart, nspec` | Dry deposition probability (−) |
| `xscav_frac1(:,:)` | `maxpart, nspec` | Scavenged fraction at receptor locations |

---

### 2.5 `interpol_mod` — Interpolation State

Thread-private arrays updated by `init_interpol` before each interpolation call.

| Variable | Description |
|---|---|
| `ix`, `jy` | Integer grid indices (lower-left of enclosing cell) |
| `ixp`, `jyp` | `ix+1`, `jy+1` (with wraparound for global grids) |
| `indz`, `indzp` | Vertical level indices |
| `p1`, `p2`, `p3`, `p4` | Bilinear interpolation weights (sum to 1) |
| `dt1`, `dt2`, `dtt` | Temporal interpolation weights (two time slots) |
| `uprof(:)`, `vprof(:)`, `wprof(:)` | Vertical wind profiles at particle position |
| `usigprof(:)`, `vsigprof(:)`, `wsigprof(:)` | Wind standard deviation profiles |
| `rhoprof(:)`, `rhogradprof(:)` | Density and vertical density gradient profiles |

---

### 2.6 `turbulence_mod` — PBL Turbulence State (Thread-Private)

| Variable | Description |
|---|---|
| `ust` | Friction velocity u* (m/s) |
| `wst` | Convective velocity scale w* (m/s) |
| `ol` | Obukhov length (m); negative = unstable |
| `h` | PBL height (m) |
| `zeta` | z/h (normalised height within PBL) |
| `sigu`, `sigv`, `sigw` | Standard deviations of u, v, w turbulent fluctuations (m/s) |
| `tlu`, `tlv`, `tlw` | Lagrangian time scales for u, v, w (s) |

---

### 2.7 `outgrid_mod` — Output Grid and Accumulated Fields

| Variable | Dimensions | Description |
|---|---|---|
| `outheight(:)` | `numzgrid` | Upper boundaries of output layers (m) |
| `outheighthalf(:)` | `numzgrid` | Mid-levels of output layers (m) |
| `oroout(:,:)` | `numxgrid,numygrid` | Topographic height at output grid points |
| `area(:,:)` | `numxgrid,numygrid` | Grid-cell horizontal area (m²) |
| `volume(:,:,:)` | `numxgrid,numygrid,numzgrid` | Grid-cell volume (m³) |
| `grid(:,:,:,:,:,:)` | `nx,ny,nz,nspec,nage,npoint` | Concentration accumulation array |
| `gridsigma` | same | Uncertainty (σ) of concentrations |
| `wetgrid`, `drygrid` | `nx,ny,nspec,nage,npoint` | Accumulated deposition arrays |
| `flux(:,:,:,:,:,:,:)` | `nx,ny,nz,6,nspec,nage,npoint` | Mass flux across 6 cell faces |
| All above with `n` suffix | — | Nested grid equivalents |

---

### 2.8 `unc_mod` — Concentration Uncertainty Grids

| Variable | Description |
|---|---|
| `gridunc(:,:,:,:,:,:,:)` | 7-D array: (nx,ny,nz,nspec,nage,npoint,**nthreads**) — each OpenMP thread writes to its own slice to avoid race conditions |
| `griduncn` | Nested equivalent |
| `drygridunc`, `wetgridunc` | Dry/wet deposition uncertainty, also thread-split |

---

### 2.9 `conv_mod` — Convection State

| Variable | Description |
|---|---|
| `pconv(:)`, `phconv(:)` | Pressure arrays for convection scheme |
| `fmass(:,:)` | Up/down mass fluxes per level |
| `fmassfrac(:,:)` | Mass flux fractions |
| `cbaseflux(nx,ny,nspec)` | Cloud-base mass flux (kg/m²/s) |
| `cbasefluxn(nx,ny,nnest,nspec)` | Nested equivalent |
| `nconvlev` | Number of convective levels |
| `nconvtop` | Index of uppermost convective level |

---

### 2.10 `getfields_mod` — Processed Intermediate Fields

These are the height-coordinate fields derived from `windfields_mod` after `verttransform_ecmwf`.

| Variable | Dimensions | Description |
|---|---|---|
| `uuh(:,:,:,:)` | `nx,ny,nz,numwfmem` | U-wind at meter-coordinate levels |
| `vvh(:,:,:,:)` | same | V-wind |
| `wwh(:,:,:,:)` | same | Vertical velocity (m/s) |
| `pvh(:,:,:,:)` | same | Potential vorticity |
| `ttlev`, `qvlev`, `ulev`, `vlev`, `zlev` | `nz` | Level arrays for verttransform calculations |
| Nested equivalents `uuhn`, `vvhn`, `pvhn`, `wwhn` | — | Nested domain fields |

---

### 2.11 `point_mod` — Release Point Geometry

| Variable | Description |
|---|---|
| `ireleasestart(:)`, `ireleaseend(:)` | Release time window (s, relative to simulation start) |
| `npart(:)` | Number of particles per release point |
| `kindz(:)` | Height type: 1=AGL, 2=ASL, 3=pressure level |
| `xpoint1/2(:)`, `ypoint1/2(:)` | Release area corner coordinates (grid units after coordtrafo) |
| `zpoint1/2(:)` | Vertical release range |
| `xmass(:,:)` | Total mass released per species (kg) |
| `dx`, `dy`, `xlon0`, `ylat0` | Met grid spacing and origin |

---

### 2.12 `drydepo_mod` — Dry Deposition Fields

| Variable | Description |
|---|---|
| `xlanduse(nx,ny,numclass)` | Fractional land-use cover at each met-grid point |
| `xlandusen` | Nested equivalent |
| `vdep(nx,ny,nspec)` | Deposition velocity at each grid point (m/s) |
| `z0_drydep(nx,ny)` | Surface roughness length for deposition (m) |

---

### 2.13 `chemistry_mod` — Chemical Loss Fields

| Variable | Description |
|---|---|
| `nxCL`, `nyCL`, `nzCL` | Chemical loss field grid dimensions |
| `lonCL(:)`, `latCL(:)`, `altCL(:)` | Coordinates of CL field |
| `CL_field(:,:,:,:)` | Chemical loss rates at native resolution |
| `clfield_cur(:,:,:)` | Current interpolated CL field |
| `memCLtime(:)`, `memCLday` | Time management for CL fields |
| `jrate_average` | Daily-averaged actinic flux |

---

### 2.14 Utility Modules

| Module | Key contents |
|---|---|
| `date_mod` | `caldate`, `juldate`: Julian ↔ Gregorian conversion |
| `random_mod` | `gasdev1` (Box-Muller), `ran3`; per-thread `rannumb` arrays |
| `erf_mod` | Error function `erf()` and complementary `erfc()` |
| `qvsat_mod` | Saturation specific humidity formula |
| `cmapf_mod` | Cartographic map-factor routines (polar stereo, lat-lon) |
| `sort_mod` | Sorting algorithms |
| `coord_ecmwf_mod` | `z_to_zeta`, `zeta_to_z`, `update_z_to_zeta` — height ↔ eta conversion |
| `class_gribfile_mod` | `detectformat`; constants `GRIBFILE_CENTRE_ECMWF`, `GRIBFILE_CENTRE_NCEP` |
| `xmass_mod` | `xmasssave` — saved initial mass for normalisation |
| `totals_mod` | `tot_mass`, `chem_loss` accumulators for LCM diagnostics |
| `mean_mod` | Running mean of particle positions/properties |
| `flux_mod` | `flux(7D)`, `calcfluxes`, `fluxoutput`, `areaeast`, `areanorth` |

---

## 3. Dependency Graph

Arrows indicate `USE` (compile-time) dependencies. `par_mod` and `com_mod` are used by almost every module and are omitted from individual edges for clarity — assume every module listed implicitly depends on them.

```
par_mod  ←────────────────────── used by all modules
com_mod  ←────────────────────── used by all non-utility modules
    │
    ├── windfields_mod
    │       └── used by: getfields_mod, timemanager_mod, advance_mod,
    │                    verttransform_mod, conv_mod, interpol_mod
    │
    ├── particle_mod
    │       └── used by: timemanager_mod, advance_mod, output_mod,
    │                    emissions_mod, initialise_mod, restart_mod
    │
    ├── interpol_mod
    │       └── used by: advance_mod, turbulence_mod, drydepo_mod,
    │                    wetdepo_mod, settling_mod
    │
    ├── turbulence_mod
    │       └── used by: advance_mod, cbl_mod
    │
    ├── getfields_mod
    │       └── used by: timemanager_mod, verttransform_mod, conv_mod
    │
    ├── outgrid_mod  ←── unc_mod
    │       └── used by: timemanager_mod, output_mod, flux_mod
    │
    ├── output_mod
    │       ├── binary_output_mod
    │       ├── netcdf_output_mod  ←── receptor_netcdf_mod
    │       └── txt_output_mod
    │
    ├── point_mod ←── used by: FLEXPART.f90, timemanager_mod, initialise_mod
    ├── readoptions_mod  ←── used by: FLEXPART.f90
    ├── initialise_mod   ←── used by: FLEXPART.f90, timemanager_mod
    ├── restart_mod      ←── used by: FLEXPART.f90, timemanager_mod
    ├── conv_mod         ←── used by: FLEXPART.f90, timemanager_mod
    ├── drydepo_mod      ←── used by: FLEXPART.f90, timemanager_mod, advance_mod
    ├── wetdepo_mod      ←── used by: advance_mod, timemanager_mod
    ├── settling_mod     ←── used by: advance_mod
    ├── chemistry_mod    ←── used by: timemanager_mod
    ├── emissions_mod    ←── used by: timemanager_mod
    ├── receptor_mod     ←── used by: timemanager_mod
    ├── flux_mod         ←── used by: timemanager_mod
    ├── plume_mod        ←── used by: FLEXPART.f90, timemanager_mod
    └── date_mod, random_mod, erf_mod, qvsat_mod, cmapf_mod,
        sort_mod, coord_ecmwf_mod, class_gribfile_mod  (utility/leaf modules)
```

### Problematic dependencies

1. **`com_mod` circular state**: Almost every module `USE`s `com_mod` directly and accesses module-level globals. This makes it impossible to test any individual physics routine in isolation.
2. **`windfields_mod` ↔ `getfields_mod`**: Both hold different representations of the same data (eta-level vs. meter-level). There is no clean interface between them.
3. **`FLEXPART.f90` contains subroutines**: `read_options_and_initialise_flexpart` and `initialise_particles` are global subroutines inside the main program file, not in their own module, preventing reuse and testing.
4. **`par_mod` as compile-time config**: Several `parameter`s that should be runtime-configurable (e.g. `numpf`, `numwfmem`, `lusekerneloutput`) are compile-time constants, forcing recompilation for different configurations.

---

## 4. Data Structures

### 4.1 Coordinate Systems

FLEXPART can be compiled in two modes selected by the preprocessor flag `-DETA`:

| Mode | Horizontal | Vertical |
|---|---|---|
| **METER** (default) | Longitude/latitude in grid index units | Height in metres above sea level |
| **ETA** | Same | ECMWF hybrid sigma-pressure (η) levels, converted to height via `coord_ecmwf_mod` |

The met-fields are always stored in eta-levels in `windfields_mod`. `verttransform_ecmwf` converts them to metre-height levels and stores the result in `getfields_mod`. In ETA mode, the particle carries both `z` (metres) and `zeta` (eta) and the lookup `z_to_zeta` / `zeta_to_z` maintains consistency.

### 4.2 Wind Field Memory Layout

```
memind(1..numwfmem)  →  slot pointers (avoid copying)

uueta(ix, jy, kz, memslot)
       ↑    ↑   ↑     ↑
       |    |   |   1 or 2 (numwfmem=2 in serial)
       |    |   eta level index (1..nuvzmax)
       |    latitude index (0..nymax-1)
       longitude index (0..nxmax-1)
```

`getfields` manages loading: it always keeps the two nearest wind-field times bracketing the current `itime`, advancing the memory slots without copying arrays.

### 4.3 Output Grid Accumulation

Concentrations are accumulated as:

```
grid(ix, jy, kz, ispec, iage, ipoint)  +=  mass × weight / volume × dt
```

where `iage` is the age-spectrum bin, `ipoint` is the release-point index (if `ioutputforeachrelease=1`), and `weight` is a temporal sampling weight (0.5 at boundaries, 1.0 in interior of the averaging window).

For uncertainty, each OpenMP thread writes to `gridunc(:,:,:,:,:,:, omp_thread_id)` to avoid atomic operations, and the slices are summed at output time.

---

## 5. Physics Subsystems

### 5.1 Particle Trajectory Integration (`advance_mod`)

Each particle is advanced by `lsynctime` seconds using the Langevin stochastic differential equation. The time step `idt` is adaptive, limited by the minimum of:
- `ctl × TL` (Lagrangian time scale / control parameter)
- `ifine × idt_wind` (CFL-like constraint for vertical wind)
- `minstep`

Inside the PBL (`z < hmix`):
- `turbulence_pbl` computes `sigu`, `sigv`, `sigw`, `tlu`, `tlv`, `tlw`
- The Langevin equation is integrated: `du = −u/TL dt + σ dW`
- In CBL conditions (`ol < 0`, `cblflag=1`), `cbl_mod` uses a bi-Gaussian skewed PDF

Above the PBL, in the free troposphere:
- Horizontal diffusion uses `d_trop` (or `d_strat` above tropopause)
- Convective parametrization is applied via `convmix`

### 5.2 Convection (`conv_mod`)

Uses Emanuel's scheme (convect43c). Called once per wind-field read. The scheme:
1. Takes column profiles of T, q, u, v
2. Computes up/down mass fluxes at each level
3. Redistributes particle masses between levels according to `fmassfrac`

### 5.3 Dry Deposition (`drydepo_mod`)

Uses Wesely (1989) resistance analogy:
```
Vd = 1 / (Ra + Rb + Rc)
```
where `Ra` is aerodynamic, `Rb` is quasi-laminar, and `Rc` is surface resistance. `Rc` depends on land-use class, season, and species properties (`reldiff`, `henry`, `f0`). The probability of deposition per time step:
```
prob = 1 - exp(-Vd × dt / href)
```

### 5.4 Wet Deposition (`wetdepo_mod`)

Scavenging coefficient Λ is parameterised as:
- **Gas**: `Λ = weta_gas × P^wetb_gas`
- **Aerosol (below-cloud)**: polynomial in `ln(dp)` and `ln(P)` (Wang et al. 2014)
- **Aerosol (in-cloud)**: `ccn_aero` and `in_aero` parameters

Deposited fraction per time step: `f = 1 - exp(-Λ × dt)`.

### 5.5 Gravitational Settling (`settling_mod`)

Iterative solution for terminal fall velocity using:
- Stokes regime: `Vs = d²(ρp−ρa)g / (18η)`
- Non-Stokes: Clift & Guavin (1971) drag coefficient
- Non-spherical shapes: `Fn`, `Fs` parameters from Bagheri & Bonadonna (2016)
- Cunningham slip correction for sub-micron particles

---

## 6. I/O Subsystems

### 6.1 Meteorological Input

```
GRIB files  →  class_gribfile_mod (format detection)
            →  read_input.f90 (gridcheck, dimension reading)
            →  readwind_ecmwf.f90 (GRIB decode → windfields_mod arrays)
            →  verttransform_mod (eta → height, cloud detection)
            →  getfields_mod (time-interpolated profiles for advance)
```

### 6.2 Output Dispatch

`output_mod` dispatches to one of two backends based on `lnetcdfout`:

```
output_conc
    ├── lnetcdfout=1  →  netcdf_output_mod::concoutput_netcdf
    └── lnetcdfout=0  →  binary_output_mod::concoutput_bin

output_particles
    ├── lnetcdfout=1  →  netcdf_output_mod::partoutput_netcdf
    └── lnetcdfout=0  →  binary_output_mod::partoutput_bin

receptoroutput
    ├── lnetcdfout=1  →  receptor_netcdf_mod
    └── lnetcdfout=0  →  binary receptor files
```

### 6.3 Restart System (`restart_mod`)

Three rotating restart files (`restart1.bin`, `restart2.bin`, `restart3.bin`) are written. Reading converts z↔zeta as needed for the coordinate system in use. Restart is triggered every `loutrestart` seconds.

---

## 7. Parallelisation Model

### OpenMP

The particle loop in `timemanager` is parallelised with `!$OMP PARALLEL DO`:

```fortran
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ip, ...) SCHEDULE(DYNAMIC,chunk)
do ip = 1, numpart
  if (part(ip)%alive) call advance(ip, itime, ...)
end do
!$OMP END PARALLEL DO
```

Thread-safety mechanisms:
- **Random numbers**: each thread has its own index into the pre-filled `rannumb` pool (or its own seed in `random_mod`)
- **Output grid**: each thread writes to a private slice of `gridunc(:,:,:,:,:,:,ithread)`, merged at output
- **Thread-private variables**: turbulence state (`ust`, `wst`, `ol`, etc.) and interpolation indices (`ix`, `jy`, `p1`…`p4`) are all `!$OMP THREADPRIVATE`
- **Deposition**: `drydepo_mod` accumulates to `drygridunc(:,:,:,:,:,:,ithread)` similarly

The number of threads is read via `OMP_GET_MAX_THREADS()` at startup and stored in `numthreads`; grid-accumulation threads are capped at `maxthreadgrid` to limit memory overhead from the uncertainty array split.

### MPI

An MPI version exists (separate code path, `numwfmem=3`). The serial/OpenMP version described here always uses `numwfmem=2`.

---

## 8. Refactoring Recommendations

The existing code is scientifically mature but carries structural debt accumulated over 25+ years. The following changes would substantially improve testability, readability, and maintainability without requiring a complete rewrite.

### 8.1 Break up `com_mod` into domain-specific modules

`com_mod` currently holds ~500 lines of disparate global state. Partition it by concern:

| Proposed module | Contents from `com_mod` |
|---|---|
| `run_config_mod` | `ibdate`, `iedate`, `ldirect`, `ideltas`, `loutstep`, `lsynctime`, and all physics switches |
| `species_mod` | `nspec`, `species`, `decay`, `weta_gas`, `density`, `dquer`, all deposition species properties |
| `grid_config_mod` | `numxgrid`, `dxout`, `outlon0`, `numbnests`, `xglobal`, etc. |
| `receptor_state_mod` | `xreceptor`, `creceptor`, `numreceptor`, satellite arrays |
| `landuse_mod` | `landinvent`, `z0` |
| `run_state_mod` | `numpart`, `memtime`, `memind`, `lwindinterv` |

This makes each module's dependencies explicit and allows unit-testing physics routines by providing only the specific state they actually need.

### 8.2 Introduce explicit interfaces for physics routines

Currently `advance`, `turbulence_pbl`, `wetdepo`, etc. read from global module state rather than through arguments. Replace with explicit argument lists:

```fortran
! Current (implicit global state):
call turbulence_pbl   ! reads ust, ol, h from com_mod implicitly

! Proposed (explicit interface):
call turbulence_pbl(ust=ust_val, ol=ol_val, h=h_val, &
                    sigu=sigu_out, sigw=sigw_out, tlu=tlu_out, tlw=tlw_out)
```

This makes dependencies visible, enables unit testing, and removes `!$OMP THREADPRIVATE` fragility.

### 8.3 Move `read_options_and_initialise_flexpart` into its own module

The two long subroutines in `FLEXPART.f90` should live in an `init_mod` module. This follows the pattern already established by `initialise_mod` and makes the main program a clean 20-line dispatcher:

```fortran
program flexpart
  use init_mod, only: read_options_and_initialise
  use timemanager_mod, only: timemanager
  call read_options_and_initialise
  call timemanager
end program
```

### 8.4 Separate compile-time from runtime configuration

Several parameters in `par_mod` should be runtime options read from `COMMAND`:

| Currently `parameter` | Should be variable |
|---|---|
| `lusekerneloutput` | Already effectively a config choice |
| `numpf` | Selects precipitation interpolation scheme (#295) |
| `d_trop`, `d_strat` | Already `real` (not `parameter`) — just expose them |
| `fturbmeso` | Already `real` — add to COMMAND namelist |
| `ispeed` | Memory/speed trade-off tunable at runtime |

### 8.5 Encapsulate the wind-field memory manager

`getfields_mod`, `windfields_mod`, and the `memtime`/`memind` machinery in `com_mod` all participate in managing the two-slot met-field cache. This should be a single opaque type:

```fortran
type :: windfield_cache
  integer :: nslots
  integer :: times(numwfmem)
  integer :: ind(numwfmem)
  real, allocatable :: uueta(:,:,:,:)
  ! ...
contains
  procedure :: load_next => windcache_load_next
  procedure :: interpolate => windcache_interpolate
end type
```

This hides the slot-pointer logic and enables straightforward unit-testing of temporal interpolation.

### 8.6 Consolidate output backends behind a polymorphic interface

The current dispatch (`if (lnetcdfout.eq.1)`) duplicates the call structure in `output_mod`. A Fortran 2003 abstract interface removes the `if`-chains and allows a third backend (HDF5, Zarr, etc.) without touching existing code:

```fortran
type, abstract :: output_backend
contains
  procedure(write_header_if), deferred :: write_header
  procedure(write_conc_if),   deferred :: write_conc
  procedure(write_parts_if),  deferred :: write_particles
end type

type, extends(output_backend) :: netcdf_backend
  ! ...
end type
type, extends(output_backend) :: binary_backend
  ! ...
end type
```

### 8.7 Replace fixed-size arrays with allocatable in `par_mod`

`landinvent(1200,600,6)` (fixed global resolution) and the OH field dimensions `maxxOH=72, maxyOH=46` should be allocated at runtime from the actual input file dimensions. This eliminates hard-coded assumptions about input resolution.

### 8.8 Add a thin testing layer

Once physics routines accept explicit arguments (8.2), add a `tests/` directory with minimal driver programs (no met-files needed) that verify:
- `date_mod`: Julian↔Gregorian round-trip
- `settling_mod`: Stokes limit for 1 µm particles at STP
- `turbulence_pbl`: σ_w = 1.3 w* at z/h = 0.5 for free-convective conditions
- `interpol_mod`: bilinear interpolation returns exact values at grid points

These require zero external data and catch regressions during refactoring.

---

*Generated from source analysis of FLEXPART v11.0 (git: 62ae0fe), src/ directory. Last updated: 2026-04-03.*
