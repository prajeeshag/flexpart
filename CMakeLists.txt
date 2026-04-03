cmake_minimum_required(VERSION 3.16)
project(FLEXPART Fortran)

# ── Build options ──────────────────────────────────────────────────────────────
option(USE_NETCDF "Enable NetCDF output"          ON)
option(USE_ETA    "Enable eta/ECMWF coordinates"  ON)
option(SERIAL     "Disable OpenMP (serial build)" OFF)
option(USE_MPI    "Enable MPI"                    OFF)
set(ARCH "native" CACHE STRING
    "CPU architecture target: native | generic (AVX2) | x86-64 (AVX) | skylake (AVX512)")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Release or Debug" FORCE)
endif()

# ── Module output directory (keeps .mod files out of source) ───────────────────
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# ── Compiler flags ─────────────────────────────────────────────────────────────
# -g is kept in all builds to match the original makefiles (always passed via FUSER)
set(_fp_flags "-m64 -mcmodel=large -g")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  string(APPEND _fp_flags " -cpp")

  if(ARCH STREQUAL "generic")
    string(APPEND _fp_flags " -march=core-avx2")
  elseif(ARCH STREQUAL "x86-64")
    string(APPEND _fp_flags " -march=x86-64")
  elseif(ARCH STREQUAL "skylake")
    string(APPEND _fp_flags " -march=skylake-avx512")
  else()  # native (default)
    string(APPEND _fp_flags " -march=native -mtune=native")
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-Og -g3 -ggdb3 -fconvert=little-endian -frecord-marker=4 \
-fmessage-length=0 -fbacktrace -Wall -fdump-core")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  string(APPEND _fp_flags " -fpp -convert little_endian -assume byterecl")

  if(ARCH STREQUAL "generic")
    string(APPEND _fp_flags " -xCORE-AVX2")
  elseif(ARCH STREQUAL "x86-64")
    string(APPEND _fp_flags " -xAVX")
  elseif(ARCH STREQUAL "skylake")
    string(APPEND _fp_flags " -xSKYLAKE-AVX512")
  else()  # native (default)
    string(APPEND _fp_flags " -xHost")
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -traceback -check bounds")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_fp_flags}")

# ── OpenMP ─────────────────────────────────────────────────────────────────────
if(NOT SERIAL)
  find_package(OpenMP REQUIRED Fortran)
endif()

# ── MPI ────────────────────────────────────────────────────────────────────────
if(USE_MPI)
  find_package(MPI REQUIRED Fortran)
endif()

# ── External libraries ─────────────────────────────────────────────────────────
find_library(ECCODES_LIB     eccodes     REQUIRED)
find_library(ECCODES_F90_LIB eccodes_f90 REQUIRED)

if(USE_NETCDF)
  find_library(NETCDFF_LIB netcdff REQUIRED)
endif()

# ── Source files ───────────────────────────────────────────────────────────────
# Order does not matter: CMake scans USE/MODULE statements and resolves
# Fortran module dependencies automatically.
set(SOURCES
  src/par_mod.f90
  src/com_mod.f90
  src/qvsat_mod.f90
  src/class_gribfile_mod.f90
  src/cmapf_mod.f90
  src/date_mod.f90
  src/erf_mod.f90
  src/mean_mod.f90
  src/sort_mod.f90
  src/random_mod.f90
  src/xmass_mod.f90
  src/particle_mod.f90
  src/point_mod.f90
  src/unc_mod.f90
  src/cbl_mod.f90
  src/pbl_profile_mod.f90
  src/turbulence_mod.f90
  src/windfields_mod.f90
  src/verttransform_mod.f90
  src/settling_mod.f90
  src/interpol_mod.f90
  src/drydepo_mod.f90
  src/advance_mod.f90
  src/plume_mod.f90
  src/outgrid_mod.f90
  src/wetdepo_mod.f90
  src/readoptions_mod.f90
  src/getfields_mod.f90
  src/flux_mod.f90
  src/txt_output_mod.f90
  src/binary_output_mod.f90
  src/receptor_mod.f90
  src/conv_mod.f90
  src/output_mod.f90
  src/restart_mod.f90
  src/initialise_mod.f90
  src/timemanager_mod.f90
  src/FLEXPART.f90
)

if(USE_NETCDF)
  list(APPEND SOURCES
    src/netcdf_output_mod.f90
    src/chemistry_mod.f90
    src/totals_mod.f90
    src/initdomain_mod.f90
    src/receptor_netcdf_mod.f90
    src/emissions_mod.f90
  )
endif()

if(USE_ETA)
  list(APPEND SOURCES src/coord_ecmwf_mod.f90)
endif()

# ── Executable name (mirrors original makefile conventions) ────────────────────
if(USE_NETCDF AND USE_ETA)
  set(EXE_NAME FLEXPART_ETA)
elseif(USE_NETCDF)
  set(EXE_NAME FLEXPART)
elseif(USE_ETA)
  set(EXE_NAME FLEXPART_ETA_BIN)
else()
  set(EXE_NAME FLEXPART_BIN)
endif()

add_executable(${EXE_NAME} ${SOURCES})

# ── Preprocessor definitions ───────────────────────────────────────────────────
if(USE_NETCDF)
  target_compile_definitions(${EXE_NAME} PRIVATE USE_NCF)
else()
  target_compile_options(${EXE_NAME} PRIVATE -UUSE_NCF)
endif()

if(USE_ETA)
  target_compile_definitions(${EXE_NAME} PRIVATE ETA)
else()
  target_compile_options(${EXE_NAME} PRIVATE -UETA)
endif()

if(USE_MPI)
  target_compile_definitions(${EXE_NAME} PRIVATE usempi)
else()
  target_compile_definitions(${EXE_NAME} PRIVATE useomp)
endif()

# ── Include / library paths from environment variables ────────────────────────
# Mirrors how the original makefiles consumed CPATH and LIBRARY_PATH
if(DEFINED ENV{CPATH})
  string(REPLACE ":" ";" _cpath_list "$ENV{CPATH}")
  target_include_directories(${EXE_NAME} PRIVATE ${_cpath_list})
endif()

if(DEFINED ENV{LIBRARY_PATH})
  string(REPLACE ":" ";" _libpath_list "$ENV{LIBRARY_PATH}")
  foreach(_lp ${_libpath_list})
    target_link_directories(${EXE_NAME} PRIVATE ${_lp})
    target_link_options(${EXE_NAME} PRIVATE "LINKER:-rpath,${_lp}")
  endforeach()
endif()

# ── Link libraries ─────────────────────────────────────────────────────────────
target_link_libraries(${EXE_NAME} PRIVATE
  ${ECCODES_LIB}
  ${ECCODES_F90_LIB}
  m
)

if(USE_NETCDF)
  target_link_libraries(${EXE_NAME} PRIVATE ${NETCDFF_LIB})
endif()

if(NOT SERIAL)
  target_link_libraries(${EXE_NAME} PRIVATE OpenMP::OpenMP_Fortran)
endif()

if(USE_MPI)
  target_link_libraries(${EXE_NAME} PRIVATE MPI::MPI_Fortran)
endif()

# ── Git version stamp ──────────────────────────────────────────────────────────
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} log --pretty=format:"%h %ad" -n 1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE _git_version
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(_git_version)
    # Write to a file in the build dir so it can be diffed on rebuild
    file(WRITE ${CMAKE_BINARY_DIR}/gitversion.txt "${_git_version}\n")
    message(STATUS "Git version: ${_git_version}")
  endif()
endif()

# ── Summary ────────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "FLEXPART build configuration:")
message(STATUS "  Compiler   : ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  Build type : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Executable : ${EXE_NAME}")
message(STATUS "  NetCDF     : ${USE_NETCDF}")
message(STATUS "  ETA coords : ${USE_ETA}")
message(STATUS "  OpenMP     : ${OPENMP_FOUND} (SERIAL=${SERIAL})")
message(STATUS "  MPI        : ${USE_MPI}")
message(STATUS "  Arch       : ${ARCH}")
message(STATUS "")
