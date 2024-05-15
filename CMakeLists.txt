cmake_minimum_required(VERSION 3.20)

project(libkmeans LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Note, flags "should" be specified in cmake presets
# Gfortran flags
set(GCC_DEBUG
     -O0                                          # 0 optimisation. Probably implicit
     -fbacktrace                                  # symbolic stack traceback
     -fcheck=all                                  # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
     -Wall                                        # Enable all warnings
 #   -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
 #   -finit-real=nan                              # All real scalars are initialised to NaN
    )
                    
set(GCC_RELEASE
    -O3
   )

# When appending existing cmake flags (which one should not do), one converts list to string
# list(JOIN GCC_BASE " " GCC_BASE)
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GCC_BASE}")  

# Interface to control scope of compiler flags and options
add_library(common_flags INTERFACE)

# Generator expressions
# $<CONFIG:Release> converts to 0 or 1
# $<bool_value:some string> expands some string "if bool_value == true"

# Build-specific flags
target_compile_options(common_flags INTERFACE 
   $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:${GCC_DEBUG}>
   $<$<AND:$<CONFIG:Release>,$<Fortran_COMPILER_ID:GNU>>:${GCC_RELEASE}>
   )

# External dependencies
# ---------------------------------
find_package(BLAS REQUIRED)    # Link with BLAS::BLAS
find_package(OpenMP REQUIRED)

# Fortuno unit testing framework
include(FetchContent)

FetchContent_Declare(
  Fortuno
  GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
  GIT_TAG "main"
)
FetchContent_MakeAvailable(Fortuno)
set(fortuno_libs Fortuno::fortuno_serial)

# MPI-dependent settings
option(MPI "Use MPI" ON)
if(MPI)
   find_package(MPI REQUIRED Fortran)

   # Fetch and assign fortuno MPI lib
   FetchContent_Declare(
     FortunoMPI
     GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno-mpi"
     GIT_TAG "main"
   )
   FetchContent_MakeAvailable(FortunoMPI)
   list(APPEND fortuno_libs Fortuno::fortuno_mpi)

   # Set preprocessor variable 
   # add_compile_definitions(USE_MPI=1)
   target_compile_definitions(common_flags INTERFACE USE_MPI)

else()
   message("-- Not using MPI")
endif()

# This library
# -----------------------------
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Define libraries
add_library(lib_kmeans "")
add_library(fortuno_wrapper OBJECT "")

# Location of source - added after defining library variables,
# as they're referenced in the source-level CMakeLists.txt
add_subdirectory(src)

set_target_properties(lib_kmeans
        PROPERTIES
        VERSION ${libkmeans_VERSION} 
        SOVERSION ${libkmeans_VERSION_MAJOR}  
        LIBRARY_OUTPUT_NAME kmeans
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        )

# Linked libraries
target_link_libraries(lib_kmeans PRIVATE 
   OpenMP::OpenMP_Fortran 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran> 
   common_flags 
   )

set_target_properties(fortuno_wrapper
        PROPERTIES
        VERSION ${libkmeans_VERSION} 
        SOVERSION ${libkmeans_VERSION_MAJOR}  
        LIBRARY_OUTPUT_NAME fortunow
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        )

target_link_libraries(fortuno_wrapper PUBLIC 
   ${fortuno_libs} 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran>
   )

# Testing
# -------------------------------

# Testing program
add_executable(run_kmeans run.f90)
target_link_libraries(run_kmeans PRIVATE 
   lib_kmeans 
   OpenMP::OpenMP_Fortran 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran>  
   common_flags) 

# Library unit tests
enable_testing()

# Grids module
add_executable(test_grids tests/test_grids_m.f90)
target_link_libraries(test_grids PRIVATE lib_kmeans fortuno_wrapper common_flags) 
add_test(NAME test_generate_real_space_grid_1D COMMAND test_grids "Real-space 1D grid")
add_test(NAME test_generate_real_space_grid_2D COMMAND test_grids "Real-space 2D grid")
add_test(NAME test_generate_real_space_grid_3D COMMAND test_grids "Real-space 3D grid")
add_test(NAME test_linspace COMMAND test_grids "Linear spaced sampling")
add_test(NAME test_linspace_to_grid2d COMMAND test_grids "2D linear grid")
add_test(NAME test_generate_gaussian_2d COMMAND test_grids "2D Gaussian")

# k-means module
add_executable(test_kmeans tests/test_kmeans_m.F90)
target_link_libraries(test_kmeans PRIVATE lib_kmeans fortuno_wrapper common_flags) 
add_test(NAME test_assign_points_to_centroids COMMAND test_kmeans "Assign points to centroids")
add_test(NAME test_update_centroids_no_movement COMMAND test_kmeans "Update centroids - no change")
add_test(NAME test_points_are_converged COMMAND test_kmeans "Difference in two sets of points")
