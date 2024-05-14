cmake_minimum_required(VERSION 3.20)

project(libkmeans LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Compiler flags for GCC only
# Note, flags "should" be specified in yaml presets, and passed to cmake as arguments
# Use c-preprocessor. Note, probably better to set target property `Fortran_PREPROCESS` instead
set(GCC_BASE 
   -cpp       
   )

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

# When using an interface target to assign compiler flags (or options), one requires a semicolon-separated string
# which is == strings appended to a list
# Convert the space-separated string into a semicolon-separated string
string(REPLACE " " ";" GCC_BASE "${GCC_BASE}")
string(REPLACE " " ";" GCC_DEBUG "${GCC_DEBUG}")
string(REPLACE " " ";" GCC_RELEASE "${GCC_RELEASE}")

# Interface to control scope of compiler flags
add_library(common_flags INTERFACE)
# Control the scope of compiler options
add_library(common_options INTERFACE)

# Build-specific flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
   target_compile_options(common_flags INTERFACE ${GCC_BASE} ${GCC_DEBUG})
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
   target_compile_options(common_flags INTERFACE ${GCC_BASE} ${GCC_RELEASE})
endif()


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
   find_package(MPI REQUIRED)

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
   target_compile_definitions(common_options INTERFACE USE_MPI)

else()
   set(MPI_Fortran_LIBRARIES "")
   message("-- Not using MPI")
   # Create an empty target with the same name
   add_library(MPI::MPI_Fortran INTERFACE IMPORTED)
endif()


# This library
# -----------------------------
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_library(lib_kmeans "")

set_target_properties(lib_kmeans
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME kmeans
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Linked libraries
target_link_libraries(lib_kmeans PRIVATE OpenMP::OpenMP_Fortran MPI::MPI_Fortran common_flags common_options)

# Source code
target_sources(lib_kmeans PUBLIC
       src/kmeans_m.F90
       src/grids_m.f90
       src/maths_m.f90
       src/mpi_m.F90
       )

# Fortuno wrapper layer
add_library(fortuno_wrapper "")
set_target_properties(fortuno_wrapper
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME fortunow
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )
target_link_libraries(fortuno_wrapper PUBLIC ${fortuno_libs} MPI::MPI_Fortran)
target_sources(fortuno_wrapper PRIVATE
       src/fortuno_interface_m.F90
       )

# Testing
# -------------------------------

# Testing program
add_executable(run_kmeans run.f90)
target_link_libraries(run_kmeans PRIVATE lib_kmeans OpenMP::OpenMP_Fortran MPI::MPI_Fortran common_flags common_options) 

# Library unit tests
enable_testing()

# Grids module
add_executable(test_grids tests/test_grids_m.f90)
target_link_libraries(test_grids PRIVATE lib_kmeans fortuno_wrapper common_flags common_options) 
add_test(NAME test_generate_real_space_grid_1D COMMAND test_grids "Real-space 1D grid")
add_test(NAME test_generate_real_space_grid_2D COMMAND test_grids "Real-space 2D grid")
add_test(NAME test_generate_real_space_grid_3D COMMAND test_grids "Real-space 3D grid")
add_test(NAME test_linspace COMMAND test_grids "Linear spaced sampling")
add_test(NAME test_linspace_to_grid2d COMMAND test_grids "2D linear grid")
add_test(NAME test_generate_gaussian_2d COMMAND test_grids "2D Gaussian")

# k-means module
add_executable(test_kmeans tests/test_kmeans_m.F90)
target_link_libraries(test_kmeans PRIVATE lib_kmeans fortuno_wrapper common_flags common_options) 
add_test(NAME test_assign_points_to_centroids COMMAND test_kmeans "Assign points to centroids")
add_test(NAME test_update_centroids_no_movement COMMAND test_kmeans "Update centroids - no change")
add_test(NAME test_points_are_converged COMMAND test_kmeans "Difference in two sets of points")


# Installation
# --------------------------------

# Populate `kmeans.pc.in` template with cmake variables, to generate pkg-config file
# configure_file(cmake/kmeans.pc.in kmeans.pc @ONLY)
# install(
#    FILES ${CMAKE_CURRENT_BINARY_DIR}/kmeans.pc
#    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
#    COMPONENT Fortuno_development
#  )