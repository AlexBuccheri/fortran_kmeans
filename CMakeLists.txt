cmake_minimum_required(VERSION 3.20)

project(libkmeans LANGUAGES Fortran)

# Build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Compiler flags for GCC only
set(GCC_DEBUG
     -g                                           # Generate symbols
     -O0                                          # 0 optimisation. Probably implicit
     -fbacktrace                                  # symbolic stack traceback
     -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
     -finit-real=nan                              # All real scalars are initialised to NaN
     -fcheck=all                                  # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
     -Wall                                        # Enable all warnings
    )

set(GCC_RELEASE
    -O3
   )

# Set flags globally
set(CMAKE_Fortran_FLAGS "-fopenmp -cpp") 
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${GCC_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${GCC_RELEASE}")

# External dependencies
# ---------------------------------
find_package(BLAS REQUIRED)
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
   # Find MPI library and headers
   find_package(MPI REQUIRED)
   add_definitions(${MPI_Fortran_COMPILE_FLAGS})
   include_directories(${MPI_Fortran_INCLUDE_PATH})
   # Remove line as these are libraries, not paths to libraries
   # link_directories(${MPI_Fortran_LIBRARIES})

   # Fetch and assign fortuno MPI lib
   FetchContent_Declare(
     FortunoMPI
     GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno-mpi"
     GIT_TAG "main"
   )
   FetchContent_MakeAvailable(FortunoMPI)
   list(APPEND fortuno_libs Fortuno::fortuno_mpi)

   # Set preprocessor variable 
   add_compile_definitions(USE_MPI=1)

else()
   set(MPI_Fortran_LIBRARIES "")
   add_compile_definitions(USE_MPI=0)
endif()


# This library
# -----------------------------
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_library(lib_kmeans "")

set_target_properties(lib_kmeans
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME intmpi
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Linked libraries
target_link_libraries(lib_kmeans OpenMP::OpenMP_Fortran ${MPI_Fortran_LIBRARIES})  # BLAS::BLAS 

# Source code
target_sources(lib_kmeans PRIVATE
       src/kmeans_m.f90
       src/grids_m.f90
       src/maths_m.f90
       src/mpi_m.f90
       )

# Fortuno wrapper layer
add_library(fortuno_wrapper "")
set_target_properties(fortuno_wrapper
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME intmpi
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )
target_link_libraries(fortuno_wrapper ${MPI_Fortran_LIBRARIES} ${fortuno_libs})
target_sources(fortuno_wrapper PRIVATE
       src/fortuno_interface_m.F90
       )

# Testing
# -------------------------------

# Testing program
add_executable(run_kmeans run.f90)
# target_compile_options(run_kmeans PRIVATE -g -fcheck=all -fbacktrace -fopenmp)
target_link_libraries(run_kmeans PRIVATE lib_kmeans OpenMP::OpenMP_Fortran ${MPI_Fortran_LIBRARIES})  # BLAS::BLAS 

# Library unit tests
enable_testing()

# Grids module
add_executable(test_grids tests/test_grids_m.f90)
target_link_libraries(test_grids PRIVATE lib_kmeans fortuno_wrapper)

add_test(NAME test_generate_real_space_grid_1D COMMAND test_grids)
add_test(NAME test_generate_real_space_grid_2D COMMAND test_grids)
add_test(NAME test_generate_real_space_grid_3D COMMAND test_grids)
add_test(NAME test_linspace COMMAND test_grids)
add_test(NAME test_linspace_to_grid2d COMMAND test_grids)
add_test(NAME test_generate_gaussian_2d COMMAND test_grids)

# k-means module
add_executable(test_kmeans tests/test_kmeans_m.f90)
target_link_libraries(test_kmeans PRIVATE lib_kmeans fortuno_wrapper) 

add_test(NAME test_assign_points_to_centroids COMMAND test_kmeans)
add_test(NAME test_update_centroids_no_movement COMMAND test_kmeans)
add_test(NAME test_points_are_converged COMMAND test_kmeans)
