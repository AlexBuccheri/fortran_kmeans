cmake_minimum_required(VERSION 3.20)

project(libkmeans LANGUAGES Fortran)

# Build directory
set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Compiler flags
set(GCC_DEBUG
     -g                                           # Generate symbols
     -fbacktrace                                  # symbolic stack traceback
     -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
     -finit-real=nan                              #  All real scalars are initialised to NaN
     -fcheck=all                                  # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
    )

# External dependencies
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

# Fortuno unit testing framework - builds the serial version
include(FetchContent)
FetchContent_Declare(
  Fortuno
  GIT_REPOSITORY "https://github.com/fortuno-repos/fortuno"
  GIT_TAG "main"
)
FetchContent_MakeAvailable(Fortuno)

# This library
add_library(lib_kmeans "")

set_target_properties(lib_kmeans
        PROPERTIES
        VERSION 0.0.1
        SOVERSION 0.0.1
        LIBRARY_OUTPUT_NAME intmpi
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_Fortran_LIB_DIRECTORY}
        )

# Flags
target_compile_options(lib_kmeans PRIVATE -fopenmp)

# Linked libraries
target_link_libraries(lib_kmeans BLAS::BLAS OpenMP::OpenMP_Fortran ${MPI_Fortran_LIBRARIES})

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Source code
target_sources(lib_kmeans PRIVATE
       src/kmeans_m.f90
       src/utils_m.f90
       )

# Library unit tests
enable_testing()

add_executable(test_utils tests/test_utils_m.f90)
target_link_libraries(test_utils PRIVATE lib_kmeans Fortuno::fortuno_serial)
add_test(NAME test_generate_real_space_grid_1D COMMAND test_utils)
add_test(NAME test_generate_real_space_grid_2D COMMAND test_utils)
add_test(NAME test_generate_real_space_grid_3D COMMAND test_utils)

add_executable(test_kmeans tests/test_kmeans_m.f90)
target_link_libraries(test_kmeans PRIVATE lib_kmeans Fortuno::fortuno_serial)
add_test(NAME test_assign_points_to_centroids COMMAND test_kmeans)
