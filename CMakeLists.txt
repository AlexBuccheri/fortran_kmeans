#[==============================================================================================[
#                            CMakeLists.txt for Kmeans Project                                  #
#                                                                                               #
#]==============================================================================================]
cmake_minimum_required(VERSION 3.20)

# Note, the project version could also be stored in version text,
# configure_file(version.txt.in version.txt) or dynamically extracted from git tags
# See: https://github.com/LecrisUT/CMakeExtraUtils/blob/main/cmake/DynamicVersion.md
project(libkmeans VERSION 1.0.0 LANGUAGES Fortran)

# Set sub-directories in the build folder
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Load and execute CMake code from provided modules
include(GNUInstallDirs)
include(FetchContent)

# Add custom modules in cmake/
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#[==============================================================================================[
#                        Compiler Flags and Options                                             #
# Note, flags "should" be specified in cmake presets rather in CMakeLists.txt                   #
#]==============================================================================================]

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

#[==============================================================================================[
#                                   External Dependencies                                      #
#]==============================================================================================]
find_package(BLAS REQUIRED)    # Link with BLAS::BLAS
find_package(OpenMP REQUIRED)

# Fortuno unit testing framework
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

#[==============================================================================================[
#                                   Core Library                                                #
#]==============================================================================================]
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Define libraries
add_library(lib_kmeans)
add_library(fortuno_wrapper OBJECT)

# Location of source - added after defining library variables,
# as they're referenced in the source-level CMakeLists.txt
add_subdirectory(src)

set_target_properties(lib_kmeans PROPERTIES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    LIBRARY_OUTPUT_NAME "kmeans"
    # ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}" Not specified as will default to the GnuInstallDir
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

# Set where to look for module files during compilation
# Specifies how downstream projects consume `lib_kmeans`
target_include_directories(lib_kmeans PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
        "$<INSTALL_INTERFACE:modules>"    # Install-time path, relative to CMAKE_INSTALL_PREFIX
        )

# Linked libraries
target_link_libraries(lib_kmeans PRIVATE 
   OpenMP::OpenMP_Fortran 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran> 
   common_flags 
   )

set_target_properties(fortuno_wrapper PROPERTIES
        VERSION "${PROJECT_VERSION}"
        SOVERSION "${PROJECT_VERSION_MAJOR}"
        LIBRARY_OUTPUT_NAME fortunow
        LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        )

target_link_libraries(fortuno_wrapper PUBLIC 
   ${fortuno_libs} 
   $<$<BOOL:${MPI}>:MPI::MPI_Fortran>
   )

#[==============================================================================================[
#                                   Testing                                                     #
#]==============================================================================================]

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

#[==============================================================================================[
#                                   Installation                                                #
#]==============================================================================================]

# Install the target library and create export-set
install(
   TARGETS lib_kmeans
   EXPORT KmeansTargets
   LIBRARY
     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
   ARCHIVE
     DESTINATION "${CMAKE_INSTALL_LIBDIR}"
 )

 # Generate and install the export file
install(EXPORT KmeansTargets
   FILE KmeansTargets.cmake
   NAMESPACE kmeans::
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/kmeans"
   )

# Install the module files
# Note, this must be consistent with cmake/kmeans.pc.in
install(
   DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
   DESTINATION "${CMAKE_INSTALL_PREFIX}"    #Copies directory ${CMAKE_Fortran_MODULE_DIRECTORY}, to this location
   )

# Generate an export file directly for the build folder
# Allows to use find_package(kmeans) and point kmeans_ROOT=${binary-dir}
export(TARGETS lib_kmeans
   FILE "${CMAKE_CURRENT_BINARY_DIR}/kmeansTargets.cmake"
   NAMESPACE kmeans::
)

# cmake export files
include(CMakePackageConfigHelpers)

configure_package_config_file(cmake/kmeansConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/kmeansConfig.cmake"
                               INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/kmeans")

write_basic_package_version_file( 
    "${CMAKE_CURRENT_BINARY_DIR}/kmeansConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/kmeansConfigVersion.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/kmeansConfig.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/kmeans"
)

# Comment from Cristian on using .pc file:
# One requires that all target_link_libraries(PRIVATE) are not included here, and all target_link_libraries(PUBLIC/INTERFACE) are included. 
# This list depends on what features you enable at the build time and how find_package determined the flags necessary.
# It's a mess to support

# pkg-config file
# Populate `kmeans.pc.in` template with cmake variables
# Would need to chain configure_file (replaces @@) -> file(GENERATE) (replaces genex $<>)
# to use Crisitian's suggestion:
# configure_file(cmake/kmeans.pc.in.in kmeans.pc.in @ONLY)
# file(GENERATE OUTPUT kmeans.pc
#     INPUT ${CMAKE_CURRENT_BINARY_DIR}/kmeans.pc.in
# )
# Else:
# configure_file(cmake/kmeans.pc.in kmeans.pc @ONLY)
# install(
#    FILES ${CMAKE_CURRENT_BINARY_DIR}/kmeans.pc
#    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
#  )
