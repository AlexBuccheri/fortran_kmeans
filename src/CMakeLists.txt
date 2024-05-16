# Source code
# Must be private if exporting a target
# source files with implementation should be listed as lib_kmeans' INTERFACE_SOURCES
# INTERFACE_SOURCES are meant to be header files or other non-compiled sources that are needed by consumers of your library. 
target_sources(lib_kmeans PRIVATE
       kmeans_m.F90
       grids_m.f90
       maths_m.f90
       mpi_m.F90
       )

target_sources(fortuno_wrapper PRIVATE
       fortuno_interface_m.F90
       )
